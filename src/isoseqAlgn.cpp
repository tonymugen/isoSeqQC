/*
 * Copyright (c) 2024-2025 Anthony J. Greenberg and Rebekah Rogers
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Read isoSeq alignments and compare to genome annotations
/** \file
 * \author Anthony J. Greenberg and Rebekah Rogers
 * \copyright Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.2
 *
 * Implementation of classes that take `.bam` files with isoSeq alignments and identify potential fused transcripts.
 *
 */

#include <cstring>
#include <cmath>
#include <algorithm>
#include <memory>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <set>
#include <unordered_set>
#include <array>
#include <string>
#include <fstream>
#include <future>
#include <thread>

#include "hts.h"
#include "sam.h"
#include "bgzf.h"

#include "isoseqAlgn.hpp"
#include "helperFunctions.hpp"

using namespace isaSpace;

//ExonGroup constants
constexpr size_t ExonGroup::strandIDidx_{6UL};
constexpr char   ExonGroup::gffDelimiter_{'\t'};
constexpr size_t ExonGroup::spanStart_{3UL};
constexpr size_t ExonGroup::spanEnd_{4UL};

//ExonGroup methods
ExonGroup::ExonGroup(std::string geneName, const char strand, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSet) :
												geneName_{std::move(geneName)}, isNegativeStrand_{strand == '-'} { // only affirmatively negative strand is marked as such
	std::copy(
		exonSet.cbegin(),
		exonSet.cend(),
		std::back_inserter(exonRanges_)
	);
}

ExonGroup::ExonGroup(std::string geneName, const char strand, const std::vector<std::string> &exonLinesFomGFF) {
	std::set< std::pair<hts_pos_t, hts_pos_t> > uniqueExons;
	std::for_each(
		exonLinesFomGFF.cbegin(),
		exonLinesFomGFF.cend(),
		[&uniqueExons](const std::string &eachGFFline) {
			std::array<std::string, nGFFfields> gffFields{parseGFFline(eachGFFline)};
			if (gffFields.at(0) != "FAIL") {
				const auto exonStart = static_cast<hts_pos_t>( stol( gffFields.at(spanStart_) ) );
				const auto exonEnd   = static_cast<hts_pos_t>( stol( gffFields.at(spanEnd_) ) );
				uniqueExons.insert({exonStart, exonEnd});
			}
		}
	);
	*this = ExonGroup(std::move(geneName), strand, uniqueExons);
}

std::pair<hts_pos_t, hts_pos_t> ExonGroup::geneSpan() const noexcept {
	return ( exonRanges_.empty() ?
				std::pair<hts_pos_t, hts_pos_t>{-1, -1} :
				std::pair<hts_pos_t, hts_pos_t>{exonRanges_.front().first, exonRanges_.back().second} 
			);
}

std::pair<hts_pos_t, hts_pos_t> ExonGroup::firstExonSpan() const noexcept {
	if ( exonRanges_.empty() ) {
		return {-1, -1};
	}
	if (isNegativeStrand_) {
		return exonRanges_.back();
	}
	return exonRanges_.front();
}

hts_pos_t ExonGroup::firstExonLength() const noexcept {
	if ( exonRanges_.empty() ) {
		return 0;
	}
	if (isNegativeStrand_) {
		return exonRanges_.back().second - exonRanges_.back().first + 1;
	}
	return exonRanges_.front().second - exonRanges_.front().first + 1;
}

uint32_t ExonGroup::firstExonAfter(const hts_pos_t &position) const noexcept {
	const auto feaIt = std::lower_bound(
		exonRanges_.cbegin(),
		exonRanges_.cend(),
		position,
		[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &position) {
			return currExonRange.first < position;
		}
	);
	const auto result = std::min(
		static_cast<uint32_t>(exonRanges_.size() - 1UL),
		static_cast<uint32_t>( std::distance(exonRanges_.cbegin(), feaIt) )
	);
	return result;
}

uint32_t ExonGroup::firstOverlappingExon(const hts_pos_t &position) const noexcept {
	const auto feaIt = std::lower_bound(
		exonRanges_.cbegin(),
		exonRanges_.cend(),
		position,
		[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &position) {
			return currExonRange.second < position;
		}
	);
	const auto result = std::min(
		static_cast<uint32_t>(exonRanges_.size() - 1UL),
		static_cast<uint32_t>( std::distance(exonRanges_.cbegin(), feaIt) )
	);
	return result;
}

uint32_t ExonGroup::lastExonBefore(const hts_pos_t &position) const noexcept {
	const auto lebIt = std::lower_bound(
		exonRanges_.crbegin(),
		exonRanges_.crend(),
		position,
		[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &position) {
			return currExonRange.second > position;
		}
	);
	const auto result = std::min(
		static_cast<uint32_t>(exonRanges_.size() - 1UL),
		static_cast<uint32_t>( std::distance( lebIt, exonRanges_.crend() ) )
	);
	return result;
}

uint32_t ExonGroup::lastOverlappingExon(const hts_pos_t &position) const noexcept {
	const auto lebIt = std::lower_bound(
		exonRanges_.crbegin(),
		exonRanges_.crend(),
		position,
		[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &position) {
			return currExonRange.first > position;
		}
	);
	const auto result = std::min(
		static_cast<uint32_t>(exonRanges_.size() - 1UL),
		static_cast<uint32_t>( std::distance( lebIt, exonRanges_.crend() ) )
	);
	return result;
}

std::pair<hts_pos_t, hts_pos_t> ExonGroup::getFirstIntronSpan() const {
	std::pair<hts_pos_t, hts_pos_t> intronSpan{-1, -1};
	if (isNegativeStrand_) {
		const hts_pos_t endOfFirstExon = exonRanges_.back().first;
		const auto trueSecondExonIt = std::lower_bound(
			std::next( exonRanges_.crbegin() ),
			exonRanges_.crend(),
			endOfFirstExon,
			[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &firstExonEnd) {
				return currExonRange.second > firstExonEnd;
			}
		);
		if ( trueSecondExonIt == exonRanges_.crend() ) {
			return intronSpan;
		}
		intronSpan.first  = trueSecondExonIt->second;
		intronSpan.second = std::prev(trueSecondExonIt)->first;
		return intronSpan;
	}
	const hts_pos_t endOfFirstExon = exonRanges_.front().second;
	const auto trueSecondExonIt = std::lower_bound(
		std::next( exonRanges_.cbegin() ),
		exonRanges_.cend(),
		endOfFirstExon,
		[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &firstExonEnd) {
			return currExonRange.first < firstExonEnd;
		}
	);
	if ( trueSecondExonIt == exonRanges_.cend() ) {
		return intronSpan;
	}
	intronSpan.second = trueSecondExonIt->first;
	intronSpan.first  = std::prev(trueSecondExonIt)->second;
	return intronSpan;
}

std::vector<float> ExonGroup::getExonCoverageQuality(const BAMrecord &alignment) const {
	// find the first overlapping exon
	const auto exonRangeIt = std::lower_bound(
		exonRanges_.cbegin(),
		exonRanges_.cend(),
		alignment.getMapStart(),
		[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &position) {
			return currExonRange.second < position;
		}
	);
	const auto nLeadingUncoveredExons = std::distance(exonRanges_.cbegin(), exonRangeIt);
	// alignment quality for all uncovered exons is set to 0.0
	std::vector<float> qualityScores(nLeadingUncoveredExons, 0.0);

	// vector that tracks reference match/mismatch status for each reference position covered by CIGAR
	const std::vector<uint32_t> cigar{alignment.getCIGARvector()};
	std::vector<float> referenceMatchStatus{getReferenceMatchStatus(cigar)};
	std::for_each(
		exonRangeIt,
		exonRanges_.cend(),
		[&referenceMatchStatus, &alignment, &qualityScores](const std::pair<hts_pos_t, hts_pos_t> &currExonSpan) {
			const hts_pos_t realExonStart{std::max(alignment.getMapStart(), currExonSpan.first)};
			const hts_pos_t realExonLength{std::max(static_cast<hts_pos_t>(0), currExonSpan.second - realExonStart + 1)};
			const auto rmsStartPosition{
				std::min(
					static_cast<std::vector<float>::difference_type>( realExonStart - alignment.getMapStart() ),
					static_cast<std::vector<float>::difference_type>(referenceMatchStatus.size() - 1)
				)
			};
			const auto rmsSpan{
				std::min(
					static_cast<std::vector<float>::difference_type>(realExonLength),
					static_cast<std::vector<float>::difference_type>(referenceMatchStatus.size() - rmsStartPosition)
				)
			};
			const auto rmsBeginIt  = referenceMatchStatus.cbegin() + rmsStartPosition;
			const float matchCount = std::accumulate(rmsBeginIt, rmsBeginIt + rmsSpan, 0.0F);
			qualityScores.push_back( matchCount / static_cast<float>(currExonSpan.second - currExonSpan.first + 1) );
		}
	);

	if (isNegativeStrand_) {
		std::reverse( qualityScores.begin(), qualityScores.end() );
	}
	return qualityScores;
}

std::vector<float> ExonGroup::getBestExonCoverageQuality(const BAMrecord &alignment) const {
	if ( !alignment.isMapped() ) {
		std::vector<float> emptyResult;
		return emptyResult;
	}
	// vector that tracks reference match/mismatch status across all secondary alignments for each reference position covered by CIGAR
	const auto referenceMatchStatus{alignment.getBestReferenceMatchStatus()};
	// find the first overlapping exon
	const auto exonRangeIt = std::lower_bound(
		exonRanges_.cbegin(),
		exonRanges_.cend(),
		referenceMatchStatus.mapStart,
		[](const std::pair<hts_pos_t, hts_pos_t> &currExonRange, const hts_pos_t &position) {
			return currExonRange.second < position;
		}
	);
	const auto nLeadingUncoveredExons = std::distance(exonRanges_.cbegin(), exonRangeIt);
	// alignment quality for all uncovered exons is set to 0.0
	std::vector<float> qualityScores(nLeadingUncoveredExons, 0.0);

	std::for_each(
		exonRangeIt,
		exonRanges_.cend(),
		[&referenceMatchStatus, &qualityScores](const std::pair<hts_pos_t, hts_pos_t> &currExonSpan) {
			const hts_pos_t realExonStart{std::max(referenceMatchStatus.mapStart, currExonSpan.first)};
			const hts_pos_t realExonLength{std::max(static_cast<hts_pos_t>(0), currExonSpan.second - realExonStart + 1)};
			const auto rmsStartPosition{
				std::min(
					static_cast<std::vector<float>::difference_type>(realExonStart - referenceMatchStatus.mapStart),
					static_cast<std::vector<float>::difference_type>(referenceMatchStatus.matchStatus.size() - 1)
				)
			};
			const auto rmsSpan{
				std::min(
					static_cast<std::vector<float>::difference_type>(realExonLength),
					static_cast<std::vector<float>::difference_type>(referenceMatchStatus.matchStatus.size() - rmsStartPosition)
				)
			};
			const auto rmsBeginIt  = referenceMatchStatus.matchStatus.cbegin() + rmsStartPosition;
			const float matchCount = std::accumulate(rmsBeginIt, rmsBeginIt + rmsSpan, 0.0F);
			qualityScores.push_back( matchCount / static_cast<float>(currExonSpan.second - currExonSpan.first + 1) );
		}
	);

	if (isNegativeStrand_) {
		std::reverse( qualityScores.begin(), qualityScores.end() );
	}
	return qualityScores;
}

// ReadMatchWindowBIC methods
ReadMatchWindowBIC::ReadMatchWindowBIC(const std::vector< std::pair<float, hts_pos_t> >::const_iterator &windowBegin, const BinomialWindowParameters &windowParameters) :
						leftProbability_{windowParameters.currentProbability}, rightProbability_{windowParameters.alternativeProbability}, nTrials_{static_cast<float>(windowParameters.windowSize)} {
	kSuccesses_ = std::accumulate(
		windowBegin,
		windowBegin + windowParameters.windowSize,
		0.0F,
		[](float currentValue, const std::pair<float, hts_pos_t> &eachElement) {
			return currentValue + eachElement.first;
		}
	);
}

float ReadMatchWindowBIC::getBICdifference() const noexcept {
	const float jFailures = nTrials_ - kSuccesses_;
	const float leftLlik  = ( kSuccesses_ * logf(leftProbability_) ) + ( jFailures * logf(1.0F - leftProbability_) );
	const float rightLlik = ( kSuccesses_ * logf(rightProbability_) ) + ( jFailures * logf(1.0F - rightProbability_) );
    return logf(nTrials_) + ( 2.0F * (leftLlik - rightLlik) );
};

// BAMrecord methods
constexpr uint16_t BAMrecord::sequenceMask_{0x00FF};
constexpr uint16_t BAMrecord::qualityShift_{8};
constexpr uint16_t BAMrecord::suppSecondaryAlgn_{BAM_FSECONDARY | BAM_FSUPPLEMENTARY};
constexpr std::array<hts_pos_t, 10> BAMrecord::queryConsumption_{
	1, 1, 0, 0, 1,
	0, 0, 1, 1, 0
};
constexpr std::array<hts_pos_t, 10> BAMrecord::referenceConsumption_{
	1, 0, 1, 1, 0,
	0, 0, 1, 1, 0
};
constexpr std::array<float, 10> BAMrecord::sequenceMatch_{
	1.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 1.0, 0.0, 0.0
};

BAMrecord::BAMrecord(const bam1_t *alignmentRecord, const sam_hdr_t *samHeader) :
			isRev_{bam_is_rev(alignmentRecord)}, mapStart_{alignmentRecord->core.pos + 1}, mapEnd_{bam_endpos(alignmentRecord) + 1},
			isMapped_{ (alignmentRecord->core.flag & BAM_FUNMAP) == 0 } {

	if ( (alignmentRecord->core.flag & suppSecondaryAlgn_) != 0 ) {
		throw std::string("ERROR: source BAM alignment record is not primary in ")
			+ std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	const auto *const refNamePtr = sam_hdr_tid2name(samHeader, alignmentRecord->core.tid); 
	if (refNamePtr != nullptr) {
		referenceName_ = std::string{refNamePtr};
	} else {
		referenceName_ = "Unknown";
	}
	readName_ = std::string{bam_get_qname(alignmentRecord)};
	cigar_    = std::vector<uint32_t>(
						bam_get_cigar(alignmentRecord),
						bam_get_cigar(alignmentRecord) + alignmentRecord->core.n_cigar
					);
	for (int32_t iSeq = 0; iSeq < alignmentRecord->core.l_qseq; ++iSeq) {
		const uint16_t qualityByte{*(bam_get_qual(alignmentRecord) + iSeq)};
		const uint16_t sequenceByte{static_cast<uint16_t>( bam_seqi(bam_get_seq(alignmentRecord), iSeq) )};
		sequenceAndQuality_.push_back( (qualityByte << qualityShift_) | sequenceByte );
	}
}

void BAMrecord::addSecondaryAlignment(const bam1_t *alignmentRecord, const sam_hdr_t *samHeader, const hts_pos_t localWindow) {
	if ( (alignmentRecord->core.flag & suppSecondaryAlgn_) == 0 ) {
		throw std::string("ERROR: source BAM alignment record is not secondary in ")
			+ std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const auto currentReadName = std::string{bam_get_qname(alignmentRecord)};
	if (currentReadName != readName_) {
		return;
	}
	const auto currentReferenceName = std::string( sam_hdr_tid2name(samHeader, alignmentRecord->core.tid) );
	if (referenceName_ != currentReferenceName) {
		++secondaryAlignmentCount_;
		return;
	}

	const bool currentIsRev{bam_is_rev(alignmentRecord)};
	BAMsecondary newSecondary;
	newSecondary.mapStart            = alignmentRecord->core.pos + 1;
	newSecondary.mapEnd              = bam_endpos(alignmentRecord) + 1;
	newSecondary.sameStrandAsPrimary = (isRev_ == currentIsRev);
    newSecondary.cigar               = std::move( std::vector<uint32_t>(bam_get_cigar(alignmentRecord), bam_get_cigar(alignmentRecord) + alignmentRecord->core.n_cigar) );

	std::pair<hts_pos_t, hts_pos_t> primaryRange(std::max(0L, mapStart_ - localWindow), mapEnd_ + localWindow);
	std::pair<hts_pos_t, hts_pos_t> secondaryRange(newSecondary.mapStart, newSecondary.mapEnd);
	if ( rangesOverlap(primaryRange, secondaryRange) ) {
		localSecondaryAlignments_.emplace_back( std::move(newSecondary) );
	}
	++secondaryAlignmentCount_;
}

uint16_t BAMrecord::localReversedSecondaryAlignmentCount() const noexcept {
	return std::count_if(
		localSecondaryAlignments_.cbegin(),
		localSecondaryAlignments_.cend(),
		[](const BAMsecondary &secondary) {
			return !secondary.sameStrandAsPrimary;
		}
	);
}

std::string BAMrecord::getCIGARstring() const {
	auto stringify = [](std::string currString, uint32_t cigarElement) {
		return std::move(currString)
			+ std::to_string( bam_cigar_oplen(cigarElement) )
			+ bam_cigar_opchr(cigarElement);
	};
	if (isRev_) {
		std::string cigar = std::accumulate(
			cigar_.crbegin(),
			cigar_.crend(),
			std::string(),
			stringify
		);
		return cigar;
	}
	std::string cigar = std::accumulate(
		cigar_.cbegin(),
		cigar_.cend(),
		std::string(),
		stringify
	);

	return cigar;
}

uint32_t BAMrecord::getFirstCIGAR() const noexcept {
	if ( cigar_.empty() ) {
		return 0;
	}
	return ( this->isRev_ ? cigar_.back() : cigar_.front() );
}

MappedReadMatchStatus BAMrecord::getBestReferenceMatchStatus() const {
	// first element of the pair will have the map start position for the read
	std::vector< std::pair< hts_pos_t, std::vector<float> > > matchVectors;
	matchVectors.emplace_back( mapStart_, getReferenceMatchStatus(cigar_) );
	std::for_each(
		localSecondaryAlignments_.cbegin(),
		localSecondaryAlignments_.cend(),
		[&matchVectors](const BAMsecondary &localSecondary) {
			if (localSecondary.sameStrandAsPrimary) {
				matchVectors.emplace_back( localSecondary.mapStart, getReferenceMatchStatus(localSecondary.cigar) );
			}
		}
	);
	std::sort(
		matchVectors.begin(),
		matchVectors.end(),
		[](const std::pair< hts_pos_t, std::vector<float> > &firstPair, const std::pair< hts_pos_t, std::vector<float> > &secondPair) {
			return firstPair.first < secondPair.first;
		}
	);
	MappedReadMatchStatus bestMatch;
	bestMatch.mapStart     = matchVectors.front().first;
	bestMatch.matchStatus  = std::move(matchVectors.front().second);
	std::for_each(
		std::next( matchVectors.cbegin() ),
		matchVectors.cend(),
		[&bestMatch](const std::pair< hts_pos_t, std::vector<float> > &eachSecondaryMatch) {
			auto bestBeginIt = bestMatch.matchStatus.begin();
			// figure out where in the best match vector the current match status starts
			// this may be past the best match end
			std::advance(bestBeginIt,
				std::min(
					eachSecondaryMatch.first - bestMatch.mapStart,
					std::distance( bestMatch.matchStatus.begin(), bestMatch.matchStatus.end() )
				)
			);
			auto bestEndIt = bestBeginIt;
			std::advance(bestEndIt,
				std::min(
					std::distance( bestBeginIt, bestMatch.matchStatus.end() ),
					std::distance( eachSecondaryMatch.second.cbegin(), eachSecondaryMatch.second.cend() )
				) 
			);
			auto currentMatchIt = eachSecondaryMatch.second.cbegin();
			std::for_each(
				bestBeginIt,
				bestEndIt,
				[&currentMatchIt](float &bestElement) {
					bestElement = std::max(bestElement, *currentMatchIt);
					++currentMatchIt;
				}
			);
			const size_t gapToSecondary = std::max(
				static_cast<hts_pos_t>(0),
				eachSecondaryMatch.first - static_cast<hts_pos_t>( bestMatch.mapStart + bestMatch.matchStatus.size() )
			);
			bestMatch.matchStatus.resize(gapToSecondary, 0.0F);

			// if the current match status vector extends beyond the current best match
			// if not, currentMatchIt will be the end iterator and nothing happens
			std::copy( currentMatchIt, eachSecondaryMatch.second.cend(), std::back_inserter(bestMatch.matchStatus) );
		}
	);
	return bestMatch;
}

std::vector< std::pair<float, hts_pos_t> > BAMrecord::getReadCentricMatchStatus() const {
	std::vector< std::pair<float, hts_pos_t> > readMatchStatus;
	hts_pos_t referencePosition{mapStart_};
	for (const auto &eachCIGAR : cigar_) {
		if (queryConsumption_.at( bam_cigar_op(eachCIGAR) ) == 0) {
			hts_pos_t iSIGLEN{0};
			const hts_pos_t rcStatus{referenceConsumption_.at( bam_cigar_op(eachCIGAR) )};
			while ( iSIGLEN < bam_cigar_oplen(eachCIGAR) ) {
				referencePosition += rcStatus;
				iSIGLEN++;
			}
			continue;
		}
		hts_pos_t iSIGLEN{0};
		const float matchStatus{sequenceMatch_.at( bam_cigar_op(eachCIGAR) )};
		const hts_pos_t rcStatus{referenceConsumption_.at( bam_cigar_op(eachCIGAR) )};
		while ( iSIGLEN < bam_cigar_oplen(eachCIGAR) ) {
			readMatchStatus.emplace_back(matchStatus, referencePosition);
			referencePosition += rcStatus;
			iSIGLEN++;
		}
	}

	return readMatchStatus;
}

std::vector<MappedReadInterval> BAMrecord::getPoorlyMappedRegions(const BinomialWindowParameters &windowParameters) const {
	std::vector<MappedReadInterval> result;
	const auto readMatchStatus{this->getReadCentricMatchStatus()};
	if ( (windowParameters.windowSize < 1) ||
			( readMatchStatus.empty() ) ||
			(readMatchStatus.size() < 2 * windowParameters.windowSize) ) {
		return result;
	}

	const float unmappedProbability{std::min(windowParameters.currentProbability, windowParameters.alternativeProbability)};
	const float mappedProbability{std::max(windowParameters.currentProbability, windowParameters.alternativeProbability)};
	const auto actualWindowSize{std::min( windowParameters.windowSize, std::distance( readMatchStatus.cbegin(), readMatchStatus.cend() ) )};

	auto windowBeginIt = readMatchStatus.cbegin();
	BinomialWindowParameters currentWindowParameters;
	currentWindowParameters.currentProbability     = unmappedProbability;
	currentWindowParameters.alternativeProbability = mappedProbability;
	currentWindowParameters.windowSize             = actualWindowSize;
	ReadMatchWindowBIC initialWindowBIC(windowBeginIt, currentWindowParameters);
	const float initialBICdiff{initialWindowBIC.getBICdifference()};

	std::vector<float> windowBICdiffs;
	windowBICdiffs.push_back(initialBICdiff);
	const auto lastWindowIt = readMatchStatus.cend() - actualWindowSize;
	while (windowBeginIt != lastWindowIt) {
		ReadMatchWindowBIC currentWindowBIC(windowBeginIt, currentWindowParameters);
		windowBICdiffs.push_back( currentWindowBIC.getBICdifference() );
		++windowBeginIt;
	}

	// Calculate BIC differences with a window size lag. Peaks and valleys will correspond to change points.
	const auto actualBICwindowSize{std::min( actualWindowSize, std::distance( windowBICdiffs.cbegin(), windowBICdiffs.cend() ) )};
	const auto lastBICwindowIt = windowBICdiffs.cend() - actualBICwindowSize;
	std::vector<float> windowDeltaBICdiffs;
	for (auto windowBICdiffsIt = windowBICdiffs.begin(); windowBICdiffsIt != lastBICwindowIt; ++windowBICdiffsIt) {
		windowDeltaBICdiffs.push_back( *windowBICdiffsIt - *std::next(windowBICdiffsIt, actualWindowSize) );
	}

	// find change points
	const auto bicDiffPeaks{getPeaks(windowDeltaBICdiffs, windowParameters.bicDifferenceThreshold)};
	const auto bicDiffValleys{getValleys(windowDeltaBICdiffs, -windowParameters.bicDifferenceThreshold)};

	// use the change points to find misaligned regions
	auto peaksIt               = bicDiffPeaks.cbegin();
	auto valleysIt             = bicDiffValleys.cbegin();
	auto currentPeakDistance   = std::distance(windowDeltaBICdiffs.cbegin(), *peaksIt);
	auto currentValleyDistance = std::distance(windowDeltaBICdiffs.cbegin(), *valleysIt);
	auto minDistance           = std::min(currentPeakDistance, currentValleyDistance);
	MappedReadInterval currentInterval;
	currentInterval.readStart      = 0;
	currentInterval.readEnd        = this->getReadLength();
	currentInterval.referenceStart = this->getMapStart();
	currentInterval.referenceEnd   = this->getMapEnd();
	while ( minDistance < windowDeltaBICdiffs.size() ) {
		if (minDistance == currentPeakDistance) {
			currentInterval.readEnd      = minDistance + actualBICwindowSize - 1;                         // subtracting 1 to not get past the end
			currentInterval.referenceEnd = readMatchStatus[minDistance + actualBICwindowSize - 1].second;
			result.push_back(currentInterval);
			std::advance(peaksIt, 1);
			currentPeakDistance       = std::distance(windowDeltaBICdiffs.cbegin(), *peaksIt);
			minDistance               = std::min(currentPeakDistance, currentValleyDistance);
			currentInterval.readStart = 0; // reset to make sure the post-loop test works
			continue;
		}
		currentInterval.readStart      = minDistance + actualBICwindowSize - 1;
		currentInterval.referenceStart = readMatchStatus[minDistance + actualBICwindowSize - 1].second;
		std::advance(valleysIt, 1);
		currentValleyDistance = std::distance(windowDeltaBICdiffs.cbegin(), *valleysIt);
		minDistance           = std::min(currentPeakDistance, currentValleyDistance);
	}

	// if the last change point goes to misalignment
	if (currentInterval.readStart != 0) {
		currentInterval.readEnd      = this->getReadLength();
		currentInterval.referenceEnd = this->getMapEnd();
		result.push_back(currentInterval);
	}

	// if no change points, return an empty vector
	return result;
}

// BAMtoGenome methods
constexpr uint16_t BAMtoGenome::suppSecondaryAlgn_{BAM_FSECONDARY | BAM_FSUPPLEMENTARY};

BAMtoGenome::BAMtoGenome(const BamAndGffFiles &bamGFFfilePairNames) {
	const auto gffExonGroups{parseGFF(bamGFFfilePairNames.gffFileName)};

	if ( gffExonGroups.empty() ) {
		throw std::string("ERROR: no mRNAs with exons found in the ")
			+ bamGFFfilePairNames.gffFileName + std::string(" GFF file in ")
			+ std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	constexpr char openMode{'r'};
	std::unique_ptr<BGZF, void(*)(BGZF *)> bamFile(
		bgzf_open(bamGFFfilePairNames.bamFileName.c_str(), &openMode),
		[](BGZF *bamFile) {
			bgzf_close(bamFile);
		}
	);
	if (bamFile == nullptr) {
		throw std::string("ERROR: failed to open the BAM file ")
			+ bamGFFfilePairNames.bamFileName + std::string(" in ")
			+ std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	// must read the header first to get to the alignments
	std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> bamHeader(
		bam_hdr_read( bamFile.get() ),
		[](sam_hdr_t *samHeader) {
			sam_hdr_destroy(samHeader);
		}
	);
	if (bamHeader == nullptr) {
		throw std::string("ERROR: failed to read the header from the BAM file ")
			+ bamGFFfilePairNames.bamFileName + std::string(" in ")
			+ std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	std::vector<std::string> outputLines; // gene information, if any, for each alignment to save
	// keeping track of the latest iterators pointing to identified exon groups, one per chromosome/reference sequence
	std::unordered_map<std::string, std::vector<ExonGroup>::const_iterator> latestExonGroupIts;
	for (const auto &eachChromosome : gffExonGroups) {
		latestExonGroupIts[eachChromosome.first] = eachChromosome.second.cbegin();
	}
	while (true) {
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		const auto nBytes = bam_read1( bamFile.get(), bamRecordPtr.get() );
		if (nBytes == -1) {
			break;
		}
		if (nBytes < -1) {
			continue;
		}
		// Is this a secondary alignment?
		if ( ( (bamRecordPtr->core.flag & suppSecondaryAlgn_) != 0 ) ) {
			if ( !readsAndExons_.empty() ) {
				readsAndExons_.back().first.addSecondaryAlignment( bamRecordPtr.get(), bamHeader.get() );
			}
			continue;
		}
		BAMrecord currentBAM( bamRecordPtr.get(), bamHeader.get() );
		std::string referenceName{currentBAM.getReferenceName()};
		const char strandID = (currentBAM.isRevComp() ? '-' : '+');
		referenceName.push_back(strandID);
		// BAM chromosome not in GFF
		if ( gffExonGroups.find(referenceName) == gffExonGroups.end() ) {
			ExonGroup empty;
			readsAndExons_.emplace_back(currentBAM, empty);
			continue;
		}
		latestExonGroupIts.at(referenceName) = findOverlappingGene_(gffExonGroups.at(referenceName), latestExonGroupIts.at(referenceName), currentBAM);
	}
}

size_t BAMtoGenome::nChromosomes() const noexcept {
	std::unordered_set<std::string> chromosomes;
	std::for_each(
		readsAndExons_.cbegin(),
		readsAndExons_.cend(),
		[&chromosomes](const std::pair<BAMrecord, ExonGroup> &currentReadEG) {
			std::string chrName = currentReadEG.first.getReadName();
			chromosomes.insert(chrName);
		}
	);

	return chromosomes.size();
}

size_t BAMtoGenome::nExonSets() const noexcept {
	return std::accumulate(
		readsAndExons_.cbegin(),
		readsAndExons_.cend(),
		0UL,
		[](const size_t &sum, const std::pair<BAMrecord, ExonGroup> &currentReadEG) {
			return sum + static_cast<size_t>(currentReadEG.second.nExons() > 0); 
		}
	);
}

void BAMtoGenome::saveReadCoverageStats(const std::string &outFileName, const size_t &nThreads) const {
	size_t actualNthreads = std::min( nThreads, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	actualNthreads        = std::min( actualNthreads, readsAndExons_.size() );
	actualNthreads        = std::max(actualNthreads, 1UL);

	const auto threadRanges{makeThreadRanges(readsAndExons_, actualNthreads)};
	std::vector<std::string> threadOutStrings(actualNthreads);
	std::vector< std::future<void> > tasks;
	tasks.reserve(actualNthreads);
	size_t iThread{0};
	std::for_each(
		threadRanges.cbegin(),
		threadRanges.cend(),
		[&iThread, &tasks, &threadOutStrings](const std::pair<bamGFFvector::const_iterator, bamGFFvector::const_iterator> &eachRange) {
			tasks.emplace_back(
				std::async(
					[iThread, eachRange, &threadOutStrings] {
						threadOutStrings.at(iThread) = stringifyAlignmentRange(eachRange.first, eachRange.second);
					}
				)
			);
			++iThread;
		}
	);
	for (const auto &eachThread : tasks) {
		eachThread.wait();
	}

	const std::string headerLine = "read_name\tchromosome\tstrand\talignment_start\talignment_end\tbest_alignment_start\tbest_alignment_end\t"
									"first_soft_clip_length\tn_secondary_alignments\tn_good_secondary_alignments\tn_local_reversed_strand\t"
									"gene_name\tn_exons\tfirst_exon_length\tfirst_exon_start\t" 
									"last_exon_end\talignment_quality_string\tbest_alignment_quality_string\n";
	std::fstream outStream;
	outStream.open(outFileName, std::ios::out | std::ios::binary | std::ios::trunc);
	outStream.write( headerLine.c_str(), static_cast<std::streamsize>( headerLine.size() ) );
	for (const auto &eachThreadString : threadOutStrings) {
		outStream.write( eachThreadString.c_str(), static_cast<std::streamsize>( eachThreadString.size() ) );
	}
	outStream.close();
};

void BAMtoGenome::saveUnmappedRegions(const std::string &outFileName, const BinomialWindowParameters &windowParameters, const size_t &nThreads) const {
	size_t actualNthreads = std::min( nThreads, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	actualNthreads        = std::min( actualNthreads, readsAndExons_.size() );
	actualNthreads        = std::max(actualNthreads, 1UL);

	const auto threadRanges{makeThreadRanges(readsAndExons_, actualNthreads)};
	std::vector<std::string> threadOutStrings(actualNthreads);
	std::vector< std::future<void> > tasks;
	tasks.reserve(actualNthreads);
	size_t iThread{0};
	std::for_each(
		threadRanges.cbegin(),
		threadRanges.cend(),
		[&iThread, &tasks, &windowParameters, &threadOutStrings](const std::pair<bamGFFvector::const_iterator, bamGFFvector::const_iterator> &eachRange) {
			tasks.emplace_back(
				std::async(
					[iThread, eachRange, &threadOutStrings, &windowParameters] {
						threadOutStrings.at(iThread) = stringifyUnmappedRegions(eachRange.first, eachRange.second, windowParameters);
					}
				)
			);
			++iThread;
		}
	);
	for (const auto &eachThread : tasks) {
		eachThread.wait();
	}

	const std::string headerLine = "read_name\tread_length\tunmapped_start\tunmapped_end\n";
	std::fstream outStream;
	outStream.open(outFileName, std::ios::out | std::ios::binary | std::ios::trunc);
	outStream.write( headerLine.c_str(), static_cast<std::streamsize>( headerLine.size() ) );
	for (const auto &eachThreadString : threadOutStrings) {
		outStream.write( eachThreadString.c_str(), static_cast<std::streamsize>( eachThreadString.size() ) );
	}
	outStream.close();
};

std::vector<ExonGroup>::const_iterator BAMtoGenome::findOverlappingGene_(const std::vector<ExonGroup> &chromosomeExonGroups,
		const std::vector<ExonGroup>::const_iterator &exonGroupSearchStart, BAMrecord &alignedRead) {
	auto searchIt = exonGroupSearchStart;

	// if the BAM file is not sorted, we may have to backtrack
	// looking for the first gene end that is after the read map start
	if (alignedRead.getMapStart() < exonGroupSearchStart->geneSpan().first) {
		auto reverseLEGI = std::make_reverse_iterator(searchIt);
		reverseLEGI      = std::lower_bound(
			reverseLEGI,
			chromosomeExonGroups.crend(),
			alignedRead.getMapStart(),
			[](const ExonGroup &currGroup, const hts_pos_t bamStart) {
				return  (currGroup.geneSpan().first > bamStart);
			}
		);
		if ( ( reverseLEGI == chromosomeExonGroups.crend() ) || ( reverseLEGI->geneSpan().second < alignedRead.getMapStart() ) ) {
			ExonGroup emptyGroup;
			readsAndExons_.emplace_back(alignedRead, emptyGroup);
			// if we are back past the first mRNA or the read does not overlap a gene, do not update the latest iterator
			return searchIt;
		}
		// convert back to the forward iterator using the reverse/forward relationship
		// safe to increment the reverse iterator here due to the test above
		searchIt = std::next(reverseLEGI).base();
		readsAndExons_.emplace_back(alignedRead, *searchIt);
		return searchIt;
	}
	searchIt = std::lower_bound(
		exonGroupSearchStart,
		chromosomeExonGroups.cend(),
		alignedRead.getMapStart(),
		[](const ExonGroup &currGroup, const hts_pos_t bamStart) {
			return  (currGroup.geneSpan().second < bamStart);
		}
	);
	if ( searchIt == chromosomeExonGroups.cend() ) {
		ExonGroup emptyGroup;
		readsAndExons_.emplace_back(alignedRead, emptyGroup);
		// return the iterator to a valid record
		std::advance(searchIt, -1);
		return searchIt;
	}
	if (alignedRead.getMapStart() < searchIt->geneSpan().first) { // no overlap with a known gene
		ExonGroup emptyGroup;
		readsAndExons_.emplace_back(alignedRead, emptyGroup);
		// do not update the iterator if no overlap
		return exonGroupSearchStart;
	}
	readsAndExons_.emplace_back(alignedRead, *searchIt);
	return searchIt;
};
