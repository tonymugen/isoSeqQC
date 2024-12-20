/*
 * Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
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
 * \version 0.1
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

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"

#include "isoseqAlgn.hpp"
#include "helperFunctions.hpp"

using namespace isaSpace;

ExonGroup::ExonGroup(std::string geneName, const char strand, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSet) :
												geneName_{std::move(geneName)}, isNegativeStrand_{strand == '-'} { // only affirmatively negative strand is marked as such
	std::copy(
		exonSet.cbegin(),
		exonSet.cend(),
		std::back_inserter(exonRanges_)
	);
	firstExonIt_ = ( isNegativeStrand_ ? std::prev( exonRanges_.end() ) : exonRanges_.begin() );
}

 std::pair<hts_pos_t, hts_pos_t> ExonGroup::geneSpan() const noexcept {
	return ( exonRanges_.empty() ?
				std::pair<hts_pos_t, hts_pos_t>{-1, -1} :
				std::pair<hts_pos_t, hts_pos_t>{exonRanges_.front().first, exonRanges_.back().second} 
			);
}

std::pair<hts_pos_t, hts_pos_t> ExonGroup::firstExonSpan() const noexcept {
	return ( exonRanges_.empty() ? std::pair<hts_pos_t, hts_pos_t>(-1, -1): *firstExonIt_);
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
	const float leftLlik  = kSuccesses_ * logf(leftProbability_)  + jFailures * logf(1.0F - leftProbability_);
	const float rightLlik = kSuccesses_ * logf(rightProbability_) + jFailures * logf(1.0F - rightProbability_);
    return logf(nTrials_) + 2.0F * (leftLlik - rightLlik);
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

	readName_      = std::string{bam_get_qname(alignmentRecord)};
	referenceName_ = std::string( sam_hdr_tid2name(samHeader, alignmentRecord->core.tid) );
	cigar_         = std::vector<uint32_t>(
						bam_get_cigar(alignmentRecord),
						bam_get_cigar(alignmentRecord) + alignmentRecord->core.n_cigar
					);
	for (int32_t iSeq = 0; iSeq < alignmentRecord->core.l_qseq; ++iSeq) {
		const uint16_t qualityByte{*(bam_get_qual(alignmentRecord) + iSeq)};
		const uint16_t sequenceByte{bam_seqi(bam_get_seq(alignmentRecord), iSeq)};
		sequenceAndQuality_.push_back( (qualityByte << qualityShift_) | sequenceByte );
	}
}

void BAMrecord::addSecondaryAlignment(const bam1_t *alignmentRecord, const sam_hdr_t *samHeader) {
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

	std::pair<hts_pos_t, hts_pos_t> primaryRange(mapStart_, mapEnd_);
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
	hts_pos_t lastMapStart = matchVectors.front().first;
	std::for_each(
		std::next( matchVectors.cbegin() ),
		matchVectors.cend(),
		[&bestMatch, &lastMapStart](const std::pair< hts_pos_t, std::vector<float> > &eachPair) {
			auto bestBeginIt = bestMatch.matchStatus.begin();
			std::advance(bestBeginIt, eachPair.first - lastMapStart); // the advance must be non-negative because the vector is sorted
			auto bestEndIt = bestMatch.matchStatus.begin();
			std::advance(bestEndIt,
				std::min(
					std::distance( bestBeginIt, bestMatch.matchStatus.end() ),
					std::distance( eachPair.second.cbegin(), eachPair.second.cend() )
				) 
			);
			auto currentIt = eachPair.second.cbegin();
			std::for_each(
				bestBeginIt,
				bestEndIt,
				[&currentIt](float &bestElement) {
					bestElement = std::max(bestElement, *currentIt);
					++currentIt;
				}
			);
			std::copy( currentIt, eachPair.second.cend(), std::back_inserter(bestMatch.matchStatus) );
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
	if ( (windowParameters.windowSize < 1) || ( readMatchStatus.empty() ) ) {
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
	std::for_each(
		gffExonGroups.cbegin(),
		gffExonGroups.cend(),
		[&latestExonGroupIts](const std::pair<std::string, std::vector<ExonGroup> > &eachChromosome){
			latestExonGroupIts[eachChromosome.first] = eachChromosome.second.cbegin();
		}
	);
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
		if ( ( (bamRecordPtr->core.flag & suppSecondaryAlgn_) != 0 ) && !readsAndExons_.empty() ) {
			readsAndExons_.back().first.addSecondaryAlignment( bamRecordPtr.get(), bamHeader.get() );
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
						threadOutStrings.at(iThread) = stringifyAlignementRange(eachRange.first, eachRange.second);
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

/*
void BAMtoGenome::saveUnmappedRegions(const std::string &outFileName, const size_t &nThreads) const {
};
*/

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
	} else {
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
	}
	readsAndExons_.emplace_back(alignedRead, *searchIt);
	return searchIt;
};
/*
void BAMtoGenome::findOverlappingGene_(const std::string &referenceName, std::vector<ExonGroup>::const_iterator &gffExonGroupStart, ReadExonCoverage &readCoverageInfo) {
	// if the BAM file is not sorted, we may have to backtrack
	// looking for the first gene end that is after the read map start
	if (readCoverageInfo.alignmentStart < gffExonGroupStart->geneSpan().first) {
		auto reverseLEGI = std::make_reverse_iterator(gffExonGroupStart);
		reverseLEGI      = std::lower_bound(
			reverseLEGI,
			gffExonGroups_[referenceName].crend(),
			readCoverageInfo.alignmentStart,
			[](const ExonGroup &currGroup, const hts_pos_t bamStart) {
				return  (currGroup.geneSpan().first > bamStart);
			}
		);
		if ( reverseLEGI == gffExonGroups_[referenceName].crend() ) {
			readCoverageInfo.geneName               = "no_overlap";
			readCoverageInfo.nExons                 =  0;
			readCoverageInfo.firstExonStart         = -1;
			readCoverageInfo.lastExonEnd            = -1;
			readCoverageInfo.firstExonLength        =  0;
			readCoverageInfo.exonCoverageScores     = std::vector<float>(1, -1.0);
			readCoverageInfo.bestExonCoverageScores = std::vector<float>(1, -1.0);
			// if we are back past the first mRNA, give up on this read but do not update the latest iterator
			return;
		}
		if (reverseLEGI->geneSpan().second < readCoverageInfo.alignmentStart) {
			readCoverageInfo.geneName               = "no_overlap";
			readCoverageInfo.nExons                 =  0;
			readCoverageInfo.firstExonStart         = -1;
			readCoverageInfo.lastExonEnd            = -1;
			readCoverageInfo.firstExonLength        =  0;
			readCoverageInfo.exonCoverageScores     = std::vector<float>(1, -1.0);
			readCoverageInfo.bestExonCoverageScores = std::vector<float>(1, -1.0);
			// the read does not overlap the gene, give up on this read but do not update the latest iterator
			return;
		}
		// convert back to the forward iterator using the reverse/forward relationship
		// safe to increment the reverse iterator here due to the test above
		gffExonGroupStart = std::next(reverseLEGI).base();
	} else {
		gffExonGroupStart = std::lower_bound(
			gffExonGroupStart,
			gffExonGroups_[referenceName].cend(),
			readCoverageInfo.alignmentStart,
			[](const ExonGroup &currGroup, const hts_pos_t bamStart) {
				return  (currGroup.geneSpan().second < bamStart);
			}
		);
	}
	if ( gffExonGroupStart == gffExonGroups_[referenceName].cend() ) {
		readCoverageInfo.geneName               = "past_last_mRNA";
		readCoverageInfo.nExons                 =  0;
		readCoverageInfo.firstExonLength        =  0;
		readCoverageInfo.firstExonStart         = -1;
		readCoverageInfo.lastExonEnd            = -1;
		readCoverageInfo.exonCoverageScores     = std::vector<float>(1, -1.0);
		readCoverageInfo.bestExonCoverageScores = std::vector<float>(1, -1.0);
		// return the iterator to a valid record
		std::advance(gffExonGroupStart, -1);
		return;
	}
	if (readCoverageInfo.alignmentStart < gffExonGroupStart->geneSpan().first) { // no overlap with a known gene
		readCoverageInfo.geneName               = "no_overlap";
		readCoverageInfo.nExons                 =  0;
		readCoverageInfo.firstExonLength        =  0;
		readCoverageInfo.firstExonStart         = -1;
		readCoverageInfo.lastExonEnd            = -1;
		readCoverageInfo.exonCoverageScores     = std::vector<float>(1, -1.0);
		readCoverageInfo.bestExonCoverageScores = std::vector<float>(1, -1.0);
		return;
	}
	readCoverageInfo.geneName        = gffExonGroupStart->geneName();
	readCoverageInfo.nExons          = gffExonGroupStart->nExons();
	readCoverageInfo.firstExonLength = gffExonGroupStart->firstExonSpan().second - gffExonGroupStart->firstExonSpan().first + 1;
	readCoverageInfo.firstExonStart  = gffExonGroupStart->geneSpan().first;
	readCoverageInfo.lastExonEnd     = gffExonGroupStart->geneSpan().second;
}
*/

/*
void BAMtoGenome::processPrimaryAlignment_(const std::string &referenceName, const BAMrecord &alignmentRecord, std::unordered_map<std::string, std::vector<ExonGroup>::const_iterator> &latestExonGroupIts) {
	const char strandID                   = (alignmentRecord.isRevComp() ? '-' : '+' );
	const std::string strandReferenceName = referenceName + strandID;
	const auto refNameIt = gffExonGroups_.find(strandReferenceName);
	if ( refNameIt == gffExonGroups_.cend() ) {  // ignore of the chromosome not in the GFF
		return;
	}
	const std::vector<uint32_t> cigarVec{alignmentRecord.getCIGARvector()};
	const uint32_t firstCIGAR = ( strandID == '+' ? cigarVec.front() : cigarVec.back() );

	ReadExonCoverage currentAlignmentInfo;
	currentAlignmentInfo.readName                 = alignmentRecord.getReadName();
	currentAlignmentInfo.chromosomeName           = referenceName;
	currentAlignmentInfo.strand                   = strandID;
	currentAlignmentInfo.alignmentStart           = alignmentRecord.getMapStart();
	currentAlignmentInfo.alignmentEnd             = alignmentRecord.getMapEnd();
	currentAlignmentInfo.bestAlignmentStart       = alignmentRecord.getMapStart();
	currentAlignmentInfo.bestAlignmentEnd         = alignmentRecord.getMapEnd();
	currentAlignmentInfo.firstSoftClipLength      = bam_cigar_oplen(firstCIGAR) * static_cast<uint32_t>(bam_cigar_opchr(firstCIGAR) == 'S');
	currentAlignmentInfo.nSecondaryAlignments     = 0;
	currentAlignmentInfo.nLocalReversedAlignments = 0;
	currentAlignmentInfo.nGoodSecondaryAlignments = 0;

	if ( latestExonGroupIts.find(strandReferenceName) == latestExonGroupIts.cend() ) {
		latestExonGroupIts[strandReferenceName] = gffExonGroups_[strandReferenceName].cbegin();
	}
	findOverlappingGene_(strandReferenceName, latestExonGroupIts[strandReferenceName], currentAlignmentInfo);
	if (currentAlignmentInfo.firstExonStart > 0) { // if an overlapping gene was found
		currentAlignmentInfo.exonCoverageScores     = latestExonGroupIts[strandReferenceName]->getExonCoverageQuality(alignmentRecord);
		currentAlignmentInfo.bestExonCoverageScores = latestExonGroupIts[strandReferenceName]->getExonCoverageQuality(alignmentRecord);
	}
	readCoverageStats_.emplace_back( std::move(currentAlignmentInfo) );
}
*/

/*
void BAMtoGenome::processSecondaryAlignment_(const std::string &referenceName, const BAMrecord &alignmentRecord, const std::unordered_map<std::string, std::vector<ExonGroup>::const_iterator> &latestExonGroupIts) {
	readCoverageStats_.back().nSecondaryAlignments++;
	const bool overlapsCurrentGene =
		(referenceName == readCoverageStats_.back().chromosomeName) &&
		rangesOverlap(readCoverageStats_.back(), alignmentRecord);
	const bool sameStrand = ( alignmentRecord.isRevComp() == (readCoverageStats_.back().strand == '-') );
	readCoverageStats_.back().nLocalReversedAlignments += static_cast<uint16_t>(overlapsCurrentGene) * static_cast<uint16_t>(!sameStrand);
	if (overlapsCurrentGene && sameStrand) {
		const std::string strandedRefName = referenceName + readCoverageStats_.back().strand;
		readCoverageStats_.back().nGoodSecondaryAlignments++;
		const std::vector<uint32_t> cigarVec{alignmentRecord.getCIGARvector()};
		const std::vector<float> exonCovQuality{latestExonGroupIts.at(strandedRefName)->getExonCoverageQuality(alignmentRecord)};
		auto ecqIt = exonCovQuality.cbegin();
		// save the element-wise max quality score
		std::transform(
			readCoverageStats_.back().bestExonCoverageScores.cbegin(),
			readCoverageStats_.back().bestExonCoverageScores.cend(),
			readCoverageStats_.back().bestExonCoverageScores.begin(),
			[&ecqIt](float score){
				const float maxValue = std::max(score, *ecqIt);
				ecqIt++;
				return maxValue;
			}
		);
		readCoverageStats_.back().bestAlignmentStart = std::min(
			readCoverageStats_.back().bestAlignmentStart,
			alignmentRecord.getMapStart()
		);
		readCoverageStats_.back().bestAlignmentEnd = std::max(
			readCoverageStats_.back().bestAlignmentEnd,
			alignmentRecord.getMapEnd()
		);
		return;
	}
}
*/
