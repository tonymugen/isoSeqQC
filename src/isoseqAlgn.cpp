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

#include <algorithm>
#include <cstring>
#include <memory>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <set>
#include <unordered_set>
#include <array>
#include <string>
#include <sstream>
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
	if ( exonSet.empty() ) {
		throw std::string("ERROR: set of exons is empty in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	std::copy(
		exonSet.cbegin(),
		exonSet.cend(),
		std::back_inserter(exonRanges_)
	);
	firstExonIt_ = ( isNegativeStrand_ ? std::prev( exonRanges_.end() ) : exonRanges_.begin() );
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
	std::vector<float> referenceMatchStatus{alignment.getReferenceMatchStatus()};
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

//BAMrecord methods
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
			isPrimary_{(alignmentRecord->core.flag & suppSecondaryAlgn_) == 0}, isMapped_{ (alignmentRecord->core.flag & BAM_FUNMAP) == 0 } {

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

std::vector<float> BAMrecord::getReferenceMatchStatus() const {
	std::vector<float> referenceMatchStatus;
	for (const auto &eachCIGAR : cigar_) {
		const std::vector<float> currCIGARfield(
			bam_cigar_oplen(eachCIGAR) * referenceConsumption_.at( bam_cigar_op(eachCIGAR) ),
			sequenceMatch_.at( bam_cigar_op(eachCIGAR) )
		);
		std::copy( currCIGARfield.cbegin(), currCIGARfield.cend(), std::back_inserter(referenceMatchStatus) );
	}
	return referenceMatchStatus;
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
			referencePosition += rcStatus;
			readMatchStatus.emplace_back(matchStatus, referencePosition);
			iSIGLEN++;
		}
	}

	return readMatchStatus;
}

// BAMtoGenome methods
constexpr char     BAMtoGenome::gffDelimiter_{'\t'};
constexpr char     BAMtoGenome::attrDelimiter_{';'};
constexpr size_t   BAMtoGenome::strandIDidx_{6UL};
constexpr size_t   BAMtoGenome::spanStart_{3UL};
constexpr size_t   BAMtoGenome::spanEnd_{4UL};
constexpr uint16_t BAMtoGenome::suppSecondaryAlgn_{BAM_FSECONDARY | BAM_FSUPPLEMENTARY};

BAMtoGenome::BAMtoGenome(const BamAndGffFiles &bamGFFfilePairNames) {
	parseGFF_(bamGFFfilePairNames.gffFileName);

	if ( gffExonGroups_.empty() ) {
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
		BAMrecord currentBAM( bamRecordPtr.get(), bamHeader.get() );
		std::string referenceName{currentBAM.getReferenceName()};
		// I want primary reads only with reverse-complement flag possibly set
		if ( currentBAM.isPrimaryMap() ) {
			processPrimaryAlignment_(referenceName, currentBAM, latestExonGroupIts);
			continue;
		}
		// This is not a primary alignment. See if it is a good secondary and process if yes.
		if ( ( !readCoverageStats_.empty() ) && currentBAM.isSecondaryMap() ) {
			processSecondaryAlignment_(referenceName, currentBAM, latestExonGroupIts);
		}
	}
}

size_t BAMtoGenome::nChromosomes() const noexcept {
	std::unordered_set<std::string> chromosomes;
	std::for_each(
		gffExonGroups_.cbegin(),
		gffExonGroups_.cend(),
		[&chromosomes](const std::pair< std::string, std::vector<ExonGroup> > &currStrandedChromosome) {
			std::string chrName = currStrandedChromosome.first;
			chrName.pop_back();
			chromosomes.insert(chrName);
		}
	);

	return chromosomes.size();
}

size_t BAMtoGenome::nExonSets() const noexcept {
	return std::accumulate(
		gffExonGroups_.cbegin(),
		gffExonGroups_.cend(),
		0UL,
		[](const size_t &sum, const std::pair< std::string, std::vector<ExonGroup> > &eachLnkGrp) {
			return sum + eachLnkGrp.second.size(); 
		}
	);
}

void BAMtoGenome::saveReadCoverageStats(const std::string &outFileName, const size_t &nThreads) const {
	size_t actualNthreads = std::min( nThreads, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	actualNthreads        = std::min( actualNthreads, readCoverageStats_.size() );
	actualNthreads        = std::max(actualNthreads, 1UL);

	const auto threadRanges{makeThreadRanges(readCoverageStats_, actualNthreads)};
	std::vector<std::string> threadOutStrings(actualNthreads);
	std::vector< std::future<void> > tasks;
	tasks.reserve(actualNthreads);
	size_t iThread{0};
	std::for_each(
		threadRanges.cbegin(),
		threadRanges.cend(),
		[&iThread, &tasks, &threadOutStrings](const std::pair<std::vector<ReadExonCoverage>::const_iterator, std::vector<ReadExonCoverage>::const_iterator> &eachRange) {
			tasks.emplace_back(
				std::async(
					[iThread, eachRange, &threadOutStrings] {
						threadOutStrings.at(iThread) = stringifyRCSrange(eachRange.first, eachRange.second);
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

void BAMtoGenome::parseGFF_(const std::string &gffFileName) {
	std::string gffLine;
	std::fstream gffStream(gffFileName, std::ios::in);
	std::set< std::pair<hts_pos_t, hts_pos_t> > exonSpans;
	std::array<std::string, nGFFfields> activeGFFfields;   // fields from the gene tracked up to now
	std::array<std::string, nGFFfields> newGFFfields;      // fields from the currently read line
	while ( std::getline(gffStream, gffLine) ) {
		if ( gffLine.empty() || (gffLine.at(0) == '#') ) {
			continue;
		}
		std::stringstream currLineStream(gffLine);
		size_t iField{0};
		while ( (iField < nGFFfields) && std::getline(currLineStream, newGFFfields.at(iField), gffDelimiter_) ) {
			++iField;
		}
		// skip incomplete lines
		if (iField < nGFFfields) {
			continue;
		}
		if (newGFFfields.at(2) == "mRNA") {
			mRNAfromGFF_(newGFFfields, activeGFFfields, exonSpans);
			continue;
		}
		if ( (newGFFfields.at(2) == "exon") && ( !activeGFFfields.back().empty() ) ) {
			const auto exonStart = static_cast<hts_pos_t>( stol( newGFFfields.at(spanStart_) ) );
			const auto exonEnd   = static_cast<hts_pos_t>( stol( newGFFfields.at(spanEnd_) ) );
			exonSpans.insert({exonStart, exonEnd});
		}
	}
	if ( !activeGFFfields.back().empty() && !exonSpans.empty() ) {
		const char strandID                  = (activeGFFfields.at(strandIDidx_) == "-" ? '-' : '+');
		const std::string strandedChromosome = activeGFFfields.front() + strandID;
		gffExonGroups_[strandedChromosome].emplace_back(activeGFFfields.back(), activeGFFfields.at(strandIDidx_).front(), exonSpans);
	}
}

void BAMtoGenome::mRNAfromGFF_(std::array<std::string, nGFFfields> &currentGFFline, std::array<std::string, nGFFfields> &previousGFFfields, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSpanSet) {
	std::stringstream attributeStream( currentGFFline.back() );
	std::string attrField;
	TokenAttibuteListPair tokenPair;
	while ( std::getline(attributeStream, attrField, attrDelimiter_) ) {
		tokenPair.attributeList.emplace_back(attrField);
	}
	tokenPair.tokenName = parentToken_;

	std::string parentName{extractAttributeName(tokenPair)};
	std::swap(currentGFFline.back(), parentName);
	if ( currentGFFline.back() == previousGFFfields.back() ) { // both could be empty
		return;
	}
	if ( currentGFFline.back().empty() ) {
		if ( !exonSpanSet.empty() ) {
			const char strandID                  = (previousGFFfields.at(strandIDidx_) == "-" ? '-' : '+');
			const std::string strandedChromosome = previousGFFfields.front() + strandID;
			gffExonGroups_[strandedChromosome].emplace_back(previousGFFfields.back(), previousGFFfields.at(strandIDidx_).front(), exonSpanSet);
			exonSpanSet.clear();
		}
		previousGFFfields.back().clear();
		previousGFFfields.front().clear();
		return;
	}
	// neither is empty and are not the same
	if ( !exonSpanSet.empty() ) {
		const char strandID                  = (previousGFFfields.at(strandIDidx_) == "-" ? '-' : '+');
		const std::string strandedChromosome = previousGFFfields.front() + strandID;
		gffExonGroups_[strandedChromosome].emplace_back(previousGFFfields.back(), previousGFFfields.at(strandIDidx_).front(), exonSpanSet);
		exonSpanSet.clear();
	}
	std::copy( currentGFFline.cbegin(), currentGFFline.cend(), previousGFFfields.begin() );
}

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
