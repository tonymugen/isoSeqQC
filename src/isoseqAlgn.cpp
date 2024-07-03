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

/// Read isoSeq alignments and save potential fusions
/** \file
 * \author Anthony J. Greenberg and Rebekah Rogers
 * \copyright Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.1
 *
 * Implementation of classes that take `.bam` files with isoSeq alignments and identify potential fused transcripts.
 *
 */

#include <algorithm>
#include <memory>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <set>
#include <array>
#include <string>
#include <sstream>
#include <fstream>

#include <iostream>

#include "hts.h"
#include "sam.h"
#include "bgzf.h"

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
	return static_cast<uint32_t>( std::distance(exonRanges_.cbegin(), feaIt) );
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
	return static_cast<uint32_t>( std::distance( lebIt, exonRanges_.crend() ) );
}

//BAMrecord methods

std::string BAMrecord::getCIGARstring() const {
	std::vector<uint32_t> cigarVec(
		bam_get_cigar( alignmentRecord_.get() ),                                    // NOLINT
		bam_get_cigar( alignmentRecord_.get() ) + alignmentRecord_->core.n_cigar    // NOLINT
	);

	auto stringify = [](std::string currString, uint32_t cigarElement){
		return std::move(currString)
			+ std::to_string( bam_cigar_oplen(cigarElement) )
			+ bam_cigar_opchr(cigarElement);
	};
	if ( bam_is_rev( alignmentRecord_.get() ) ) {
		std::string cigar = std::accumulate(
			cigarVec.crbegin(),
			cigarVec.crend(),
			std::string(),
			stringify
		);
		return cigar;
	}
	std::string cigar = std::accumulate(
		cigarVec.cbegin(),
		cigarVec.cend(),
		std::string(),
		stringify
	);

	return cigar;
}

// BAMtoGenome methods
constexpr char   BAMtoGenome::gffDelimiter_{'\t'};
constexpr char   BAMtoGenome::attrDelimiter_{';'};
constexpr size_t BAMtoGenome::strandIDidx_{6UL};
constexpr size_t BAMtoGenome::spanStart_{3UL};
constexpr size_t BAMtoGenome::spanEnd_{4UL};

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
		CbamRecordDeleter localDeleter;
		std::unique_ptr<bam1_t, CbamRecordDeleter> bamRecordPtr(bam_init1(), localDeleter);
		const auto nBytes = bam_read1( bamFile.get(), bamRecordPtr.get() );
		if (nBytes == -1) {
			break;
		}
		if (nBytes < -1) {
			continue;
		}
		// I want reads only with reverse-complement flag possibly set
		const bool isUselessMapping = (bamRecordPtr->core.flag & ~BAM_FREVERSE) != 0;
		if (isUselessMapping) {
			continue;
		}
		std::string referenceName( sam_hdr_tid2name(bamHeader.get(), bamRecordPtr->core.tid) );
		const auto refNameIt = gffExonGroups_.find(referenceName);
		if ( refNameIt == gffExonGroups_.cend() ) {  // ignore of the chromosome not in the GFF
			continue;
		}
		BAMrecord currentBAM( std::move(bamRecordPtr) );

		ReadExonCoverage currentAlignmentInfo;
		currentAlignmentInfo.readName       = currentBAM.getReadName();
		currentAlignmentInfo.chromosomeName = referenceName;
		currentAlignmentInfo.strand         = (currentBAM.isRevComp() ? '-' : '+');
		currentAlignmentInfo.cigarString    = currentBAM.getCIGARstring();
		currentAlignmentInfo.alignmentStart = currentBAM.getMapStart();
		currentAlignmentInfo.alignmentEnd   = currentBAM.getMapEnd();

		if ( latestExonGroupIts.find(referenceName) == latestExonGroupIts.cend() ) {
			latestExonGroupIts[referenceName] = gffExonGroups_[referenceName].cbegin();
		}
		// if the BAM file is not sorted, we may have to backtrack
		// looking for the first gene end that is after the read map start
		if (currentAlignmentInfo.alignmentStart < latestExonGroupIts[referenceName]->firstExonSpan().first) {
			// TODO: fill in the actual treatment of this case
			continue;
		}
		latestExonGroupIts[referenceName] = std::lower_bound(
			latestExonGroupIts[referenceName],
			gffExonGroups_[referenceName].cend(),
			currentAlignmentInfo.alignmentStart,
			[&currentBAM](const ExonGroup &currGroup, const hts_pos_t bamStart) {
				// must be the same strand as well as closest
				return ( currentBAM.isRevComp() != currGroup.isNegativeStrand() ) && (currGroup.geneSpan().second < bamStart);
			}
		);
		if ( latestExonGroupIts[referenceName] == gffExonGroups_[referenceName].cend() ) {
			// not break because other chromosomes may still be in play
			continue;
		}
		if (currentAlignmentInfo.alignmentStart < latestExonGroupIts[referenceName]->geneSpan().first) { // no overlap with a known gene
			currentAlignmentInfo.geneName            = "NA";
			currentAlignmentInfo.nExons              =  0;
			currentAlignmentInfo.firstExonStart      = -1;
			currentAlignmentInfo.lastExonEnd         = -1;
			currentAlignmentInfo.firstCoveredExonIdx =  0;
			currentAlignmentInfo.lastCoveredExonIdx  =  0;
			std::cout << stringify(currentAlignmentInfo) << "\n";
			readCoverageStats_[referenceName].emplace_back( std::move(currentAlignmentInfo) );
			continue;
		}
		currentAlignmentInfo.geneName            = latestExonGroupIts[referenceName]->geneName();
		currentAlignmentInfo.nExons              = latestExonGroupIts[referenceName]->nExons();
		currentAlignmentInfo.firstExonStart      = latestExonGroupIts[referenceName]->geneSpan().first;
		currentAlignmentInfo.lastExonEnd         = latestExonGroupIts[referenceName]->geneSpan().second;
		currentAlignmentInfo.firstCoveredExonIdx = latestExonGroupIts[referenceName]->firstExonAfter(currentAlignmentInfo.alignmentStart);
		currentAlignmentInfo.lastCoveredExonIdx  = std::min(
			static_cast<uint32_t>(currentAlignmentInfo.nExons - 1),
			latestExonGroupIts[referenceName]->lastExonBefore(currentAlignmentInfo.alignmentEnd)
		);
		std::cout << stringify(currentAlignmentInfo) << "\n";
		readCoverageStats_[referenceName].emplace_back( std::move(currentAlignmentInfo) );
	}
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
		gffExonGroups_[activeGFFfields.front()].emplace_back(activeGFFfields.back(), activeGFFfields.at(strandIDidx_).front(), exonSpans);
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
			gffExonGroups_[previousGFFfields.front()].emplace_back(previousGFFfields.back(), previousGFFfields.at(strandIDidx_).front(), exonSpanSet);
			exonSpanSet.clear();
		}
		previousGFFfields.back().clear();
		previousGFFfields.front().clear();
		return;
	}
	// neither is empty and are not the same
	if ( !exonSpanSet.empty() ) {
		gffExonGroups_[previousGFFfields.front()].emplace_back(previousGFFfields.back(), previousGFFfields.at(strandIDidx_).front(), exonSpanSet);
		exonSpanSet.clear();
	}
	std::copy( currentGFFline.cbegin(), currentGFFline.cend(), previousGFFfields.begin() );
}

