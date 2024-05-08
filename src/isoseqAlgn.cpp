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

#include <cstdint>
#include <algorithm>
#include <memory>
#include <iterator>
#include <numeric>
#include <vector>
#include <set>
#include <array>
#include <string>
#include <sstream>
#include <fstream>

#include "hts.h"
#include "sam.h"

#include "isoseqAlgn.hpp"
#include "helperFunctions.hpp"

using namespace isaSpace;

ExonGroup::ExonGroup(std::string geneName, const char strand, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSet) : geneName_{std::move(geneName)} {
	if ( exonSet.empty() ) {
		throw std::string("ERROR: set of exons is empty in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	std::copy(
		exonSet.cbegin(),
		exonSet.cend(),
		std::back_inserter(exonRanges_)
	);
	firstExonIt_ = ( strand == '-' ? std::prev( exonRanges_.end() ) : exonRanges_.begin() );
}

//SAMrecord methods
SAMrecord::SAMrecord(const std::unique_ptr<bam1_t, CbamRecordDeleter> &alignmentRecord) : 
			readName_{bam_get_qname( alignmentRecord.get() )}, // NOLINT
			mappingQuality_{alignmentRecord->core.qual} {
	this->appendSecondary(alignmentRecord);
}

void SAMrecord::appendSecondary(const std::unique_ptr<bam1_t, CbamRecordDeleter> &alignmentRecord) {
	constexpr std::array<char, 2> alignmentScoreToken{'A', 'S'};
	cigar_.emplace_back(bam_get_cigar(alignmentRecord), bam_get_cigar(alignmentRecord) + alignmentRecord->core.n_cigar); // NOLINT
	positionOnReference_.emplace_back(alignmentRecord->core.pos);
	auto *const alignmentScoreRecord{bam_aux_get( alignmentRecord.get(), alignmentScoreToken.data() )};
	if (alignmentScoreRecord == nullptr) {
		alignmentScore_.emplace_back(0);
		return;
	}
	alignmentScore_.emplace_back( static_cast<uint16_t>( bam_aux2i(alignmentScoreRecord) ) );
}

// FirstExonRemap methods
constexpr size_t FirstExonRemap::nGFFfields_{9};
constexpr char   FirstExonRemap::gffDelimiter_{'\t'};
constexpr char   FirstExonRemap::attrDelimiter_{';'};
constexpr size_t FirstExonRemap::strandIDidx_{6};
constexpr size_t FirstExonRemap::spanStart_{3};
constexpr size_t FirstExonRemap::spanEnd_{4};

FirstExonRemap::FirstExonRemap(const BamAndGffFiles &bamGFFfilePairNames) {
	parseGFF_(bamGFFfilePairNames.gffFileName);
}

size_t FirstExonRemap::nExonSets() const noexcept {
	return std::accumulate(
		gffExonGroups_.cbegin(),
		gffExonGroups_.cend(),
		0UL,
		[](const size_t &sum, const std::pair< std::string, std::vector<ExonGroup> > &eachLnkGrp) {
			return sum + eachLnkGrp.second.size(); 
		}
	);
}

void FirstExonRemap::parseGFF_(const std::string &gffFileName) {
	std::string gffLine;
	std::fstream gffStream(gffFileName, std::ios::in);
	std::string latestGeneName;
	std::set< std::pair<hts_pos_t, hts_pos_t> > exonSpans;
	std::array<std::string, nGFFfields_> gffFields;
	while ( std::getline(gffStream, gffLine) ) {
		if ( gffLine.empty() || (gffLine.at(0) == '#') ) {
			continue;
		}
		std::stringstream currLineStream(gffLine);
		size_t iField{0};
		while ( (iField < nGFFfields_) && std::getline(currLineStream, gffFields.at(iField), gffDelimiter_) ) {
			++iField;
		}
		if (iField < nGFFfields_) {
			continue;
		}
		if (gffFields.at(2) == "mRNA") {
			mRNAfromGFF_(gffFields.front(), gffFields.at(strandIDidx_).front(), gffFields.back(), latestGeneName, exonSpans);
			continue;
		}
		if ( (gffFields.at(2) == "exon") && ( !latestGeneName.empty() ) ) {
			const auto exonStart = static_cast<hts_pos_t>( stol( gffFields.at(spanStart_) ) );
			const auto exonEnd   = static_cast<hts_pos_t>( stol( gffFields.at(spanEnd_) ) );
			exonSpans.insert({exonStart, exonEnd});
		}
	}
	if ( !latestGeneName.empty() && !exonSpans.empty() ) {
		gffExonGroups_[gffFields.front()].emplace_back(latestGeneName, gffFields.at(strandIDidx_).front(), exonSpans);
	}
}

void FirstExonRemap::mRNAfromGFF_(const std::string &lgField, char strandID, const std::string &attributeField, std::string &latestGeneName, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSpanSet) {
	std::stringstream attributeStream(attributeField);
	std::string attrField;
	TokenAttibuteListPair tokenPair;
	while ( std::getline(attributeStream, attrField, attrDelimiter_) ) {
		tokenPair.attributeList.emplace_back(attrField);
	}
	tokenPair.tokenName = parentToken_;

	std::string parentName{extractAttributeName(tokenPair)};
	if ( latestGeneName.empty() || parentName.empty() ) {
		std::swap(latestGeneName, parentName);
		return;
	}
	if ( (latestGeneName != parentName) && ( !exonSpanSet.empty() ) ) {
		gffExonGroups_[lgField].emplace_back(latestGeneName, strandID, exonSpanSet);
		std::swap(latestGeneName, parentName);
		exonSpanSet.clear();
	}
}

