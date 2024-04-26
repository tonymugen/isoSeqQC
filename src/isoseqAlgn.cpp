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
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.1
 *
 * Implementation of classes that take `.bam` files with isoSeq alignments and identify potential fused transcripts.
 *
 */

#include <cstdint>
#include <iterator>
#include <memory>
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include <sstream>

#include "sam.h"

#include "isoseqAlgn.hpp"

using namespace isaSpace;

// ExonGroup methods
constexpr std::string::difference_type parentTokenSize_{7};

ExonGroup::ExonGroup(const std::string &geneName, std::vector< std::stringstream > &exonGFFlines) {
	if ( exonGFFlines.empty() ) {
		throw std::string("ERROR: vector of exons is empty in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	geneName_ = geneName;

	for (auto &eachLine : exonGFFlines) {
		std::string field;
		eachLine >> field;
		const auto exonStart{static_cast<hts_pos_t>( std::stoi(field) )};
		eachLine >> field;
		const auto exonEnd{static_cast<hts_pos_t>( std::stoi(field) )};

		// make sure the range is increasing as per GFF spec
		exonRanges_.emplace_back( std::min(exonStart, exonEnd), std::max(exonStart, exonEnd) );
	}
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
