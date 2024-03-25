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
 * Interface definitions of classes that take `.bam` files with isoSeq alignments and identify potential fused transcripts.
 *
 */

#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <string>

#include "sam.h"

namespace isaSpace {
	class SAMrecord;

	/** \brief Summary of a SAM record
	 *
	 * Stores relevant information from a SAM format alignment record.
	 * Includes primary and secondary alignments of a read.
	 *
	 */
	class SAMrecord {
	public:
		/** \brief Default constructor */
		SAMrecord() = default;
		/** \brief Constructor with data 
		 *
		 * Constructs an object from alignment records.
		 * The records group the primary and all the secondary alignments of a given read.
		 *
		 * \param[in] alignmentGroup alignment records of a read
		 */
		SAMrecord(const std::vector< std::unique_ptr<bam1_t> > &alignmentGroup);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		SAMrecord(const SAMrecord &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `SAMrecord` object
		 */
		SAMrecord& operator=(const SAMrecord &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		SAMrecord(SAMrecord &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `SAMrecord` object
		 */
		SAMrecord& operator=(SAMrecord &&toMove) noexcept = default;
		/** \brief Destructor */
		~SAMrecord() = default;
	private:
		/** \brief Read (i.e. query) name */
		std::string readName_;
		/** \brief Mapping quality (one for all alignments) */
		uint8_t mappingQuality_ = 0;
		/** \brief Start of the alignment position on the reference (one per alignment) */
		std::vector<uint32_t> positionOnReference_;
		/** \brief `minimap2` alignment scores (one per alignment) */
		std::vector<uint16_t> alignmentScore_;
		/** \brief CIGAR strings (one per alignment) */
		std::vector< std::vector<uint32_t> > cigar_;
	};
}
