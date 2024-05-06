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
 * Interface definitions of classes that take `.bam` files with isoSeq alignments and identify potential fused transcripts.
 *
 */

#pragma once

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>
#include <set>
#include <string>

#include "sam.h"

namespace isaSpace {
	struct CbamRecordDeleter;
	struct BamAndGffFiles;
	struct TokenAttibuteListPair;
	class ExonGroup;
	class SAMrecord;
	class FirstExonRemap;

	/** \brief Deleter of the C BAM record */
	struct CbamRecordDeleter {
		void operator()(bam1_t *bamRecordPtr) const {
			bam_destroy1(bamRecordPtr);
		}
	};
	/** \brief BAM and GFF file name pair */
	struct BamAndGffFiles {
		std::string bamFileName;
		std::string gffFileName;
	};
	/** \brief Token name and GFF attribute list pair */
	struct TokenAttibuteListPair {
		std::string tokenName;
		std::vector<std::string> attributeList;
	};

	/** \brief Group of exons from the same gene
	 *
	 * Gathers exons belonging to all transcripts of a gene.
	 */
	class ExonGroup {
	public:
		/** \brief Default constructor */
		ExonGroup() = default;
		/** \brief Constructor with lines from a GFF file 
		 *
		 * The exon set is ordered by the exon starts.
		 * The first exon can be the first of the last, depending on the strand.
		 * The strand is assumed positive, unless explicitly specified as negative by passing the `-` character.
		 *
		 * \param[in] geneName gene name
		 * \param[in] strand mRNA strand
		 * \param[in] exonSet set of exons from the same gene
		 */
		ExonGroup(std::string geneName, const char strand, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSet);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		ExonGroup(const ExonGroup &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `ExonGroup` object
		 */
		ExonGroup& operator=(const ExonGroup &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		ExonGroup(ExonGroup &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `ExonGroup` object
		 */
		ExonGroup& operator=(ExonGroup &&toMove) noexcept = default;
		/** \brief Destructor */
		~ExonGroup() = default;
	private:
		/** \brief Gene name */
		std::string geneName_;
		/** \brief Start and end positions of each exon in order
		 *
		 * `hts_pos_t` is `int64_t`
		 */
		std::vector< std::pair<hts_pos_t, hts_pos_t> > exonRanges_;
		/** \brief iterator pointing to the first exon */
		std::vector< std::pair<hts_pos_t, hts_pos_t> >::iterator firstExonIt_;
	};

	/** \brief Summary of a SAM record set
	 *
	 * Stores relevant information from SAM format alignment records.
	 * Includes primary and secondary alignments of a read.
	 *
	 */
	class SAMrecord {
	public:
		/** \brief Default constructor */
		SAMrecord() = default;
		/** \brief Constructor with data 
		 *
		 * Constructs an object from an alignment record.
		 * This should be the first record of a read encountered in a `.bam` file,
		 * typically the primary alignment.
		 *
		 * \param[in] alignmentRecord alignment record of a read
		 */
		SAMrecord(const std::unique_ptr<bam1_t, CbamRecordDeleter> &alignmentRecord);
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

		/** \brief Output read name 
		 *
		 * \return read name
		 */
		[[gnu::warn_unused_result]] std::string getReadName() const {return readName_; };
		/** \brief Append a secondary read
		 *
		 * Append a secondary alignment record to the current.
		 *
		 * \param[in] alignmentRecord alignment record of a read
		 */
		void appendSecondary(const std::unique_ptr<bam1_t, CbamRecordDeleter> &alignmentRecord);
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
	/** \brief Candidate first exon re-map alignments
	 *
	 * Collects alignments that are candidates for first exon re-mapping.
	 *
	 */
	class FirstExonRemap {
	public:
		/** \brief Default constructor */
		FirstExonRemap() = default;
		/** \brief Constructor intersecting iso-Seq alignments and exons
		 *
		 * Uses exon positions from the provided GFF file to find iso-Seq alignments from the
		 * provided BAM file that may have mis-mapped first exons.
		 *
		 * \param[in] BamAndGffFiles BAM and GFF file name pair
		 *
		 */
		FirstExonRemap(const BamAndGffFiles &bamGFFfilePairNames);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		FirstExonRemap(const FirstExonRemap &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `FirstExonRemap` object
		 */
		FirstExonRemap& operator=(const FirstExonRemap &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		FirstExonRemap(FirstExonRemap &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `FirstExonRemap` object
		 */
		FirstExonRemap& operator=(FirstExonRemap &&toMove) noexcept = default;
		/** \brief Destructor */
		~FirstExonRemap() = default;

		/** \brief Number of exon sets (genes with exons)
		 *
		 * \return number of exon sets
		 */
		[[gnu::warn_unused_result]] size_t nExonSets() const noexcept { return gffExonGroups_.size(); };
	private:
		/** \brief Number of fields in a GFF file */
		static const size_t nGFFfields_;
		/** \brief GFF file column delimiter */
		static const char gffDelimiter_;
		/** \brief GFF file attribute list delimiter */
		static const char attrDelimiter_;
		/** \brief Strand ID position */
		static const size_t strandIDidx_;
		/** \brief Span start position */
		static const size_t spanStart_;
		/** \brief Span end position */
		static const size_t spanEnd_;

		/** \brief GFF file parent record identifier token */
		std::string parentToken_{"Parent="};
		/** \brief GFF file ID token */
		std::string idToken_{"ID="};

		/** \brief Records of failed GFF parsing events 
		 *
		 * Line number/failure description pairs.
		 */
		std::vector< std::pair<uint64_t, std::string> > failedGFFparsingRecords_;

		/** \brief Vector of exon groups (one group per gene) */
		std::vector<ExonGroup> gffExonGroups_;
		/** \brief Vector of abridged SAM/BAM records */
		std::vector<SAMrecord> candidateAlignments_;
	};
}
