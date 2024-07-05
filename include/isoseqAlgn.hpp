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

#include <memory>
#include <utility>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <string>

#include "sam.h"

namespace isaSpace {
	constexpr size_t nGFFfields{9UL};

	struct CbamRecordDeleter;
	struct BamAndGffFiles;
	struct TokenAttibuteListPair;
	struct ReadExonCoverage;
	class  ExonGroup;
	class  BAMrecord;
	class  BAMtoGenome;

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
	/** \brief Exons covered by a read
	 *
	 * For a given alignment, stores information on exons covered by the read.
	 * All positions are 1-based, indexes are 0-based.
	 */
	struct ReadExonCoverage {
		std::string chromosomeName;
		std::string readName;
		std::string cigarString;
		// smaller value first for the negative and positive strand
		hts_pos_t   alignmentStart;
		hts_pos_t   alignmentEnd;
		std::string geneName;       // NA if no known gene
		char        strand;         // '+' or '-'
		uint16_t    nExons;
		hts_pos_t   firstExonStart;
		hts_pos_t   lastExonEnd;
		// completely covered exons
		uint16_t    firstCoveredExonIdx;
		uint16_t    lastCoveredExonIdx;
	};

	/** \brief Group of exons from the same gene
	 *
	 * Gathers exons belonging to all transcripts of a gene.
	 */
	class ExonGroup {
		//friend class BAMtoGenome;
	public:
		/** \brief Default constructor */
		ExonGroup() = default;
		/** \brief Constructor with lines from a GFF file 
		 *
		 * The exon set is ordered by the exon starts.
		 * The first exon can be the first or the last, depending on the strand.
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

		/** \brief Report the gene name 
		 *
		 * \return gene name
		 */
		[[gnu::warn_unused_result]] std::string geneName() const { return geneName_; };
		/** \brief Number of exons in the gene 
		 *
		 * \return number of exons
		 */
		[[gnu::warn_unused_result]] size_t nExons() const noexcept { return exonRanges_.size(); };
		/** \brief Strand ID
		 *
		 * \return strand ID (`+` or `-`)
		 */
		[[gnu::warn_unused_result]] char strand() const noexcept { return isNegativeStrand_ ? '-' : '+'; };
		/** \brief Gene span 
		 *
		 * Returns the position span of the gene.
		 * First element of the pair is always smaller than the last
		 * regardless of the strand.
		 *
		 * \return first and last nucleotide position (1-based) of the gene
		 */
		[[gnu::warn_unused_result]] std::pair<hts_pos_t, hts_pos_t> geneSpan() const {
			return std::pair<hts_pos_t, hts_pos_t>{exonRanges_.front().first, exonRanges_.back().second};
		};
		/** \brief First exon span 
		 *
		 * Returns the position of the first exon, depends on the strand.
		 * First position is always smaller than the second regardless of strand,
		 * in keeping with the GFF3 specification.
		 *
		 * \return first exon nucleotide position pair (1-based)
		 */
		[[gnu::warn_unused_result]] std::pair<hts_pos_t, hts_pos_t> firstExonSpan() const { return *firstExonIt_; };
		/** \brief Index of the first exon after a given position
		 * 
		 * 0-based index of the first exon found entirely after the given position.
		 * Equal to the total number of exons if the position is after the last exon.
		 *
		 * \param[in] position genome position to test
		 * \return index of the first exon after the given position
		 */
		[[gnu::warn_unused_result]] uint32_t firstExonAfter(const hts_pos_t &position) const noexcept;
		/** \brief Index of the last exon before a given position
		 * 
		 * 0-based index of the last exon found entirely before the given position.
		 * Equal to the total number of exons if the position is after the last exon.
		 *
		 * \param[in] position genome position to test
		 * \return index of the last exon before the given position
		 */
		[[gnu::warn_unused_result]] uint32_t lastExonBefore(const hts_pos_t &position) const noexcept;
	private:
		/** \brief Gene name */
		std::string geneName_;
		/** \brief Is the mRNA on the negative strand? */
		bool isNegativeStrand_{false};
		/** \brief Start and end positions of each exon in order
		 *
		 * `hts_pos_t` is `int64_t`
		 */
		std::vector< std::pair<hts_pos_t, hts_pos_t> > exonRanges_;
		/** \brief Iterator pointing to the first exon */
		std::vector< std::pair<hts_pos_t, hts_pos_t> >::iterator firstExonIt_;
	};

	/** \brief Summary of a BAM record set
	 *
	 * Stores relevant information from BAM format alignment records.
	 */
	class BAMrecord {
	public:
		/** \brief Default constructor */
		BAMrecord() = default;
		/** \brief Constructor with data 
		 *
		 * Constructs an object from an HTSLIB alignment record.
		 *
		 * \param[in] alignmentRecordPointer pointer to a read alignment record
		 */
		BAMrecord(std::unique_ptr<bam1_t, CbamRecordDeleter> &&alignmentRecordPointer) : alignmentRecord_{std::move(alignmentRecordPointer)} {};
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		BAMrecord(const BAMrecord &toCopy) = delete;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `BAMrecord` object
		 */
		BAMrecord& operator=(const BAMrecord &toCopy) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		BAMrecord(BAMrecord &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `BAMrecord` object
		 */
		BAMrecord& operator=(BAMrecord &&toMove) noexcept = default;
		/** \brief Destructor */
		~BAMrecord() = default;

		/** \brief Output read name 
		 *
		 * \return read name
		 */
		[[gnu::warn_unused_result]] std::string getReadName() const {return std::string{bam_get_qname( alignmentRecord_.get() )}; }; // NOLINT
		/** \brief Map start position
		 *
		 * Position of the first mapped nucleotide.
		 *
		 * \return 1-based read map start position
		 */
		[[gnu::warn_unused_result]] hts_pos_t getMapStart() const noexcept { return alignmentRecord_->core.pos + 1; };
		/** \brief Map end position
		 *
		 * Position of the first past the mapped region of the reference.
		 *
		 * \return 1-based read map end position
		 */
		[[gnu::warn_unused_result]] hts_pos_t getMapEnd() const noexcept { return bam_endpos( alignmentRecord_.get() ) + 1; };
		/** \brief mRNA start position
		 *
		 * Position of the first mRNA read nucleotide, taking into account possible reverse-complement.
		 *
		 * \return 1-based read nRNA start position
		 */
		[[gnu::warn_unused_result]] hts_pos_t getmRNAstart() const noexcept {
			return bam_is_rev( alignmentRecord_.get() ) ? bam_endpos( alignmentRecord_.get() ) + 1 : alignmentRecord_->core.pos + 1 ;
		};
		/** \brief Is the read reverse-complemented?
		 *
		 * \return true if the read is reverse-complemented
		 */
		[[gnu::warn_unused_result]] bool isRevComp() const noexcept { return bam_is_rev( alignmentRecord_.get() ); };
		/** \brief CIGAR string 
		 *
		 * Reversed if the read is reverse-complemented.
		 *
		 * \return CIGAR string
		 */
		[[gnu::warn_unused_result]] std::string getCIGARstring() const;
	private:
		/** \brief Pointer to the BAM record */
		std::unique_ptr<bam1_t, CbamRecordDeleter> alignmentRecord_;
	};

	/** \brief Relate BAM alignments to exons
	 *
	 * Relates each BAM alignment to a gene and check which exons are covered by the read.
	 *
	 */
	class BAMtoGenome {
	public:
		/** \brief Default constructor */
		BAMtoGenome() = default;
		/** \brief Constructor intersecting iso-Seq alignments and exons
		 *
		 * Uses exon positions from the provided GFF file to find iso-Seq alignments from the
		 * provided BAM file that may have mis-mapped first exons.
		 *
		 * \param[in] BamAndGffFiles BAM and GFF file name pair
		 *
		 */
		BAMtoGenome(const BamAndGffFiles &bamGFFfilePairNames);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		BAMtoGenome(const BAMtoGenome &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `BAMtoGenome` object
		 */
		BAMtoGenome& operator=(const BAMtoGenome &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		BAMtoGenome(BAMtoGenome &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `BAMtoGenome` object
		 */
		BAMtoGenome& operator=(BAMtoGenome &&toMove) noexcept = default;
		/** \brief Destructor */
		~BAMtoGenome() = default;

		/** \brief Number of chromosomes/scaffolds/linkage groups */
		[[gnu::warn_unused_result]] size_t nChromosomes() const noexcept { return gffExonGroups_.size(); };
		/** \brief Number of exon sets (genes with exons)
		 *
		 * \return number of exon sets
		 */
		[[gnu::warn_unused_result]] size_t nExonSets() const noexcept;
	private:
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

		/** \brief Vector of exon groups (one group per gene)
		 *
		 * The map keys are linkage groups, scaffolds, or chromosomes.
		 */
		std::unordered_map< std::string, std::vector<ExonGroup> > gffExonGroups_;
		/** \brief Vector of abridged SAM/BAM records
		 *
		 * The map keys are linkage groups, scaffolds, or chromosomes.
		 */
		std::unordered_map< std::string, std::vector<ReadExonCoverage> > readCoverageStats_;

		/** \brief Parse a GFF file
		 *
		 * Extract exons from a GFF file and group them by gene.
		 *
		 * \param[in] gffFileName GFF file name
		 */
		void parseGFF_(const std::string &gffFileName);
		/** \brief Extract mRNA information from GFF 
		 *
		 * Extract information from a GFF line that has mRNA data.
		 * Changes the latest gene name if the parent of the current mRNA is different and updates the exon groups vector if necessary.
		 *
		 * \param[in,out] currentGFFline fields from the GFF file line just read
		 * \param[in,out] previousGFFfields GFF fields from previous lines; the attribute field only has the current gene ID
		 * \param[in,out] exonSpanSet set of unique exon start/end pairs
		 */
		void mRNAfromGFF_(std::array<std::string, nGFFfields> &currentGFFline, std::array<std::string, nGFFfields> &previousGFFfields, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSpanSet);
	};
}
