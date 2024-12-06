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
 * Interface definitions of classes that take `.bam` files with isoSeq alignments and identify potential fused transcripts.
 *
 */

#pragma once

#include <utility>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <string>

#include "htslib/sam.h"

namespace isaSpace {
	constexpr size_t nGFFfields{9UL};

	struct BamAndGffFiles;
	struct TokenAttibuteListPair;
	struct MappedReadInterval;
	struct BinomialWindowParameters;
	struct ReadExonCoverage;
	class  ExonGroup;
	class  ReadMatchWindowBIC;
	class  BAMrecord;
	class  BAMtoGenome;

	/** \brief BAM and GFF file name pair */
	struct BamAndGffFiles {
		/** \brief BAM file name */
		std::string bamFileName;
		/** \brief GFF file name */
		std::string gffFileName;
	};
	/** \brief Token name and GFF attribute list pair */
	struct TokenAttibuteListPair {
		/** \brief GFF token name */
		std::string tokenName;
		/** \brief Attribute list */
		std::vector<std::string> attributeList;
	};
	/** \brief Delineates an interval in a mapped read 
	 *
	 * Provides a start and end of a read interval and corresponding
	 * positions on the reference.
	 */
	struct MappedReadInterval {
		/** \brief Base-0 start position on the read */
		hts_pos_t readStart{-1};
		/** \brief Base-0 end position on the read */
		hts_pos_t readEnd{-1};
		/** \brief Base-1 start position on the reference */
		hts_pos_t referenceStart{-1};
		/** \brief Base-1 end position on the reference */
		hts_pos_t referenceEnd{-1};
	};
	/** \brief Binomial window parameters */
	struct BinomialWindowParameters {
		/** \brief Probability of success up to this window */
		float currentProbability{0.0F};
		/** \brief Probability of success to consider as an alternative 
		 *
		 * Test for a switch to this probability when assessing evidence for a change point 
		 * between mapped and unmapped regions of a read.
		 */
		float alternativeProbability{0.0F};
		/** \brief Minimum difference in BIC between windows with current and alternative probabilities */
		float bicDifferenceThreshold{0.0F};
		/** \brief Window size */
		std::vector< std::pair<float, hts_pos_t> >::difference_type windowSize{0};
	};
	/** \brief Exons covered by a read
	 *
	 * For a given alignment, stores information on exons covered by the read.
	 * All positions are 1-based, indexes are 0-based.
	 */
	struct ReadExonCoverage {
		/** \brief Chromosome name 
		 * 
		 * This can also be a linkage group or scaffold name,
		 * as listed in the GFF file in the `seqid` column.
		 */
		std::string chromosomeName;
		/** \brief Read name */
		std::string readName;
		/** \brief Alignment start
		 *
		 * Base-1 position of the read start from the primary alignment.
		 * Never larger than `alignmentEnd` regardless of strand.
		 */
		hts_pos_t alignmentStart{0};
		/** \brief Alignment end
		 *
		 * Base-1 position of the read end from the primary alignment.
		 * Never smaller than `alignmentStart` regardless of strand.
		 */
		hts_pos_t alignmentEnd{0};
		/** \brief Best alignment start
		 *
		 * Base-1 position of the earliest read start from the primary and all good secondary alignments.
		 * Never larger than `bestAlignmentEnd` regardless of strand.
		 */
		hts_pos_t bestAlignmentStart{0};
		/** \brief Best alignment end
		 *
		 * Base-1 position of the latest read end from the primary and all good secondary alignments.
		 * Never smaller than `bestAlignmentStart` regardless of strand.
		 */
		hts_pos_t bestAlignmentEnd{0};
		/** \brief Length of the soft clip at read start
		 *
		 * End of the CIGAR string if the read is reverse-complemented.
		 * Set to 0 if there is no soft clip.
		 */
		uint32_t firstSoftClipLength{0};
		/** \brief Number of secondary alignments */
		uint16_t nSecondaryAlignments;
		/** \brief Number of good secondary alignments 
		 *
		 * Secondary alignments that are on the same strand as
		 * the primary and overlap the same gene.
		 */
		uint16_t nGoodSecondaryAlignments{0};
		/** \brief Number of locally mapping reverse-complemented reads 
		 *
		 * Number of alignments on the opposite strand from the 
		 * primary that overlap the same gene as the primary.
		 */
		uint16_t nLocalReversedAlignments{0};
		/** \brief Number of exons */
		uint16_t nExons;
		/** \brief Gene name
		 *
		 * Set to `no_overlap` or `past_last_mRNA` if there is no
		 * gene in the GFF file that overlaps the alignment.
		 */
		std::string geneName;
		/** \brief Strand ID
		 *
		 * Must be `+` or `-`.
		 */
		char strand{'\0'};
		/** \brief First exon length 
		 *
		 * Actual first exon, last in the sequence if the strand is negative.
		 */
		hts_pos_t firstExonLength{0};
		/** \brief First exon start
		 *
		 * Base-1 position of the first exon start from the GFF file.
		 * End of the last exon if the strand negative.
		 */
		hts_pos_t firstExonStart{0};
		/** \brief Last exon end
		 *
		 * Base-1 position of the last exon end from the GFF file.
		 * Start of the first exon if the strand negative.
		 */
		hts_pos_t lastExonEnd{0};
		/** \brief Exon coverage scores
		 *
		 * Fraction of reference bases in each exon covered by a matching base
		 * in the read from the primary alignment.
		 */
		std::vector<float> exonCoverageScores;
		/** \brief Best exon coverage scores
		 *
		 * Fraction of reference bases in each exon covered by a matching base
		 * in the read from the best alignment for that exon among all alignments for the read.
		 */
		std::vector<float> bestExonCoverageScores;
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
		/** \brief Range covered by a given exon
		 *
		 * \param[in] idx exon index
		 * \return the exon start and end nucleotide position pair
		 */
		[[gnu::warn_unused_result]] std::pair<hts_pos_t, hts_pos_t> operator[](const size_t &idx) const {return exonRanges_.at(idx);};
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
		 * Equal to the index of the last exon if the position is after the last exon.
		 *
		 * \param[in] position genome position to test
		 * \return index of the first exon after the given position
		 */
		[[gnu::warn_unused_result]] uint32_t firstExonAfter(const hts_pos_t &position) const noexcept;
		/** \brief Index of the first exon overlapping a given position
		 * 
		 * 0-based index of the first exon found at least partially after the given position.
		 * Equal to the index of the last exon if the position is after the last exon.
		 *
		 * \param[in] position genome position to test
		 * \return index of the first exon overlapping the given position
		 */
		[[gnu::warn_unused_result]] uint32_t firstOverlappingExon(const hts_pos_t &position) const noexcept;
		/** \brief Index of the last exon before a given position
		 * 
		 * 0-based index of the last exon found entirely before the given position.
		 * Equal to the index of the last exon if the position is after the last exon.
		 *
		 * \param[in] position genome position to test
		 * \return index of the last exon before the given position
		 */
		[[gnu::warn_unused_result]] uint32_t lastExonBefore(const hts_pos_t &position) const noexcept;
		/** \brief Index of the last exon overlapping a given position
		 * 
		 * 0-based index of the last exon found at least before after the given position.
		 * Equal to the index of the last exon if the position is after the last exon.
		 *
		 * \param[in] position genome position to test
		 * \return index of the last exon overlapping the given position
		 */
		[[gnu::warn_unused_result]] uint32_t lastOverlappingExon(const hts_pos_t &position) const noexcept;
		/** \brief First intron span
		 *
		 * Smaller value first regardless of strand, but the intron is always the first in the gene,
		 * last in the sequence if the strand is negative.
		 * If there is only one exon, the start and end are reported as `-1`.
		 *
		 * \return first intron start and end positions
		 */
		[[gnu::warn_unused_result]] std::pair<hts_pos_t, hts_pos_t> getFirstIntronSpan() const;
		/** \brief Get read coverage quality per exon 
		 *
		 * Use CIGAR information to extract alignment quality for each exon.
		 * Quality is the fraction of reference nucleotides covered by matching read bases.
		 *
		 * \param[in] alignment BAM alignment record
		 * \return vector of alignment qualities, one per exon
		 *
		 */
		[[gnu::warn_unused_result]] std::vector<float> getExonCoverageQuality(const BAMrecord &alignment) const;
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
		// TODO: to implement tracking transcripts, additional vector of iterators pointing to the next exon
		// Keep the number the same and equal to the number of transcripts listed in the GFF
	};

	/** \brief Binomial BIC for a read window
	 *
	 * Estimates of the Bayesian information criterion for a window of the read-centric match status vector.
	 * The BIC compares the model with the same success probability as the previous window to a change.
	 * The success probabilities reflect identity between unrelated or matching sequences.
	 */
	class ReadMatchWindowBIC {
	public:
		/** \brief Default constructor */
		ReadMatchWindowBIC() = default;
		/** \brief Constructor with a read match vector window 
		 *
		 * \param[in] windowBegin start of a read match status vector
		 * \param[in] windowParameters binomial probability parameters of the window
		 */
		ReadMatchWindowBIC(const std::vector< std::pair<float, hts_pos_t> >::const_iterator &windowBegin, const BinomialWindowParameters &windowParameters);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		ReadMatchWindowBIC(const ReadMatchWindowBIC &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `ReadMatchWindowBIC` object
		 */
		ReadMatchWindowBIC& operator=(const ReadMatchWindowBIC &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		ReadMatchWindowBIC(ReadMatchWindowBIC &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `ReadMatchWindowBIC` object
		 */
		ReadMatchWindowBIC& operator=(ReadMatchWindowBIC &&toMove) noexcept = default;
		/** \brief Destructor */
		~ReadMatchWindowBIC() = default;

		/** \brief Get BIC difference 
		 *
		 * Reports the BIC difference between the model with the same mapping quality
		 * as the previous window and a model with a switch to a different quality.
		 *
		 * \return evidence for a switch in mapping quality
		 */
	 	[[gnu::warn_unused_result]] float getBICdifference() const noexcept;
	private:
		/** \brief Success probability for the previous window */
		float leftProbability_{0.0F};
		/** \brief Potential success probability for the current window 
		 *
		 * If the `leftProbability_` corresponds to an unmapped region,
		 * takes the value of the mapped region probability and vice versa.
		 */
		float rightProbability_{0.0F};
		/** \brief Read to reference match count */
		float kSuccesses_{0.0F};
		/** \brief Window size */
		float nTrials_{0.0F};
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
		 * Constructs an object from an HTSLIB alignment record and the corresponding header.
		 *
		 * \param[in] alignmentRecord pointer to a read alignment record
		 * \param[in] samHeader pointer to the corresponding BAM/SAM header
		 */
		BAMrecord(const bam1_t *alignmentRecord, const sam_hdr_t *samHeader);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		BAMrecord(const BAMrecord &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `BAMrecord` object
		 */
		BAMrecord& operator=(const BAMrecord &toCopy) = default;
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
		[[gnu::warn_unused_result]] std::string getReadName() const {return readName_; };
		/** \brief Map start position
		 *
		 * Position of the first mapped nucleotide.
		 *
		 * \return 1-based read map start position
		 */
		[[gnu::warn_unused_result]] hts_pos_t getMapStart() const noexcept { return mapStart_; };
		/** \brief Map end position
		 *
		 * Position of the first past the mapped region of the reference.
		 *
		 * \return 1-based read map end position
		 */
		[[gnu::warn_unused_result]] hts_pos_t getMapEnd() const noexcept { return mapEnd_; };
		/** \brief mRNA start position
		 *
		 * Position of the first mRNA read nucleotide, taking into account possible reverse-complement.
		 *
		 * \return 1-based read mRNA start position
		 */
		[[gnu::warn_unused_result]] hts_pos_t getmRNAstart() const noexcept {
			return isRev_ ? mapEnd_ : mapStart_;
		};
		/** \brief Is the read reverse-complemented?
		 *
		 * \return true if the read is reverse-complemented
		 */
		[[gnu::warn_unused_result]] bool isRevComp() const noexcept { return isRev_; };
		/** \brief Is this the primary alignment?
		 *
		 * \return true if the alignment is primary
		 */
		[[gnu::warn_unused_result]] bool isPrimaryMap() const noexcept { return isPrimary_ && isMapped_; };
		/** \brief Is this a valid secondary alignment?
		 *
		 * \return true if the read is mapped and the alignment is secondary
		 */
		[[gnu::warn_unused_result]] bool isSecondaryMap() const noexcept { return (!isPrimary_) && isMapped_; };
		/** \brief Read length 
		 *
		 * \return read length in bases
		 */
		[[gnu::warn_unused_result]] hts_pos_t getReadLength() const noexcept { return static_cast<hts_pos_t>( sequenceAndQuality_.size() ); };
		/** \brief CIGAR vector 
		 *
		 * Orientation independent of strand
		 *
		 * \return CIGAR vector
		 */
		[[gnu::warn_unused_result]] std::vector<uint32_t> getCIGARvector() const { return cigar_; };
		/** \brief CIGAR string 
		 *
		 * Reversed if the read is reverse-complemented.
		 *
		 * \return CIGAR string
		 */
		[[gnu::warn_unused_result]] std::string getCIGARstring() const;
		/** \brief Get reference name 
		 *
		 * Reference sequence (e.g., chromosome) name. If absent, returns `*` like samtools.
		 *
		 * \return reference name
		 */ 
		[[gnu::warn_unused_result]] std::string getReferenceName() const {return referenceName_; };
		/** \brief Read match status along the reference
		 *
		 * Parses CIGAR to track read (query) match/mismatch (1.0 for match, 0.0 for mismatch) status along the reference.
		 * This means that insertions in the read are ignored.
		 * The vector start begins at the position closest to the first exon start regardless of strand.
		 *
		 * \return vector of match status
		 */
		[[gnu::warn_unused_result]] std::vector<float> getReferenceMatchStatus() const;
		/** \brief Reference match status along the read 
		 *
		 * Parses CIGAR to track reference match/mismatch (1.0 for match, 0.0 for mismatch) status along the read,
		 * relating each read position to the corresponding reference base-1 nucleotide position.
		 * The vector start begins at the position closest to the first exon start regardless of strand.
		 *
		 * \return vector of match status/reference position pairs
		 */
		[[gnu::warn_unused_result]] std::vector< std::pair<float, hts_pos_t> > getReadCentricMatchStatus() const;
		/** \brief Identify unmapped portions of the read 
		 *
		 * Returns a vector of poorly mapped portions of the read. Vector is empty if the read is mapped.
		 *
		 * \param[in] windowParameters window parameters: size and the mapped/unmapped region binomial probabilities 
		 *
		 * \return vector of poorly mapped region coordinates
		 */
		[[gnu::warn_unused_result]] std::vector<MappedReadInterval> getPoorlyMappedRegions(const BinomialWindowParameters &windowParameters) const;
	private:
		/** \brief Mask isolating the sequence byte */
		static const uint16_t sequenceMask_;
		/** \brief Size of shift to get the quality byte */
		static const uint16_t qualityShift_;
		/** \brief Mask to test for secondary alignment */
		static const uint16_t suppSecondaryAlgn_;
		/** \brief Query (read) consumption status array 
		 *
		 * Can be indexed into using the CIGAR operation bit field. 
		 */
		static const std::array<hts_pos_t, 10> queryConsumption_;
		/** \brief Reference consumption status array 
		 *
		 * Can be indexed into using the CIGAR operation bit field. 
		 */
		static const std::array<hts_pos_t, 10> referenceConsumption_;
		/** \brief Sequence match status array 
		 *
		 * Can be indexed into using the CIGAR operation bit field. 
		 */
		static const std::array<float, 10> sequenceMatch_;
		/** \brief Is this the primary alignment? */
		bool isPrimary_{true};
		/** \brief Is the read mapped? */
		bool isMapped_{false};
		/** \brief Is the read reverse-complemented? */
		bool isRev_{false};
		/** \brief Map start 
		 *
		 * Base-1 position on the reference where the read alignment starts.
		 * As in the BAM record, always precedes map end regardless of strand.
		 */
		hts_pos_t mapStart_{0};
		/** \brief Map end 
		 *
		 * Base-1 position on the reference where the read alignment ends.
		 * As in the BAM record, always after the map start regardless of strand.
		 */
		hts_pos_t mapEnd_{0};
		/** \brief Read name */
		std::string readName_;
		/** \brief Reference name */
		std::string referenceName_;
		/** \brief CIGAR vector */
		std::vector<uint32_t> cigar_;
		/** \brief Sequence its quality
		 *
		 * Lower byte is the sequence in BAM format,
		 * higher byte the BAM (no +33 adjustment) quality score.
		 */
		std::vector<uint16_t> sequenceAndQuality_;
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

		/** \brief Number of chromosomes/scaffolds/linkage groups
		 *
		 * Counts the number of unique GFF `seqid` elements.
		 * These are typically chromosomes, scaffolds, or linkage groups.
		 *
		 * \return number of chromosomes
		 */
		[[gnu::warn_unused_result]] size_t nChromosomes() const noexcept;
		/** \brief Number of exon sets (genes with exons)
		 *
		 * \return number of exon sets
		 */
		[[gnu::warn_unused_result]] size_t nExonSets() const noexcept;
		/** \brief Save read coverage to file
		 *
		 * Saves the read coverage statistics to a file.
		 * If a file with the same name exists it is overwritten.
		 *
		 * \param[in] outFileName output file name
		 * \param[in] nThreads number of concurrent threads
		 */
		void saveReadCoverageStats(const std::string &outFileName, const size_t &nThreads) const;
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
		/** \brief Flag testing the two possible secondary alignment markers */
		static const uint16_t suppSecondaryAlgn_;

		/** \brief GFF file parent record identifier token */
		std::string parentToken_{"Parent="};

		/** \brief Vector of exon groups (one group per gene)
		 *
		 * The map keys are linkage groups, scaffolds, or chromosomes plus strand ID.
		 */
		std::unordered_map< std::string, std::vector<ExonGroup> > gffExonGroups_;
		/** \brief Vector of abridged SAM/BAM records */
		std::vector<ReadExonCoverage> readCoverageStats_;

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
		/** \brief Find the gene overlapping the read
		 *
		 * Finds the collection of exons belonging to a gene that is covered by the current isoSeq read.
		 * Updates the GFF iterator and the coverage statistics object.
		 *
		 * \param[in] referenceName the name of the reference sequence (chromosome, linkage group, or contig)
		 * \param[in,out] gffExonGroupStart the iterator to start the search
		 * \param[in,out] readCoverageInfo the object with read coverage information
		 */
		void findOverlappingGene_(const std::string &referenceName, std::vector<ExonGroup>::const_iterator &gffExonGroupStart, ReadExonCoverage &readCoverageInfo);
		/** \brief Process a primary alignment 
		 *
		 * \param[in] referenceName reference sequence name
		 * \param[in] alignmentRecord current alignment record
		 * \param[in,out] latestExonGroupIts iterators to the latest exon groups for each reference and strand
		 */
		void processPrimaryAlignment_(const std::string &referenceName, const BAMrecord &alignmentRecord, std::unordered_map<std::string, std::vector<ExonGroup>::const_iterator> &latestExonGroupIts);
		/** \brief Process a secondary alignment 
		 *
		 * \param[in] referenceName reference sequence name
		 * \param[in] alignmentRecord current alignment record
		 * \param[in] latestExonGroupIts iterators to the latest exon groups for each reference and strand
		 */
		void processSecondaryAlignment_(const std::string &referenceName, const BAMrecord &alignmentRecord, const std::unordered_map<std::string, std::vector<ExonGroup>::const_iterator> &latestExonGroupIts);
	};
}
