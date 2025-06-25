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

/// Genomic analyses helper functions
/** \file
 * \author Anthony J. Greenberg and Rebekah Rogers
 * \copyright Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.2
 *
 * Definitions of class-external functions needed by genomic analyses.
 *
 */

#pragma once

#include <memory>
#include <array>
#include <string>
#include <utility> // for std::pair
#include <vector>

#include "sam.h"

#include "isoseqAlgn.hpp"

namespace isaSpace {
	constexpr size_t nGFFfields{9UL};
	using bamGFFvector = std::vector< std::pair<BAMrecord, ExonGroup> >;

	/** \brief Extract a name 
	 *
	 * Extract a token name from GFF attributes.
	 *
	 * \param[in] tokenAndAttrList field token and the list of attributes
	 */
	[[gnu::warn_unused_result]] std::string extractAttributeName(const TokenAttibuteListPair &tokenAndAttrList);

	/** \brief Extract parent name
	 *
	 * Extract the value of the `Parent=` attribute from the provided GFF attribute string.
	 *
	 * \param[in] attributeString GFF attribute string
	 * \return parent name
	 */
	[[gnu::warn_unused_result]] std::string extractParentName(const std::string &attributeString);

	/** \brief Test for range overlap
	 *
	 * The first position in each range need not be before the second in the same range.
	 *
	 * \param[in] range1 first range
	 * \param[in] range2 second range
	 * \return `true` if there is overlap
	 */
	[[gnu::warn_unused_result]] bool rangesOverlap(const std::pair<hts_pos_t, hts_pos_t> &range1, const std::pair<hts_pos_t, hts_pos_t> &range2) noexcept;

	/** \brief Parse a GFF line into fields
	 *
	 *
	 * Places `FAIL` in the first element if the number of fields is not `nGFFfields` as required by the GFF specification.
	 *
	 * \param[in] gffLine one line of a GFF file
	 * \return each GFF field in a separate element
	 */
	[[gnu::warn_unused_result]] std::array<std::string, nGFFfields> parseGFFline(const std::string &gffLine);

	/** \brief Parse a GFF file
	 *
	 * Extract exons from a GFF file and group them by gene and chromosome/linkage group.
	 * The map keys are linkage groups, scaffolds, or chromosomes plus strand ID.
	 *
	 * \param[in] gffFileName GFF file name
	 * \return collection of exon group vectors by chromosome and strand
	 */
	[[gnu::warn_unused_result]] std::unordered_map< std::string, std::vector<ExonGroup> > parseGFF(const std::string &gffFileName);

	/** \brief Identify peaks in numerical data 
	 *
	 * Returns iterators to the elements in a vector that correspond to peaks above the provided threshold.
	 * Last element is the `const_iterator` to the end of the vector unless it is empty, and is the only element if there are no peaks.
	 *
	 * \param[in] values vector of values
	 * \param[in] threshold value that must be exceeded for a peak call
	 * \return vector of iterators to peak elements
	 *
	 */
	[[gnu::warn_unused_result]] std::vector<std::vector<float>::const_iterator> getPeaks(const std::vector<float> &values, const float &threshold);

	/** \brief Identify valleys in numerical data 
	 *
	 * Returns iterators to the elements in a vector that correspond to valleys below the provided threshold.
	 * Last element is the `const_iterator` to the end of the vector unless it is empty, and is the only element if there are no valleys.
	 *
	 * \param[in] values vector of values
	 * \param[in] threshold value that must exceed the valley values
	 * \return vector of iterators to valley elements
	 *
	 */
	[[gnu::warn_unused_result]] std::vector<std::vector<float>::const_iterator> getValleys(const std::vector<float> &values, const float &threshold);

	/** \brief Read match status along the reference
	 *
	 * Parses CIGAR to track read (query) match/mismatch (1.0 for match, 0.0 for mismatch) status along the reference.
	 * This means that insertions in the read are ignored.
	 * The vector start begins at the position closest to the first exon start regardless of strand.
	 *
	 * \param[in] cigar CIGAR vector
	 * \return vector of match status
	 */
	[[gnu::warn_unused_result]] std::vector<float> getReferenceMatchStatus(const std::vector<uint32_t> &cigar);

	/** \brief Extract exon coverage statistics for a read 
	 *
	 * \param[in] readAndExons read alignment with the corresponding exon group
	 * \return exon coverage object
	 */
	[[gnu::warn_unused_result]] ReadExonCoverage getExonCoverageStats(const std::pair<BAMrecord, ExonGroup> &readAndExons);

	/** \brief Convert `ReadExonCoverage` to string
	 *
	 * \param[in] readRecord individual read record
	 * \param[in] separator field separator
	 * \return `std::string` with the read record elements, without a new line at the end
	 */
	[[gnu::warn_unused_result]] std::string stringifyExonCoverage(const ReadExonCoverage &readRecord, char separator = '\t');

	/** \brief Produce a string from a range of read alignments
	 *
	 * \param[in] begin start iterator
	 * \param[in] end end iterator
	 * \return string with coverage information
	 */
	[[gnu::warn_unused_result]] std::string stringifyAlignmentRange(const bamGFFvector::const_iterator &begin, const bamGFFvector::const_iterator &end);
	/** \brief Produce a string of poorly aligned region statistics from an alignment range 
	 *
	 * Only saves information from reads that have poorly mapped regions, potentially multiple per read.
	 *
	 * \param[in] begin start iterator
	 * \param[in] end end iterator
	 * \param[in] windowParameters sliding window parameters
	 * \return string with coverage information
	 */
	[[gnu::warn_unused_result]] std::string stringifyUnmappedRegions(const bamGFFvector::const_iterator &begin, const bamGFFvector::const_iterator &end, const BinomialWindowParameters &windowParameters);
	/** \brief Produce a string of poorly aligned region statistics and corresponding FASTQ records from an alignment range 
	 *
	 * Only saves information from reads that have poorly mapped regions, potentially multiple per read.
	 *
	 * \param[in] begin start iterator
	 * \param[in] end end iterator
	 * \param[in] windowParameters sliding window parameters
	 * \return strings with coverage information (`.first`) and FASTQ (`.second`)
	 */
	[[gnu::warn_unused_result]] std::pair<std::string, std::string> getUnmappedRegionsAndFASTQ(const bamGFFvector::const_iterator &begin, const bamGFFvector::const_iterator &end, const BinomialWindowParameters &windowParameters);

	/** \brief Extract original read name and coordinates
	 *
	 * The remapped read names are original names, with `_`-separated start and end coordinates.
	 * Underscores in the original read name are allowed.
	 * If the coordinates are absent, returns an empty object.
	 *
	 * \param[in] remappedReadName re-mapped read portion name
	 * \return original read name and segment coordinates
	 */
	[[gnu::warn_unused_result]] ReadPortion parseRemappedReadName(const std::string &remappedReadName);
	/** \brief Modify the BAM CIGAR string to erase alignment in a range 
	 *
	 * Substitute operations in a BAM record within the given range with non-matching operations.
	 *
	 * \param[in] modRange modification range
	 * \param[in] bamRecord BAM record to be modified
	 * \return BAM record with the CIGAR vector replaced
	 */
	[[gnu::warn_unused_result]] std::unique_ptr<bam1_t, BAMrecordDeleter> modifyCIGAR(const ReadPortion &modRange, const std::unique_ptr<bam1_t, BAMrecordDeleter> &bamRecord);
	/** \brief Add a re-mapped secondary alignment
	 *
	 * Add a read portion remap as a secondary alignment to a vector of BAM records.
	 * Only reads that pass the identity threshold are added.
	 *
	 * \param[in] newRecordHeader header corresponding to the re-mapped record
	 * \param[in] newRecord remapped BAM record
	 * \param[in] remapInfo original name and read segment range
	 * \param[in] originalHeader header corresponding to the original record
	 * \param[in] remapIdentityCutoff fraction of sites in the remapped read that are identical to the reference
	 * \param[in,out] readMapVector vector of alignments of a read, first element is the primary alignment
	 */
	void addRemappedSecondaryAlignment(
			const std::unique_ptr<sam_hdr_t, BAMheaderDeleter> &newRecordHeader, const std::unique_ptr<bam1_t, BAMrecordDeleter> &newRecord, const ReadPortion &remapInfo,
			const std::unique_ptr<sam_hdr_t, BAMheaderDeleter> &originalHeader, const float &remapIdentityCutoff, std::vector< std::unique_ptr<bam1_t, BAMrecordDeleter> > &readMapVector);

	/** \brief Open a BGZF file handle for appending
	 *
	 * Opens a handle to the BAM file for appending, deleting the original if it exists.
	 *
	 * \param[in] bamFileName BAM file name
	 */
	std::unique_ptr<BGZF, BGZFhandleDeleter> openBGZFtoAppend(const std::string &bamFileName);

	/** \brief Make per-thread alignment record/annotation vector ranges
	 *
	 * Constructs a vector of iterator pairs bracketing chunks of a vector to be processed in parallel.
	 *
	 * \param[in] targetVector the vector to be processed
	 * \param[in] threadCount number of threads
	 *
	 * \return vector of iterator pairs for each thread
	 */
	[[gnu::warn_unused_result]] std::vector< std::pair<bamGFFvector::const_iterator, bamGFFvector::const_iterator> > 
		makeThreadRanges(const bamGFFvector &targetVector, const size_t &threadCount);

	/** \brief Command line parser
	 *
	 * Maps flags to values. Flags assumed to be of the form `--flag-name value`.
	 *
	 * \param[in] argc size of the `argv` array
	 * \param[in] argv command line input array
	 * \return map of tags to values
	 */
	[[gnu::warn_unused_result]] std::unordered_map<std::string, std::string> parseCL(int &argc, char **argv);

	/** \brief Extract parameters from parsed command line interface flags
	 *
	 * Extracts needed variable values, indexed by `std::string` encoded variable names.
	 *
	 * \param[in] parsedCLI flag values parsed from the command line
	 * \param[out] intVariables indexed `int` variables for use by `main()`
	 * \param[out] floatVariables indexed `float` variables for use by `main()`
	 * \param[out] stringVariables indexed `std::string` variables for use by `main()`
	 */
	void extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI,
			std::unordered_map<std::string, int> &intVariables,
			std::unordered_map<std::string, float> &floatVariables,
			std::unordered_map<std::string, std::string> &stringVariables);
}
