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

/// Genomic analyses helper functions
/** \file
 * \author Anthony J. Greenberg and Rebekah Rogers
 * \copyright Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.1
 *
 * Definitions of class-external functions needed by genomic analyses.
 *
 */

#pragma once

#include <string>
#include <utility> // for std::pair
#include <vector>

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

	/** \brief Test for range overlap
	 *
	 * The first position in each range need not be before the second in the same range.
	 *
	 * \param[in] range1 first range
	 * \param[in] range2 second range
	 * \return `true` if there is overlap
	 */
	[[gnu::warn_unused_result]] bool rangesOverlap(const std::pair<hts_pos_t, hts_pos_t> &range1, const std::pair<hts_pos_t, hts_pos_t> &range2) noexcept;

	/** \brief Extract mRNA information from GFF 
	 *
	 * Extract information from a GFF line that has mRNA data.
	 * Changes the latest gene name if the parent of the current mRNA is different and updates the exon groups vector if necessary.
	 *
	 * \param[in,out] currentGFFline fields from the GFF file line just read
	 * \param[in,out] previousGFFfields GFF fields from previous lines; the attribute field only has the current gene ID
	 * \param[in,out] exonSpanSet set of unique exon start/end pairs
	 * \return `ExonGroup` object
	 */
	[[gnu::warn_unused_result]] ExonGroup mRNAfromGFF(std::array<std::string, nGFFfields> &currentGFFline, std::array<std::string, nGFFfields> &previousGFFfields, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSpanSet);

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
	[[gnu::warn_unused_result]] std::string stringifyAlignementRange(
			const std::vector< std::pair<BAMrecord, ExonGroup> >::const_iterator &begin,
			const std::vector< std::pair<BAMrecord, ExonGroup> >::const_iterator &end
		);

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
	 * \param[out] stringVariables indexed `std::string` variables for use by `main()`
	 */
	void extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI,
			std::unordered_map<std::string, int> &intVariables, std::unordered_map<std::string, std::string> &stringVariables);
}
