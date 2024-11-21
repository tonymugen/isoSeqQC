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
	/** \brief Extract a name 
	 *
	 * Extract a token name from GFF attributes.
	 *
	 * \param[in] tokenAndAttrList field token and the list of attributes
	 */
	[[gnu::warn_unused_result]] std::string extractAttributeName(const TokenAttibuteListPair &tokenAndAttrList);

	/** \brief Test for range overlap
	 *
	 * Checks if the BAM record overlaps the nucleotide range covered by a gene.
	 *
	 * \param[in] geneInfo gene information
	 * \param[in] candidateBAM BAM record
	 * \return `true` if there is overlap
	 */
	[[gnu::warn_unused_result]] bool rangesOverlap(const ReadExonCoverage &geneInfo, const BAMrecord &candidateBAM) noexcept;

	/** \brief Convert `ReadExonCoverage` to string
	 *
	 * \param[in] readRecord individual read record
	 * \param[in] separator field separator
	 * \return `std::string` with the read record elements, without a new line at the end
	 */
	[[gnu::warn_unused_result]] std::string stringify(const ReadExonCoverage &readRecord, char separator = '\t');

	/** \brief Produce a string from a range of `ReadExonCoverage` elements
	 *
	 * \param[in] begin start iterator
	 * \param[in] end end iterator
	 * \return string with coverage information
	 */
	[[gnu::warn_unused_result]] std::string stringifyRCSrange(const std::vector<ReadExonCoverage>::const_iterator &begin, const std::vector<ReadExonCoverage>::const_iterator &end);

	/** \brief Make per-thread `ReadExonCoverage` vector ranges
	 *
	 * Constructs a vector of iterator pairs bracketing chunks of a vector to be processed in parallel.
	 *
	 * \param[in] targetVector the vector to be processed
	 * \param[in] threadCount number of threads
	 *
	 * \return vector of iterator pairs for each thread
	 */
	[[gnu::warn_unused_result]] std::vector< std::pair<std::vector<ReadExonCoverage>::const_iterator, std::vector<ReadExonCoverage>::const_iterator> > 
		makeThreadRanges(const std::vector<ReadExonCoverage> &targetVector, const size_t &threadCount);

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
