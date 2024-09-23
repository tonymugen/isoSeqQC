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
 * Implementation of class-external functions needed by genomic analyses.
 *
 */

#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>
#include <utility>

#include "helperFunctions.hpp"
#include "isoseqAlgn.hpp"

using namespace isaSpace;

std::string isaSpace::extractAttributeName(const TokenAttibuteListPair &tokenAndAttrList) {
	constexpr char attrDelimEscape{'\\'};
	constexpr char attrDelimiter{';'};
	auto tokenIt = std::find_if(
		tokenAndAttrList.attributeList.cbegin(),
		tokenAndAttrList.attributeList.cend(),
		[tokenAndAttrList](const std::string &eachAttr) {
			return std::equal( tokenAndAttrList.tokenName.cbegin(), tokenAndAttrList.tokenName.cend(), eachAttr.cbegin() );
		}
	);
	if ( tokenIt == tokenAndAttrList.attributeList.cend() ) {
		return "";
	}
	std::string attributeField;
	const auto tokenNameDistance = static_cast<std::string::difference_type>( tokenAndAttrList.tokenName.size() );
	std::copy(
		tokenIt->cbegin() + tokenNameDistance,
		tokenIt->cend(),
		std::back_inserter(attributeField)
	);
	// deal with any escaped ';' delimiters
	std::advance(tokenIt, 1);
	while ( (attributeField.back() == attrDelimEscape) && ( tokenIt != tokenAndAttrList.attributeList.cend() ) ) {
		attributeField.push_back(attrDelimiter);
		std::copy(
			tokenIt->cbegin(),
			tokenIt->cend(),
			std::back_inserter(attributeField)
		);
	}
	return attributeField;
}

std::string isaSpace::stringify(const ReadExonCoverage &readRecord, char separator) {
	std::string coverages = std::accumulate(
		readRecord.exonCoverageScores.cbegin(),
		readRecord.exonCoverageScores.cend(),
		std::string("{"),
		[](std::string strVal, float val) {
			return std::move(strVal) + std::to_string(val) + ',';
		}
	);
	coverages.back() = '}';
	std::string result =  readRecord.readName                                + separator
						+ readRecord.chromosomeName                          + separator
						+ readRecord.strand                                  + separator
						+ std::to_string(readRecord.alignmentStart)          + separator
						+ std::to_string(readRecord.alignmentEnd)            + separator
						+ readRecord.geneName                                + separator
						+ std::to_string(readRecord.nExons)                  + separator
						+ std::to_string(readRecord.firstExonStart)          + separator
						+ std::to_string(readRecord.lastExonEnd)             + separator
						+ coverages;

	return result;
}

std::string isaSpace::stringifyRCSrange(const std::vector<ReadExonCoverage>::const_iterator &begin, const std::vector<ReadExonCoverage>::const_iterator &end) {
	std::string outString;
	std::for_each(
		begin,
		end,
		[&outString](const ReadExonCoverage &currStat) {
			outString += stringify(currStat) + "\n";
		}
	);
	return outString;
}

std::vector< std::pair<std::vector<ReadExonCoverage>::const_iterator, std::vector<ReadExonCoverage>::const_iterator> > 
											isaSpace::makeThreadRanges(const std::vector<ReadExonCoverage> &targetVector, const size_t &threadCount) {
	std::vector<std::vector<ReadExonCoverage>::difference_type> chunkSizes(
		threadCount,
		static_cast<std::vector<ReadExonCoverage>::difference_type>(targetVector.size() / threadCount)
	);
	// spread the left over elements among chunks
	std::for_each(
		chunkSizes.begin(),
		chunkSizes.begin() +
			static_cast<std::vector<ReadExonCoverage>::difference_type >(targetVector.size() % threadCount),
		[](std::vector<ReadExonCoverage>::difference_type &currSize) {return ++currSize;}
	);
	std::vector<
		std::pair<
			std::vector<ReadExonCoverage>::const_iterator,
			std::vector<ReadExonCoverage>::const_iterator
		>
	> threadRanges;
	auto chunkBeginIt = targetVector.cbegin();

	std::for_each(
		chunkSizes.cbegin(),
		chunkSizes.cend(),
		[&chunkBeginIt, &threadRanges](std::vector<ReadExonCoverage>::difference_type currDiff) {
			std::pair<
				std::vector<ReadExonCoverage>::const_iterator,
				std::vector<ReadExonCoverage>::const_iterator
			> currItPair{chunkBeginIt, chunkBeginIt + currDiff};
			chunkBeginIt = currItPair.second;
			threadRanges.emplace_back( std::move(currItPair) );
		}
	);
	return threadRanges;
}

std::unordered_map<std::string, std::string> isaSpace::parseCL(int &argc, char **argv) {
	std::unordered_map<std::string, std::string> cliResult;
	// set to true after encountering a flag token (the characters after the dash)
	bool val = false;
	// store the token value here
	std::string curFlag;

	for (int iArg = 1; iArg < argc; iArg++) {
		const char *pchar = argv[iArg];
		if ( (pchar[0] == '-') && (pchar[1] == '-') ) { // encountered the double dash, look for the token after it
			if (val) { // A previous flag had no value
				cliResult[curFlag] = "set";
			}
			// what follows the dash?
			val     = true;
			curFlag = pchar + 2;
			continue;
		}
		if (val) {
			val                = false;
			cliResult[curFlag] = pchar;
		}
	}
	return cliResult;
}

void isaSpace::extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI,
		std::unordered_map<std::string, int> &intVariables, std::unordered_map<std::string, std::string> &stringVariables) {
	intVariables.clear();
	stringVariables.clear();
	const std::array<std::string, 3> requiredStringVariables{"input-bam", "input-gff", "out"};
	const std::array<std::string, 1> optionalIntVariables{"threads"};

	const std::unordered_map<std::string, int> defaultIntValues{ {"threads", -1} };

	if ( parsedCLI.empty() ) {
		throw std::string("No command line flags specified;");
	}
	for (const auto &eachFlag : optionalIntVariables) {
		try {
			intVariables[eachFlag] = stoi( parsedCLI.at(eachFlag));
		} catch(const std::exception &problem) {
			intVariables[eachFlag] = defaultIntValues.at(eachFlag);
		}
	}
	for (const auto &eachFlag : requiredStringVariables) {
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			throw std::string("ERROR: ") + eachFlag + std::string(" specification is required");
		}
	}
}
