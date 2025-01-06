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
#include <cmath>
#include <iterator>
#include <numeric>
#include <string>
#include <sstream>
#include <fstream>
#include <utility>

#include "sam.h"

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

bool isaSpace::rangesOverlap(const std::pair<hts_pos_t, hts_pos_t> &range1, const std::pair<hts_pos_t, hts_pos_t> &range2) noexcept {
	const auto local1 = std::minmax(range1.first, range1.second);
	const auto local2 = std::minmax(range2.first, range2.second);
	return (local2.first <= local1.second) && (local2.second >= local1.first);
}

ExonGroup isaSpace::mRNAfromGFF(std::array<std::string, nGFFfields> &currentGFFline, std::array<std::string, nGFFfields> &previousGFFfields, std::set< std::pair<hts_pos_t, hts_pos_t> > &exonSpanSet) {
	const std::string parentToken{"Parent="};
	constexpr char attrDelimiter{';'};
	constexpr size_t strandIDidx{6UL};
	ExonGroup result;
	std::stringstream attributeStream( currentGFFline.back() );
	std::string attrField;
	TokenAttibuteListPair tokenPair;
	while ( std::getline(attributeStream, attrField, attrDelimiter) ) {
		tokenPair.attributeList.emplace_back(attrField);
	}
	tokenPair.tokenName = parentToken;

	std::string parentName{extractAttributeName(tokenPair)};
	std::swap(currentGFFline.back(), parentName);
	if ( currentGFFline.back() == previousGFFfields.back() ) { // both could be empty
		return result;
	}
	if ( currentGFFline.back().empty() ) {
		if ( !exonSpanSet.empty() ) {
			const char strandID                  = (previousGFFfields.at(strandIDidx) == "-" ? '-' : '+');
			const std::string strandedChromosome = previousGFFfields.front() + strandID;
			result = ExonGroup(previousGFFfields.back(), previousGFFfields.at(strandIDidx).front(), exonSpanSet);
			exonSpanSet.clear();
		}
		previousGFFfields.back().clear();
		previousGFFfields.front().clear();
		return result;
	}
	// neither is empty and are not the same
	if ( !exonSpanSet.empty() ) {
		const char strandID                  = (previousGFFfields.at(strandIDidx) == "-" ? '-' : '+');
		const std::string strandedChromosome = previousGFFfields.front() + strandID;
		result = ExonGroup(previousGFFfields.back(), previousGFFfields.at(strandIDidx).front(), exonSpanSet);
		exonSpanSet.clear();
	}
	std::copy( currentGFFline.cbegin(), currentGFFline.cend(), previousGFFfields.begin() );
	return result;
}

std::unordered_map< std::string, std::vector<ExonGroup> > isaSpace::parseGFF(const std::string &gffFileName) {
	constexpr size_t strandIDidx{6UL};
	constexpr char   gffDelimiter{'\t'};
	constexpr size_t spanStart{3UL};
	constexpr size_t spanEnd{4UL};

	std::unordered_map< std::string, std::vector<ExonGroup> > result;
	std::string gffLine;
	std::fstream gffStream(gffFileName, std::ios::in);
	std::set< std::pair<hts_pos_t, hts_pos_t> > exonSpans;
	std::array<std::string, nGFFfields> activeGFFfields;   // fields from the gene tracked up to now
	std::array<std::string, nGFFfields> newGFFfields;      // fields from the currently read line
	while ( std::getline(gffStream, gffLine) ) {
		if ( gffLine.empty() || (gffLine.at(0) == '#') ) {
			continue;
		}
		std::stringstream currLineStream(gffLine);
		size_t iField{0};
		while ( (iField < nGFFfields) && std::getline(currLineStream, newGFFfields.at(iField), gffDelimiter) ) {
			++iField;
		}
		// skip incomplete lines
		if (iField < nGFFfields) {
			continue;
		}
		if (newGFFfields.at(2) == "mRNA") {
			if ( activeGFFfields.front().empty() ) { // this is the first mRNA in the file
				activeGFFfields = newGFFfields;
				continue;
			}
			const char strandID                  = (activeGFFfields.at(strandIDidx) == "-" ? '-' : '+');
			const std::string strandedChromosome = activeGFFfields.front() + strandID;
			result[strandedChromosome].emplace_back( mRNAfromGFF(newGFFfields, activeGFFfields, exonSpans) );
			continue;
		}
		if ( (newGFFfields.at(2) == "exon") && ( !activeGFFfields.back().empty() ) ) {
			const auto exonStart = static_cast<hts_pos_t>( stol( newGFFfields.at(spanStart) ) );
			const auto exonEnd   = static_cast<hts_pos_t>( stol( newGFFfields.at(spanEnd) ) );
			exonSpans.insert({exonStart, exonEnd});
		}
	}
	if ( !activeGFFfields.back().empty() && !exonSpans.empty() ) {
		const char strandID                  = (activeGFFfields.at(strandIDidx) == "-" ? '-' : '+');
		const std::string strandedChromosome = activeGFFfields.front() + strandID;
		result[strandedChromosome].emplace_back(activeGFFfields.back(), activeGFFfields.at(strandIDidx).front(), exonSpans);
	}
	return result;
}

std::vector<std::vector<float>::const_iterator> isaSpace::getPeaks(const std::vector<float> &values, const float &threshold) {
	std::vector<std::vector<float>::const_iterator> result;
	auto peakIt = values.cbegin();
	while ( peakIt != values.cend() ) {
		const auto peakBeginIt = std::find_if(peakIt,      values.cend(), [&threshold](float value){return value >= threshold;});
		const auto peakEndIt   = std::find_if(peakBeginIt, values.cend(), [&threshold](float value){return value < threshold;});
		auto localPeakIt       = std::max_element(peakBeginIt, peakEndIt);
		result.push_back(localPeakIt);
		peakIt = peakEndIt;
	}

	return result;
}

std::vector<std::vector<float>::const_iterator> isaSpace::getValleys(const std::vector<float> &values, const float &threshold) {
	std::vector<std::vector<float>::const_iterator> result;
	auto peakIt = values.cbegin();
	while ( peakIt != values.cend() ) {
		const auto peakBeginIt = std::find_if(peakIt,      values.cend(), [&threshold](float value){return value <= threshold;});
		const auto peakEndIt   = std::find_if(peakBeginIt, values.cend(), [&threshold](float value){return value > threshold;});
		auto localPeakIt       = std::min_element(peakBeginIt, peakEndIt);
		result.push_back(localPeakIt);
		peakIt = peakEndIt;
	}

	return result;
}

std::vector<float> isaSpace::getReferenceMatchStatus(const std::vector<uint32_t> &cigar) {
	constexpr std::array<hts_pos_t, 10> referenceConsumption{
		1, 0, 1, 1, 0,
		0, 0, 1, 1, 0
	};
	constexpr std::array<float, 10> sequenceMatch{
		1.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0, 0.0
	};
	std::vector<float> referenceMatchStatus;
	for (const auto &eachCIGAR : cigar) {
		const std::vector<float> currCIGARfield(
			bam_cigar_oplen(eachCIGAR) * referenceConsumption.at( bam_cigar_op(eachCIGAR) ),
			sequenceMatch.at( bam_cigar_op(eachCIGAR) )
		);
		std::copy( currCIGARfield.cbegin(), currCIGARfield.cend(), std::back_inserter(referenceMatchStatus) );
	}
	return referenceMatchStatus;
}

ReadExonCoverage isaSpace::getExonCoverageStats(const std::pair<BAMrecord, ExonGroup> &readAndExons) {
	const char strandID = (readAndExons.first.isRevComp() ? '-' : '+' );
	const auto firstCIGAR{readAndExons.first.getFirstCIGAR()};
	ReadExonCoverage currentAlignmentInfo;
	if ( readAndExons.first.isMapped() && (firstCIGAR > 0) ) {
		currentAlignmentInfo.readName                 = readAndExons.first.getReadName();
		currentAlignmentInfo.chromosomeName           = readAndExons.first.getReferenceName();
		currentAlignmentInfo.strand                   = strandID;
		currentAlignmentInfo.alignmentStart           = readAndExons.first.getMapStart();
		currentAlignmentInfo.alignmentEnd             = readAndExons.first.getMapEnd();
		currentAlignmentInfo.bestAlignmentStart       = readAndExons.first.getMapStart();
		currentAlignmentInfo.bestAlignmentEnd         = readAndExons.first.getMapEnd();
		currentAlignmentInfo.firstSoftClipLength      = bam_cigar_oplen(firstCIGAR) * static_cast<uint32_t>(bam_cigar_opchr(firstCIGAR) == 'S');
		currentAlignmentInfo.nSecondaryAlignments     = readAndExons.first.secondaryAlignmentCount();
		currentAlignmentInfo.nLocalReversedAlignments = readAndExons.first.localReversedSecondaryAlignmentCount();
		currentAlignmentInfo.nGoodSecondaryAlignments = readAndExons.first.localSecondaryAlignmentCount() - currentAlignmentInfo.nLocalReversedAlignments;
		currentAlignmentInfo.exonCoverageScores       = readAndExons.second.getExonCoverageQuality(readAndExons.first);
		currentAlignmentInfo.bestExonCoverageScores   = readAndExons.second.getBestExonCoverageQuality(readAndExons.first);
	}

	return currentAlignmentInfo;
}

std::string isaSpace::stringifyExonCoverage(const ReadExonCoverage &readRecord, char separator) {
	std::string coverages = std::accumulate(
		readRecord.exonCoverageScores.cbegin(),
		readRecord.exonCoverageScores.cend(),
		std::string("{"),
		[](std::string strVal, float val) {
			return std::move(strVal) + std::to_string(val) + ',';
		}
	);
	coverages.resize(std::max( 2UL, coverages.size() ) - 1);
	coverages += "}";
	std::string bestCoverages = std::accumulate(
		readRecord.bestExonCoverageScores.cbegin(),
		readRecord.bestExonCoverageScores.cend(),
		std::string("{"),
		[](std::string strVal, float val) {
			return std::move(strVal) + std::to_string(val) + ',';
		}
	);
	bestCoverages.resize(std::max( 2UL, bestCoverages.size() ) - 1);
	bestCoverages += "}";
	std::string result =  readRecord.readName                                 + separator
						+ readRecord.chromosomeName                           + separator
						+ readRecord.strand                                   + separator
						+ std::to_string(readRecord.alignmentStart)           + separator
						+ std::to_string(readRecord.alignmentEnd)             + separator
						+ std::to_string(readRecord.bestAlignmentStart)       + separator
						+ std::to_string(readRecord.bestAlignmentEnd)         + separator
						+ std::to_string(readRecord.firstSoftClipLength)      + separator
						+ std::to_string(readRecord.nSecondaryAlignments)     + separator
						+ std::to_string(readRecord.nGoodSecondaryAlignments) + separator
						+ std::to_string(readRecord.nLocalReversedAlignments) + separator
						+ readRecord.geneName                                 + separator
						+ std::to_string(readRecord.nExons)                   + separator
						+ std::to_string(readRecord.firstExonLength)          + separator
						+ std::to_string(readRecord.firstExonStart)           + separator
						+ std::to_string(readRecord.lastExonEnd)              + separator
						+ coverages + separator + bestCoverages;

	return result;
}

std::string isaSpace::stringifyAlignementRange(
		const std::vector< std::pair<BAMrecord, ExonGroup> >::const_iterator &begin,
		const std::vector< std::pair<BAMrecord, ExonGroup> >::const_iterator &end) {
	std::string outString;
	std::for_each(
		begin,
		end,
		[&outString](const  std::pair<BAMrecord, ExonGroup>  &currentAlignment) {
			const ReadExonCoverage exonCoverage = getExonCoverageStats(currentAlignment);
			outString += stringifyExonCoverage(exonCoverage) + "\n";
		}
	);
	return outString;
}

std::vector< std::pair<bamGFFvector::const_iterator, bamGFFvector::const_iterator> > 
											isaSpace::makeThreadRanges(const bamGFFvector &targetVector, const size_t &threadCount) {
	std::vector<bamGFFvector::difference_type> chunkSizes(
		threadCount,
		static_cast<bamGFFvector::difference_type>(targetVector.size() / threadCount)
	);
	// spread the left over elements among chunks
	std::for_each(
		chunkSizes.begin(),
		chunkSizes.begin() +
			static_cast<bamGFFvector::difference_type >(targetVector.size() % threadCount),
		[](bamGFFvector::difference_type &currSize) {return ++currSize;}
	);
	std::vector<
		std::pair<
			bamGFFvector::const_iterator,
			bamGFFvector::const_iterator
		>
	> threadRanges;
	auto chunkBeginIt = targetVector.cbegin();

	std::for_each(
		chunkSizes.cbegin(),
		chunkSizes.cend(),
		[&chunkBeginIt, &threadRanges](bamGFFvector::difference_type currDiff) {
			std::pair<
				bamGFFvector::const_iterator,
				bamGFFvector::const_iterator
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
