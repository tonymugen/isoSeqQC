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
 * Implementation of class-external functions needed by genomic analyses.
 *
 */

#include <cassert>
#include <algorithm>
#include <array>
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

std::string isaSpace::extractParentName(const std::string &attributeString) {
	const std::string parentToken{"Parent="};
	constexpr char attrDelimiter{';'};
	std::stringstream attributeStream(attributeString);
	std::string attrField;
	TokenAttibuteListPair tokenPair;
	while ( std::getline(attributeStream, attrField, attrDelimiter) ) {
		tokenPair.attributeList.emplace_back(attrField);
	}
	tokenPair.tokenName = parentToken;

	return extractAttributeName(tokenPair);
}

bool isaSpace::rangesOverlap(const std::pair<hts_pos_t, hts_pos_t> &range1, const std::pair<hts_pos_t, hts_pos_t> &range2) noexcept {
	const auto local1 = std::minmax(range1.first, range1.second);
	const auto local2 = std::minmax(range2.first, range2.second);
	return (local2.first <= local1.second) && (local2.second >= local1.first);
}

std::array<std::string, nGFFfields> isaSpace::parseGFFline(const std::string &gffLine) {
	constexpr char gffDelimiter{'\t'};
	std::array<std::string, nGFFfields> result;
	std::stringstream lineStream(gffLine);
	size_t iField{0};
	while ( (iField < nGFFfields) && std::getline(lineStream, result.at(iField), gffDelimiter) ) {
		++iField;
	}
	// skip incomplete lines
	if (iField < nGFFfields) {
		result.at(0) = "FAIL";
	}

	return result;
}

std::unordered_map< std::string, std::vector<ExonGroup> > isaSpace::parseGFF(const std::string &gffFileName) {
	constexpr size_t strandIDidx{6UL};
	constexpr size_t typeIDidx{2UL};

	std::unordered_map< std::string, std::vector<ExonGroup> > result;
	std::string gffLine;
	std::fstream gffStream(gffFileName, std::ios::in);
	std::string currentGeneName;
	std::string currentStrandedChromosome;
	char currentStrandID{'+'};
	std::vector<std::string> exonLines;
	while ( std::getline(gffStream, gffLine) ) {
		if ( gffLine.empty() || (gffLine.at(0) == '#') ) {
			continue;
		}
		std::array<std::string, nGFFfields> gffFields{parseGFFline(gffLine)};
		// skip incomplete lines
		if (gffFields.at(0) == "FAIL") {
			continue;
		}
		if (gffFields.at(typeIDidx) == "exon") {
			exonLines.emplace_back(gffLine);
			continue;
		}
		if (gffFields.at(2) == "mRNA") {
			std::string thisGeneName{extractParentName( gffFields.back() )};
			if (thisGeneName == currentGeneName) {
				continue;
			}
			if ( currentGeneName.empty() ) { // first mRNA in the file
				currentGeneName           = thisGeneName;
				currentStrandID           = (gffFields.at(strandIDidx) == "-" ? '-' : '+');
				currentStrandedChromosome = gffFields.front() + currentStrandID;
				continue;
			}
			result[currentStrandedChromosome].emplace_back(currentGeneName, currentStrandID, exonLines);
			currentStrandID           = (gffFields.at(strandIDidx) == "-" ? '-' : '+');
			currentStrandedChromosome = gffFields.front() + currentStrandID;
			currentGeneName           = thisGeneName;
			exonLines.clear();
		}
	}
	if ( ( !exonLines.empty() ) && ( !currentGeneName.empty() ) ) {
		result[currentStrandedChromosome].emplace_back(currentGeneName, currentStrandID, exonLines);
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
	currentAlignmentInfo.readName       = readAndExons.first.getReadName();
	currentAlignmentInfo.strand         = strandID;
	currentAlignmentInfo.chromosomeName = "not_mapped";
	if ( readAndExons.first.isMapped() && (firstCIGAR > 0) ) {
		const MappedReadMatchStatus bestMatchStats{readAndExons.first.getBestReferenceMatchStatus()};
		const std::pair<hts_pos_t, hts_pos_t> geneSpan{readAndExons.second.geneSpan()};
		currentAlignmentInfo.chromosomeName           = readAndExons.first.getReferenceName();
		currentAlignmentInfo.alignmentStart           = readAndExons.first.getMapStart();
		currentAlignmentInfo.alignmentEnd             = readAndExons.first.getMapEnd();
		currentAlignmentInfo.bestAlignmentStart       = bestMatchStats.mapStart;
		currentAlignmentInfo.bestAlignmentEnd         = bestMatchStats.mapStart + static_cast<hts_pos_t>( bestMatchStats.matchStatus.size() );
		currentAlignmentInfo.firstSoftClipLength      = bam_cigar_oplen(firstCIGAR) * static_cast<uint32_t>(bam_cigar_opchr(firstCIGAR) == 'S');
		currentAlignmentInfo.nSecondaryAlignments     = readAndExons.first.secondaryAlignmentCount();
		currentAlignmentInfo.nLocalReversedAlignments = readAndExons.first.localReversedSecondaryAlignmentCount();
		currentAlignmentInfo.nGoodSecondaryAlignments = readAndExons.first.localSecondaryAlignmentCount() - currentAlignmentInfo.nLocalReversedAlignments;
		currentAlignmentInfo.geneName                 = readAndExons.second.geneName();
		currentAlignmentInfo.nExons                   = readAndExons.second.nExons();
		currentAlignmentInfo.firstExonStart           = geneSpan.first;
		currentAlignmentInfo.lastExonEnd              = geneSpan.second;
		currentAlignmentInfo.firstExonLength          = readAndExons.second.firstExonLength();
		currentAlignmentInfo.exonCoverageScores       = readAndExons.second.getExonCoverageQuality(readAndExons.first);
		currentAlignmentInfo.bestExonCoverageScores   = readAndExons.second.getBestExonCoverageQuality(readAndExons.first);
	}

	if ( currentAlignmentInfo.geneName.empty() ) {
		currentAlignmentInfo.geneName = "no_overlap";
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

std::string isaSpace::stringifyAlignmentRange(const bamGFFvector::const_iterator &begin, const bamGFFvector::const_iterator &end) {
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

std::string isaSpace::stringifyUnmappedRegions(const bamGFFvector::const_iterator &begin, const bamGFFvector::const_iterator &end, const BinomialWindowParameters &windowParameters) {
	std::string badAlignmentString;
	std::for_each(
		begin,
		end,
		[&badAlignmentString, &windowParameters](const std::pair<BAMrecord, ExonGroup> &currentRAG) {
			const std::vector<MappedReadInterval> badRegions{currentRAG.first.getPoorlyMappedRegions(windowParameters)};
			std::for_each(
				badRegions.cbegin(),
				badRegions.cend(),
				[&badAlignmentString, &windowParameters, &currentRAG](const MappedReadInterval &eachInterval) {
					std::string regionStats =
						currentRAG.first.getReadName() + "\t" +
						std::to_string( currentRAG.first.getReadLength() ) + "\t" +
						std::to_string(eachInterval.readStart) + "\t" +
						std::to_string(eachInterval.readEnd) + "\t" +
						std::to_string(windowParameters.windowSize) + "\n";
					badAlignmentString += regionStats;
				}
			);
		}
	);
	return badAlignmentString;
}

std::pair<std::string, std::string> isaSpace::getUnmappedRegionsAndFASTQ(const bamGFFvector::const_iterator &begin, const bamGFFvector::const_iterator &end, const BinomialWindowParameters &windowParameters) {
	std::pair<std::string, std::string> badAlignmentInfoFQ;
	std::for_each(
		begin,
		end,
		[&badAlignmentInfoFQ, &windowParameters](const std::pair<BAMrecord, ExonGroup> &currentRAG) {
			const std::vector<MappedReadInterval> badRegions{currentRAG.first.getPoorlyMappedRegions(windowParameters)};
			std::for_each(
				badRegions.cbegin(),
				badRegions.cend(),
				[&badAlignmentInfoFQ, &windowParameters, &currentRAG](const MappedReadInterval &eachInterval) {
					std::string regionStats =
						currentRAG.first.getReadName() + "\t" +
						std::to_string( currentRAG.first.getReadLength() ) + "\t" +
						std::to_string(eachInterval.readStart) + "\t" +
						std::to_string(eachInterval.readEnd) + "\t" +
						std::to_string(windowParameters.windowSize) + "\n";
					badAlignmentInfoFQ.first += regionStats;
					std::string fastqString{currentRAG.first.getSequenceAndQuality(eachInterval)};
					fastqString =
						"@"  + currentRAG.first.getReadName() +
						"_"  + std::to_string(eachInterval.readStart) +
						"_"  + std::to_string(eachInterval.readEnd) +
						"\n" + fastqString;
					badAlignmentInfoFQ.second += fastqString;
				}
			);
		}
	);
	return badAlignmentInfoFQ;
}

ReadPortion isaSpace::parseRemappedReadName(const std::string &remappedReadName) {
	ReadPortion result;
	if ( remappedReadName.empty() ) {
		return result;
	}
	constexpr size_t nStartEnd{2};
	std::array<std::string, nStartEnd> startAndEndStrings;
	auto trailerBeginIt = std::prev( remappedReadName.end() );
	size_t iRange{0};
	while ( (iRange < nStartEnd) && ( trailerBeginIt != remappedReadName.begin() ) ) {
		if (*trailerBeginIt == '_') {
			++iRange;
			--trailerBeginIt;
			continue;
		}
		startAndEndStrings.at(iRange) = *trailerBeginIt + startAndEndStrings.at(iRange); // reversing the order since we are moving from the back
		--trailerBeginIt;
	}

	result.originalName = std::string( remappedReadName.begin(), std::next(trailerBeginIt) ); // next because we advance back on '_'
	try {
		// start is in the back of the array because we are reading back to front
		result.start = std::stoul( startAndEndStrings.at(1) );
		result.end   = std::stoul( startAndEndStrings.at(0) );
	} catch(const std::exception &invalid) {
		result.originalName.clear();
		result.start = 0;
		result.end   = 0;
		return result;
	}
	return result;
}

void isaSpace::modifyCIGAR(const ReadPortion &modRange, std::unique_ptr<bam1_t, BAMrecordDeleter> &bamRecord) {
}

void isaSpace::addRemappedSecondaryAlignment(
		const std::unique_ptr<sam_hdr_t, BAMheaderDeleter> &newRecordHeader, const std::unique_ptr<bam1_t, BAMrecordDeleter> &newRecord, const ReadPortion &remapInfo,
		const std::unique_ptr<sam_hdr_t, BAMheaderDeleter> &originalHeader, const float &remapIdentityCutoff, std::vector< std::unique_ptr<bam1_t, BAMrecordDeleter> > &readMapVector) {
	BAMrecordDeleter newSecondaryDeleter;
	std::unique_ptr<bam1_t, BAMrecordDeleter> secondaryFromNew(bam_init1(), newSecondaryDeleter);

	// reference indexes may not be the same in the original and new headers
	// retrieve the original index from the reference name
	const auto *const newRefNamePtr = sam_hdr_tid2name(newRecordHeader.get(), newRecord->core.tid); 
	if (newRefNamePtr == nullptr) {
		return;
	}
	const int32_t originalTID = bam_name2id(originalHeader.get(), newRefNamePtr);
	if (originalTID < 0) {
		return;
	}

	constexpr std::array<float, 10> sequenceMatch{
		1.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0, 0.0
	};
	constexpr std::array<float, 10> readConsumption{
		1.0, 1.0, 0.0, 0.0, 1.0,
		0.0, 0.0, 1.0, 1.0, 0.0
	};
	// Add soft clips to the CIGAR string of the remapped portion if necessary
	auto softClip = static_cast<uint32_t>(remapInfo.start);
	std::vector<uint32_t> remapCIGAR(
		static_cast<uint32_t>(remapInfo.start > 0),
		bam_cigar_gen(softClip, BAM_CSOFT_CLIP)
	);
	// only the remapped portion of the read counts towards the identity fraction calculation
	float matchCount{0.0};
	float readLength{0.0};
	for (uint32_t iCIGAR = 0; iCIGAR < newRecord->core.n_cigar; ++iCIGAR) {
		remapCIGAR.push_back( *(bam_get_cigar( newRecord.get() ) + iCIGAR) );
		const auto cigarOpLength = static_cast<float>( bam_cigar_oplen( remapCIGAR.back() ) );
		matchCount += cigarOpLength * sequenceMatch.at( bam_cigar_op( remapCIGAR.back() ) );
		readLength += cigarOpLength * readConsumption.at( bam_cigar_op( remapCIGAR.back() ) );
	}
	if (matchCount/readLength < remapIdentityCutoff) {
		return;
	}
	if (remapInfo.end <= readMapVector.front()->core.l_qseq) {
		softClip = static_cast<uint32_t>(readMapVector.front()->core.l_qseq - remapInfo.end);
		remapCIGAR.push_back( bam_cigar_gen(softClip, BAM_CSOFT_CLIP) );
	}
	assert(
		bam_cigar2qlen( static_cast<int32_t>( remapCIGAR.size() ), remapCIGAR.data() ) == 
			bam_cigar2qlen(
				static_cast<int32_t>(readMapVector.front()->core.n_cigar), bam_get_cigar( readMapVector.front().get() )
			) &&
		"ERROR: original and new CIGAR strings imply different read lengths"
	);
	const int32_t success = bam_set1(
		secondaryFromNew.get(),
		remapInfo.originalName.size(),
		remapInfo.originalName.c_str(),
		newRecord->core.flag | BAM_FSECONDARY,
		originalTID,
		newRecord->core.pos,
		newRecord->core.qual,
		remapCIGAR.size(),
		remapCIGAR.data(),
		0, 0, newRecord->core.isize, 0, nullptr, nullptr, 0
	);

	if (success >= 0) {
		readMapVector.emplace_back( std::move(secondaryFromNew) );
	}
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
	const std::array<std::string, 1> optionalStringVariables{"out-fastq"};
	const std::array<std::string, 2> optionalIntVariables{"threads", "window-size"};

	const std::unordered_map<std::string, int> defaultIntValues{ {"threads", -1}, {"window-size", 75} };
	const std::unordered_map<std::string, std::string> defaultStringValues{ {"out-fastq", "NULL"} };

	if ( parsedCLI.empty() ) {
		throw std::string("No command line flags specified;");
	}
	for (const auto &eachFlag : optionalIntVariables) {
		try {
			intVariables[eachFlag] = stoi( parsedCLI.at(eachFlag) );
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
	for (const auto &eachFlag : optionalStringVariables) {
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			stringVariables[eachFlag] = defaultStringValues.at(eachFlag);
		}
	}
}
