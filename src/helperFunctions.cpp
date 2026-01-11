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

std::unique_ptr<bam1_t, BAMrecordDeleter> isaSpace::modifyCIGAR(const ReadPortion &modRange, const std::unique_ptr<bam1_t, BAMrecordDeleter> &bamRecord) {
	assert( (modRange.end >= modRange.start)
		&& "ERROR: end of the range must be no smaller than the start");
	assert( (bamRecord->core.n_cigar > 0)
		&& "ERROR: CIGAR string is empty");

	constexpr std::array<uint32_t, 2>  mismatchKind{BAM_CDIFF, BAM_CSOFT_CLIP};
	constexpr std::array<uint32_t, 10> readConsumption{1, 1, 0, 0, 1, 0, 0, 1, 1, 0};
	constexpr std::array<uint32_t, 10> referenceConsumption{1, 0, 1, 1, 0, 0, 0, 1, 1, 0};
	constexpr std::array<char, 16>     seqNT16str{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

	std::vector<uint32_t> newCIGAR;
	auto *oldCIGARptr = bam_get_cigar( bamRecord.get() );
	const auto readLength{
		static_cast<size_t>( bam_cigar2qlen(static_cast<int32_t>(bamRecord->core.n_cigar), oldCIGARptr) )
	};
	const uint32_t actualMismatchKind{mismatchKind.at(
		static_cast<size_t>( (modRange.start == 0) || (modRange.end == readLength) ) 
	)};

	uint32_t iCIGAR{0};
	uint32_t readConsumptionCounter{0};
	// accumulate any elements before the replacement range start
	while ( (readConsumptionCounter < modRange.start) && (iCIGAR < bamRecord->core.n_cigar) ) {
		readConsumptionCounter += bam_cigar_oplen(oldCIGARptr[iCIGAR]) * readConsumption.at( bam_cigar_op(oldCIGARptr[iCIGAR]) );
		newCIGAR.push_back(oldCIGARptr[iCIGAR]);
		++iCIGAR;
	}
	// deal with over-shoot, if any
	if (iCIGAR > 0) {
		auto lastCIGARlen   = bam_cigar_oplen( newCIGAR.back() );
		const auto overhang = readConsumptionCounter  -  std::min(static_cast<uint32_t>(modRange.start), readConsumptionCounter);

		assert( (overhang < lastCIGARlen)
			&& "ERROR: CIGAR field overhang is greater than the current field length");

		lastCIGARlen   -= overhang;
		lastCIGARlen    = bam_cigar_gen( lastCIGARlen, bam_cigar_op( newCIGAR.back() ) );
		newCIGAR.back() = lastCIGARlen;
		// I originally thought that I need to put the reference consumption counter back to modRange.start
		// but that is not correct because I would have to add back the overhang to track reference consumption
		// within the substitution range
	}
	uint32_t referenceConsumptionCounter{0};
	while ( (readConsumptionCounter < modRange.end) && (iCIGAR < bamRecord->core.n_cigar) ) {
		// the read consumption counter is necessary to calculate the remainder of a possible overshoot
		readConsumptionCounter      += bam_cigar_oplen(oldCIGARptr[iCIGAR]) * readConsumption.at( bam_cigar_op(oldCIGARptr[iCIGAR]) );
		referenceConsumptionCounter += bam_cigar_oplen(oldCIGARptr[iCIGAR]) * referenceConsumption.at( bam_cigar_op(oldCIGARptr[iCIGAR]) );
		++iCIGAR;
	}

	// advance the reference start position if the unmapped region is at the beginning
	hts_pos_t newRefPos{bamRecord->core.pos};
	if (modRange.start == 0) {
		newRefPos += referenceConsumptionCounter;
	}
	newCIGAR.push_back( bam_cigar_gen(modRange.end - modRange.start, actualMismatchKind) );

	// only deal with the remainder if we are not at the read start;
	// cannot be more than one past the end
	if (iCIGAR > 0) {
		const uint32_t currentCIGARlen{bam_cigar_oplen(oldCIGARptr[iCIGAR - 1])};
		const uint32_t remainderCIGARlen =  readConsumptionCounter - std::min(static_cast<uint32_t>(modRange.end), readConsumptionCounter);
		if (remainderCIGARlen > 0) {
			newCIGAR.push_back( bam_cigar_gen( remainderCIGARlen, bam_cigar_op(oldCIGARptr[iCIGAR - 1]) ) );
		}
	}

	while (iCIGAR < bamRecord->core.n_cigar) { // this iCIGAR is already past the previous one dealt with above
		newCIGAR.push_back(oldCIGARptr[iCIGAR]);
		++iCIGAR;
	}

	BAMrecordDeleter localDeleter;
	std::unique_ptr<bam1_t, BAMrecordDeleter> modifiedBAM(bam_init1(), localDeleter);

	// Casting from uint8_t*, should be safe
	// Necessary to match the function signature
	//auto *const seqPtr    = reinterpret_cast<char*>( bam_get_seq( bamRecord.get() ) );
	auto *const seqPtr    = bam_get_seq( bamRecord.get() );
	auto *const qualPtr   = reinterpret_cast<char*>( bam_get_qual( bamRecord.get() ) );

	// Must unpack the sequence before presenting it to bam_set1
	std::string unpackedSequence;
	for (size_t iSeq = 0; iSeq < bamRecord->core.l_qseq; ++iSeq) {
		unpackedSequence += seqNT16str.at( static_cast<char>( bam_seqi(seqPtr, iSeq) ) );
	}
	const int32_t success = bam_set1(
		modifiedBAM.get(),
		modRange.originalName.size(),
		modRange.originalName.c_str(),
		bamRecord->core.flag,
		bamRecord->core.tid,
		newRefPos,
		bamRecord->core.qual,
		newCIGAR.size(),
		newCIGAR.data(),
		bamRecord->core.mtid,
		bamRecord->core.mpos,
		bamRecord->core.isize,
		bamRecord->core.l_qseq,
		//seqPtr,
		unpackedSequence.c_str(),
		qualPtr,
		0 // no need for the aux field
	);

	assert( (success >= 0)
		&& "ERROR: failed to set new BAM record");

	return modifiedBAM;
}

void isaSpace::addRemappedSecondaryAlignment(
		const std::unique_ptr<sam_hdr_t, BAMheaderDeleter> &newRecordHeader, const std::unique_ptr<bam1_t, BAMrecordDeleter> &newRecord, const ReadPortion &remapInfo,
		const std::unique_ptr<sam_hdr_t, BAMheaderDeleter> &originalHeader, const float &remapIdentityCutoff, std::vector< std::unique_ptr<bam1_t, BAMrecordDeleter> > &readMapVector) {
	BAMrecordDeleter newSecondaryDeleter;
	std::unique_ptr<bam1_t, BAMrecordDeleter> secondaryFromNew(bam_init1(), newSecondaryDeleter);

	// reference indexes may not be the same in the original and new headers
	// retrieve the original index from the reference name
	const auto *const newRefNamePtr = sam_hdr_tid2name(newRecordHeader.get(), newRecord->core.tid); 
	if ( (newRefNamePtr == nullptr) || (newRecord->core.n_cigar == 0) ) {
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
	// we made sure there is at least one CIGAR field above
	// must test if the first CIGAR is a soft clip because we can only have one
	const uint32_t firstCIGAR{*bam_get_cigar( newRecord.get() )};
	if ( (bam_cigar_op(firstCIGAR) == BAM_CSOFT_CLIP) && !remapCIGAR.empty() ) {
		const uint32_t summedFirstSoftClips{bam_cigar_oplen(remapCIGAR.back() + bam_cigar_oplen(firstCIGAR))};
		remapCIGAR.back() = bam_cigar_gen(summedFirstSoftClips, BAM_CSOFT_CLIP);
		const auto cigarOpLength{static_cast<float>( bam_cigar_oplen(firstCIGAR) )};
		// match count is 0 for a soft clip
		readLength += cigarOpLength * readConsumption.at( bam_cigar_op( remapCIGAR.back() ) );
	} else {
		remapCIGAR.push_back( *( bam_get_cigar( newRecord.get() ) ) );
		const auto cigarOpLength{static_cast<float>( bam_cigar_oplen( remapCIGAR.back() ) )};
		matchCount += cigarOpLength * sequenceMatch.at( bam_cigar_op( remapCIGAR.back() ) );
		readLength += cigarOpLength * readConsumption.at( bam_cigar_op( remapCIGAR.back() ) );
	}
	for (uint32_t iCIGAR = 1; iCIGAR < newRecord->core.n_cigar; ++iCIGAR) {
		remapCIGAR.push_back( *(bam_get_cigar( newRecord.get() ) + iCIGAR) );
		const auto cigarOpLength{static_cast<float>( bam_cigar_oplen( remapCIGAR.back() ) )};
		matchCount += cigarOpLength * sequenceMatch.at( bam_cigar_op( remapCIGAR.back() ) );
		readLength += cigarOpLength * readConsumption.at( bam_cigar_op( remapCIGAR.back() ) );
	}
	if (matchCount/readLength < remapIdentityCutoff) {
		return;
	}
	if (remapInfo.end <= readMapVector.front()->core.l_qseq) {
		softClip = static_cast<uint32_t>(readMapVector.front()->core.l_qseq - remapInfo.end);
		if (bam_cigar_op( remapCIGAR.back() ) == BAM_CSOFT_CLIP) {
			softClip         += bam_cigar_oplen( remapCIGAR.back() );
			remapCIGAR.back() = bam_cigar_gen(softClip, BAM_CSOFT_CLIP);
		} else {
			remapCIGAR.push_back( bam_cigar_gen(softClip, BAM_CSOFT_CLIP) );
		}
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
		readMapVector.front() = modifyCIGAR( remapInfo, readMapVector.front() );
	}
}

std::unique_ptr<BGZF, BGZFhandleDeleter> isaSpace::openBGZFtoAppend(const std::string &bamFileName) {
	// we will be appending, so must delete this file if it exists
	const auto rmvSuccess = std::remove( bamFileName.c_str() );

	constexpr char openMode{'a'};
	BGZFhandleDeleter handleDeleter;
	std::unique_ptr<BGZF, BGZFhandleDeleter> outputBAMfile(
		bgzf_open(bamFileName.c_str(), &openMode),
		handleDeleter
	);
	if (outputBAMfile == nullptr) {
		throw std::string("ERROR: failed to open the BAM file ")
			+ bamFileName + " for writing in "
			+ std::string( static_cast<const char*>(__PRETTY_FUNCTION__) )
			+ std::string( strerror(errno) ); // NOLINT
	}

	return outputBAMfile;
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
	if (val) { // The last flag had no value
		cliResult[curFlag] = "set";
	}
	return cliResult;
}

void isaSpace::extractCLinfo(
		const std::unordered_map<std::string, std::string> &parsedCLI,
		std::unordered_map<std::string, int> &intVariables,
		std::unordered_map<std::string, float> &floatVariables,
		std::unordered_map<std::string, std::string> &stringVariables) {
	intVariables.clear();
	floatVariables.clear();
	stringVariables.clear();
	const std::array<std::string, 4> requiredStringVariables{"input-bam", "input-gff", "out", "remapped-bam"};
	const std::array<std::string, 2> optionalStringVariables{"out-fastq", "unsorted-output"};
	const std::array<std::string, 2> optionalIntVariables{"threads", "window-size"};
	const std::array<std::string, 1> optionalFloatVariables{"remap-cutoff"};

	const std::unordered_map<std::string, int> defaultIntValues{ {"threads", -1}, {"window-size", 75} };
	const std::unordered_map<std::string, float> defaultFloatValues{ {"remap-cutoff", 0.99F} };
	const std::unordered_map<std::string, std::string> defaultStringValues{ {"out-fastq", "NULL"}, {"unsorted-output", "unset"} };

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
	for (const auto &eachFlag : optionalFloatVariables) {
		try {
			floatVariables[eachFlag] = stof( parsedCLI.at(eachFlag) );
		} catch(const std::exception &problem) {
			floatVariables[eachFlag] = defaultFloatValues.at(eachFlag);
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
