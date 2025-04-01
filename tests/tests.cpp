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

#include <cstddef>
#include <memory>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <utility>
#include <set>
#include <unordered_set>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "bgzf.h"
#include "sam.h"

#include "isoseqAlgn.hpp"
#include "helperFunctions.hpp"

#include "catch2/catch_test_macros.hpp"
//#include "catch2/matchers/catch_matchers.hpp"
//#include "catch2/matchers/catch_matchers_string.hpp"
#include "catch2/catch_approx.hpp"

TEST_CASE("Helper functions work") {
	SECTION("Extract attribute name") {
		constexpr size_t nAttr{4};
		std::array<std::string, nAttr> attributeArr{
			"ID=rna-XM_043206943.1",
			"Parent=gene-LOC6\\",
			"526243",
			"Dbxref=GeneID:6526243,Genbank:XM_043206943.1"
		};
		isaSpace::TokenAttibuteListPair testTokAttr;
		std::copy(
			attributeArr.cbegin(),
			attributeArr.cend(),
			std::back_inserter(testTokAttr.attributeList)
		);
		testTokAttr.tokenName = "Parent=";
		const std::string parentResult{isaSpace::extractAttributeName(testTokAttr)};
		REQUIRE(parentResult == "gene-LOC6\\;526243");
		testTokAttr.tokenName = "Random=";
		const std::string absentAttrResult{isaSpace::extractAttributeName(testTokAttr)};
		REQUIRE( absentAttrResult.empty() );
	}

	SECTION("Test range overlap") {
		constexpr hts_pos_t earlyES{30};
		constexpr hts_pos_t midES{90};
		constexpr hts_pos_t lateES{245};
		constexpr hts_pos_t earlyEE{45};
		constexpr hts_pos_t midEE{110};
		constexpr hts_pos_t lateEE{345};

		constexpr std::pair<hts_pos_t, hts_pos_t> first1(earlyES, lateES);
		constexpr std::pair<hts_pos_t, hts_pos_t> second1(earlyEE, lateEE);
		REQUIRE( isaSpace::rangesOverlap(first1, second1) );
		constexpr std::pair<hts_pos_t, hts_pos_t> first2(earlyES, midES);
		constexpr std::pair<hts_pos_t, hts_pos_t> second2(midEE, lateEE);
		REQUIRE( !isaSpace::rangesOverlap(first2, second2) );
		constexpr std::pair<hts_pos_t, hts_pos_t> second3(earlyEE, midEE);
		REQUIRE( isaSpace::rangesOverlap(first1, second3) );
		constexpr std::pair<hts_pos_t, hts_pos_t> first3(midES, lateES);
		REQUIRE( isaSpace::rangesOverlap(first3, second1) );
		constexpr std::pair<hts_pos_t, hts_pos_t> firstRev(lateES, earlyES);
		constexpr std::pair<hts_pos_t, hts_pos_t> secondRev(lateEE, earlyEE);
		REQUIRE( isaSpace::rangesOverlap(firstRev, secondRev) );
	}

	SECTION("Peak and valley functions") {
		// Peaks function
		// one peak
        const std::vector<float> values1 = {1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 4.0F, 3.0F, 2.0F, 1.0F};
        constexpr float threshold1{3.0F};
		constexpr float correctPeakValue1{5.0F};
        std::vector<std::vector<float>::const_iterator> peaks1{isaSpace::getPeaks(values1, threshold1)};
        REQUIRE(peaks1.size()   == 2);
        REQUIRE(*peaks1.front() == correctPeakValue1);
		REQUIRE( peaks1.back()  == values1.cend() );

		// multiple peaks
        const std::vector<float> values2 = {1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 4.0F, 3.0F, 6.0F, 7.0F, 6.0F, 5.0F, 4.0F, 3.0F, 2.0F, 1.0F};
        constexpr float threshold2{4.0F};
		constexpr std::array<float, 2> correctPeakValues2{5.0F, 7.0F};
        std::vector<std::vector<float>::const_iterator> peaks2 = isaSpace::getPeaks(values2, threshold2);
        REQUIRE(peaks2.size()  == correctPeakValues2.size() + 1);
        REQUIRE( *peaks2.at(0) == correctPeakValues2.at(0) );
        REQUIRE( *peaks2.at(1) == correctPeakValues2.at(1) );
		REQUIRE( peaks2.back() == values2.cend() );

		// no peaks above the threshold
		constexpr float highThreshold{10.0F};
        std::vector<std::vector<float>::const_iterator> peaks3{isaSpace::getPeaks(values1, highThreshold)};
		REQUIRE(peaks3.size() == 1);
		REQUIRE( peaks3.back() == values1.cend() );

		// empty input vector
		const std::vector<float> emptyValues;
        std::vector<std::vector<float>::const_iterator> peaks4{isaSpace::getPeaks(emptyValues, threshold1)};
		REQUIRE( peaks4.empty() );

		// Valleys function
		// one peak
        const std::vector<float> vValues1 = {-1.0F, -2.0F, -3.0F, -4.0F, -5.0F, -4.0F, -3.0F, -2.0F, -1.0F};
        constexpr float vThreshold1{-3.0F};
		constexpr float correctValleyValue1{-5.0F};
        std::vector<std::vector<float>::const_iterator> valleys1{isaSpace::getValleys(vValues1, vThreshold1)};
        REQUIRE(valleys1.size()   == 2);
        REQUIRE(*valleys1.front() == correctValleyValue1);
        REQUIRE( valleys1.back()  == vValues1.cend() );

		// multiple peaks
        const std::vector<float> vValues2 = {-1.0F, -2.0F, -3.0F, -4.0F, -5.0F, -4.0F, -3.0F, -6.0F, -7.0F, -6.0F, -5.0F, -4.0F, -3.0F, -2.0F, -1.0F};
        constexpr float vThreshold2{-4.0F};
		constexpr std::array<float, 2> correctValleyValues2{-5.0F, -7.0F};
        std::vector<std::vector<float>::const_iterator> valleys2 = isaSpace::getValleys(vValues2, vThreshold2);
        REQUIRE(valleys2.size()  == correctValleyValues2.size() + 1);
        REQUIRE( *valleys2.at(0) == correctValleyValues2.at(0) );
        REQUIRE( *valleys2.at(1) == correctValleyValues2.at(1) );
        REQUIRE( valleys2.back()  == vValues2.cend() );

		// no peaks above the threshold
		constexpr float lowThreshold{-10.0F};
        std::vector<std::vector<float>::const_iterator> valleys3{isaSpace::getValleys(vValues1, lowThreshold)};
		REQUIRE(valleys3.size() == 1);
		REQUIRE( valleys3.back() == vValues1.cend() );

		// empty input vector
        std::vector<std::vector<float>::const_iterator> valleys4{isaSpace::getValleys(emptyValues, vThreshold1)};
		REQUIRE( valleys4.empty() );
	}

	/*
	 * Reference match status function tests are in the BAM record section
	 */

	SECTION("Exon coverage and stringify functions") {
		// Exon coverage
		const std::pair<isaSpace::BAMrecord, isaSpace::ExonGroup> emptyPair;
		const auto emptyExonCoverage{isaSpace::getExonCoverageStats(emptyPair)};
		REQUIRE( emptyExonCoverage.geneName == "no_overlap" );
		REQUIRE( emptyExonCoverage.chromosomeName == "not_mapped" );
		REQUIRE( emptyExonCoverage.readName.empty() );
		REQUIRE( emptyExonCoverage.exonCoverageScores.empty() );
		REQUIRE( emptyExonCoverage.bestExonCoverageScores.empty() );
		REQUIRE(emptyExonCoverage.strand == '+');

		// Single stringify function
		const std::string correctStatsLine(
			"m54312U_201215_225530/460252/ccs	NC_052529.2	+	2983523	2988359	2983523	2988359	0	2	2	0	"
			"gene-LOC6532627	4	101	2980697	2988366	{0.000000,0.930000,1.000000,0.500000}	{0.500000,0.930000,1.000000,0.500000}"
		);
		isaSpace::ReadExonCoverage coverageStats;
		coverageStats.readName                 = "m54312U_201215_225530/460252/ccs";
		coverageStats.chromosomeName           = "NC_052529.2";
		coverageStats.strand                   = '+';
		coverageStats.alignmentStart           = 2983523; // NOLINT
		coverageStats.alignmentEnd             = 2988359; // NOLINT
		coverageStats.bestAlignmentStart       = 2983523; // NOLINT
		coverageStats.bestAlignmentEnd         = 2988359; // NOLINT
		coverageStats.firstSoftClipLength      = 0;
		coverageStats.nSecondaryAlignments     = 2;
		coverageStats.nGoodSecondaryAlignments = 2;
		coverageStats.nLocalReversedAlignments = 0;
		coverageStats.geneName                 = "gene-LOC6532627";
		coverageStats.nExons                   = 4;       // NOLINT
		coverageStats.firstExonLength          = 101;     // NOLINT
		coverageStats.firstExonStart           = 2980697; // NOLINT
		coverageStats.lastExonEnd              = 2988366; // NOLINT
		coverageStats.exonCoverageScores       = std::vector<float>{0.0, 0.93, 1.0, 0.5}; // NOLINT
		coverageStats.bestExonCoverageScores   = std::vector<float>{0.5, 0.93, 1.0, 0.5}; // NOLINT
		REQUIRE(isaSpace::stringifyExonCoverage(coverageStats) == correctStatsLine);

		// thread ranges function
		constexpr size_t nThreads{4};
		const isaSpace::bamGFFvector testRECvec(11);
		const auto threadRanges{isaSpace::makeThreadRanges(testRECvec, nThreads)};
		REQUIRE(threadRanges.size() == nThreads);
		REQUIRE(
			std::all_of(
				threadRanges.cbegin(),
				threadRanges.cend(),
				[](const std::pair<isaSpace::bamGFFvector::const_iterator, isaSpace::bamGFFvector::const_iterator> &currPair) {
					return std::distance(currPair.first, currPair.second) >= 0;
				}
			)
		);
		REQUIRE(
			std::all_of(
				threadRanges.cbegin(),
				threadRanges.cend(),
				[&testRECvec, &nThreads](const std::pair<isaSpace::bamGFFvector::const_iterator, isaSpace::bamGFFvector::const_iterator> &currPair) {
					return std::distance(currPair.first, currPair.second) >= testRECvec.size() / nThreads;
				}
			)
		);
		REQUIRE(
			std::distance(threadRanges.front().first, threadRanges.front().second) >
			std::distance(threadRanges.back().first, threadRanges.back().second)
		);
		REQUIRE( threadRanges.front().first == testRECvec.cbegin() );
		REQUIRE( threadRanges.back().second == testRECvec.cend() );

		// Range stringify function
		const std::vector< std::pair<isaSpace::BAMrecord, isaSpace::ExonGroup> > emptyPairVector{emptyPair, emptyPair};
		const auto emptyMultipleRecords{isaSpace::stringifyAlignmentRange( emptyPairVector.cbegin(), emptyPairVector.cend() )};
		REQUIRE(
			std::count_if(
				emptyMultipleRecords.cbegin(),
				emptyMultipleRecords.cend(),
				[](char eachChar) {
					return eachChar == '\n';
				}
			) == emptyPairVector.size()
		);
	}

	SECTION("GFF parsing") {
		const std::string gffLine(
			"NC_052526.2	Gnomon	exon	50812	50970	.	+	.	"
			"ID=exon-XM_043206943.1-1;Parent=rna-XM_043206943.1;Dbxref=GeneID:6526243,Genbank:XM_043206943.1;experiment=COORDINATES: "
			"polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC6526243;product=active breakpoint cluster region-related protein%2C transcript variant X1;transcript_id=XM_043206943.1"
		);
		const std::array<std::string, isaSpace::nGFFfields> goodGFFfields{isaSpace::parseGFFline(gffLine)};
		REQUIRE(goodGFFfields.at(0) == "NC_052526.2");
		REQUIRE(goodGFFfields.at(1) == "Gnomon");
		REQUIRE(goodGFFfields.at(2) == "exon");
		REQUIRE(goodGFFfields.at(3) == "50812");
		REQUIRE(goodGFFfields.at(4) == "50970");
		REQUIRE(goodGFFfields.at(6) == "+");

		const std::string badGFFline(
			"NC_052526.2	exon	50812	50970	.	+	.	"
			"ID=exon-XM_043206943.1-1;Parent=rna-XM_043206943.1;Dbxref=GeneID:6526243,Genbank:XM_043206943.1;experiment=COORDINATES: "
			"polyA evidence [ECO:0006239];gbkey=mRNA;gene=LOC6526243;product=active breakpoint cluster region-related protein%2C transcript variant X1;transcript_id=XM_043206943.1"
		);
		const std::array<std::string, isaSpace::nGFFfields> badGFFfields{isaSpace::parseGFFline(badGFFline)};
		REQUIRE(badGFFfields.at(0) == "FAIL");

		constexpr size_t nReferences{3};
		const std::string goodGFFname("../tests/goodGFF.gff");
		const auto parsedGFF{isaSpace::parseGFF(goodGFFname)};
		REQUIRE(parsedGFF.size() == nReferences);
		REQUIRE(
			std::count_if(
				parsedGFF.cbegin(),
				parsedGFF.cend(),
				[](const std::pair<std::string, std::vector<isaSpace::ExonGroup> > &eachReference){
					return eachReference.second.size() == 1;
				}
			) == 2
		);
		REQUIRE(
			std::count_if(
				parsedGFF.cbegin(),
				parsedGFF.cend(),
				[](const std::pair<std::string, std::vector<isaSpace::ExonGroup> > &eachReference){
					return eachReference.second.size() == 2;
				}
			) == 1
		);
	}

	SECTION("Re-alignment parsing") {
		const std::string goodSubreadName("test_1234_567");
		const isaSpace::ReadPortion goodRP{isaSpace::parseRemappedReadName(goodSubreadName)};

		REQUIRE(goodRP.originalName == "test");
		REQUIRE(goodRP.start        == 1234);
		REQUIRE(goodRP.end          == 567);

		const std::string goodSubreadNameUS("test_test_1234_567");
		const isaSpace::ReadPortion goodRPUS{isaSpace::parseRemappedReadName(goodSubreadNameUS)};

		REQUIRE(goodRPUS.originalName == "test_test");
		REQUIRE(goodRPUS.start        == 1234);
		REQUIRE(goodRPUS.end          == 567);

		const std::string badSubreadName("test");
		const isaSpace::ReadPortion badRP{isaSpace::parseRemappedReadName(badSubreadName)};

		REQUIRE( badRP.originalName.empty() );
		REQUIRE(badRP.start == 0);
		REQUIRE(badRP.end   == 0);

		const std::string badSubreadNameUS("test_test");
		const isaSpace::ReadPortion badRPUS{isaSpace::parseRemappedReadName(badSubreadNameUS)};

		REQUIRE( badRPUS.originalName.empty() );
		REQUIRE(badRPUS.start == 0);
		REQUIRE(badRPUS.end   == 0);

		const std::string emptySubreadName;
		const isaSpace::ReadPortion emptyRP{isaSpace::parseRemappedReadName(emptySubreadName)};

		REQUIRE( emptyRP.originalName.empty() );
		REQUIRE(emptyRP.start == 0);
		REQUIRE(emptyRP.end   == 0);
	}
}

TEST_CASE("Exon range extraction works") {
	constexpr size_t nExons{5};
	constexpr std::array<std::pair<hts_pos_t, hts_pos_t>, nExons> testExonSpans{
		std::pair<hts_pos_t, hts_pos_t>{50812, 50970},
		std::pair<hts_pos_t, hts_pos_t>{52164, 52649},
		std::pair<hts_pos_t, hts_pos_t>{52164, 52649},
		std::pair<hts_pos_t, hts_pos_t>{55522, 56102},
		std::pair<hts_pos_t, hts_pos_t>{56205, 56835}
	};
	const std::string testGeneName("testGene");
	// negative strand
	char strand{'-'};
	std::set< std::pair<hts_pos_t, hts_pos_t> > testSet;
	std::copy( testExonSpans.cbegin(), testExonSpans.cend(), std::inserter( testSet, testSet.end() ) );
	isaSpace::ExonGroup testExonGroupNeg(testGeneName, strand, testSet);
	constexpr size_t correctNexons{4};
	REQUIRE(testExonGroupNeg.nExons()  == correctNexons);
	REQUIRE(testExonGroupNeg.strand()  == strand);
	REQUIRE(testExonGroupNeg.at(1).first  == testExonSpans.at(1).first);
	REQUIRE(testExonGroupNeg.at(1).second == testExonSpans.at(1).second);
	auto fullSpan{testExonGroupNeg.geneSpan()};
	REQUIRE(fullSpan.first  == testExonSpans.front().first);
	REQUIRE(fullSpan.second == testExonSpans.back().second);
	auto firstExon{testExonGroupNeg.firstExonSpan()};
	REQUIRE(firstExon.first  == testExonSpans.back().first);
	REQUIRE(firstExon.second == testExonSpans.back().second);
	auto firstIntron{testExonGroupNeg.getFirstIntronSpan()};
	REQUIRE(firstIntron.first  == testExonSpans.at(3).second);
	REQUIRE(firstIntron.second == testExonSpans.at(4).first);
	// positive strand
	strand = '+';
	isaSpace::ExonGroup testExonGroupPos(testGeneName, strand, testSet);
	REQUIRE(testExonGroupPos.nExons()  == correctNexons);
	REQUIRE(testExonGroupPos.strand()  == strand);
	REQUIRE(testExonGroupNeg.at(1).first  == testExonSpans.at(1).first);
	REQUIRE(testExonGroupNeg.at(1).second == testExonSpans.at(1).second);
	fullSpan = testExonGroupPos.geneSpan();
	REQUIRE(fullSpan.first  == testExonSpans.front().first);
	REQUIRE(fullSpan.second == testExonSpans.back().second);
	firstExon = testExonGroupPos.firstExonSpan();
	REQUIRE(firstExon.first  == testExonSpans.front().first);
	REQUIRE(firstExon.second == testExonSpans.front().second);
	firstIntron = testExonGroupPos.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == testExonSpans.at(0).second);
	REQUIRE(firstIntron.second == testExonSpans.at(1).first);
	// undetermined strand, assumed positive
	strand = '.';
	isaSpace::ExonGroup testExonGroupUnd(testGeneName, strand, testSet);
	REQUIRE(testExonGroupUnd.nExons()  == correctNexons);
	REQUIRE(testExonGroupUnd.strand()  == '+');
	REQUIRE(testExonGroupNeg.at(1).first  == testExonSpans.at(1).first);
	REQUIRE(testExonGroupNeg.at(1).second == testExonSpans.at(1).second);
	fullSpan = testExonGroupUnd.geneSpan();
	REQUIRE(fullSpan.first  == testExonSpans.front().first);
	REQUIRE(fullSpan.second == testExonSpans.back().second);
	firstExon = testExonGroupUnd.firstExonSpan();
	REQUIRE(firstExon.first  == testExonSpans.front().first);
	REQUIRE(firstExon.second == testExonSpans.front().second);
	firstIntron = testExonGroupUnd.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == testExonSpans.at(0).second);
	REQUIRE(firstIntron.second == testExonSpans.at(1).first);

	// Overlapping first exons and other anomalies
	constexpr std::array<std::pair<hts_pos_t, hts_pos_t>, 3> overExonSpans{
		std::pair<hts_pos_t, hts_pos_t>{50812, 50970},
		std::pair<hts_pos_t, hts_pos_t>{50812, 51000},
		std::pair<hts_pos_t, hts_pos_t>{52164, 52649}
	};
	testSet.clear();
	std::copy( overExonSpans.cbegin(), overExonSpans.cend(), std::inserter( testSet, testSet.end() ) );
	strand = '+';
	isaSpace::ExonGroup testExonGroupOvr(testGeneName, strand, testSet);
	firstIntron = testExonGroupOvr.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == overExonSpans.at(1).second);
	REQUIRE(firstIntron.second == overExonSpans.at(2).first);
	strand = '-';
	isaSpace::ExonGroup testExonGroupOvrNeg(testGeneName, strand, testSet);
	firstIntron = testExonGroupOvrNeg.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == overExonSpans.at(1).second);
	REQUIRE(firstIntron.second == overExonSpans.at(2).first);
	strand = '+';
	testSet.clear();
	std::copy( overExonSpans.cbegin(), overExonSpans.cbegin() + 1, std::inserter( testSet, testSet.end() ) );
	isaSpace::ExonGroup testExonGroupOE(testGeneName, strand, testSet);
	firstIntron = testExonGroupOE.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == -1);
	REQUIRE(firstIntron.second == -1);
	strand = '-';
	isaSpace::ExonGroup testExonGroupOEneg(testGeneName, strand, testSet);
	firstIntron = testExonGroupOEneg.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == -1);
	REQUIRE(firstIntron.second == -1);
	strand = '+';
	testSet.clear();
	std::copy( overExonSpans.cbegin(), overExonSpans.cbegin() + 2, std::inserter( testSet, testSet.end() ) );
	isaSpace::ExonGroup testExonGroupTE(testGeneName, strand, testSet);
	firstIntron = testExonGroupTE.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == -1);
	REQUIRE(firstIntron.second == -1);
	strand = '-';
	isaSpace::ExonGroup testExonGroupTEneg(testGeneName, strand, testSet);
	firstIntron = testExonGroupTEneg.getFirstIntronSpan();
	REQUIRE(firstIntron.first  == -1);
	REQUIRE(firstIntron.second == -1);
	
	// exon index functions only tested once: they do not depend on strand
	constexpr hts_pos_t positionBefore{49000};
	constexpr hts_pos_t positionInMiddle{52500};
	constexpr hts_pos_t positionAfter{59000};
	constexpr uint16_t  correctPosBefore{0};
	constexpr uint16_t  correctPosMidF{2};
	constexpr uint16_t  correctPosMidL{1};
	constexpr uint16_t  correctPosAfter{3};
	REQUIRE(testExonGroupPos.firstExonAfter(positionBefore)   == correctPosBefore);
	REQUIRE(testExonGroupPos.firstExonAfter(positionInMiddle) == correctPosMidF);
	REQUIRE(testExonGroupPos.firstExonAfter(positionAfter)    == correctPosAfter);
	REQUIRE(testExonGroupPos.firstOverlappingExon(positionBefore)   == correctPosBefore);
	REQUIRE(testExonGroupPos.firstOverlappingExon(positionInMiddle) == correctPosMidL);
	REQUIRE(testExonGroupPos.firstOverlappingExon(positionAfter)    == correctPosAfter);
	REQUIRE(testExonGroupPos.lastExonBefore(positionBefore)   == correctPosBefore);
	REQUIRE(testExonGroupPos.lastExonBefore(positionInMiddle) == correctPosMidL);
	REQUIRE(testExonGroupPos.lastExonBefore(positionAfter)    == correctPosAfter);
	REQUIRE(testExonGroupPos.lastOverlappingExon(positionBefore)   == correctPosBefore);
	REQUIRE(testExonGroupPos.lastOverlappingExon(positionInMiddle) == correctPosMidF);
	REQUIRE(testExonGroupPos.lastOverlappingExon(positionAfter)    == correctPosAfter);

	// Per-exon alignment quality tests
	std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> commonBAMheader(
		sam_hdr_init(),
		[](sam_hdr_t *samHeader) {
			sam_hdr_destroy(samHeader);
		}
	);
	int headRes = sam_hdr_add_line(commonBAMheader.get(), "HD", "VN", "1.6", "SO", "unknown", NULL); // NOLINT
	headRes     = sam_hdr_add_line(commonBAMheader.get(), "SQ", "SN", "test", "LN", "100000", NULL); // NOLINT
	const std::string tstQname("testQueryName");
	constexpr uint16_t  flag{0};
	constexpr int32_t   tid{0};
	constexpr uint8_t   mapq{60};

	// Simple no mismatch alignments
	constexpr std::array<uint32_t, 7> simpleCIGAR{
		bam_cigar_gen(159,  BAM_CMATCH),
		bam_cigar_gen(1193, BAM_CREF_SKIP),
		bam_cigar_gen(486,  BAM_CMATCH),
		bam_cigar_gen(2872, BAM_CREF_SKIP),
		bam_cigar_gen(581,  BAM_CMATCH),
		bam_cigar_gen(102,  BAM_CREF_SKIP),
		bam_cigar_gen(631,  BAM_CMATCH)
	};
	const std::string tstSeq1(bam_cigar2qlen( simpleCIGAR.size(), simpleCIGAR.data() ), 'A');
	const std::string tstQual1(tstSeq1.size(), '~');
	std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamRecordPtr1(
		bam_init1(),
		[](bam1_t *bamRecord){
			bam_destroy1(bamRecord);
		}
	);
	int bamSetRes = bam_set1(
		bamRecordPtr1.get(),
		tstQname.size(),
		tstQname.c_str(),
		flag,
		tid,
		testExonSpans.front().first - 1,
		mapq,
		simpleCIGAR.size(),
		simpleCIGAR.data(),
		0,
		0,
		static_cast<hts_pos_t>( tstSeq1.size() ),
		tstSeq1.size(),
		tstSeq1.c_str(),
		tstQual1.c_str(),
		0
	);
	isaSpace::BAMrecord bamRecord1( bamRecordPtr1.get(), commonBAMheader.get() );
	std::vector<float> exonCoverage{testExonGroupPos.getExonCoverageQuality(bamRecord1)};
	REQUIRE( exonCoverage.size() == testExonGroupPos.nExons() );
	REQUIRE(
		std::all_of(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val == 1.0F;})
	);
	std::vector<float> bestExonCoverage{testExonGroupPos.getBestExonCoverageQuality(bamRecord1)};
	REQUIRE(
		std::equal(
			exonCoverage.cbegin(),
			exonCoverage.cend(),
			bestExonCoverage.cbegin()
		)
	);

	// Read starts early
	constexpr hts_pos_t earlyStartPos{49999}; // must be base-0
	constexpr float qsCutOff{0.95};
	constexpr std::array<uint32_t, 14> earlyCigarFields{
		bam_cigar_gen(100,  BAM_CSOFT_CLIP),
		bam_cigar_gen(821,  BAM_CMATCH),
		bam_cigar_gen(2,    BAM_CINS),
		bam_cigar_gen(50,   BAM_CMATCH),
		bam_cigar_gen(5,    BAM_CDEL),
		bam_cigar_gen(95,   BAM_CMATCH),
		bam_cigar_gen(1193, BAM_CREF_SKIP),
		bam_cigar_gen(486,  BAM_CMATCH),
		bam_cigar_gen(2872, BAM_CDEL),
		bam_cigar_gen(500,  BAM_CMATCH),
		bam_cigar_gen(1,    BAM_CDIFF),
		bam_cigar_gen(80,   BAM_CMATCH),
		bam_cigar_gen(102,  BAM_CREF_SKIP),
		bam_cigar_gen(604,  BAM_CMATCH)
	};
	const std::string tstSeq2(bam_cigar2qlen( earlyCigarFields.size(), earlyCigarFields.data() ), 'A');
	const std::string tstQual2(tstSeq2.size(), '~');
	std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamRecordPtr2(
		bam_init1(),
		[](bam1_t *bamRecord){
			bam_destroy1(bamRecord);
		}
	);
	bamSetRes = bam_set1(
		bamRecordPtr2.get(),
		tstQname.size(),
		tstQname.c_str(),
		flag,
		tid,
		earlyStartPos,
		mapq,
		earlyCigarFields.size(),
		earlyCigarFields.data(),
		0,
		0,
		static_cast<hts_pos_t>( tstSeq2.size() ),
		tstSeq2.size(),
		tstSeq2.c_str(),
		tstQual2.c_str(),
		0
	);
	isaSpace::BAMrecord bamRecord2( bamRecordPtr2.get(), commonBAMheader.get() );
	exonCoverage.clear();
	exonCoverage = testExonGroupPos.getExonCoverageQuality(bamRecord2);
	REQUIRE( exonCoverage.size() == testExonGroupPos.nExons() );
	REQUIRE(
		std::count_if(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val == 1.0F;}) == 1
	);
	REQUIRE(
		std::all_of(exonCoverage.cbegin(), exonCoverage.cend(), [&qsCutOff](float val) {return val > qsCutOff;})
	);
	bestExonCoverage.clear();
	bestExonCoverage = testExonGroupPos.getBestExonCoverageQuality(bamRecord2);
	REQUIRE(
		std::equal(
			exonCoverage.cbegin(),
			exonCoverage.cend(),
			bestExonCoverage.cbegin()
		)
	);

	// Read starts in the middle of an exon
	constexpr hts_pos_t lateStartPos{52299};
	constexpr std::array<uint32_t, 5> lateCigarFields{
		bam_cigar_gen(350,  BAM_CMATCH),
		bam_cigar_gen(2872, BAM_CREF_SKIP),
		bam_cigar_gen(581,  BAM_CMATCH),
		bam_cigar_gen(102,  BAM_CREF_SKIP),
		bam_cigar_gen(701,  BAM_CMATCH)
	};
	const std::string tstSeq3(bam_cigar2qlen( lateCigarFields.size(), lateCigarFields.data() ), 'A');
	const std::string tstQual3(tstSeq3.size(), '~');
	std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamRecordPtr3(
		bam_init1(),
		[](bam1_t *bamRecord){
			bam_destroy1(bamRecord);
		}
	);
	bamSetRes = bam_set1(
		bamRecordPtr3.get(),
		tstQname.size(),
		tstQname.c_str(),
		flag,
		tid,
		lateStartPos,
		mapq,
		lateCigarFields.size(),
		lateCigarFields.data(),
		0,
		0,
		static_cast<hts_pos_t>( tstSeq3.size() ),
		tstSeq3.size(),
		tstSeq3.c_str(),
		tstQual3.c_str(),
		0
	);
	isaSpace::BAMrecord bamRecord3( bamRecordPtr3.get(), commonBAMheader.get() );
	exonCoverage.clear();
	exonCoverage = testExonGroupPos.getExonCoverageQuality(bamRecord3);
	REQUIRE( exonCoverage.size() == testExonGroupPos.nExons() );
	REQUIRE(
		std::count_if(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val == 1.0F;}) == 2
	);
	REQUIRE(
		std::count_if(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val == 0.0F;}) == 1
	);
	REQUIRE(
		std::count_if(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val > 0.0F;}) == 3
	);
	bestExonCoverage.clear();
	bestExonCoverage = testExonGroupPos.getBestExonCoverageQuality(bamRecord3);
	REQUIRE(
		std::equal(
			exonCoverage.cbegin(),
			exonCoverage.cend(),
			bestExonCoverage.cbegin()
		)
	);

	// Negative strand tests
	constexpr std::array<uint32_t, 14> negativeCigarFields{
		bam_cigar_gen(821,  BAM_CMATCH),
		bam_cigar_gen(2,    BAM_CINS),
		bam_cigar_gen(50,   BAM_CMATCH),
		bam_cigar_gen(5,    BAM_CDEL),
		bam_cigar_gen(95,   BAM_CMATCH),
		bam_cigar_gen(1193, BAM_CREF_SKIP),
		bam_cigar_gen(486,  BAM_CMATCH),
		bam_cigar_gen(2872, BAM_CDEL),
		bam_cigar_gen(500,  BAM_CMATCH),
		bam_cigar_gen(1,    BAM_CDIFF),
		bam_cigar_gen(80,   BAM_CMATCH),
		bam_cigar_gen(102,  BAM_CREF_SKIP),
		bam_cigar_gen(500,  BAM_CMATCH),
		bam_cigar_gen(100,  BAM_CSOFT_CLIP)
	};
	const std::string tstSeq4(bam_cigar2qlen( negativeCigarFields.size(), negativeCigarFields.data() ), 'A');
	const std::string tstQual4(tstSeq4.size(), '~');
	std::unique_ptr<bam1_t, void(*)(bam1_t *bamRecord)> bamRecordPtr4(
		bam_init1(),
		[](bam1_t *bamRecord){
			bam_destroy1(bamRecord);
		}
	);
	bamSetRes = bam_set1(
		bamRecordPtr4.get(),
		tstQname.size(),
		tstQname.c_str(),
		flag,
		tid,
		earlyStartPos,
		mapq,
		negativeCigarFields.size(),
		negativeCigarFields.data(),
		0,
		0,
		static_cast<hts_pos_t>( tstSeq4.size() ),
		tstSeq4.size(),
		tstSeq4.c_str(),
		tstQual4.c_str(),
		0
	);
	isaSpace::BAMrecord bamRecord4( bamRecordPtr4.get(), commonBAMheader.get() );
	exonCoverage.clear();
	exonCoverage = testExonGroupNeg.getExonCoverageQuality(bamRecord4);
	REQUIRE( exonCoverage.size() == testExonGroupPos.nExons() );
	REQUIRE(
		std::count_if(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val == 1.0F;}) == 1
	);
	REQUIRE(
		std::count_if(exonCoverage.cbegin(), exonCoverage.cend(), [&qsCutOff](float val) {return val > qsCutOff;}) == 3
	);
	bestExonCoverage.clear();
	bestExonCoverage = testExonGroupNeg.getBestExonCoverageQuality(bamRecord4);
	REQUIRE(
		std::equal(
			exonCoverage.cbegin(),
			exonCoverage.cend(),
			bestExonCoverage.cbegin()
		)
	);
}

TEST_CASE("Read match window statistics work") {
	constexpr size_t vecLength{20};
	constexpr std::vector< std::pair<float, hts_pos_t> >::difference_type subWindow{5};
	constexpr float hiProb{0.99F};
	constexpr float loProb{0.25F};
	constexpr float bicDiffCutOff{50.0F};
	isaSpace::BinomialWindowParameters windowParameters;
	windowParameters.currentProbability     = loProb;
	windowParameters.alternativeProbability = hiProb;
	windowParameters.windowSize             = vecLength;
	windowParameters.bicDifferenceThreshold = bicDiffCutOff;

	constexpr float correctBIChi{118.758};
	constexpr float correctBIClo{-112.767};
	std::vector< std::pair<float, hts_pos_t> > window(vecLength, {0.0F, 0});
	std::for_each(
		window.begin(),
		window.begin() + subWindow,
		[](std::pair<float, hts_pos_t> &eachPair) {
			eachPair.first = 1.0F;
		}
	);

	isaSpace::ReadMatchWindowBIC hiProbSwitch(window.cbegin(), windowParameters);
	REQUIRE( hiProbSwitch.getBICdifference() == Catch::Approx(correctBIChi) );

	windowParameters.currentProbability     = hiProb;
	windowParameters.alternativeProbability = loProb;
	isaSpace::ReadMatchWindowBIC loProbSwitch(window.cbegin(), windowParameters);
	REQUIRE( loProbSwitch.getBICdifference() == Catch::Approx(correctBIClo) );
}

TEST_CASE("Reading individual BAM records works") {
	constexpr char openMode{'r'};
	isaSpace::BinomialWindowParameters binomialParams;
	binomialParams.windowSize             = 50;     // NOLINT
	binomialParams.currentProbability     = 0.25F;  // NOLINT
	binomialParams.alternativeProbability = 0.99F;  // NOLINT
	binomialParams.bicDifferenceThreshold = 200.0F; // NOLINT

	SECTION("Reading a single straight record") {
		const std::string oneRecordBAMname("../tests/oneRecord.bam");
		const std::string correctReadName("m54312U_201215_225530/657093/ccs");
		const std::string correctReferenceName("NC_052530.2");
		const std::string correctCIGAR("22M1D667M1552N264M485N193M361N75M81N196M53N1706M62N126M236N170M");
		constexpr std::array<uint32_t, 17> correctCIGV{
			352,  18,   10672, 24835, 4224,
			7763, 3088, 5779,  1200,  1299,
			3136, 851,  27296, 995,   2016,
			3779, 2720
		};
		constexpr hts_pos_t correctMapPosition{2468434};
		constexpr hts_pos_t correctMapEndPosition{2474684};
		constexpr hts_pos_t correctReadLength{3419};
		constexpr size_t    correctRefMatchLength{6250};
		constexpr std::vector< std::pair<float, hts_pos_t> >::difference_type correctNjumps{8};

		// read one record from the BAM file
		std::unique_ptr<BGZF, void(*)(BGZF *)> orBAMfile(
			bgzf_open(oneRecordBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);

		// must read the header first to get to the alignments
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> orBAMheader(
			bam_hdr_read( orBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		auto nBytes = bam_read1( orBAMfile.get(), bamRecordPtr.get() );
		isaSpace::BAMrecord bamRecord( bamRecordPtr.get(), orBAMheader.get() );
		REQUIRE(nBytes > 0);
		REQUIRE( !bamRecord.isRevComp() );
		REQUIRE(bamRecord.getReadName()    == correctReadName);
		REQUIRE(bamRecord.getmRNAstart()   == correctMapPosition);
		REQUIRE(bamRecord.getMapStart()    == correctMapPosition);
		REQUIRE(bamRecord.getMapEnd()      == correctMapEndPosition);
		REQUIRE(bamRecord.getReadLength()  == correctReadLength);
		REQUIRE(bamRecord.getCIGARstring() == correctCIGAR);

		const auto bamRecordBadAln{bamRecord.getPoorlyMappedRegions(binomialParams)};
		REQUIRE( bamRecordBadAln.empty() );

		const std::vector<uint32_t> cigarVec{bamRecord.getCIGARvector()};
		REQUIRE(
			std::equal( cigarVec.cbegin(), cigarVec.cend(), correctCIGV.cbegin() )
		);
		REQUIRE(bamRecord.getReferenceName() == correctReferenceName);

		std::vector<float> referenceMatches{isaSpace::getReferenceMatchStatus( bamRecord.getCIGARvector() )};
		REQUIRE(referenceMatches.size() == correctRefMatchLength);
		float matchSum = std::accumulate(referenceMatches.cbegin(), referenceMatches.cend(), 0.0F);
		REQUIRE( matchSum == static_cast<float>(correctReadLength) );

		std::vector< std::pair<float, hts_pos_t> > readMatches{bamRecord.getReadCentricMatchStatus()};
		REQUIRE(readMatches.front().second == correctMapPosition);
		REQUIRE(
			std::all_of(
				readMatches.cbegin(),
				readMatches.cend(),
				[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.first == 1.0;}
			)
		);
		std::vector< std::pair<float, hts_pos_t> > posDiffs;
		std::adjacent_difference(
			readMatches.cbegin(),
			readMatches.cend(),
			std::back_inserter(posDiffs),
			[](std::pair<float, hts_pos_t> prevVal, const std::pair<float, hts_pos_t> &eachPair){
				prevVal.second -= eachPair.second;
				return prevVal;
			}
		);
		REQUIRE(
			std::count_if(
				std::next( posDiffs.cbegin() ),
				posDiffs.cend(),
				[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.second > 1;}
			) == correctNjumps
		);

		constexpr auto startTooLate{correctReadLength + 2};
		constexpr auto endTooLate{correctReadLength + 10};
		isaSpace::MappedReadInterval chunkInterval;
		chunkInterval.readStart = startTooLate;
		chunkInterval.readEnd   = endTooLate;
		REQUIRE(bamRecord.getSequenceAndQuality(chunkInterval).size() == 4);
		constexpr hts_pos_t goodStart{2};
		constexpr hts_pos_t goodEnd{10};
		constexpr hts_pos_t goodLength{goodEnd - goodStart};
		chunkInterval.readStart = goodStart;
		chunkInterval.readEnd   = goodEnd;
		const std::string goodSAQ{bamRecord.getSequenceAndQuality(chunkInterval)};
		REQUIRE(goodSAQ.size() == (2 * goodLength) + 4);
		REQUIRE(
			std::count_if(
				goodSAQ.cbegin(),
				goodSAQ.cend(),
				[](char eachChar) {
					return eachChar == '\n';
				}
			) == 3
		);
		REQUIRE(
			std::all_of(
				goodSAQ.cbegin() + goodLength + 4,
				goodSAQ.cend() - 1,
				[](char eachChar) {
					return eachChar == '~';
				}
			) 
		);
		chunkInterval.readEnd = endTooLate;
		const std::string tooLateSAQ{bamRecord.getSequenceAndQuality(chunkInterval)};
		REQUIRE(tooLateSAQ.size() == ( 2 * (correctReadLength - goodStart) ) + 4);
		// read too short for the window
		const std::string shortBAMname("../tests/shortRead.bam");
		std::unique_ptr<BGZF, void(*)(BGZF *)> shortBAMfile(
			bgzf_open(shortBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);

		// must read the header first to get to the alignments
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> shortBAMheader(
			bam_hdr_read( shortBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> shortBamRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		nBytes = bam_read1( shortBAMfile.get(), shortBamRecordPtr.get() );
		isaSpace::BAMrecord shortBamRecord( shortBamRecordPtr.get(), shortBAMheader.get() );
		binomialParams.windowSize = 80;//NOLINT
		const auto regions{shortBamRecord.getPoorlyMappedRegions(binomialParams)};
		REQUIRE( regions.empty() );
	}

	SECTION("Reading a single reverse-complemented record") {
		const std::string oneRecordRevBAMname("../tests/oneRecordRev.bam");
		const std::string correctRevReadName("m54312U_201215_225530/38470652/ccs");
		const std::string correctRevReferenceName("NC_052528.2");
		const std::string correctCIGARrev("454M68N701M61N133M61N1652M78N106M56N280M");
		constexpr std::array<uint32_t, 11> correctCIGVrev{
			4480, 899,  1696, 1251,  26432,
			979,  2128, 979,  11216, 1091, 7264
		};
		constexpr hts_pos_t correctRevMapPosition{22550967};
		constexpr hts_pos_t correctRevmRNAposition{22554617};
		constexpr hts_pos_t correctRevReadLength{3326};
		constexpr hts_pos_t correctRevMapEndPosition = correctRevmRNAposition;
		constexpr size_t    correctRevRefMatchLength{3650};
		constexpr std::vector< std::pair<float, hts_pos_t> >::difference_type correctNjumpsRev{5};

		std::unique_ptr<BGZF, void(*)(BGZF *)> orRevBAMfile(
			bgzf_open(oneRecordRevBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);

		// must read the header first to get to the alignments
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> orRevBAMheader(
			bam_hdr_read( orRevBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamRevRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		const auto nBytes = bam_read1( orRevBAMfile.get(), bamRevRecordPtr.get() );
		isaSpace::BAMrecord bamRecordRev( bamRevRecordPtr.get(), orRevBAMheader.get() );
		REQUIRE(nBytes > 0);
		REQUIRE( bamRecordRev.isRevComp() );
		REQUIRE(bamRecordRev.getReadName()      == correctRevReadName);
		REQUIRE(bamRecordRev.getmRNAstart()     == correctRevmRNAposition);
		REQUIRE(bamRecordRev.getMapStart()      == correctRevMapPosition);
		REQUIRE(bamRecordRev.getMapEnd()        == correctRevMapEndPosition);
		REQUIRE(bamRecordRev.getReadLength()    == correctRevReadLength);
		REQUIRE(bamRecordRev.getCIGARstring()   == correctCIGARrev);
		REQUIRE(bamRecordRev.getReferenceName() == correctRevReferenceName);

		const auto bamRecordBadAlnRev{bamRecordRev.getPoorlyMappedRegions(binomialParams)};
		REQUIRE( bamRecordBadAlnRev.empty() );

		const std::vector<uint32_t> cigarVecRev{bamRecordRev.getCIGARvector()};
		REQUIRE(
			std::equal( cigarVecRev.cbegin(), cigarVecRev.cend(), correctCIGVrev.cbegin() )
		);

		std::vector<float> referenceMatches = isaSpace::getReferenceMatchStatus( bamRecordRev.getCIGARvector() );
		REQUIRE(referenceMatches.size() == correctRevRefMatchLength);
		float matchSum = std::accumulate(referenceMatches.cbegin(), referenceMatches.cend(), 0.0F);
		REQUIRE( matchSum == static_cast<float>(correctRevReadLength) );

		std::vector< std::pair<float, hts_pos_t> > readMatches = bamRecordRev.getReadCentricMatchStatus();
		REQUIRE(readMatches.front().second == correctRevMapPosition);
		REQUIRE(
			std::all_of(
				readMatches.cbegin(),
				readMatches.cend(),
				[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.first == 1.0;}
			)
		);
		std::vector< std::pair<float, hts_pos_t> > posDiffs;
		std::adjacent_difference(
			readMatches.cbegin(),
			readMatches.cend(),
			std::back_inserter(posDiffs),
			[](std::pair<float, hts_pos_t> prevVal, const std::pair<float, hts_pos_t> &eachPair){
				prevVal.second -= eachPair.second;
				return prevVal;
			}
		);
		REQUIRE(
			std::count_if(
				std::next( posDiffs.cbegin() ),
				posDiffs.cend(),
				[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.second > 1;}
			) == correctNjumpsRev
		);
		constexpr auto startTooLate{correctRevReadLength + 2};
		constexpr auto endTooLate{correctRevReadLength + 10};
		isaSpace::MappedReadInterval chunkInterval;
		chunkInterval.readStart = startTooLate;
		chunkInterval.readEnd   = endTooLate;
		REQUIRE(bamRecordRev.getSequenceAndQuality(chunkInterval).size() == 4);
		constexpr hts_pos_t goodStart{2};
		constexpr hts_pos_t goodEnd{10};
		constexpr hts_pos_t goodLength{goodEnd - goodStart};
		chunkInterval.readStart = goodStart;
		chunkInterval.readEnd   = goodEnd;
		const std::string goodSAQ{bamRecordRev.getSequenceAndQuality(chunkInterval)};
		REQUIRE(goodSAQ.size() == (2 * goodLength) + 4);
		REQUIRE(
			std::count_if(
				goodSAQ.cbegin(),
				goodSAQ.cend(),
				[](char eachChar) {
					return eachChar == '\n';
				}
			) == 3
		);
		REQUIRE(
			std::all_of(
				goodSAQ.cbegin() + goodLength + 4,
				goodSAQ.cend() - 1,
				[](char eachChar) {
					return eachChar == '~';
				}
			) 
		);
		chunkInterval.readEnd = endTooLate;
		const std::string tooLateSAQ{bamRecordRev.getSequenceAndQuality(chunkInterval)};
		REQUIRE(tooLateSAQ.size() == ( 2 * (correctRevReadLength - goodStart) ) + 4);
	}

	SECTION("Reading a single soft clipped record") {
		const std::string oneRecordSoftBAMname("../tests/oneRecordSC.bam");
		const std::string correctSoftReadName("m54312U_201215_225530/1115205/ccs");
		const std::string correctSoftReferenceName("NC_052528.2");
		const std::string correctCIGARsoft("13S119M2003N264M1538N210M221N887M168N79M765N1176M67N906M701N64M440N164M");
		constexpr std::array<uint32_t, 18> correctCIGVsoft{
			2624,  7043,  1024, 11219, 14496, 1075,
			18816, 12243, 1264, 2691,  14192, 3539,
			3360,  24611, 4224, 32051, 1904,  212
		};
		constexpr hts_pos_t correctSoftMapPosition{6894852};
		constexpr hts_pos_t correctSoftmRNAposition{6904624};
		constexpr hts_pos_t correctSoftReadLength{3882};
		constexpr hts_pos_t correctSoftMapEndPosition{6904624};
		constexpr hts_pos_t correctSoftMisMatchCount{13};
		constexpr size_t    correctSoftRefMatchLength{9772};
		constexpr float     correctSoftMatchCount{3869};
		constexpr std::vector< std::pair<float, hts_pos_t> >::difference_type correctNjumpsSoft{8};

		std::unique_ptr<BGZF, void(*)(BGZF *)> orSoftBAMfile(
			bgzf_open(oneRecordSoftBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);

		// must read the header first to get to the alignments
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> orSoftBAMheader(
			bam_hdr_read( orSoftBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamSoftRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		const auto nBytes = bam_read1( orSoftBAMfile.get(), bamSoftRecordPtr.get() );
		isaSpace::BAMrecord bamRecordSoft( bamSoftRecordPtr.get(), orSoftBAMheader.get() );
		REQUIRE(nBytes > 0);
		REQUIRE( bamRecordSoft.isRevComp() );
		REQUIRE(bamRecordSoft.getReadName()      == correctSoftReadName);
		REQUIRE(bamRecordSoft.getMapStart()      == correctSoftMapPosition);
		REQUIRE(bamRecordSoft.getmRNAstart()     == correctSoftmRNAposition);
		REQUIRE(bamRecordSoft.getMapEnd()        == correctSoftMapEndPosition);
		REQUIRE(bamRecordSoft.getReadLength()    == correctSoftReadLength);
		REQUIRE(bamRecordSoft.getCIGARstring()   == correctCIGARsoft);
		REQUIRE(bamRecordSoft.getReferenceName() == correctSoftReferenceName);

		const auto bamRecordBadAlnSoft{bamRecordSoft.getPoorlyMappedRegions(binomialParams)};
		REQUIRE( bamRecordBadAlnSoft.empty() );

		const std::vector<uint32_t> cigarVecSoft{bamRecordSoft.getCIGARvector()};
		REQUIRE(
			std::equal( cigarVecSoft.cbegin(), cigarVecSoft.cend(), correctCIGVsoft.cbegin() )
		);

		std::vector<float> referenceMatches = isaSpace::getReferenceMatchStatus( bamRecordSoft.getCIGARvector() );
		REQUIRE(referenceMatches.size() == correctSoftRefMatchLength);
		float matchSum = std::accumulate(referenceMatches.cbegin(), referenceMatches.cend(), 0.0F);
		REQUIRE(matchSum == correctSoftMatchCount);

		std::vector< std::pair<float, hts_pos_t> > readMatches = bamRecordSoft.getReadCentricMatchStatus();
		REQUIRE(readMatches.front().second == correctSoftMapPosition);
		REQUIRE(
			std::count_if(
				readMatches.cbegin(),
				readMatches.cend(),
				[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.first == 0.0;}
			) == correctSoftMisMatchCount
		);
		std::vector< std::pair<float, hts_pos_t> > posDiffs;
		std::adjacent_difference(
			readMatches.cbegin(),
			readMatches.cend(),
			std::back_inserter(posDiffs),
			[](std::pair<float, hts_pos_t> prevVal, const std::pair<float, hts_pos_t> &eachPair){
				prevVal.second -= eachPair.second;
				return prevVal;
			}
		);
		REQUIRE(
			std::count_if(
				std::next( posDiffs.cbegin() ),
				posDiffs.cend(),
				[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.second > 1;}
			) == correctNjumpsSoft
		);
	}

	SECTION("Reading a single record with a poorly mapped region in front") {
		const std::string longSoftClipBAMname("../tests/longSC.bam");
		std::unique_ptr<BGZF, void(*)(BGZF *)> longSoftClipBAMfile(
			bgzf_open(longSoftClipBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> longSoftClipBAMheader(
			bam_hdr_read( longSoftClipBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> longSoftClipRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		const auto nBytes = bam_read1( longSoftClipBAMfile.get(), longSoftClipRecordPtr.get() );
		isaSpace::BAMrecord longSoftClipBAM( longSoftClipRecordPtr.get(), longSoftClipBAMheader.get() );

		constexpr hts_pos_t correctSoftClipSize{258};
		const auto poorlyMappedRegion{longSoftClipBAM.getPoorlyMappedRegions(binomialParams)};
		REQUIRE(poorlyMappedRegion.size() == 1);
		REQUIRE(poorlyMappedRegion.front().referenceStart == poorlyMappedRegion.front().referenceEnd);
		REQUIRE(poorlyMappedRegion.front().readStart == 0);
		REQUIRE(poorlyMappedRegion.front().readEnd == correctSoftClipSize);
	}

	SECTION("Reading a single record with a poorly mapped region at the end") {
		const std::string longSoftClipEndBAMname("../tests/longSCrc.bam");
		std::unique_ptr<BGZF, void(*)(BGZF *)> longSoftClipEndBAMfile(
			bgzf_open(longSoftClipEndBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> longSoftClipEndBAMheader(
			bam_hdr_read( longSoftClipEndBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> longSoftClipEndRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		const auto nBytes = bam_read1( longSoftClipEndBAMfile.get(), longSoftClipEndRecordPtr.get() );
		isaSpace::BAMrecord longSoftClipEndBAM( longSoftClipEndRecordPtr.get(), longSoftClipEndBAMheader.get() );

		constexpr hts_pos_t correctSoftClipEndSize{215};
		const auto poorlyMappedEndRegion{longSoftClipEndBAM.getPoorlyMappedRegions(binomialParams)};
		REQUIRE(poorlyMappedEndRegion.size() == 1);
		REQUIRE(poorlyMappedEndRegion.front().referenceStart == poorlyMappedEndRegion.front().referenceEnd);
		REQUIRE( (poorlyMappedEndRegion.front().readEnd - poorlyMappedEndRegion.front().readStart) == correctSoftClipEndSize );
		REQUIRE( poorlyMappedEndRegion.front().readEnd == longSoftClipEndBAM.getReadLength() );
	}

	SECTION("Reading a single record with a poorly mapped region in the middle") {
		const std::string longSoftClipMidBAMname("../tests/longSCmid.bam");
		std::unique_ptr<BGZF, void(*)(BGZF *)> longSoftClipMidBAMfile(
			bgzf_open(longSoftClipMidBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> longSoftClipMidBAMheader(
			bam_hdr_read( longSoftClipMidBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> longSoftClipMidRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		const auto nBytes = bam_read1( longSoftClipMidBAMfile.get(), longSoftClipMidRecordPtr.get() );
		isaSpace::BAMrecord longSoftClipMidBAM( longSoftClipMidRecordPtr.get(), longSoftClipMidBAMheader.get() );

		constexpr hts_pos_t correctSoftClipMidSize{258};
		constexpr hts_pos_t correctSoftClipMidStart{741};
		constexpr hts_pos_t correctSoftClipMidEnd{999};
		const auto poorlyMappedMidRegion{longSoftClipMidBAM.getPoorlyMappedRegions(binomialParams)};
		REQUIRE(poorlyMappedMidRegion.size() == 1);
		REQUIRE(poorlyMappedMidRegion.front().referenceStart == poorlyMappedMidRegion.front().referenceEnd);
		REQUIRE( (poorlyMappedMidRegion.front().readEnd - poorlyMappedMidRegion.front().readStart) == correctSoftClipMidSize );
		REQUIRE(poorlyMappedMidRegion.front().readStart == correctSoftClipMidStart);
		REQUIRE(poorlyMappedMidRegion.front().readEnd   == correctSoftClipMidEnd);
	}

	SECTION("Reading a record with secondary alignments") {
		constexpr uint16_t correctNsecondary{5};
		constexpr uint16_t correctNlocalSecondary{4};
		constexpr uint16_t correctNlocalSecondaryRev{1};
		const std::string secondaryBAMname("../tests/oneRecordWithSecondary.bam");
		std::unique_ptr<BGZF, void(*)(BGZF *)> secondaryBAMfile(
			bgzf_open(secondaryBAMname.c_str(), &openMode),
			[](BGZF *bamFile) {
				bgzf_close(bamFile);
			}
		);
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> secondaryBAMheader(
			bam_hdr_read( secondaryBAMfile.get() ),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> primaryRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		auto nBytes = bam_read1( secondaryBAMfile.get(), primaryRecordPtr.get() );
		isaSpace::BAMrecord primaryWithSecondaryBAM( primaryRecordPtr.get(), secondaryBAMheader.get() );
		while (true) {
			std::unique_ptr<bam1_t, void(*)(bam1_t *)> secondaryRecordPtr(
				bam_init1(),
				[](bam1_t *bamRecord){
					bam_destroy1(bamRecord);
				}
			);
			nBytes = bam_read1( secondaryBAMfile.get(), secondaryRecordPtr.get() );
			if (nBytes < 0) {
				break;
			}
			primaryWithSecondaryBAM.addSecondaryAlignment( secondaryRecordPtr.get(), secondaryBAMheader.get() );
		}
		REQUIRE( primaryWithSecondaryBAM.hasSecondaryAlignments() );
		REQUIRE(primaryWithSecondaryBAM.secondaryAlignmentCount() == correctNsecondary);
		REQUIRE(primaryWithSecondaryBAM.localSecondaryAlignmentCount() == correctNlocalSecondary);
		REQUIRE(primaryWithSecondaryBAM.localReversedSecondaryAlignmentCount() == correctNlocalSecondaryRev);
		const isaSpace::MappedReadMatchStatus primaryBAQ = primaryWithSecondaryBAM.getBestReferenceMatchStatus();
		const std::vector<float> primaryAQ               = isaSpace::getReferenceMatchStatus(primaryWithSecondaryBAM.getCIGARvector());
		REQUIRE(primaryWithSecondaryBAM.getMapStart() == primaryBAQ.mapStart);
		REQUIRE( primaryAQ.size() < primaryBAQ.matchStatus.size() );
		const auto aqMatchCount = std::count_if(
			primaryAQ.cbegin(),
			primaryAQ.cend(),
			[](float matchStatus) {
				return matchStatus == 1.0F;
			}
		);
		const auto baqMatchCount = std::count_if(
			primaryBAQ.matchStatus.cbegin(),
			primaryBAQ.matchStatus.cend(),
			[](float matchStatus) {
				return matchStatus == 1.0F;
			}
		);
		REQUIRE(aqMatchCount < baqMatchCount);
	}
}

/*
 Since HTSLIB is not in my control, I cannot test everything
 For example, testing the throw on wrong file name is iffy
 because sometimes HTSLIB seems to create the file if there is none, but not consistently
 There also are occasional unexplained crashes, perhaps HTSLIB is not thread safe?
 I have commented out the bad BAM tests for now
TEST_CASE("Catching bad GFF and BAM files works") {
	const std::string goodGFFname("../tests/goodGFF.gff");
	const std::string nomrnaGFFname("../tests/nomRNA.gff");
	const std::string goodBAMname("../tests/testAlgn.bam");
	const std::string randomNoiseBAMname("../tests/randomNoise.bam"); // completely random noise, no magic bytes, but valid EOF sequence
	const std::string headerlessBAMname("../tests/headerless.bam");   // random noise with magic bytes and EOF but no header
	isaSpace::BamAndGffFiles gffPair;
	gffPair.gffFileName = goodGFFname;
	gffPair.bamFileName = randomNoiseBAMname;
	REQUIRE_THROWS_WITH(
		isaSpace::BAMtoGenome(gffPair),
		Catch::Matchers::StartsWith("ERROR: failed to read the header from the BAM file")
	);
	gffPair.bamFileName = headerlessBAMname;
	REQUIRE_THROWS_WITH(
		isaSpace::BAMtoGenome(gffPair),
		Catch::Matchers::StartsWith("ERROR: failed to read the header from the BAM file")
	);
	gffPair.gffFileName = nomrnaGFFname;
	gffPair.bamFileName = goodBAMname;
	REQUIRE_THROWS_WITH(
		isaSpace::BAMtoGenome(gffPair),
		Catch::Matchers::StartsWith("ERROR: no mRNAs with exons found in the")
	);
}
*/

TEST_CASE("GFF and BAM parsing works") {
	const std::string gffName("../tests/posNegYak.gff");
	const std::string testAlignmentBAMname("../tests/testAlignment.bam");
	const std::string outFileName("../tests/testResults.tsv");
	constexpr size_t nThreads{2};
	isaSpace::BamAndGffFiles gffPair;
	gffPair.gffFileName = gffName;
	gffPair.bamFileName = testAlignmentBAMname;
	isaSpace::BAMtoGenome testBTG(gffPair);
	testBTG.saveReadCoverageStats(outFileName, nThreads);

	std::fstream saveResultFile(outFileName, std::ios::in);
	std::string line;
	std::vector<char> strands;
	std::vector<int32_t> fscLengths;                // first soft clip length
	std::vector<int32_t> nSecondary;
	std::vector<std::string> geneNames;
	std::vector<int32_t> nExons;
	std::vector<int32_t> firstEL;
	std::vector<int32_t> fesValues;                 // first exon starts
	std::vector<int32_t> felValues;                 // first exon lengths
	std::getline(saveResultFile, line);             // get rid of the header
	while ( std::getline(saveResultFile, line) ) {
		std::stringstream lineStream;
		lineStream.str(line);
		std::string field;
		lineStream >> field;
		lineStream >> field;
		lineStream >> field;
		strands.push_back( field.at(0) );
		lineStream >> field;
		lineStream >> field;
		lineStream >> field;
		lineStream >> field;
		lineStream >> field;
		fscLengths.push_back( stoi(field) );
		lineStream >> field;
		nSecondary.push_back( stoi(field) );
		lineStream >> field;
		lineStream >> field;
		lineStream >> field;
		geneNames.emplace_back(field);
		lineStream >> field;
		nExons.push_back( stoi(field) );
		lineStream >> field;
		felValues.push_back( stoi(field) );
		lineStream >> field;
		fesValues.push_back( stoi(field) );
	}
	saveResultFile.close();
	constexpr size_t correctResSize{12};
	constexpr int32_t correctNpos{6};
	constexpr int32_t correctNneg{6};
	constexpr int32_t correctNgenes{8};
	constexpr int32_t correctNfailed{4};
	constexpr int32_t correctNsoftClip{4};
	constexpr int32_t correctNsecondary{3};

	REQUIRE(strands.size() == correctResSize);
	REQUIRE(std::count(strands.cbegin(), strands.cend(), '+') == correctNpos);
	REQUIRE(std::count(strands.cbegin(), strands.cend(), '-') == correctNneg);

	REQUIRE(
		std::count_if(
			geneNames.cbegin(),
			geneNames.cend(),
			[](const std::string &name) {
				return name == "no_overlap";
			}
		) == correctNfailed
	);
	REQUIRE(
		std::count_if(
			geneNames.cbegin(),
			geneNames.cend(),
			[](const std::string &name) {
				return name.substr(0, 4) == "gene";
			}
		) == correctNgenes
	);

	REQUIRE(std::count(nExons.cbegin(), nExons.cend(), 0) == correctNfailed);
	REQUIRE(std::count(fesValues.cbegin(), fesValues.cend(), -1) == correctNfailed);

	REQUIRE(
		std::count_if(
			fscLengths.cbegin(),
			fscLengths.cend(),
			[](int32_t val){return val > 0;}
		) == correctNsoftClip
	);
	REQUIRE(
		std::count_if(
			felValues.cbegin(),
			felValues.cend(),
			[](int32_t val){return val > 0;}
		) == correctNgenes
	);
	REQUIRE(
		std::count_if(
			nSecondary.cbegin(),
			nSecondary.cend(),
			[](int32_t val){return val > 0;}
		) == correctNsecondary
	);
	// poorly aligned portions test
	const std::string badRegionFileName("../tests/testBadRegions.tsv");
	constexpr float hiProb{0.99F};
	constexpr float loProb{0.25F};
	constexpr float bicDiffCutOff{100.0F};
	constexpr hts_pos_t windowSize{80};
	isaSpace::BinomialWindowParameters bwParams;
	bwParams.currentProbability     = loProb;
	bwParams.alternativeProbability = hiProb;
	bwParams.windowSize             = windowSize;
	bwParams.bicDifferenceThreshold = bicDiffCutOff;
	testBTG.saveUnmappedRegions(badRegionFileName, bwParams, nThreads);

	std::fstream badRgionResultFile(badRegionFileName, std::ios::in);
	constexpr size_t nIntFields{4};
	std::unordered_set<std::string> uniqueReadNames;
	std::vector< std::array<int32_t, nIntFields> > intFields;
	std::getline(badRgionResultFile, line);             // get rid of the header
	while ( std::getline(badRgionResultFile, line) ) {
		std::stringstream lineStream;
		std::array<int32_t, nIntFields> currentIntFields{0, 0, 0};
		lineStream.str(line);
		std::string field;
		lineStream >> field;
		uniqueReadNames.insert(field);
		lineStream >> field;
		currentIntFields.at(0) = stoi(field);
		lineStream >> field;
		currentIntFields.at(1) = stoi(field);
		lineStream >> field;
		currentIntFields.at(2) = stoi(field);
		lineStream >> field;
		currentIntFields.at(3) = stoi(field);
		intFields.emplace_back(currentIntFields);
	}
	badRgionResultFile.close();

	constexpr size_t correctBadRegionResSize{3};
	constexpr size_t correctNuniqueReads{2};
	REQUIRE(uniqueReadNames.size() == correctNuniqueReads);
    REQUIRE(intFields.size() == correctBadRegionResSize);
	REQUIRE(
		std::all_of(
			intFields.cbegin(),
			intFields.cend(),
			[](const std::array<int, nIntFields> &eachLine) {
				return eachLine.at(0) >= eachLine.at(2);
			}
		)
	);
	REQUIRE(
		std::all_of(
			intFields.cbegin(),
			intFields.cend(),
			[&bwParams](const std::array<int, nIntFields> &eachLine) {
				return ( eachLine.at(2) - eachLine.at(1) ) >= bwParams.windowSize;
			}
		)
	);
	REQUIRE(
		std::all_of(
			intFields.cbegin(),
			intFields.cend(),
			[&bwParams](const std::array<int, nIntFields> &eachLine) {
				return eachLine.at(3)== bwParams.windowSize;
			}
		)
	);

	isaSpace::StatsAndFastqFiles badRegionFiles;
	badRegionFiles.statsFileName = "../tests/testBadRegionsWFQ.tsv";
	badRegionFiles.fastqFileName = "../tests/testBadRegionsWFQ.fq";
	testBTG.saveUnmappedRegions(badRegionFiles, bwParams, nThreads);
	std::fstream badRgionFASTQfile(badRegionFiles.fastqFileName, std::ios::in);
	uint16_t nPluses{0};
	uint16_t nHeaders{0};
	while ( std::getline(badRgionFASTQfile, line) ) {
		if (line == "+") {
			nPluses++;
		}
		if (line.front() == '@') {
			nHeaders++;
		}
	}
	badRgionResultFile.close();
	REQUIRE(nPluses  == correctBadRegionResSize);
    REQUIRE(nHeaders == correctBadRegionResSize);
}

TEST_CASE("Test adding and saving unmapped regions from a BAM file") {
	const std::string testInputAlgnBAMname("../tests/testRealign.bam");
	constexpr size_t correctNprimary{12};
	isaSpace::BAMfile testInputAlgnBAM(testInputAlgnBAMname);
	REQUIRE(testInputAlgnBAM.getPrimaryAlignmentCount() == correctNprimary);
}
