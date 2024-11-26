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

#include <cstddef>
#include <memory>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <utility>
#include <set>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/sam.h"

#include "isoseqAlgn.hpp"
#include "helperFunctions.hpp"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"
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
		std::unique_ptr<bam1_t, void(*)(bam1_t *)> bamRecordPtr(
			bam_init1(),
			[](bam1_t *bamRecord){
				bam_destroy1(bamRecord);
			}
		);
		std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> tstBAMheader(
			sam_hdr_init(),
			[](sam_hdr_t *samHeader) {
				sam_hdr_destroy(samHeader);
			}
		);
		const std::string tstSeq("GCCTACTGCAGTCCACAAGGAGGCCACATCCACAGCCACGACAACGGCGACATACGCCAATGGCAATCCCAATTCTAACGCAAATCCTAGCCAGAGTCAG");
		const std::string tstQual(tstSeq.size(), '~');
		const std::string tstQname("testQueryName");
		constexpr uint16_t  flag{0};
		constexpr int32_t   tid{0};
		constexpr hts_pos_t position{49};
		constexpr uint8_t   mapq{60};
		const std::array<uint32_t, 1> tstCigar{bam_cigar_gen(tstSeq.size(), BAM_CMATCH)};

		const int bamSetRes = bam_set1(
			bamRecordPtr.get(),
			tstQname.size(),
			tstQname.c_str(),
			flag,
			tid,
			position,
			mapq,
			tstCigar.size(),
			tstCigar.data(),
			0,
			0,
			static_cast<hts_pos_t>( tstSeq.size() ),
			tstSeq.size(),
			tstSeq.c_str(),
			tstQual.c_str(),
			0
		);
		int headRes = sam_hdr_add_line(tstBAMheader.get(), "HD", "VN", "1.6", "SO", "unknown", NULL); // NOLINT
		headRes     = sam_hdr_add_line(tstBAMheader.get(), "SQ", "SN", "test", "LN", "100000", NULL); // NOLINT
		isaSpace::BAMrecord tstBAMrecord( bamRecordPtr.get(), tstBAMheader.get() );
		constexpr hts_pos_t earlyES{30};
		constexpr hts_pos_t midES{90};
		constexpr hts_pos_t lateES{245};
		constexpr hts_pos_t earlyEE{45};
		constexpr hts_pos_t midEE{110};
		constexpr hts_pos_t lateEE{345};

		isaSpace::ReadExonCoverage tstREC;
		tstREC.readName               = "m54312U_201215_225530/460252/ccs";
		tstREC.chromosomeName         = "NC_052529.2";
		tstREC.strand                 = '+';
		tstREC.geneName               = "gene-LOC6532627";
		tstREC.firstExonStart         = earlyES;
		tstREC.lastExonEnd            = lateEE;
		tstREC.exonCoverageScores     = std::vector<float>{0.0, 0.93, 1.0, 0.5}; // NOLINT
		tstREC.bestExonCoverageScores = std::vector<float>{0.5, 0.93, 1.0, 0.5}; // NOLINT

		REQUIRE( isaSpace::rangesOverlap(tstREC, tstBAMrecord) );
		tstREC.firstExonStart = earlyES;
		tstREC.lastExonEnd    = midEE;
		REQUIRE( isaSpace::rangesOverlap(tstREC, tstBAMrecord) );
		tstREC.firstExonStart = earlyES;
		tstREC.lastExonEnd    = earlyEE;
		REQUIRE( !isaSpace::rangesOverlap(tstREC, tstBAMrecord) );
		tstREC.firstExonStart = midES;
		tstREC.lastExonEnd    = lateEE;
		REQUIRE( isaSpace::rangesOverlap(tstREC, tstBAMrecord) );
		tstREC.firstExonStart = lateES;
		tstREC.lastExonEnd    = lateEE;
		REQUIRE( !isaSpace::rangesOverlap(tstREC, tstBAMrecord) );
		tstREC.firstExonStart = midES;
		tstREC.lastExonEnd    = midEE;
		REQUIRE( isaSpace::rangesOverlap(tstREC, tstBAMrecord) );
	}

	SECTION("Peak and valley functions") {
		// Peaks function
		// one peak
        const std::vector<float> values = {1.0F, 2.0F, 3.0F, 4.0F, 5.0F, 4.0F, 3.0F, 2.0F, 1.0F};
        constexpr float threshold{3.0F};
		constexpr float correctPeakValue1{5.0F};
        std::vector<std::vector<float>::const_iterator> peaks1{isaSpace::getPeaks(values, threshold)};
        REQUIRE(peaks1.size() == 1);
        REQUIRE(*peaks1.front() == correctPeakValue1);
	}

	SECTION("Stringify functions") {
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
		REQUIRE(isaSpace::stringify(coverageStats) == correctStatsLine);

		constexpr size_t recVecSize{7};
		const std::vector<isaSpace::ReadExonCoverage> recVec(recVecSize, coverageStats);
		const auto testVecResult{isaSpace::stringifyRCSrange(recVec.cbegin() + 1, recVec.cbegin() + 4)};
		REQUIRE(testVecResult.size() == 3 * correctStatsLine.size() + 3); // the newlines are the extra three characters
		REQUIRE(
			std::includes(
				testVecResult.cbegin(),
				testVecResult.cend(),
				correctStatsLine.cbegin(),
				correctStatsLine.cend()
			)
		);
		// thread ranges function
		constexpr size_t nThreads{4};
		const std::vector<isaSpace::ReadExonCoverage> testRECvec(11);
		const auto threadRanges{isaSpace::makeThreadRanges(testRECvec, nThreads)};
		REQUIRE(threadRanges.size() == nThreads);
		REQUIRE(
			std::all_of(
				threadRanges.cbegin(),
				threadRanges.cend(),
				[](const std::pair<std::vector<isaSpace::ReadExonCoverage>::const_iterator, std::vector<isaSpace::ReadExonCoverage>::const_iterator> &currPair) {
					return std::distance(currPair.first, currPair.second) >= 0;
				}
			)
		);
		REQUIRE(
			std::all_of(
				threadRanges.cbegin(),
				threadRanges.cend(),
				[&testRECvec, &nThreads](const std::pair<std::vector<isaSpace::ReadExonCoverage>::const_iterator, std::vector<isaSpace::ReadExonCoverage>::const_iterator> &currPair) {
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
	REQUIRE(testExonGroupNeg[1].first  == testExonSpans.at(1).first);
	REQUIRE(testExonGroupNeg[1].second == testExonSpans.at(1).second);
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
	REQUIRE(testExonGroupNeg[1].first  == testExonSpans.at(1).first);
	REQUIRE(testExonGroupNeg[1].second == testExonSpans.at(1).second);
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
	REQUIRE(testExonGroupNeg[1].first  == testExonSpans.at(1).first);
	REQUIRE(testExonGroupNeg[1].second == testExonSpans.at(1).second);
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

	// throwing on empty set
	std::set< std::pair<hts_pos_t, hts_pos_t> > emptyExonSet;
	REQUIRE_THROWS_WITH(
		isaSpace::ExonGroup(testGeneName, strand, emptyExonSet),
		Catch::Matchers::StartsWith("ERROR: set of exons is empty in")
	);
}

TEST_CASE("Read match window statistics work") {
	constexpr size_t vecLength{20};
	constexpr std::vector< std::pair<float, hts_pos_t> >::difference_type subWindow{5};
	constexpr float hiProb{0.99F};
	constexpr float loProb{0.25F};
	isaSpace::BinomialWindowParameters windowParameters;
	windowParameters.alternativeProbability = loProb;
	windowParameters.currentProbability     = hiProb;
	windowParameters.windowSize             = vecLength;

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

	windowParameters.currentProbability     = loProb;
	windowParameters.alternativeProbability = hiProb;
	isaSpace::ReadMatchWindowBIC loProbSwitch(window.cbegin(), windowParameters);
	REQUIRE( loProbSwitch.getBICdifference() == Catch::Approx(correctBIClo) );
}

TEST_CASE("Reading individual BAM records works") {
	// straight read
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
	constexpr char openMode{'r'};
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
	REQUIRE( bamRecord.isPrimaryMap() );
	REQUIRE( !bamRecord.isSecondaryMap() );
	REQUIRE(bamRecord.getReadName()    == correctReadName);
	REQUIRE(bamRecord.getmRNAstart()   == correctMapPosition);
	REQUIRE(bamRecord.getMapStart()    == correctMapPosition);
	REQUIRE(bamRecord.getMapEnd()      == correctMapEndPosition);
	REQUIRE(bamRecord.getReadLength()  == correctReadLength);
	REQUIRE(bamRecord.getCIGARstring() == correctCIGAR);

	const std::vector<uint32_t> cigarVec{bamRecord.getCIGARvector()};
	REQUIRE(
		std::equal( cigarVec.cbegin(), cigarVec.cend(), correctCIGV.cbegin() )
	);
	REQUIRE(bamRecord.getReferenceName() == correctReferenceName);

	std::vector<float> referenceMatches{bamRecord.getReferenceMatchStatus()};
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

	// reverse-complemented read
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
	nBytes = bam_read1( orRevBAMfile.get(), bamRevRecordPtr.get() );
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

	const std::vector<uint32_t> cigarVecRev{bamRecordRev.getCIGARvector()};
	REQUIRE(
		std::equal( cigarVecRev.cbegin(), cigarVecRev.cend(), correctCIGVrev.cbegin() )
	);

	referenceMatches = bamRecordRev.getReferenceMatchStatus();
	REQUIRE(referenceMatches.size() == correctRevRefMatchLength);
	matchSum = std::accumulate(referenceMatches.cbegin(), referenceMatches.cend(), 0.0F);
	REQUIRE( matchSum == static_cast<float>(correctRevReadLength) );

	readMatches = bamRecordRev.getReadCentricMatchStatus();
	REQUIRE(readMatches.front().second == correctRevMapPosition);
	REQUIRE(
		std::all_of(
			readMatches.cbegin(),
			readMatches.cend(),
			[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.first == 1.0;}
		)
	);
	posDiffs.clear();
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

	// Read alignment with a soft clip
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
	nBytes = bam_read1( orSoftBAMfile.get(), bamSoftRecordPtr.get() );
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

	const std::vector<uint32_t> cigarVecSoft{bamRecordSoft.getCIGARvector()};
	REQUIRE(
		std::equal( cigarVecSoft.cbegin(), cigarVecSoft.cend(), correctCIGVsoft.cbegin() )
	);

	referenceMatches = bamRecordSoft.getReferenceMatchStatus();
	REQUIRE(referenceMatches.size() == correctSoftRefMatchLength);
	matchSum = std::accumulate(referenceMatches.cbegin(), referenceMatches.cend(), 0.0F);
	REQUIRE(matchSum == correctSoftMatchCount);

	readMatches = bamRecordSoft.getReadCentricMatchStatus();
	REQUIRE(readMatches.front().second == correctSoftMapPosition);
	REQUIRE(
		std::count_if(
			readMatches.cbegin(),
			readMatches.cend(),
			[](const std::pair<float, hts_pos_t> &eachPair){return eachPair.first == 0.0;}
		) == correctSoftMisMatchCount
	);
	posDiffs.clear();
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
	isaSpace::BAMtoGenome testBAM(gffPair);
	testBAM.saveReadCoverageStats(outFileName, nThreads);

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
	constexpr size_t correctResSize{11};
	constexpr int32_t correctNpos{6};
	constexpr int32_t correctNneg{5};
	constexpr int32_t correctNgenes{8};
	constexpr int32_t correctNfailed{3};
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
				return name == "past_last_mRNA";
			}
		) == 1
	);
	REQUIRE(
		std::count_if(
			geneNames.cbegin(),
			geneNames.cend(),
			[](const std::string &name) {
				return name == "no_overlap";
			}
		) == 2
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
}
