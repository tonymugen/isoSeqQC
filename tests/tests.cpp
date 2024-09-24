/*
 * Copyright (c) 2023 Anthony J. Greenberg and Rebekah Rogers
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
#include <htslib/sam.h>
#include <iterator>
#include <algorithm>
#include <utility>
#include <set>
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
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

TEST_CASE("Helper functions work") {
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

	// stringify tests
	const std::string correctStatsLine("m54312U_201215_225530/460252/ccs	NC_052529.2	+	2983523	2988359	0	gene-LOC6532627	4	2980697	2988366	{0.000000,0.930000,1.000000,0.500000}");
	isaSpace::ReadExonCoverage coverageStats;
	coverageStats.readName            = "m54312U_201215_225530/460252/ccs";
	coverageStats.chromosomeName      = "NC_052529.2";
	coverageStats.strand              = '+';
	coverageStats.alignmentStart      = 2983523; // NOLINT
	coverageStats.alignmentEnd        = 2988359; // NOLINT
	coverageStats.firstSoftClipLength = 0;
	coverageStats.geneName            = "gene-LOC6532627";
	coverageStats.nExons              = 4;       // NOLINT
	coverageStats.firstExonStart      = 2980697; // NOLINT
	coverageStats.lastExonEnd         = 2988366; // NOLINT
	coverageStats.exonCoverageScores  = std::vector<float>{0.0, 0.93, 1.0, 0.5}; // NOLINT
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
	std::vector<uint32_t> cigarVec;
	std::copy( simpleCIGAR.cbegin(), simpleCIGAR.cend(), std::back_inserter(cigarVec) );
	std::vector<float> exonCoverage{testExonGroupPos.getExonCoverageQuality(cigarVec, testExonSpans.front().first)};
	REQUIRE( exonCoverage.size() == testExonGroupPos.nExons() );
	REQUIRE(
		std::all_of(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val == 1.0F;})
	);

	// Read starts early
	constexpr hts_pos_t earlyStartPos{50000};
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
	cigarVec.clear();
	std::copy( earlyCigarFields.cbegin(), earlyCigarFields.cend(), std::back_inserter(cigarVec) );
	exonCoverage.clear();
	exonCoverage = testExonGroupPos.getExonCoverageQuality(cigarVec, earlyStartPos);
	REQUIRE( exonCoverage.size() == testExonGroupPos.nExons() );
	REQUIRE(
		std::count_if(exonCoverage.cbegin(), exonCoverage.cend(), [](float val) {return val == 1.0F;}) == 1
	);
	REQUIRE(
		std::all_of(exonCoverage.cbegin(), exonCoverage.cend(), [&qsCutOff](float val) {return val > qsCutOff;})
	);

	// Read starts in the middle of an exon
	constexpr hts_pos_t lateStartPos{52300};
	constexpr std::array<uint32_t, 5> lateCigarFields{
		bam_cigar_gen(350,  BAM_CMATCH),
		bam_cigar_gen(2872, BAM_CREF_SKIP),
		bam_cigar_gen(581,  BAM_CMATCH),
		bam_cigar_gen(102,  BAM_CREF_SKIP),
		bam_cigar_gen(701,  BAM_CMATCH)
	};
	cigarVec.clear();
	std::copy( lateCigarFields.cbegin(), lateCigarFields.cend(), std::back_inserter(cigarVec) );
	exonCoverage.clear();
	exonCoverage = testExonGroupPos.getExonCoverageQuality(cigarVec, lateStartPos);
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
	cigarVec.clear();
	std::copy( negativeCigarFields.cbegin(), negativeCigarFields.cend(), std::back_inserter(cigarVec) );
	exonCoverage.clear();
	exonCoverage = testExonGroupNeg.getExonCoverageQuality(cigarVec, earlyStartPos);
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

TEST_CASE("Reading individual BAM records works") {
	// straight read
	const std::string oneRecordBAMname("../tests/oneRecord.bam");
	const std::string correctReadName("m54312U_201215_225530/657093/ccs");
	const std::string correctCIGAR("22M1D667M1552N264M485N193M361N75M81N196M53N1706M62N126M236N170M");
	constexpr std::array<uint32_t, 17> correctCIGV{
		352,  18,   10672, 24835, 4224,
		7763, 3088, 5779,  1200,  1299,
		3136, 851,  27296, 995,   2016,
		3779, 2720
	};
	constexpr hts_pos_t correctMapPosition{2468434};
	constexpr hts_pos_t correctMapEndPosition{2474684};
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
	isaSpace::CbamRecordDeleter localDeleter;
	std::unique_ptr<bam1_t, isaSpace::CbamRecordDeleter> bamRecordPtr(bam_init1(), localDeleter);
	auto nBytes = bam_read1( orBAMfile.get(), bamRecordPtr.get() );
	isaSpace::BAMrecord bamRecord( std::move(bamRecordPtr) );
	REQUIRE(nBytes > 0);
	REQUIRE( !bamRecord.isRevComp() );
	REQUIRE(bamRecord.getReadName()    == correctReadName);
	REQUIRE(bamRecord.getmRNAstart()   == correctMapPosition);
	REQUIRE(bamRecord.getMapStart()    == correctMapPosition);
	REQUIRE(bamRecord.getMapEnd()      == correctMapEndPosition);
	REQUIRE(bamRecord.getCIGARstring() == correctCIGAR);

	const std::vector<uint32_t> cigarVec{bamRecord.getCIGARvector()};
	REQUIRE(
		std::equal( cigarVec.cbegin(), cigarVec.cend(), correctCIGV.cbegin() )
	);

	// reverse-complemented read
	const std::string oneRecordRevBAMname("../tests/oneRecordRev.bam");
	const std::string correctRevReadName("m54312U_201215_225530/38470652/ccs");
	const std::string correctCIGARrev("454M68N701M61N133M61N1652M78N106M56N280M");
	constexpr std::array<uint32_t, 11> correctCIGVrev{
		4480, 899,  1696, 1251,  26432,
		979,  2128, 979,  11216, 1091, 7264
	};
	constexpr hts_pos_t correctRevMapPosition{22550967};
	constexpr hts_pos_t correctRevmRNAposition{22554617};
	constexpr hts_pos_t correctRevMapEndPosition = correctRevmRNAposition;
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
	std::unique_ptr<bam1_t, isaSpace::CbamRecordDeleter> bamRevRecordPtr(bam_init1(), localDeleter);
	nBytes = bam_read1( orRevBAMfile.get(), bamRevRecordPtr.get() );
	isaSpace::BAMrecord bamRecordRev( std::move(bamRevRecordPtr) );
	REQUIRE(nBytes > 0);
	REQUIRE( bamRecordRev.isRevComp() );
	REQUIRE(bamRecordRev.getReadName()    == correctRevReadName);
	REQUIRE(bamRecordRev.getmRNAstart()   == correctRevmRNAposition);
	REQUIRE(bamRecordRev.getMapStart()    == correctRevMapPosition);
	REQUIRE(bamRecordRev.getMapEnd()      == correctRevMapEndPosition);
	REQUIRE(bamRecordRev.getCIGARstring() == correctCIGARrev);
	const std::vector<uint32_t> cigarVecRev{bamRecordRev.getCIGARvector()};
	REQUIRE(
		std::equal( cigarVecRev.cbegin(), cigarVecRev.cend(), correctCIGVrev.cbegin() )
	);
}

TEST_CASE("Catching bad GFF and BAM files works") {
	// Since HTSLIB is not in my control, I cannot test everything
	// For example, testing the throw on wrong file name is iffy
	// because sometimes HTSLIB seems to create the file if there is none,
	// but not consistently
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
	std::vector<std::string> geneNames;
	std::vector<int32_t> nExons;
	std::vector<int32_t> fesValues;                 // first exon starts
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
		geneNames.emplace_back(field);
		lineStream >> field;
		nExons.push_back( stoi(field) );
		lineStream >> field;
		fesValues.push_back( stoi(field) );
	}
	saveResultFile.close();
	constexpr size_t correctResSize{9};
	constexpr int32_t correctNpos{5};
	constexpr int32_t correctNneg{4};
	constexpr int32_t correctNgenes{6};
	constexpr int32_t correctNfailed{3};

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
}
