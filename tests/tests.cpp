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
#include <iterator>
#include <algorithm>
#include <utility>
#include <set>
#include <array>
#include <string>

//#include "bgzf.h"
//#include "sam.h"

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
}

TEST_CASE("Record classes work") {
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
	REQUIRE(testExonGroupNeg.nExons() == correctNexons);
	auto fullSpan{testExonGroupNeg.geneSpan()};
	REQUIRE(fullSpan.first  == testExonSpans.front().first);
	REQUIRE(fullSpan.second == testExonSpans.back().second);
	auto firstExon{testExonGroupNeg.firsExonSpan()};
	REQUIRE(firstExon.first  == testExonSpans.back().first);
	REQUIRE(firstExon.second == testExonSpans.back().second);
	// positive strand
	strand = '+';
	isaSpace::ExonGroup testExonGroupPos(testGeneName, strand, testSet);
	REQUIRE(testExonGroupPos.nExons() == correctNexons);
	fullSpan = testExonGroupPos.geneSpan();
	REQUIRE(fullSpan.first  == testExonSpans.front().first);
	REQUIRE(fullSpan.second == testExonSpans.back().second);
	firstExon = testExonGroupPos.firsExonSpan();
	REQUIRE(firstExon.first  == testExonSpans.front().first);
	REQUIRE(firstExon.second == testExonSpans.front().second);
	// undetermined strand, assumed positive
	strand = '.';
	isaSpace::ExonGroup testExonGroupUnd(testGeneName, strand, testSet);
	REQUIRE(testExonGroupUnd.nExons() == correctNexons);
	fullSpan = testExonGroupUnd.geneSpan();
	REQUIRE(fullSpan.first  == testExonSpans.front().first);
	REQUIRE(fullSpan.second == testExonSpans.back().second);
	firstExon = testExonGroupUnd.firsExonSpan();
	REQUIRE(firstExon.first  == testExonSpans.front().first);
	REQUIRE(firstExon.second == testExonSpans.front().second);

	// throwing on empty set
	std::set< std::pair<hts_pos_t, hts_pos_t> > emptyExonSet;
	REQUIRE_THROWS_WITH(
		isaSpace::ExonGroup(testGeneName, strand, emptyExonSet),
		Catch::Matchers::StartsWith("ERROR: set of exons is empty in")
	);
}

TEST_CASE("Catching bad GFF and BAM files works") {
	const std::string goodGFFname("../tests/goodGFF.gff");
	const std::string nomrnaGFFname("../tests/nomRNA.gff");
	const std::string goodBAMname("../tests/testAlgn.bam");
	const std::string wrongBAMname("../tests/wrong.bam");
	const std::string randomNoiseBAMname("../tests/randomNoise.bam"); // completely random noise, no magic bytes, but valid EOF sequence
	const std::string headerlessBAMname("../tests/headerless.bam");   // random noise with magic bytes and EOF but no header
	isaSpace::BamAndGffFiles gffPair;
	gffPair.gffFileName = goodGFFname;
	gffPair.bamFileName = wrongBAMname;
	REQUIRE_THROWS_WITH(
		isaSpace::FirstExonRemap(gffPair),
		Catch::Matchers::StartsWith("ERROR: failed to open the BAM file")
	);
	gffPair.bamFileName = randomNoiseBAMname;
	REQUIRE_THROWS_WITH(
		isaSpace::FirstExonRemap(gffPair),
		Catch::Matchers::StartsWith("ERROR: failed to read the header from the BAM file")
	);
	gffPair.bamFileName = headerlessBAMname;
	REQUIRE_THROWS_WITH(
		isaSpace::FirstExonRemap(gffPair),
		Catch::Matchers::StartsWith("ERROR: failed to read the header from the BAM file")
	);
	gffPair.gffFileName = nomrnaGFFname;
	gffPair.bamFileName = goodBAMname;
	REQUIRE_THROWS_WITH(
		isaSpace::FirstExonRemap(gffPair),
		Catch::Matchers::StartsWith("ERROR: no mRNAs with exons found in the")
	);
}

TEST_CASE("GFF and BAM parsing works") {
	const std::string goodGFFname("../tests/goodGFF.gff");
	const std::string goodBAMname("../tests/testAlgn.bam");
	isaSpace::BamAndGffFiles gffPair;
	gffPair.gffFileName = goodGFFname;
	gffPair.bamFileName = goodBAMname;
	constexpr size_t correctNsets{4};
	constexpr size_t correctNchrom{2};

	isaSpace::FirstExonRemap parsedGoodGFF(gffPair);
	REQUIRE(parsedGoodGFF.nChromosomes() == correctNchrom);
	REQUIRE(parsedGoodGFF.nExonSets()    == correctNsets);

	/*
	const std::string messyGFFname("../tests/messyGFF.gff");
	gffPair.gffFileName = messyGFFname;
	isaSpace::FirstExonRemap parsedMessyGFF(gffPair);
	REQUIRE(parsedMessyGFF.nChromosomes() == correctNchrom);
	REQUIRE(parsedMessyGFF.nExonSets()    == correctNsets);
	*/
}
/*
TEST_CASE("HTSLIB doodles") {
	const std::string testBAMname("../tests/testAlgn.bam");
	const char openMode{'r'};
	std::unique_ptr<BGZF, void(*)(BGZF *)> testBAMfile(
		bgzf_open(testBAMname.c_str(), &openMode),
		[](BGZF *bamFile) {
			bgzf_close(bamFile);
		}
	);

	// must read the header first to get to the alignments
	std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> testBAMheader(
		bam_hdr_read( testBAMfile.get() ),
		[](sam_hdr_t *samHeader) {
			sam_hdr_destroy(samHeader);
		}
	);

	const std::string testHeaderText{sam_hdr_str( testBAMheader.get() )};
	std::cout << "=================\n";
	constexpr uint16_t isNotPrimary{BAM_FSECONDARY | BAM_FSUPPLEMENTARY};
	std::vector<isaSpace::SAMrecord> primaryRecords;
	uint32_t iRecord{1};
	while (true) {
		isaSpace::CbamRecordDeleter localDeleter;
		std::unique_ptr<bam1_t, isaSpace::CbamRecordDeleter> bamRecordPtr(bam_init1(), localDeleter);
		int32_t nBytes = bam_read1( testBAMfile.get(), bamRecordPtr.get() );
		if (nBytes == -1) {
			break;
		}
		if (nBytes < -1) {
			continue;
		}
		++iRecord;
		const std::string currQname( bam_get_qname( bamRecordPtr.get() ) ); // NOLINT
		auto *const currCIGARptr = bam_get_cigar( bamRecordPtr.get() );     // NOLINT
		// is it a primary alignment with a soft clip?
		if ( (bam_cigar_opchr(*currCIGARptr) == 'S') && ( (bamRecordPtr->core.flag & isNotPrimary) == 0 ) ) { // NOLINT
			primaryRecords.emplace_back(bamRecordPtr);
			continue;
		}
		if ( ( !primaryRecords.empty() ) && ( currQname == primaryRecords.back().getReadName() ) ) {
			primaryRecords.back().appendSecondary(bamRecordPtr);
		}
		//std::vector<uint32_t> cigar(bam_get_cigar( bamRecords.back().get() ), bam_get_cigar( bamRecords.back().get() ) +  bamRecords.back()->core.n_cigar); // NOLINT
		//std::cout << "first CIGAR: " << bam_cigar_oplen( cigar.front() ) << bam_cigar_opchr( cigar.front() ) << "\n"; // NOLINT
		//std::cout << std::bitset<16>(0x900) << "\n" << std::bitset<16>(bamRecordPtr->core.flag) << "\n\n";//NOLINT
	}
	std::cout << "Number saved: " << primaryRecords.size() << "; number considered: " << iRecord << "\n";
	std::cout << "+++++++++++++++++\n";

}
*/
