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
#include <cstdint>
#include <iterator>
#include <memory>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>

#include <iostream>

#include "bgzf.h"
#include "sam.h"

#include "isoseqAlgn.hpp"

#include "catch2/catch_test_macros.hpp"
//#include "catch2/matchers/catch_matchers.hpp"
//#include "catch2/matchers/catch_matchers_string.hpp"

TEST_CASE("GFF parsing works") {
	constexpr size_t nGFFfields{9};
	constexpr char gffDelimiter{'\t'};
	constexpr char attrDelimiter{';'};
	constexpr char attrDelimEscape{'\\'};
	const std::string parentToken("Parent=");
	const auto parentTokenSize = std::distance( parentToken.cbegin(), parentToken.cend() );
	const std::string idToken("ID=");
	const auto idTokenSize = std::distance( idToken.cbegin(), idToken.cend() );
	std::unordered_map<std::string, std::string> mRNAtoGeneName;
	const std::string goodGFFname("../tests/goodGFF.gff");
	std::vector<isaSpace::ExonGroup> exonGroups;
	std::string gffLine;
	std::fstream goodGFF(goodGFFname, std::ios::in);
	while ( std::getline(goodGFF, gffLine) ) {
		std::array<std::string, nGFFfields> gffFields;
		std::stringstream currLineStream(gffLine);
		size_t iField{0};
		while ( (iField < nGFFfields) && std::getline(currLineStream, gffFields.at(iField), gffDelimiter)) {
			++iField;
		}
		if (iField < nGFFfields) {
			// in the actual implementation, save the line number where there is an incorrect number of fields
			continue;
		}
		if (gffFields.at(2) == "mRNA") {
			std::stringstream attributeStream( gffFields.back() );
			std::vector<std::string> attributes;
			std::string attrField;
			while ( std::getline(attributeStream, attrField, attrDelimiter) ) {
				attributes.emplace_back(attrField);
			}
			auto mRNAidIt = std::find_if(
				attributes.cbegin(),
				attributes.cend(),
				[&idToken](const std::string &eachAttr) {
					return std::equal( idToken.cbegin(), idToken.cend(), eachAttr.cbegin() );
				}
			);
			if ( mRNAidIt == attributes.cend() ) {
				// also save the line where there is a missing ID for an mRNA
				continue;
			}
			std::string mRNAid;
			std::copy(
				mRNAidIt->cbegin() + idTokenSize,
				mRNAidIt->cend(),
				std::back_inserter(mRNAid)
			);
			// deal with any escaped ';' delimiters
			std::advance(mRNAidIt, 1);
			while ( (mRNAid.back() == attrDelimEscape) && ( mRNAidIt != attributes.cend() ) ) {
				mRNAid.push_back(attrDelimiter);
				std::copy(
					mRNAidIt->cbegin(),
					mRNAidIt->cend(),
					std::back_inserter(mRNAid)
				);
			}
			auto parentIt = std::find_if(
				attributes.cbegin(),
				attributes.cend(),
				[&parentToken](const std::string &eachAttr) {
					return std::equal( parentToken.cbegin(), parentToken.cend(), eachAttr.cbegin() );
				}
			);
			if ( parentIt == attributes.cend() ) {
				// also save the line where there is a missing parent for an mRNA
				continue;
			}
			std::string parentName;
			std::copy(
				parentIt->cbegin() + parentTokenSize,
				parentIt->cend(),
				std::back_inserter(parentName)
			);
			// deal with any escaped ';' delimiters
			std::advance(parentIt, 1);
			while ( (parentName.back() == attrDelimEscape) && ( parentIt != attributes.cend() ) ) {
				parentName.push_back(attrDelimiter);
				std::copy(
					parentIt->cbegin(),
					parentIt->cend(),
					std::back_inserter(parentName)
				);
			}
			std::cout << "mRNA ID: " << mRNAid << "; parnet ID: " << parentName << "\n";
			mRNAtoGeneName[mRNAid] = parentName;
			continue;
		}
		//break;
	}
	goodGFF.close();
	std::cout << "mRNA to gene ID:\n";
	for (const auto &eachPair : mRNAtoGeneName) {
		std::cout << eachPair.first << " <--> " << eachPair.second << "\n";
	}
}

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
