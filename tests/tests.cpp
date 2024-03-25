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

#include <cstdint>
#include <memory>
#include <string>
#include <iostream>

#include "bgzf.h"
#include "sam.h"

#include "catch2/catch_test_macros.hpp"
//#include "catch2/matchers/catch_matchers.hpp"
//#include "catch2/matchers/catch_matchers_string.hpp"

TEST_CASE("HTSLIB doodles") {
	const std::string testBAMname("../tests/testAlgn.bam");
	const char openMode{'r'};
	std::unique_ptr<BGZF, void(*)(BGZF *)> testBAMfile(
		bgzf_open(testBAMname.c_str(), &openMode),
		[](BGZF *bamFile) {
			bgzf_close(bamFile);
		}
	);
	std::unique_ptr<sam_hdr_t, void(*)(sam_hdr_t *)> testBAMheader(
		bam_hdr_read( testBAMfile.get() ),
		[](sam_hdr_t *samHeader) {
			sam_hdr_destroy(samHeader);
		}
	);

	const std::string testHeaderText{sam_hdr_str( testBAMheader.get() )};
	std::cout << testHeaderText << "\n";

	std::unique_ptr<bam1_t, void(*)(bam1_t *)> firstBAMalnPtr(
		bam_init1(),
		[](bam1_t *bamPtr){
			bam_destroy1(bamPtr);
		}
	);
	int32_t nBytes = bam_read1( testBAMfile.get(), firstBAMalnPtr.get() );

	const std::string queryName{bam_get_qname( firstBAMalnPtr.get() )}; // NOLINT
	std::cout << "QNAME: " << queryName << "; number of bytes: " << nBytes << "\n";

	auto *const alignmentScoreRecord{bam_aux_get(firstBAMalnPtr.get(), "AS")};
	std::cout << "=================\n";
	if (alignmentScoreRecord != nullptr) {
		std::cout << "value = " << bam_aux2i(alignmentScoreRecord) << "\n";
	}
	std::cout << "+++++++++++++++++\n";

}
