/*
 * Copyright (c) 2025 Anthony J. Greenberg and Rebekah Rogers
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


/// Add remapped read portions as secondary alignments
/** \file
 * \author Anthony J. Greenberg and Rebekah Rogers
 * \copyright Copyright (c) 2025 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.2
 *
 * Goes through read portion realignments and add the successfully re-mapped regions as secondary alignments.
 *
 */

#include <string>
#include <unordered_map>
#include <iostream>

#include "isoseqAlgn.hpp"
#include "helperFunctions.hpp"

int main(int argc, char *argv[]) {
	// set usage message
	const std::string cliHelp = "Available command line flags (in any order):\n" 
		"  --input-bam        bam_file_name (input BAM file name; required).\n"
		"  --remapped-bam     remapped_bam_file_name (BAM file with remapped read portions; required).\n"
		"  --out              out_file_name (output BAM file name; required).\n"
		"  --remap-cutoff     identity cutoff for read portion remapping (defaults to 0.99).\n"
		"  --unsorted-output  if set (with no value), the output BAM file is unsorted (otherwise, sorted).\n";
	try {
		std::unordered_map <std::string, std::string> stringVariables;
		std::unordered_map <std::string, int>         intVariables;
		std::unordered_map <std::string, float>       floatVariables;
		auto clInfo{isaSpace::parseCL(argc, argv)};
		clInfo["input-gff"] = "NULL";
		isaSpace::extractCLinfo(clInfo, intVariables, floatVariables, stringVariables);

		isaSpace::BAMfile primaryBAM( stringVariables.at("input-bam") );
		primaryBAM.addRemaps( stringVariables.at("remapped-bam"), floatVariables.at("remap-cutoff") );
		if (stringVariables.at("unsorted-output") == "unset") {
			const std::vector<std::string> failedReads{primaryBAM.saveSortedRemappedBAM( stringVariables.at("out") )};
			if ( !failedReads.empty() ) {
                std::cerr << "WARNING: failed to save " << failedReads.size() << " reads" << "\n";
			}
			return 0;
		}
		const std::vector<std::string> failedReads{primaryBAM.saveRemappedBAM( stringVariables.at("out") )};
		if ( !failedReads.empty() ) {
			std::cerr << "WARNING: failed to save " << failedReads.size() << " reads" << "\n";
		}
		return 0;
	} catch(std::string &problem) {
		std::cerr << problem << "\n";
		std::cerr << cliHelp;
		return 1;
	}
}
