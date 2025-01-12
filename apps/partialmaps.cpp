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


/// Find reads with unmapped regions
/** \file
 * \author Anthony J. Greenberg and Rebekah Rogers
 * \copyright Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.1
 *
 * Goes through isoSeq read alignments and looks for reads that have unmapped regions.
 *
 */

#include <string>
#include <unordered_map>
#include <iostream>
#include <thread>

#include "isoseqAlgn.hpp"
#include "helperFunctions.hpp"

int main(int argc, char *argv[]) {
	// set usage message
	const std::string cliHelp = "Available command line flags (in any order):\n" 
		"  --input-bam   bam_file_name (input BAM file name; required).\n"
		"  --input-gff   gff_file_name (input GFF file name; required).\n"
		"  --out         out_file_name (output file name; required).\n"
		"  --threads     number_of_threads (maximal number of threads to use; defaults to maximal available).\n";
	try {
		std::unordered_map <std::string, std::string> stringVariables;
		std::unordered_map <std::string, int>         intVariables;
		const auto clInfo{isaSpace::parseCL(argc, argv)};
		isaSpace::extractCLinfo(clInfo, intVariables, stringVariables);
		size_t nThreads{0};
		if (intVariables.at("threads") < 1) {
			nThreads = static_cast<size_t>( std::thread::hardware_concurrency() );
		} else {
			nThreads = static_cast<size_t>( intVariables.at("threads") );
		}
		isaSpace::BamAndGffFiles bamAndGFF;
		bamAndGFF.bamFileName = stringVariables.at("input-bam");
		bamAndGFF.gffFileName = stringVariables.at("input-gff");
		constexpr float hiProb{0.99F};
		constexpr float loProb{0.25F};
		constexpr int32_t windowSize{80};
		isaSpace::BinomialWindowParameters windowParameters;
		windowParameters.currentProbability     = loProb;
		windowParameters.alternativeProbability = hiProb;
		windowParameters.windowSize             = windowSize;
		isaSpace::BAMtoGenome bam2genome(bamAndGFF);
		bam2genome.saveUnmappedRegions(stringVariables.at("out"), windowParameters, nThreads);
		return 0;
	} catch(std::string &problem) {
		std::cerr << problem << "\n";
		std::cerr << cliHelp;
		return 1;
	}
}
