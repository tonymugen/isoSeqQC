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

/// Genomic analyses helper functions
/** \file
 * \author Anthony J. Greenberg and Rebekah Rogers
 * \copyright Copyright (c) 2024 Anthony J. Greenberg and Rebekah Rogers
 * \version 0.1
 *
 * Definitions of class-external functions needed by genomic analyses.
 *
 */

#pragma once

#include <string>

#include "isoseqAlgn.hpp"

namespace isaSpace {
	/** \brief Extract a name 
	 *
	 * Extract a token name from GFF attributes.
	 *
	 * \param[in] tokenAndAttrList field token and the list of attributes
	 */
	[[gnu::warn_unused_result]] std::string extractAttributeName(const TokenAttibuteListPair &tokenAndAttrList);
	/** \brief Convert `ReadExonCoverage` to string
	 *
	 * \param[in] readRecord individual read record
	 * \param[in] separator field separator
	 * \return `std::string` with the read record elements, without a new line at the end
	 */
	[[gnu::warn_unused_result]] std::string stringify(const ReadExonCoverage &readRecord, char separator = '\t');
}
