# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```sh
# Configure and build (Release)
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .

# Install (may require root)
cmake --install .

# Build and run tests (uses address+leak sanitizers automatically)
mkdir buildTest && cd buildTest
cmake -DCMAKE_BUILD_TYPE=Test -DBUILD_TESTS=ON ..
cmake --build .
./tests

# Run a single test by name
./tests "Safe BAM reading works"

# Build documentation (requires doxygen and pdflatex)
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCS=ON ..
cmake --build .
```

Build types: `Release`, `Debug`, `Profile`, `Test`. The `Test` build enables Catch2 (fetched automatically) and address/leak sanitizers.

## Architecture

This is a C++17 library (`isoseqAln`) and three CLI tools that assess quality of [Iso-Seq](https://www.pacb.com/wp-content/uploads/2018-10-NA-UGM-Iso-Seq-Method.pdf) long-read RNA alignments. htslib is fetched and built automatically by CMake via `ExternalProject_Add`.

**Source layout:**
- `include/` — public headers (`isoseqAlgn.hpp`, `helperFunctions.hpp`); API documentation lives here as Doxygen comments
- `src/` — library implementation (`isoseqAlgn.cpp`, `helperFunctions.cpp`)
- `apps/` — CLI entry points: `exoncoverage.cpp`, `partialmaps.cpp`, `addremaps.cpp`
- `tests/` — Catch2 test suite (`tests.cpp`) and test data (BAM/GFF/SAM files)

All code lives in `namespace isaSpace`.

**Core data flow:**

1. `BAMsafeReader` — wraps htslib BGZF file I/O; temporarily sets BAM file read-only during reads to guard against an observed htslib bug that can delete files
2. `BAMrecord` — parses a primary BAM record plus its local secondary alignments from raw `bam1_t`; exposes CIGAR-based match status vectors and poorly-mapped region detection via a sliding-window BIC change-point algorithm (`ReadMatchWindowBIC`)
3. `ExonGroup` — represents all exons of a gene parsed from GFF3; provides exon coverage quality calculations against a `BAMrecord`
4. `BAMtoGenome` — pairs `BAMrecord` objects with their overlapping `ExonGroup` by streaming both the BAM and GFF together; drives `saveReadCoverageStats()` and `saveUnmappedRegions()` used by `exoncoverage` and `partialmaps`
5. `BAMfile` — loads all primary BAM records for `addremaps`, which merges re-mapped read portions as secondary alignments and rewrites CIGAR strings

**Helper functions** (`helperFunctions.hpp/cpp`) handle GFF parsing (`parseGFF`, `parseGFFline`), CIGAR-to-match-vector conversion, thread range partitioning (`makeThreadRanges`), CLI parsing (`parseCL`, `extractCLinfo`), and BAM record manipulation (`modifyCIGAR`, `addRemappedSecondaryAlignment`).

**Parallelism:** `saveReadCoverageStats` and `saveUnmappedRegions` split work using `makeThreadRanges` and `std::thread`.

## Code Style

- Compiler warnings are treated strictly (Wall, Wextra, Wconversion, Wpedantic, Wshadow, etc. — see `CMakeLists.txt`)
- `clang-tidy` is configured via `.clang-tidy` (most checks enabled; excludes abseil, boost, google, llvm, and a few others)
- Use `[[nodiscard]]` on all functions returning values that should not be ignored
- Custom deleters (`BAMheaderDeleter`, `BAMrecordDeleter`, `BGZFhandleDeleter`) wrap htslib raw pointers into `std::unique_ptr`
- Positions: read positions are 0-based, reference/GFF positions are 1-based throughout
- All variable and function names are in camelCase except for `main`. Class and struct names are always capitalized.
- when a closing parenthesis is immediately followed by another, insert a space before and after the internal function call
