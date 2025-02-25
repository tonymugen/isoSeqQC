# Overview

A C++14 library and software to assess the quality of [iso-Seq](https://www.pacb.com/wp-content/uploads/2018-10-NA-UGM-Iso-Seq-Method.pdf) read alignments. The project is still in development. Planned functionality includes:

    1. Assessing mapping quality and relating mapped regions to annotated genes.
    2. Identification of mistakes in mapping and read regions with poor mapping quality.
    3. Re-alignment of poorly mapped read portions and re-integration of the results into alignment BAM files.

Goals (1) and (2) have been completed. Goal (3) is under development.

# Dependencies

The project requires a C++14 compiler and depends on [htslib](https://github.com/samtools/htslib) for handling BAM files. A `cmake` script is included in the repository for easy building and installation. It requires `cmake` version 3.21 or later and automatically downloads `htslib`. If one wants to optionally build tests, Catch2 is also required and is automatically downloaded using the `cmake` script.

# Download and install

Clone the repository:

```{sh}
git clone https://github.com/tonymugen/isoSeqQC
```
Next, create a build directory

```{sh}
cd isoSeqQC
mkdir build
```
Finally, run `cmake` to build and install the software

```{sh}
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .
```
Installation may require root privileges.

# Use `exoncoverage` 

`exoncoverage` is a command line tool that takes a BAM alignment and a GFF annotation file and reports exon coverage statistics for each BAM read record, taking into account secondary alignments. Running it without any command line flags prints usage information. The flags are

```{sh}
  --input-bam   bam_file_name (input BAM file name; required).
  --input-gff   gff_file_name (input GFF file name; required).
  --out         out_file_name (output file name; required).
  --threads     number_of_threads (maximal number of threads to use; defaults to maximal available).
```

and can be specified in any order.

The output is a tab-delimited text file with largely self-explanatory column names. Columns that require special explanation are:

    best_alignment_start          -- start position of best alignment, incorporating secondary alignments.
    best_alignment_end            -- the same for the end position of the alignment.
	n_good_secondary_alignments   -- number of secondary alignments that are cover the same region as the primary and are on the same strand.
    n_local_reversed_strand       -- number of secondary alignments that are on the opposite strand but still in vicinity of the primary.
    alignment_quality_string      -- a string of alignment match fractions for each exon, by position on the reference regardless of transcipt arrangement or strand.
    best_alignment_quality_string -- the same, but taking into account secondary alignments.

# Use `partialmaps` 

`partialmaps` is a command line tool that takes a BAM alignment file and a GFF annotation file, identifies poorly mapped regions, exports data on these alignment gaps, and optionally saves the unmapped sequences to a FATSQ file for subsequent realignment. Alignment gap identification relies on a change point detection algorithm that uses sliding windows. Sliding window size is a user-controlled parameter. Running `paritalmaps` without any command line flags prints usage information. The flags are

```{sh}
  --input-bam    bam_file_name (input BAM file name; required).
  --input-gff    gff_file_name (input GFF file name; required).
  --out          out_file_name (output file name; required).
  --out-fastq    out_fastq_file_name (output FASTQ file name; optional, no FASTQ file created if omitted).
  --window-size  sliding window size for poorly aligned region identification (defaults to 75).
  --threads      number_of_threads (maximal number of threads to use; defaults to maximal available).
```

and can be specified in any order.

The output is a tab-delimited text file with the following columns:

    read_name      -- read name.
    read_length    -- read length.
	unmapped_start -- 0-based index of the unmapped region start in the read.
	unmapped_end   -- 0-based index of the unmapped region end in the read.
    window_size    -- sliding window size.

# Tests

Unit tests can be optionally built, without installing the software on the system:

```{sh}
mkdir buildTest
cd buildTest
cmake -DCMAKE_BUILD_TYPE=Test ..
cmake --build .
./tests
```

# Funding

This project is funded in part by an NIH NIGMS MIRA R35 GM133376 grant to Rebekah L. Rogers.
