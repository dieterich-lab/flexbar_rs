# flexbar_rs

## Overview

This tool is a rewrite of Flexbar in Rust. The original Fexbar can be found here:

* https://github.com/seqan/flexbar

It is heavily based on Rust-bio and uses the rust-htslib library to handle BAM files.

## Code Maturity

This code is BETA. It is working and feature complete, but lacks extensive use in the field.

## Quickstart for Flexbar Veterans

### Why should you use it?

* It complies in 2026 on up-to-date operating systems with modern dependencies. 
* It handles ZSTD compression out-of-the-box (see tuning hints below!).
* It supports BAM files as input (-r/--reads option) and as output. With
  the `--barcode-no-split` option it is possible  to populate the BC auxiliary
  tag with the barcode ID without splitting the file.

### What changed?

Most options work as for Flexbar, but:

* Multi-letter options including 2-letter shortcuts (e.g. `-qt` for `--trim-treshold`) must be written now
  with double-dash: `--qt`. You may now cluster single-letter options (e.g. `-df` instead of `-d -f`).
* The gap scores have been split into two parameters, e.g. `--barcode-gap` was replaced by `--barcode-gap-open`
  and `--barcode-gap-extend`.
* the `--fasta-output` was removed. Use `--output-format=FASTA` instead.

Dropped features we will not implement at all:

* Support for interleaved FASTQ files
* Trimming of homopolymers
* A bunch of adapter-related options
  * --adapter-relaxed
  * --adapter-pair-overlap
  * --adapter-min-poverlap
  * --adapter-trimmed-out
  * --adapter-tail-length
* --barcode-tail-length
* --number-tags

### Are there known differences in the output?

Yes, differences have been observed with similar parameters:

* In some cases multiple optimal alignments exist with the same score.
  In this case flexbar and flexbar_rs may select different alignments
  and come to different results. Both flexbar and flexbar_rs do a
  validation step based on the X-error-rate and X-overlap parameters
  after selecting the optimal alignment. As a result, any of the two
  may miss a valid alignment.

* flexbar_rs may detect slightly more alignments in the case of short
  overlaps with trim end RIGHT, LEFT, and ANY and low values for the
  related overlap parameters. That seems even to happen if the above
  issue does not apply.

With real-life examples, such effects have been observed in about 0.1%
of the reads. Presently, no case is known in which flexbar_rs will use
a wrong alignment.

## Tuning

### Optimizing for Execution Time

**Stop using gzip – use ZSTD compression.** In typical HPC setups (>10 vCores)
it is essential to use compression – not only to save disk space but also
to reduce the I/O load. With typical flexbar workloads the time needed to
decompress and compress the data will dominate the overall runtime. When
using ZSTD compression, expect the wall-clock time to be cut by 50% with
similar compression rates and energy usage.

Of course, the number of threads should not exceed the number of vCPUs.
The benefit of adding more threads diminishes due to the overhead above
about 20 threads.

### Optimizing for Energy Consumption

While more threads may reduce the wall-clock time, expect to get
better results per Watt with thread counts below 20. 

### Optimizing for Memory Usage

If you run out of memory, consider using fewer threads. You may
also reduce the size of the IO buffers (--write-buffer, --read-buffer).
If you create lots of output files (=many barcodes). The write buffer
size is especially important because a separate one is kept for each
output file.

## Building

Before building make sure that you have checked out the version
you want, including the submodule. So don't forget `--recursive`
option with clone command:

```bash
git clone --recursive https://github.com/dieterich-lab/flexbar_rs.git
```

To build the htslib submodule a bunch of compile time depedenicies are
needed. On Debian-like systems most fo them can be installed like this:

```bash
sudo apt-get install build-essential libclang-dev
```

To create an optimized production binary, you need to first get some support software:

```bash
cargo install --locked cargo-multivers
cargo install cargo-pgo
rustup component add llvm-tools-preview
```

For details see:
* https://crates.io/crates/cargo-multivers
* https://github.com/Kobzol/cargo-pgo

To create binary run:

```bash
./scripts/do_prod_build.sh
```

The binary will be created in `target/cargo-multivers/x86_64-unknown-linux-gnu/release/flexbar_rs`.

## References

### Flexbar

Flexbar - flexible barcode and adapter removal, version 3.5.0
Developed with SeqAn, the library for sequence analysis

Johannes T. Roehr, Christoph Dieterich, Knut Reinert:
Flexbar 3.0 - SIMD and multicore parallelization. Bioinformatics 2017.

Matthias Dodt, Johannes T. Roehr, Rina Ahmed, Christoph Dieterich:
Flexbar - flexible barcode and adapter processing for next-generation
sequencing platforms. Biology 2012.

Available on https://github.com/seqan/flexbar

### Rust-bio

Köster, J. (2016). Rust-Bio: a fast and safe bioinformatics library. Bioinformatics, 32(3), 444-446.

DOI: [29821195](https://zenodo.org/badge/latestdoi/29821195)

### rust-htslib

The rust-htslib library is used to handle BAM files and is
integrated as a git submodule.

Available on https://github.com/rust-bio/rust-htslib
