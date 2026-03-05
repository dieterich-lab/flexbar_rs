# flexbar_rs

## Overview

This tool is a rewrite of Flexbar in Rust. The original Fexbar can be found here:

* https://github.com/seqan/flexbar

It is heavily based on Rust-bio.

## Quickstart for Flexbar Veterans

### Why should you use it?

* It complies in 2026 on up-to-date operating systems with modern dependencies. 
* It handles ZSTD compression out-of-the-box (see tuning hints below).

### What changed?

Most options work as for Flexbar, but:

* Multi-letter options including 2-letter shortcuts (e.g. `-qt` for `--trim-treshold`) must be written now
  with double-dash: `--qt`. You may now cluster single-letter options (e.g. `-df` instead of `-d -f`).
* The gap scores have been split into two parameters, e.g. `--barcode-gap` was replaced by `--barcode-gap-open`
  and `--barcode-gap-extend`.

Not yet implemented:

* Paired reads
* Trimming of homopolymers
* pre-trim-left/right

Dropped features we will not implement at all:

* Support for interleaved FASTQ files

## Tuning

### Optimizing for Execution Time

**Stop using gzip - use ZSTD compression.** In typical HPC setups (>20 vCPUs)
it is essential to use compression - not only to save disk space but also
to reduce the I/O load. With typical flexbar workloads the time needed to
decompress and compress the data will dominate the overall runtime. When
using ZSTD compression, expect the wall-clock time to be cut by 50% with
similar compression rates and energy usage.

Of course, the number of threads should not exceed the number of vCPUs.
The benefit of adding more threads diminishes due to the overhead above
about 20 threads.

### Optimizing for Energy Consumption

While more threads may reduce the wall-clock time, expect that to get
better results per Watt with thread counts below 20. 

### Optimizing for Memory Usage

If you run out of memory, consider using fewer threads. You may
also reduce the size of the IO buffers (--write-buffer, --read-buffer).
If you create lots of output files (=many barcodes). The write buffer
size is especially important because a separate one is kept for each
output file.

## Building

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