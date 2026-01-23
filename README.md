# flexbar_rs

## Overview

This tool is a rewrite of Flexbar in Rust. The original Fexbar can be found here:

* https://github.com/seqan/flexbar

It is heavily based on Rust-bio.

## Quickstart for Flexbar Veterans

Most options work as for Flexbar, but:

* Multi-letter options including 2-letter shortcuts (e.g. `-qt` for `--trim-treshold`) must be written now
  with double-dash: `--qt`. You may now cluster single-letter options (e.g. `-df` instead of `-d -f`).
* The gap scores have been split into two parameters, e.g. `--barcode-gap` was replaced by `--barcode-gap-open`
  and `--barcode-gap-extend`.

Not yet implemented:

* Paired reads
* Trimming of homopolymers
* pre-trim-left/right
* min-read-length

Dropped features:

* Support for interleaved FASTQ files

## Building

To create an optimised production binary, run:
```bash
cargo +nightly multivers
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