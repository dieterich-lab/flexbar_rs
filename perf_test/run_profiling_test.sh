#!/bin/bash

# To run this script you need to point the input_file variable to an input
# file not included in the repo. That should be a fastq.gz file - 10,000 reads
# seem to be a reasonable size. The file should use barcodes as listed in barcodes.fa
# (barcode_file). You may have to replace them too. If you are not interested in
# the thread-join overhead best run with "threads=1". You will need to install flamegrapth:
#
# cargo install flamegraph
#
# Also add this temporarily to your cargo.toml:
#
#   [target.debug]
#   linker = "/usr/bin/clang"
#   rustflags = ["-Clink-arg=-fuse-ld=lld", "-Clink-arg=-Wl,--no-rosegment"]
#
# You also need to have build a binary with:
#
#   cargo build
#
# Do that again after reverting your cargo.toml.
#
# In addition you may have to install perf:
#
#   sudo apt-get install perf-tools-unstable
#
# And relax your kernel's security checks:
#
#   echo 2 > /proc/sys/kernel/perf_event_paranoid
#

set -eu

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
input_file="$script_dir/EV19002.R1_10000.fastq.gz"
output_prefix="$script_dir/out/flexbar"
barcode_file="$script_dir/barcodes.fa"
flexbar_bin="$script_dir/../target/debug/flexbar_rs"
threads=1

rm -rf "$script_dir/out"
flamegraph --freq 99 -- $flexbar_bin -r "$input_file" -t "$output_prefix" -n $threads \
  -z GZ -u 2 -q TAIL --qt 10 -k 248 -b "$barcode_file" \
  --bo 11 --bt RIGHT --qf sanger -m 16 -g -d
