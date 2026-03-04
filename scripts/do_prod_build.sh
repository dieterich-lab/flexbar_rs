#!/bin/bash

set -eu

#!/bin/bash

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
perf_test_dir="$script_dir/../perf_test"
input_file="$perf_test_dir/EV19002.R1_10000.fastq.gz"
output_prefix="$perf_test_dir/out/flexbar"
barcode_file="$perf_test_dir/barcodes.fa"
flexbar_bin="$perf_test_dir/../target/debug/flexbar_rs"
threads=1

echo '###'
echo '### Instrumenting binary ...'

cargo pgo instrument build

echo '###'
echo '### Collecting perfdata ...'
echo '###'

rm -rf "$perf_test_dir/out"
cargo pgo run -- -- -r "$input_file" -t "$output_prefix" -n $threads \
  -z GZ -u 2 -q TAIL --qt 10 -k 248 -b "$barcode_file" \
  --bo 11 --bt RIGHT --qf sanger -m 16 -g -d


echo '###'
echo '### Doing production build ...'
echo '###'

cargo pgo optimize

echo '###'
echo '### Doing multivers production build ...'
echo '###'

pgo_profile_dir=$(dirname "$script_dir")/target/pgo-profiles
pref_data="$(echo $pgo_profile_dir/merged-*.profdata)"
RUSTFLAGS="-Cprofile-use=$pref_data" cargo +nightly multivers


