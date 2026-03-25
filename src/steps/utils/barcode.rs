use crate::args::Args;
use crate::multi_aligner::AlignerConfig;

pub fn get_barcode_aligner_config(args: &Args) -> AlignerConfig {
    AlignerConfig {
        match_score: args.barcode_match,
        mismatch_score: args.barcode_mismatch,
        gap_open: args.barcode_gap_open,
        gap_extend: args.barcode_gap_extend,
        error_rate: args.barcode_error_rate,
        min_overlap: args.barcode_min_overlap,
        cycles: if args.barcode_keep { 1 } else { args.barcode_cycles },
    }
}
