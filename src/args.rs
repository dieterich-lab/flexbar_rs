use crate::types::*;
use clap::Parser;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::PathBuf;

#[derive(Parser, Debug, Serialize, Deserialize, Clone)]
#[serde(default)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// Number of threads to employ
    #[arg(short = 'n', long, default_value_t = 1)]
    pub threads: usize,

    /// Prefix for output file names or paths
    #[arg(short = 't', long, default_value = "flexbarOut")]
    pub target: String,

    /// Fasta/q file or stdin (-) with reads that may contain barcodes
    /// Files compressed with Gzip (.gz), bzip2 (.bz2), or zstd (.zst) are supported
    /// and will be piped through a suitable decompression tool automatically.
    #[arg(short = 'r', long)]
    pub reads: Option<PathBuf>,

    /// Second input file of paired reads
    #[arg(short = 'p', long)]
    pub reads2: Option<PathBuf>,

    /// Accept iupac symbols in reads and convert to N if not ATCG
    #[arg(short = 'I', long)]
    pub iupac: bool,

    // Barcode detection
    /// Fasta file with barcodes for demultiplexing, may contain N
    #[arg(short = 'b', long)]
    pub barcodes: Option<PathBuf>,

    /// Additional barcodes file for second read set in paired mode
    #[arg(long = "barcodes2", alias = "b2")]
    pub barcodes2: Option<PathBuf>,

    /// Fasta/q file containing separate barcode reads for detection
    #[arg(long = "barcode-reads", alias = "br")]
    pub barcode_reads: Option<PathBuf>,

    /// Minimum overlap of barcode and read
    #[arg(long = "barcode-min-overlap", alias = "bo")]
    pub barcode_min_overlap: Option<usize>,

    /// Error rate threshold for mismatches and gaps
    #[arg(long = "barcode-error-rate", alias = "be", default_value_t = 0.0)]
    pub barcode_error_rate: f32,

    /// Type of detection: ANY, LEFT, RIGHT, LTAIL, RTAIL
    #[arg(long = "barcode-trim-end", alias = "bt", default_value = "LTAIL")]
    pub barcode_trim_end: TrimEnd,

    /// Region size in tail trim-end modes
    #[arg(long = "barcode-tail-length", alias = "bn", default_value_t = 0)]
    pub barcode_tail_length: usize,

    /// Keep barcodes within reads instead of removal
    #[arg(long = "barcode-keep", alias = "bk")]
    pub barcode_keep: bool,

    /// Include unassigned reads in output generation
    #[arg(long = "barcode-unassigned", alias = "bu")]
    pub barcode_unassigned: bool,

    /// Alignment match score
    #[arg(long = "barcode-match", alias = "bm", default_value_t = 1)]
    pub barcode_match: i32,

    /// Alignment mismatch score
    #[arg(long = "barcode-mismatch", alias = "bi", default_value_t = -1)]
    pub barcode_mismatch: i32,

    /// Alignment gap score
    #[arg(long = "barcode-gap", alias = "bg", default_value_t = -9)]
    pub barcode_gap: i32,

    // Adapter removal
    /// Fasta file with adapters for removal that may contain N
    #[arg(short = 'a', long)]
    pub adapters: Option<PathBuf>,

    /// File with extra adapters for the second read set in paired mode
    #[arg(long = "adapters2", alias = "a2")]
    pub adapters2: Option<PathBuf>,

    /// Single adapter sequence as an alternative to adapters option
    #[arg(long = "adapter-seq", alias = "as")]
    pub adapter_seq: Option<String>,

    /// Adapter preset: TruSeq, SmallRNA, Methyl, Ribo, Nextera, NexteraMP
    #[arg(long = "adapter-preset", alias = "aa")]
    pub adapter_preset: Option<AdapterPreset>,

    /// Minimum overlap for removal without pair overlap
    #[arg(long = "adapter-min-overlap", alias = "ao", default_value_t = 3)]
    pub adapter_min_overlap: usize,

    /// Error rate threshold for mismatches and gaps
    #[arg(long = "adapter-error-rate", alias = "ae", default_value_t = 0.1)]
    pub adapter_error_rate: f32,

    /// Type of removal: ANY, LEFT, RIGHT, LTAIL, RTAIL
    #[arg(long = "adapter-trim-end", alias = "at", default_value = "RIGHT")]
    pub adapter_trim_end: TrimEnd,

    /// Region size for tail trim-end modes
    #[arg(long = "adapter-tail-length", alias = "an", default_value_t = 0)]
    pub adapter_tail_length: usize,

    /// Batch size for record processing
    #[arg(long = "batch-size", default_value_t = 256)]
    pub batch_size: usize,

    /// Skip restriction to pass read ends in right and left modes
    #[arg(long = "adapter-relaxed", alias = "ax")]
    pub adapter_relaxed: bool,

    /// Overlap detection of paired reads: OFF (default), ON, SHORT, ONLY
    #[arg(long = "adapter-pair-overlap", alias = "ap", default_value = "OFF")]
    pub adapter_pair_overlap: PairOverlapMode,

    /// Minimum overlap of paired reads for detection
    #[arg(long = "adapter-min-poverlap", alias = "av", default_value_t = 40)]
    pub adapter_min_poverlap: usize,

    /// Include reverse complements of adapters: OFF (default) ON, ONLY
    #[arg(long = "adapter-revcomp", alias = "ac", default_value = "OFF")]
    pub adapter_revcomp: RevCompMode,

    /// Use different trim-end for reverse complements of adapters
    #[arg(long = "adapter-revcomp-end", alias = "ad")]
    pub adapter_revcomp_end: Option<TrimEnd>,

    /// Add reverse complement of detected barcode to adapters
    #[arg(long = "adapter-add-barcode", alias = "ab")]
    pub adapter_add_barcode: bool,

    /// Consider only a single read set for adapters: 1, 2
    #[arg(long = "adapter-read-set", alias = "ar", default_value = "None")]
    pub adapter_read_set: AdapterReadset,

    /// Modify that trimmed reads are kept: OFF, ONLY
    #[arg(long = "adapter-trimmed-out", alias = "ak")]
    pub adapter_trimmed_out: Option<AdapterTrimmed>,

    /// Number of adapter removal cycles
    #[arg(long = "adapter-cycles", alias = "ay", default_value_t = 1)]
    pub adapter_cycles: usize,

    /// Alignment match score
    #[arg(long = "adapter-match", alias = "am", default_value_t = 1)]
    pub adapter_match: i32,

    /// Alignment mismatch score
    #[arg(long = "adapter-mismatch", alias = "ai", default_value_t = -1)]
    pub adapter_mismatch: i32,

    /// Alignment gap score
    #[arg(long = "adapter-gap", alias = "ag", default_value_t = -6)]
    pub adapter_gap: i32,

    // Filtering and trimming
    /// Allowed uncalled bases N for each read
    #[arg(short = 'u', long = "max-uncalled", default_value_t = 0)]
    pub max_uncalled: usize,

    /// Trim given number of bases on 5' read end before detection, default 0
    #[arg(short = 'x', long = "pre-trim-left", default_value_t = 0)]
    pub pre_trim_left: usize,

    /// Trim a specified number of bases on 3' end prior to detection, default 0
    #[arg(short = 'y', long = "pre-trim-right", default_value_t = 0)]
    pub pre_trim_right: usize,

    /// Trim to specified read length from 3' end after removal, default 0
    #[arg(short = 'k', long = "post-trim-length", default_value_t = 0)]
    pub post_trim_length: usize,

    /// Minimum read length to remain after removal, default 18
    #[arg(short = 'm', long = "min-read-length", default_value_t = 18)]
    pub min_read_length: usize,

    // Quality-based trimming
    /// Quality-based trimming mode: OFF (default), TAIL, WIN, BWA
    #[arg(short = 'q', long, default_value = "OFF")]
    pub qtrim: QualTrimType,

    /// Quality format: sanger (default), solexa, i1.3, i1.5, i1.8
    #[arg(long = "qtrim-format", alias = "qf", default_value = "sanger")]
    pub qtrim_format: QualityType,

    /// Minimum quality as threshold for trimming, default 20
    #[arg(long = "qtrim-threshold", alias = "qt", default_value_t = 20)]
    pub qtrim_threshold: usize,

    /// Region size for sliding window approach
    #[arg(long = "qtrim-win-size", alias = "qw", default_value_t = 5)]
    pub qtrim_win_size: usize,

    /// Perform quality-based trimming after removal steps
    #[arg(long = "qtrim-post-removal", alias = "qa")]
    pub qtrim_post_removal: bool,

    // Trimming of homopolymers
    /// Trim specific homopolymers on left read end after removal
    #[arg(long = "htrim-left", alias = "hl")]
    pub htrim_left: Option<String>,

    /// Trim certain homopolymers on right read end after removal
    #[arg(long = "htrim-right", alias = "hr")]
    pub htrim_right: Option<String>,

    /// Minimum length of homopolymers at read ends
    #[arg(long = "htrim-min-length", alias = "hi", default_value_t = 3)]
    pub htrim_min_length: usize,

    /// Minimum length for homopolymers specified after the first one, default 0
    #[arg(long = "htrim-min-length2", alias = "h2",  default_value_t = 0)]
    pub htrim_min_length2: usize,

    /// Maximum length of homopolymers on the left and right read end, default 0
    #[arg(long = "htrim-max-length", alias = "hx",  default_value_t = 0)]
    pub htrim_max_length: usize,

    /// Apply maximum length of homopolymers only for first one
    #[arg(long = "htrim-max-first", alias = "hf")]
    pub htrim_max_first: bool,

    /// Error rate threshold for mismatches
    #[arg(long = "htrim-error-rate", alias = "he", default_value_t = 0.1)]
    pub htrim_error_rate: f32,

    /// Trim only in case of adapter removal on same side
    #[arg(long = "htrim-adapter", alias = "ha")]
    pub htrim_adapter: bool,

    // Output selection
    /// Prefer non-quality format fasta for output
    #[arg(short = 'f', long)]
    pub fasta_output: bool,

    /// Direct compression of output files: Uncompressed, GZ, BZ2, ZSTD, KEEP (default,
    /// use same compression as input)
    #[arg(short = 'z', long, default_value = "KEEP")]
    pub output_compression: Option<CompressionType>,

    /// Write reads to stdout, tagged and interleaved if needed
    #[arg(long = "stdout-reads")]
    pub stdout_reads: bool,

    /// Output file for reads instead of target prefix usage
    #[arg(short = 'R', long = "output-reads")]
    pub output_reads: Option<PathBuf>,

    /// Output file for reads2 instead of target prefix usage
    #[arg(short = 'P', long = "output-reads2")]
    pub output_reads2: Option<PathBuf>,

    /// Generate length distribution for read output files
    #[arg(short = 'j', long = "length-dist")]
    pub length_dist: bool,

    /// Write single reads for too short counterparts in pairs
    #[arg(short = 's', long = "single-reads")]
    pub single_reads: bool,

    /// Write paired single reads with N for short counterparts
    #[arg(short = 'S', long = "single-reads-paired")]
    pub single_reads_paired: bool,

    // Logging and tagging
    /// Print chosen read alignments: ALL, MOD, TAB
    #[arg(short = 'l', long)]
    pub align_log: Option<LogAlign>,

    /// Write statistics to stdout instead of target log file
    #[arg(short = 'o', long = "stdout-log")]
    pub stdout_log: bool,

    /// Output file for logging instead of target prefix usage
    #[arg(short = 'O', long = "output-log")]
    pub output_log: Option<PathBuf>,

    /// Tag reads that are subject to adapter or barcode removal
    #[arg(short = 'g', long = "removal-tags")]
    pub removal_tags: bool,

    /// Replace read tags by ascending number to save space
    #[arg(short = 'e', long = "number-tags")]
    pub number_tags: bool,

    /// Capture UMIs in reads at barcode or adapter N positions
    #[arg(short = 'd', long = "umi-tags")]
    pub umi_tags: bool,

    /// Print debug information
    #[arg(long)]
    pub debug: bool,

    /// Path to a parameter file in TOML format
    #[arg(long = "parameter-file")]
    pub parameter_file: Option<PathBuf>,
}

impl Default for Args {
    fn default() -> Self {
        Self::parse_from(["flexbar"])
    }
}

impl Args {
    pub fn new() -> Self {
        Self::parse()
    }

    pub fn get_args() -> Self {
        let mut args = Self::new();

        if let Some(ref param_file) = args.parameter_file {
            let file = fs::File::open(param_file).expect("Could not open parameter file");
            let file_args = Self::load_parameter_file(file);
            args.merge_from_file(file_args);
        }
        Self::validate_args(&args);

        args
    }

    pub fn load_parameter_file<R: std::io::Read>(mut reader: R) -> Self {
        let mut content = String::new();
        reader.read_to_string(&mut content).expect("Could not read from parameter source");
        toml::from_str(&content).expect("Could not parse parameter file as TOML")
    }

    fn merge_from_file(&mut self, file_args: Self) {
        use clap::CommandFactory;
        let matches = Self::command().get_matches();
        self.merge_with_matches(file_args, &matches);
    }

    fn merge_with_matches(&mut self, file_args: Self, matches: &clap::ArgMatches) {
        // Helper to check if an argument was provided via CLI
        let is_present = |name: &str| {
            matches.value_source(name) == Some(clap::parser::ValueSource::CommandLine)
        };

        // We only want to overwrite values that were NOT provided via CLI.
        // clap's `get_matches()` knows what was actually on the command line.

        if !is_present("threads") { self.threads = file_args.threads; }
        if !is_present("target") { self.target = file_args.target; }
        if !is_present("reads") { self.reads = file_args.reads; }
        if !is_present("reads2") { self.reads2 = file_args.reads2; }
        if !is_present("iupac") { self.iupac = file_args.iupac; }

        // Barcode detection
        if !is_present("barcodes") { self.barcodes = file_args.barcodes; }
        if !is_present("barcodes2") { self.barcodes2 = file_args.barcodes2; }
        if !is_present("barcode_reads") { self.barcode_reads = file_args.barcode_reads; }
        if !is_present("barcode_min_overlap") { self.barcode_min_overlap = file_args.barcode_min_overlap; }
        if !is_present("barcode_error_rate") { self.barcode_error_rate = file_args.barcode_error_rate; }
        if !is_present("barcode_trim_end") { self.barcode_trim_end = file_args.barcode_trim_end; }
        if !is_present("barcode_tail_length") { self.barcode_tail_length = file_args.barcode_tail_length; }
        if !is_present("barcode_keep") { self.barcode_keep = file_args.barcode_keep; }
        if !is_present("barcode_unassigned") { self.barcode_unassigned = file_args.barcode_unassigned; }
        if !is_present("barcode_match") { self.barcode_match = file_args.barcode_match; }
        if !is_present("barcode_mismatch") { self.barcode_mismatch = file_args.barcode_mismatch; }
        if !is_present("barcode_gap") { self.barcode_gap = file_args.barcode_gap; }

        // Adapter removal
        if !is_present("adapters") { self.adapters = file_args.adapters; }
        if !is_present("adapters2") { self.adapters2 = file_args.adapters2; }
        if !is_present("adapter_seq") { self.adapter_seq = file_args.adapter_seq; }
        if !is_present("adapter_preset") { self.adapter_preset = file_args.adapter_preset; }
        if !is_present("adapter_min_overlap") { self.adapter_min_overlap = file_args.adapter_min_overlap; }
        if !is_present("adapter_error_rate") { self.adapter_error_rate = file_args.adapter_error_rate; }
        if !is_present("adapter_trim_end") { self.adapter_trim_end = file_args.adapter_trim_end; }
        if !is_present("adapter_tail_length") { self.adapter_tail_length = file_args.adapter_tail_length; }
        if !is_present("adapter_relaxed") { self.adapter_relaxed = file_args.adapter_relaxed; }
        if !is_present("adapter_pair_overlap") { self.adapter_pair_overlap = file_args.adapter_pair_overlap; }
        if !is_present("adapter_min_poverlap") { self.adapter_min_poverlap = file_args.adapter_min_poverlap; }
        if !is_present("adapter_revcomp") { self.adapter_revcomp = file_args.adapter_revcomp; }
        if !is_present("adapter_revcomp_end") { self.adapter_revcomp_end = file_args.adapter_revcomp_end; }
        if !is_present("adapter_add_barcode") { self.adapter_add_barcode = file_args.adapter_add_barcode; }
        if !is_present("adapter_read_set") { self.adapter_read_set = file_args.adapter_read_set; }
        if !is_present("adapter_trimmed_out") { self.adapter_trimmed_out = file_args.adapter_trimmed_out; }
        if !is_present("adapter_cycles") { self.adapter_cycles = file_args.adapter_cycles; }
        if !is_present("adapter_match") { self.adapter_match = file_args.adapter_match; }
        if !is_present("adapter_mismatch") { self.adapter_mismatch = file_args.adapter_mismatch; }
        if !is_present("adapter_gap") { self.adapter_gap = file_args.adapter_gap; }

        // Filtering and trimming
        if !is_present("max_uncalled") { self.max_uncalled = file_args.max_uncalled; }
        if !is_present("pre_trim_left") { self.pre_trim_left = file_args.pre_trim_left; }
        if !is_present("pre_trim_right") { self.pre_trim_right = file_args.pre_trim_right; }
        if !is_present("post_trim_length") { self.post_trim_length = file_args.post_trim_length; }
        if !is_present("min_read_length") { self.min_read_length = file_args.min_read_length; }

        // Quality-based trimming
        if !is_present("qtrim") { self.qtrim = file_args.qtrim; }
        if !is_present("qtrim_format") { self.qtrim_format = file_args.qtrim_format; }
        if !is_present("qtrim_threshold") { self.qtrim_threshold = file_args.qtrim_threshold; }
        if !is_present("qtrim_win_size") { self.qtrim_win_size = file_args.qtrim_win_size; }
        if !is_present("qtrim_post_removal") { self.qtrim_post_removal = file_args.qtrim_post_removal; }

        // Trimming of homopolymers
        if !is_present("htrim_left") { self.htrim_left = file_args.htrim_left; }
        if !is_present("htrim_right") { self.htrim_right = file_args.htrim_right; }
        if !is_present("htrim_min_length") { self.htrim_min_length = file_args.htrim_min_length; }
        if !is_present("htrim_min_length2") { self.htrim_min_length2 = file_args.htrim_min_length2; }
        if !is_present("htrim_max_length") { self.htrim_max_length = file_args.htrim_max_length; }
        if !is_present("htrim_max_first") { self.htrim_max_first = file_args.htrim_max_first; }
        if !is_present("htrim_error_rate") { self.htrim_error_rate = file_args.htrim_error_rate; }
        if !is_present("htrim_adapter") { self.htrim_adapter = file_args.htrim_adapter; }

        // Output selection
        if !is_present("fasta_output") { self.fasta_output = file_args.fasta_output; }
        if !is_present("output_compression") { self.output_compression = file_args.output_compression; }
        if !is_present("stdout_reads") { self.stdout_reads = file_args.stdout_reads; }
        if !is_present("output_reads") { self.output_reads = file_args.output_reads; }
        if !is_present("output_reads2") { self.output_reads2 = file_args.output_reads2; }
        if !is_present("length_dist") { self.length_dist = file_args.length_dist; }
        if !is_present("single_reads") { self.single_reads = file_args.single_reads; }
        if !is_present("single_reads_paired") { self.single_reads_paired = file_args.single_reads_paired; }

        // Logging and tagging
        if !is_present("align_log") { self.align_log = file_args.align_log; }
        if !is_present("stdout_log") { self.stdout_log = file_args.stdout_log; }
        if !is_present("output_log") { self.output_log = file_args.output_log; }
        if !is_present("removal_tags") { self.removal_tags = file_args.removal_tags; }
        if !is_present("number_tags") { self.number_tags = file_args.number_tags; }
        if !is_present("umi_tags") { self.umi_tags = file_args.umi_tags; }
        if !is_present("debug") { self.debug = file_args.debug; }
    }

    fn validate_args(&self) {
        let count_adapters_specified =
            if self.adapters.is_some() { 1 } else { 0 }
            + if self.adapter_seq.is_some() { 1 } else { 0 }
            + if self.adapter_preset.is_some() { 1 } else { 0 };
        if count_adapters_specified > 1 {
            panic!("Only one of --adapters, --adapter-seq or --adapter-preset can be specified.");
        }

        if self.adapter_cycles < 1 {
            panic!("--adapter-cycles must be greater than 0.");
        }

        let is_paired = self.reads2.is_some();
        match self.adapter_preset {
            Some(AdapterPreset::NexteraMp) => {
                if !is_paired {
                    panic!("--adapter-preset NexteraMp requires paired reads.");
                }
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_from_file() {
        use clap::CommandFactory;
        
        // CLI arguments: --threads 4 --target cli_target
        let args_vec = vec!["flexbar", "--threads", "4", "--target", "cli_target"];
        let matches = Args::command().get_matches_from(args_vec.clone());
        let mut cli_args = Args::parse_from(args_vec);
        
        let file_args = Args {
            threads: 8,                       // Should NOT overwrite CLI
            target: "file_target".to_string(), // Should NOT overwrite CLI
            ..Args::parse_from(["flexbar"])   // defaults for others
        };

        cli_args.merge_with_matches(file_args, &matches);

        assert_eq!(cli_args.threads, 4);
        assert_eq!(cli_args.target, "cli_target");
    }

    #[test]
    fn test_load_parameter_file() {
        use std::io::Cursor;
        let content = "threads = 10\ntarget = 'toml_target'\n";
        let reader = Cursor::new(content);
        
        let loaded_args = Args::load_parameter_file(reader);
        
        assert_eq!(loaded_args.threads, 10);
        assert_eq!(loaded_args.target, "toml_target");
    }
}
