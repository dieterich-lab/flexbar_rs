use log::debug;
use crate::args::Args;
use crate::types::*;

#[derive(Debug, Clone)]
pub struct Config {
    pub is_paired: bool,
    pub use_adapter_file: bool,
    pub log_stdout: bool,
    pub use_stdin: bool,
    pub use_rc_trim_end: bool,
    pub add_barcode_adapter: bool,
    pub interleaved_input: bool,
    pub cut_len_begin: usize,
    pub cut_len_end: usize,
    pub cut_len_read: usize,
    pub a_tail_len: usize,
    pub b_tail_len: usize,
    pub qtrim_thresh: usize,
    pub htrim_min_length2: usize,
    pub htrim_max_length: usize,
    pub n_bundles: usize,
    pub a_cycles: usize,

    pub arc_end: TrimEnd,
    pub format: FileFormat,
    pub qual: QualityType,
    pub q_trim: QualTrimType,
    pub log_align: LogAlign,
    pub cmprs_type: CompressionType,
    pub run_type: RunType,
    pub bar_detect: BarcodeDetect,
    pub adap_rm: AdapterRemoval,
    pub rc_mode: RevCompMode,
    pub po_mode: PairOverlap,
    pub a_preset: AdapterPreset,
    pub a_trimmed: AdapterTrimmed,

    pub args: Args,
}

impl Config {
    pub fn new(args: Args) -> Self {
        let is_paired = args.reads2.is_some() || args.interleaved;
        let interleaved_input = args.interleaved;
        let (run_type, bar_detect) = Self::determine_run_and_barcode_type(&args, is_paired);
        
        let use_stdin = args.reads.as_ref().map_or(false, |p| p.to_string_lossy() == "-");
        let format = Self::determine_file_format(&args);

        let (adap_rm, use_adapter_file, a_preset) = Self::determine_adapter_settings(&args, is_paired);
        let (q_trim, qual, qtrim_thresh) = Self::determine_quality_settings(&args, format);

        let cmprs_type = args.zip_output.unwrap_or(CompressionType::Uncompressed);
        let log_align = args.align_log.unwrap_or(LogAlign::None);

        let po_mode = Self::determine_pair_overlap_mode(&args, is_paired, adap_rm);

        let (rc_mode, arc_end, use_rc_trim_end) = Self::determine_revcomp_settings(&args, args.adapter_trim_end);

        let add_barcode_adapter = args.adapter_add_barcode && is_paired && args.adapter_trim_end == TrimEnd::Right && rc_mode != RevCompMode::On &&
                                  bar_detect != BarcodeDetect::BarcodeRead && bar_detect != BarcodeDetect::Off && args.barcode_trim_end == TrimEnd::LTail;

        let a_trimmed = args.adapter_trimmed_out.unwrap_or(AdapterTrimmed::On);

        let final_adap_rm = if is_paired && adap_rm != AdapterRemoval::Normal2 {
            args.adapter_read_set.unwrap_or(adap_rm)
        } else {
            adap_rm
        };

        let mut a_cycles = args.adapter_cycles;
        if a_preset == AdapterPreset::NexteraMp && a_cycles < 3 {
            a_cycles = 3;
        }

        let config = Config {
            args: args.clone(),
            is_paired,
            use_adapter_file,
            log_stdout: args.stdout_log && !args.stdout_reads,
            use_stdin,
            use_rc_trim_end,
            add_barcode_adapter,
            interleaved_input,
            cut_len_begin: args.pre_trim_left.unwrap_or(0),
            cut_len_end: args.pre_trim_right.unwrap_or(0),
            cut_len_read: args.post_trim_length.unwrap_or(0),
            a_tail_len: args.adapter_tail_length.unwrap_or(0),
            b_tail_len: args.barcode_tail_length.unwrap_or(0),
            qtrim_thresh,
            htrim_min_length2: args.htrim_min_length2.unwrap_or(0),
            htrim_max_length: args.htrim_max_length.unwrap_or(0),
            n_bundles: args.bundles.unwrap_or(0),
            a_cycles,
            arc_end,
            format,
            qual,
            q_trim,
            log_align,
            cmprs_type,
            run_type,
            bar_detect,
            adap_rm: final_adap_rm,
            rc_mode,
            po_mode,
            a_preset,
            a_trimmed,
        };
        debug!("{:#?}", config);
        config
    }

    fn determine_run_and_barcode_type(args: &Args, is_paired: bool) -> (RunType, BarcodeDetect) {
        let mut run_type = if is_paired { RunType::Paired } else { RunType::Single };
        let mut bar_detect = BarcodeDetect::Off;

        if args.barcodes.is_some() {
            bar_detect = if args.barcode_reads.is_some() {
                BarcodeDetect::BarcodeRead
            } else if args.barcode_keep {
                BarcodeDetect::WithinRead
            } else {
                BarcodeDetect::WithinReadRemoval
            };

            run_type = match run_type {
                RunType::Single => RunType::SingleBarcoded,
                RunType::Paired => RunType::PairedBarcoded,
                _ => run_type,
            };

            if args.barcodes2.is_some() && bar_detect != BarcodeDetect::BarcodeRead && is_paired {
                bar_detect = match bar_detect {
                    BarcodeDetect::WithinReadRemoval => BarcodeDetect::WithinReadRemoval2,
                    BarcodeDetect::WithinRead => BarcodeDetect::WithinRead2,
                    _ => bar_detect,
                };
            }
        }
        (run_type, bar_detect)
    }

    fn determine_file_format(args: &Args) -> FileFormat {
        args.reads.as_ref().map_or(FileFormat::Fasta, |path| {
            let ext = path.extension().and_then(|s| s.to_str()).map(|s| s.to_lowercase());
            if let Some(e) = ext {
                if e == "fasta" || e == "fa" {
                    return FileFormat::Fasta;
                }
            }
            FileFormat::Fastq
        })
    }

    fn determine_adapter_settings(args: &Args, is_paired: bool) -> (AdapterRemoval, bool, AdapterPreset) {
        let mut adap_rm = AdapterRemoval::Off;
        let mut use_adapter_file = false;
        let mut a_preset = AdapterPreset::Off;

        if args.adapters.is_some() {
            adap_rm = AdapterRemoval::Normal;
            use_adapter_file = true;
        } else if args.adapter_seq.is_some() {
            adap_rm = AdapterRemoval::Normal;
        }

        if args.adapters2.is_some() && adap_rm == AdapterRemoval::Normal && is_paired {
            adap_rm = AdapterRemoval::Normal2;
        }

        if let Some(preset) = args.adapter_preset {
            a_preset = preset;
            adap_rm = AdapterRemoval::Normal;
            if is_paired && matches!(a_preset, AdapterPreset::TruSeq | AdapterPreset::Methyl | AdapterPreset::NexteraMp) {
                adap_rm = AdapterRemoval::Normal2;
            }
        }
        (adap_rm, use_adapter_file, a_preset)
    }

    fn determine_quality_settings(args: &Args, format: FileFormat) -> (QualTrimType, QualityType, usize) {
        let mut q_trim = QualTrimType::Off;
        let mut qual = QualityType::Sanger;
        let mut qtrim_thresh = args.qtrim_threshold;

        if let (Some(qt), FileFormat::Fastq) = (args.qtrim, format) {
            q_trim = qt;

            if let Some(qf) = args.qtrim_format {
                qual = qf;
            }

            qtrim_thresh += match qual {
                QualityType::Sanger => 33,
                QualityType::Solexa => 59,
                QualityType::Illumina => 64,
            };
        }
        (q_trim, qual, qtrim_thresh)
    }

    fn determine_pair_overlap_mode(args: &Args, is_paired: bool, adap_rm: AdapterRemoval) -> PairOverlap {
        if !is_paired {
            return PairOverlap::Off;
        }
        if let Some(po_mode) = args.adapter_pair_overlap {
            if adap_rm == AdapterRemoval::Off && matches!(po_mode, PairOverlap::On | PairOverlap::Short) {
                PairOverlap::Off
            } else {
                po_mode
            }
        } else {
            PairOverlap::Off
        }
    }

    fn determine_revcomp_settings(args: &Args, a_end: TrimEnd) -> (RevCompMode, TrimEnd, bool) {
        let mut rc_mode = RevCompMode::Off;
        let mut arc_end = TrimEnd::Right;
        let mut use_rc_trim_end = false;

        if let Some(rc) = args.adapter_revcomp {
            rc_mode = rc;

            if let (RevCompMode::On, Some(rce)) = (rc_mode, args.adapter_revcomp_end) {
                arc_end = rce;
                if arc_end != a_end {
                    use_rc_trim_end = true;
                }
            }
        }
        (rc_mode, arc_end, use_rc_trim_end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::path::PathBuf;

    #[test]
    fn test_config_basic() {
        let args = Args {
            threads: 4,
            target: "test_target".to_string(),
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert_eq!(config.args.threads, 4);
        assert_eq!(config.is_paired, false);
        assert_eq!(config.run_type, RunType::Single);
    }

    #[test]
    fn test_config_paired() {
        let args = Args {
            reads: Some(PathBuf::from("reads1.fastq")),
            reads2: Some(PathBuf::from("reads2.fastq")),
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert!(config.is_paired);
        assert_eq!(config.run_type, RunType::Paired);
    }

    #[test]
    fn test_config_qtrim_thresh() {
        let args = Args {
            reads: Some(PathBuf::from("reads.fastq")),
            qtrim: Some(QualTrimType::Tail),
            qtrim_format: Some(QualityType::Sanger),
            qtrim_threshold: 20,
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert_eq!(config.q_trim, QualTrimType::Tail);
        assert_eq!(config.qtrim_thresh, 20 + 33);
    }

    #[test]
    fn test_config_adapter_preset() {
        let args = Args {
            adapter_preset: Some(AdapterPreset::TruSeq),
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert_eq!(config.a_preset, AdapterPreset::TruSeq);
        assert_eq!(config.adap_rm, AdapterRemoval::Normal);
    }

    #[test]
    fn test_config_adapter_preset_paired() {
        let args = Args {
            reads2: Some(PathBuf::from("r2.fq")),
            adapter_preset: Some(AdapterPreset::TruSeq),
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert_eq!(config.a_preset, AdapterPreset::TruSeq);
        assert_eq!(config.adap_rm, AdapterRemoval::Normal2);
    }
}
