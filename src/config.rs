use log::{debug, info, warn};
use crate::args::Args;
use crate::types::*;

#[derive(Debug, Clone)]
pub struct Config {
    pub log_stdout: bool,
    pub use_rc_trim_end: bool,
    pub corrected_qtrim_threshold: usize,

    pub format: FileFormat,
    pub run_type: RunType,
    pub bar_detect: BarcodeDetect,
    pub adapter_removal: AdapterRemoval,
    pub a_trimmed: AdapterTrimmed,

    pub batch_size: usize,

    pub args: Args,
}

impl Config {
    pub fn new(args: Args) -> Self {
        let has_adapters = Self::determine_has_adapters(&args);
        let format = Self::determine_file_format(&args);
        let bar_detect = Self::determine_barcode_detection(&args);
        let a_trimmed = args.adapter_trimmed_out.unwrap_or(AdapterTrimmed::On);

        let raw_config = Config {
            args: args.clone(),
            log_stdout: args.stdout_log && !args.stdout_reads,
            use_rc_trim_end: Self::determine_use_rc_trim_end(&args),
            corrected_qtrim_threshold: Self::determine_qtrim_threshold(&args, format),
            format,
            run_type: Self::determine_run_type(&args),
            bar_detect,
            adapter_removal: Self::determine_adapter_removal(&args, has_adapters),
            a_trimmed,
            batch_size: args.batch_size,
        };
        let config = Self::sanitize_config(raw_config, has_adapters);
        debug!("{:#?}", config);
        config
    }

    fn determine_has_adapters(args: &Args) -> bool {
        args.adapter_preset.is_some() || args.adapters.is_some() || args.adapter_seq.is_some()
    }

    fn determine_run_type(args: &Args) -> RunType {
        match args.barcodes {
            Some(_) => if args.reads2.is_some() { RunType::PairedBarcoded } 
                        else { RunType::SingleBarcoded },
            None => if args.reads2.is_some() { RunType::Paired } else { RunType::Single },
        }
    }

    fn determine_barcode_detection(args: &Args) -> BarcodeDetect {
        match args.barcodes {
            Some(_) => {
                let raw = if args.barcode_reads.is_some() {
                    BarcodeDetect::BarcodeRead
                } else if args.barcode_keep {
                    BarcodeDetect::WithinRead
                } else {
                    BarcodeDetect::WithinReadRemoval
                };
                if args.reads2.is_some() && args.barcodes2.is_some() {
                    match raw {
                        BarcodeDetect::WithinReadRemoval => BarcodeDetect::WithinReadRemoval2,
                        BarcodeDetect::WithinRead => BarcodeDetect::WithinRead2,
                        _ => raw,
                    }
                } else {
                    raw
                }
            },
            None => BarcodeDetect::Off,
        }
    }

    fn determine_use_rc_trim_end(args: &Args) -> bool {
        match args.adapter_revcomp_end {
            Some(x) => {
                match args.adapter_revcomp {
                    RevCompMode::On => x != args.adapter_trim_end,
                    _ => false,
                }
            },
            _ => false,
        }
    }

    fn determine_file_format(args: &Args) -> FileFormat {
        args.reads.as_ref().map_or(FileFormat::Fasta, |path| {
            if path.to_string_lossy() == "-" {
                return FileFormat::Fastq;
            }
            let ext = path.extension().and_then(|s| s.to_str()).map(|s| s.to_lowercase());
            if let Some(e) = ext {
                if e == "fasta" || e == "fa" {
                    return FileFormat::Fasta;
                }
            }
            FileFormat::Fastq
        })
    }

    fn determine_adapter_removal(args: &Args, has_adapters: bool) -> AdapterRemoval {
        if !has_adapters {
            return AdapterRemoval::Off;
        }
        let paired_adapter_removal = Self::is_paired_adapter_removal(&args);
        match args.adapter_read_set {
            AdapterReadset::One => {
                if paired_adapter_removal {AdapterRemoval::Normal2} else {AdapterRemoval::One}
            },
            AdapterReadset::Two => {
                if paired_adapter_removal {AdapterRemoval::Normal2} else {AdapterRemoval::Two}
            }
            AdapterReadset::None => {
                if paired_adapter_removal {AdapterRemoval::Normal2} else {AdapterRemoval::Normal}
            }
        }
    }

    fn is_paired_adapter_removal(args: &Args) -> bool {
        if args.reads2.is_none() { return false; }
        match args.adapter_preset {
            Some(AdapterPreset::TruSeq) => return true,
            Some(AdapterPreset::NexteraMp) => return true,
            Some(AdapterPreset::Methyl) => return true,
            Some(_) => return false,
            _ => ()
        }
        if args.adapters2.is_some() && args.reads2.is_some() {
            return true;
        }
        false
    }

    fn determine_qtrim_threshold(args: &Args, format: FileFormat) -> usize {
        if format == FileFormat::Fastq {
            args.qtrim_threshold + match args.qtrim_format {
                QualityType::Sanger => 33,
                QualityType::Solexa => 59,
                QualityType::Illumina => 64,
            }
        } else {
            args.qtrim_threshold
        }
    }

    fn sanitize_config(raw_config: Config, has_adapters: bool) -> Config {
        let mut config = raw_config.clone();

        if config.batch_size == 0 {
            config.batch_size = 256;
        }

        let is_paired = raw_config.args.reads2.is_some();

        match raw_config.args.adapter_preset {
            Some(AdapterPreset::NexteraMp) => {
                if config.args.adapter_cycles < 3 {
                    info!("Setting adapter_cycles to 3 (minimum for NexteraMP preset).");
                    config.args.adapter_cycles = 3;
                }
            }
            _ => ()
        }

        if !is_paired && raw_config.args.adapter_pair_overlap != PairOverlapMode::Off {
            warn!("Adapter pair overlap mode is not supported for unpaired reads - disabling");
            config.args.adapter_pair_overlap = PairOverlapMode::Off;
        }
        if !has_adapters && matches!(
            raw_config.args.adapter_pair_overlap, PairOverlapMode::On | PairOverlapMode::Short
        ) {
            warn!("Adapter pair overlap mode is not supported for no adapters - disabling");
            config.args.adapter_pair_overlap = PairOverlapMode::Off;
        }
        if raw_config.args.adapter_add_barcode {
            if !is_paired {
                warn!("Adapter barcode addition is not supported for unpaired reads - disabling");
                config.args.adapter_add_barcode = false;
            }
            if raw_config.args.adapter_trim_end != TrimEnd::Right {
                warn!("Adapter barcode addition is only supported for right trim end - disabling");
                config.args.adapter_add_barcode = false;
            }
            if raw_config.args.adapter_revcomp == RevCompMode::On {
                warn!("Adapter barcode addition is not supported for reverse complement - disabling");
                config.args.adapter_add_barcode = false;
            }
            if raw_config.bar_detect == BarcodeDetect::BarcodeRead
                || raw_config.bar_detect == BarcodeDetect::Off {
                warn!("Adapter barcode addition is not supported with this settings - disabling");
                config.args.adapter_add_barcode = false;
            }
        }
        config
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
        assert_eq!(config.run_type, RunType::Paired);
    }

    #[test]
    fn test_config_qtrim_thresh() {
        let args = Args {
            reads: Some(PathBuf::from("reads.fastq")),
            qtrim: QualTrimType::Tail,
            qtrim_format: QualityType::Sanger,
            qtrim_threshold: 20,
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert_eq!(config.corrected_qtrim_threshold, 20 + 33);
    }

    #[test]
    fn test_config_adapter_preset() {
        let args = Args {
            adapter_preset: Some(AdapterPreset::TruSeq),
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert_eq!(config.adapter_removal, AdapterRemoval::Normal);
    }

    #[test]
    fn test_config_adapter_preset_paired() {
        let args = Args {
            reads2: Some(PathBuf::from("r2.fq")),
            adapter_preset: Some(AdapterPreset::TruSeq),
            ..Args::parse_from(["flexbar"])
        };
        let config = Config::new(args);
        assert_eq!(config.adapter_removal, AdapterRemoval::Normal2);
    }
}
