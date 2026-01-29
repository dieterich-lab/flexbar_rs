use crate::args::Args;
use crate::types::{AdapterPreset, AdapterReadset, AdapterRemoval, BarcodeDetect, FileFormat, PairOverlapMode, QualTrimType, QualityType, RevCompMode, RunType, TrimEnd};
use crate::utils::io::get_file_type;
use bio::io::fasta;
use bio::alphabets::dna::revcomp;
use log::{debug, info, warn};
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use crate::logging::log_fatal;


#[derive(Clone)]
pub struct SequenceWithId {
    pub id: String,
    pub sequence: String,
}

pub struct Context {
    pub adapters: Vec<String>,
    pub adapters2: Vec<String>,
    pub barcodes: Vec<SequenceWithId>,
    pub barcodes2: Vec<SequenceWithId>,
    pub use_rc_trim_end: bool,
    pub corrected_qtrim_threshold: usize,
    pub format: FileFormat,
    pub run_type: RunType,
    pub bar_detect: BarcodeDetect,
    pub adapter_removal: AdapterRemoval,
}

impl Context {
    pub fn new(args: &Args) -> Self {
        let has_adapters = Self::determine_has_adapters(args);
        let format = get_file_type(args.reads.as_ref());
        let bar_detect = Self::determine_barcode_detection(args);

        let mut context = Context {
            adapters: Vec::new(),
            adapters2: Vec::new(),
            barcodes: Vec::new(),
            barcodes2: Vec::new(),
            use_rc_trim_end: Self::determine_use_rc_trim_end(args),
            corrected_qtrim_threshold: Self::determine_qtrim_threshold(args, format),
            format,
            run_type: Self::determine_run_type(args),
            bar_detect,
            adapter_removal: Self::determine_adapter_removal(args, has_adapters),
        };
        context.sanitize(args, has_adapters);
        debug!("{:#?}", context.run_type);
        context.load(args);
        context
    }

    fn determine_has_adapters(args: &Args) -> bool {
        args.adapter_preset.is_some() || args.adapters.is_some() || args.adapter_seq.is_some()
    }

    fn determine_run_type(args: &Args) -> RunType {
        match args.barcodes {
            Some(_) => {
                if args.reads2.is_some() {
                    RunType::PairedBarcoded
                } else {
                    RunType::SingleBarcoded
                }
            }
            None => {
                if args.reads2.is_some() {
                    RunType::Paired
                } else {
                    RunType::Single
                }
            }
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
            }
            None => BarcodeDetect::Off,
        }
    }

    fn determine_use_rc_trim_end(args: &Args) -> bool {
        match args.adapter_revcomp_end {
            Some(x) => match args.adapter_revcomp {
                RevCompMode::On => x != args.adapter_trim_end,
                _ => false,
            },
            _ => false,
        }
    }

    fn determine_adapter_removal(args: &Args, has_adapters: bool) -> AdapterRemoval {
        if !has_adapters {
            return AdapterRemoval::Off;
        }
        let paired_adapter_removal = Self::is_paired_adapter_removal(args);
        match args.adapter_read_set {
            AdapterReadset::One => {
                if paired_adapter_removal {
                    AdapterRemoval::Normal2
                } else {
                    AdapterRemoval::One
                }
            }
            AdapterReadset::Two => {
                if paired_adapter_removal {
                    AdapterRemoval::Normal2
                } else {
                    AdapterRemoval::Two
                }
            }
            AdapterReadset::None => {
                if paired_adapter_removal {
                    AdapterRemoval::Normal2
                } else {
                    AdapterRemoval::Normal
                }
            }
        }
    }

    fn is_paired_adapter_removal(args: &Args) -> bool {
        if args.reads2.is_none() {
            return false;
        }
        match args.adapter_preset {
            Some(AdapterPreset::TruSeq) => return true,
            Some(AdapterPreset::NexteraMp) => return true,
            Some(AdapterPreset::Methyl) => return true,
            Some(_) => return false,
            _ => (),
        }
        if args.adapters2.is_some() && args.reads2.is_some() {
            return true;
        }
        false
    }

    fn determine_qtrim_threshold(args: &Args, format: FileFormat) -> usize {
        if format == FileFormat::Fastq {
            args.qtrim_threshold
                + match args.qtrim_format {
                    QualityType::Sanger => 33,
                    QualityType::Solexa => 59,
                    QualityType::Illumina => 64,
                }
        } else {
            args.qtrim_threshold
        }
    }

    fn sanitize(&mut self, args: &Args, has_adapters: bool) {
        let is_paired = args.reads2.is_some();
        match args.adapter_preset {
            Some(AdapterPreset::NexteraMp) => {
                if args.adapter_cycles < 3 {
                    info!("Setting adapter_cycles to 3 (minimum for NexteraMP preset).");
                    // Note: args is immutable here. In config.rs it was mutable.
                }
            }
            _ => ()
        }

        if !is_paired && args.adapter_pair_overlap != PairOverlapMode::Off {
            warn!("Adapter pair overlap mode is not supported for unpaired reads - disabling");
        }
        if !has_adapters && matches!(
            args.adapter_pair_overlap, PairOverlapMode::On | PairOverlapMode::Short
        ) {
            warn!("Adapter pair overlap mode is not supported for no adapters - disabling");
        }
        if args.adapter_add_barcode {
            if !is_paired {
                warn!("Adapter barcode addition is not supported for unpaired reads - disabling");
            }
            if args.adapter_trim_end != TrimEnd::Right {
                warn!("Adapter barcode addition is only supported for right trim end - disabling");
            }
            if args.adapter_revcomp == RevCompMode::On {
                warn!("Adapter barcode addition is not supported for reverse complement - disabling");
            }
            if self.bar_detect == BarcodeDetect::BarcodeRead
                || self.bar_detect == BarcodeDetect::Off {
                warn!("Adapter barcode addition is not supported with this settings - disabling");
            }
        }
    }

    pub fn load(&mut self, args: &Args) {
        let adapter_removal = self.adapter_removal;
        if let Some(preset) = &args.adapter_preset {
            if matches!(preset, AdapterPreset::Off) {
                return;
            }
            self.load_preset(preset);
        }

        if let Some(adapter_seq) = &args.adapter_seq {
            self.adapters.push(adapter_seq.clone());
        }

        if let Some(adapter_file) = &args.adapters {
            let mut new_items = self.load_fasta_file(adapter_file)
                .into_iter()
                .map(|item| item.sequence)
                .collect();
            self.adapters.append(&mut new_items);
        }

        if adapter_removal == AdapterRemoval::Normal2 {
            if let Some(adapter_file2) = &args.adapters2 {
                self.adapters2 = self.load_fasta_file(adapter_file2)
                    .into_iter()
                    .map(|item| item.sequence)
                    .collect();
            }
        }

        if let Some(barcode_file) = &args.barcodes {
            self.barcodes = self.load_fasta_file(barcode_file);
        }

        if let Some(barcode_file2) = &args.barcodes2 {
            self.barcodes2 = self.load_fasta_file(barcode_file2);
        }

        self.apply_revcomp(args.adapter_revcomp);
        
        for (i, adapter) in self.adapters.iter().enumerate() {
            info!("Adapter 1-{}: {}", i + 1, adapter);
        }
        for (i, adapter) in self.adapters2.iter().enumerate() {
            info!("Adapter 2-{}: {}", i + 1, adapter);
        }
        for (i, b) in self.barcodes.iter().enumerate() {
            info!("Barcode-1-{}: {} {}", i, b.id, b.sequence);
        }
        for (i, b) in self.barcodes2.iter().enumerate() {
            info!("Barcode 2-{}: {} {}", i, b.id, b.sequence);
        }
    }

    fn load_preset(&mut self, preset: &AdapterPreset) {
        match preset {
            AdapterPreset::TruSeq => {
                self.adapters.push("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC".to_string());
                self.adapters2.push("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string());
            }
            AdapterPreset::SmallRna => {
                self.adapters.push("TGGAATTCTCGGGTGCCAAGG".to_string());
            }
            AdapterPreset::Methyl => {
                self.adapters.push("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC".to_string());
                self.adapters2.push("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string());
            }
            AdapterPreset::Ribo => {
                self.adapters.push("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC".to_string());
                self.adapters2.push("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string());
            }
            AdapterPreset::Nextera => {
                self.adapters.push("CTGTCTCTTATACACATCT".to_string());
                self.adapters2.push("CTGTCTCTTATACACATCT".to_string());
            }
            AdapterPreset::NexteraMp => {
                self.adapters.push("CTGTCTCTTATACACATCT".to_string());
                self.adapters2.push("CTGTCTCTTATACACATCT".to_string());
            }
            AdapterPreset::Off => (),
        }
    }

    fn load_fasta_file(&mut self, path: &std::path::Path) -> Vec<SequenceWithId> {
        let file = File::open(path).unwrap_or_else(|e| {
            log_fatal(&format!("Failed to open barcode file {}: {}", path.display(), e))
        });
        let reader = fasta::Reader::new(BufReader::new(file));
        let mut last_successful: Option<SequenceWithId> = None;
        let mut barcodes: Vec<SequenceWithId> = Vec::new();
        for result in reader.records() {
            match result {
                Ok(record) => {
                    let seq = String::from_utf8(record.seq().to_vec()).unwrap_or_else(|e| {
                        let context = self.format_last_successful(last_successful);
                        log_fatal(&format!(
                            "Failed to parse sequence in {} as UTF-8: {}{}",
                            path.display(),
                            e,
                            context
                        ))
                    });
                    let new_item = SequenceWithId {
                        id: record.id().to_string(),
                        sequence: seq.clone()
                    };
                    last_successful = Some(new_item.clone());
                    barcodes.push(new_item);
                }
                Err(e) => {
                    let context = self.format_last_successful(last_successful);
                    log_fatal(&format!(
                        "Error reading FASTA record from {}: {}{}",
                        path.display(),
                        e,
                        context
                    ))
                }
            }
        }
        barcodes
    }

    fn format_last_successful(&self, last_successful: Option<SequenceWithId>) -> String {
        if let Some(b) = last_successful {
            let mut info = format!(" (last successful: {} {})", b.id, b.sequence);
            if info.len() > 80 {
                info.truncate(77);
                info.push_str("...");
            }
            info
        } else {
            " (no records successfully loaded yet)".to_string()
        }
    }

    fn apply_revcomp(&mut self, mode: RevCompMode) {
        match mode {
            RevCompMode::Off => (),
            RevCompMode::On => {
                let rc1: Vec<String> = self.adapters.iter().map(|s| String::from_utf8(revcomp(s.as_bytes())).unwrap()).collect();
                let rc2: Vec<String> = self.adapters2.iter().map(|s| String::from_utf8(revcomp(s.as_bytes())).unwrap()).collect();
                self.adapters.extend(rc1);
                self.adapters2.extend(rc2);
            }
            RevCompMode::Only => {
                self.adapters = self.adapters.iter().map(|s| String::from_utf8(revcomp(s.as_bytes())).unwrap()).collect();
                self.adapters2 = self.adapters2.iter().map(|s| String::from_utf8(revcomp(s.as_bytes())).unwrap()).collect();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use std::io::Write;
    use clap::Parser;
    use tempfile::NamedTempFile;

    #[test]
    fn test_context_load_preset() {
        let mut args = Args::default();
        args.adapter_preset = Some(AdapterPreset::TruSeq);
        let context = Context::new(&args);
        assert!(!context.adapters.is_empty());
        assert!(!context.adapters2.is_empty());
        assert_eq!(context.adapters[0], "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
    }

    #[test]
    fn test_context_load_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">id1\nACGT").unwrap();
        let path = file.path().to_path_buf();

        let mut args = Args::default();
        args.adapters = Some(path);
        let context = Context::new(&args);
        assert_eq!(context.adapters.len(), 1);
        assert_eq!(context.adapters[0], "ACGT");
    }

    #[test]
    fn test_context_load_barcodes() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">bc1\nGGGG").unwrap();
        let path = file.path().to_path_buf();

        let mut args = Args::default();
        args.barcodes = Some(path);
        let context = Context::new(&args);
        assert_eq!(context.barcodes.len(), 1);
        assert_eq!(context.barcodes[0].id, "bc1");
        assert_eq!(context.barcodes[0].sequence, "GGGG");
    }

    #[test]
    fn test_context_revcomp() {
        let mut args = Args::default();
        args.adapter_seq = Some("ACGT".to_string());
        args.adapter_revcomp = RevCompMode::On;
        let context = Context::new(&args);
        assert_eq!(context.adapters.len(), 2);
        assert_eq!(context.adapters[0], "ACGT");
        assert_eq!(context.adapters[1], "ACGT"); // RC of ACGT is ACGT
        
        let mut args = Args::default();
        args.adapter_seq = Some("AAAA".to_string());
        args.adapter_revcomp = RevCompMode::Only;
        let context = Context::new(&args);
        assert_eq!(context.adapters.len(), 1);
        assert_eq!(context.adapters[0], "TTTT");
    }

    #[test]
    fn test_context_basic() {
        let args = Args {
            threads: 4,
            target: "test_target".to_string(),
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(args.threads, 4);
        assert_eq!(context.run_type, RunType::Single);
    }

    #[test]
    fn test_context_paired() {
        let args = Args {
            reads: Some(PathBuf::from("reads1.fastq")),
            reads2: Some(PathBuf::from("reads2.fastq")),
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(context.run_type, RunType::Paired);
    }

    #[test]
    fn test_context_qtrim_thresh() {
        let args = Args {
            reads: Some(PathBuf::from("reads.fastq")),
            qtrim: QualTrimType::Tail,
            qtrim_format: QualityType::Sanger,
            qtrim_threshold: 20,
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(context.corrected_qtrim_threshold, 20 + 33);
    }

    #[test]
    fn test_context_adapter_preset() {
        let args = Args {
            adapter_preset: Some(AdapterPreset::TruSeq),
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(context.adapter_removal, AdapterRemoval::Normal);
    }

    #[test]
    fn test_context_adapter_preset_paired() {
        let args = Args {
            reads2: Some(PathBuf::from("r2.fq")),
            adapter_preset: Some(AdapterPreset::TruSeq),
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(context.adapter_removal, AdapterRemoval::Normal2);
    }
}
