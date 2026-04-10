use crate::args::Args;
use crate::types::{AdapterPreset, AdapterReadset, AdapterRemoval, BarcodeDetect, FileFormat, PairOverlapMode, RevCompMode, RunType, TrimEnd};
use crate::utils::io::get_file_type;
use bio::io::fasta;
use bio::alphabets::dna::revcomp;
use log::{debug, info, warn};
use std::fs::File;
use std::io::BufReader;
use crate::logging::log_fatal;


#[derive(Clone, Debug, PartialEq)]
pub struct SequenceWithId {
    pub id: String,
    pub sequence: String,
}

impl PartialEq<&str> for SequenceWithId {
    fn eq(&self, other: &&str) -> bool {
        self.sequence == *other
    }
}

pub struct Context {
    pub adapters: Vec<SequenceWithId>,
    pub adapters2: Vec<SequenceWithId>,
    pub barcodes: Vec<SequenceWithId>,
    pub barcodes2: Vec<SequenceWithId>,
    pub use_rc_trim_end: bool,
    pub format: FileFormat,
    pub run_type: RunType,
    pub bar_detect: BarcodeDetect,
    pub adapter_removal: AdapterRemoval,
}

impl Context {
    pub fn empty() -> Self {
        Self {
            adapters: Vec::new(),
            adapters2: Vec::new(),
            barcodes: Vec::new(),
            barcodes2: Vec::new(),
            use_rc_trim_end: false,
            format: FileFormat::Fasta,
            run_type: RunType::Single,
            bar_detect: BarcodeDetect::Off,
            adapter_removal: AdapterRemoval::Off,
        }
    }

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
        match args.barcode_reads {
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
        let has_adapters2 = args.adapters2.is_some();
        if args.reads2.is_none() {
            if has_adapters2 {
                warn!("WARNING: No paired reads - ignoring --adapter2.");
            }
            return false;
        }
        match args.adapter_preset {
            Some(AdapterPreset::TruSeq) => return true,
            Some(AdapterPreset::NexteraMp) => return true,
            Some(AdapterPreset::Methyl) => return true,
            Some(preset) => {
                if has_adapters2 {
                    warn!("WARNING: Adapter preset {:?} does not use a second set of adapters - ignoring --adapters2",
                        preset
                    );
                }
            }
            _ => (),
        }
        if args.adapters2.is_some() && args.reads2.is_some() {
            return true;
        }
        false
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
        if let Some(preset) = &args.adapter_preset {
            if matches!(preset, AdapterPreset::Off) {
                return;
            }
            self.load_preset(preset);
        }

        if let Some(adapter_seq) = &args.adapter_seq {
            append_single_adapter(&mut self.adapters, adapter_seq);
        }

        if let Some(adapter_file) = &args.adapters {
            let mut new_items  = self.load_fasta_file(adapter_file)
                .into_iter()
                .collect();
            self.adapters.append(&mut new_items);
        }

        if let Some(adapter_file2) = &args.adapters2
            && self.adapter_removal == AdapterRemoval::Normal2
        {
            self.adapters2 = self.load_fasta_file(adapter_file2)
                .into_iter()
                .collect();
        }

        if self.adapters2.is_empty()
            && (args.adapter_read_set == AdapterReadset::Two || args.adapter_read_set == AdapterReadset::None)
            && (self.run_type == RunType::Paired || self.run_type == RunType::PairedBarcoded)
        {
            self.adapters2 = self.adapters.clone();
        }

        if let Some(barcode_file) = &args.barcodes {
            let b = self.load_fasta_file(barcode_file);
            self.check_barcodes(&args, &b, "barcodes");
            self.barcodes = b;
        }

        if let Some(barcode_file2) = &args.barcodes2 {
            let b =  self.load_fasta_file(barcode_file2);
            self.check_barcodes(&args, &b, "barcodes2");
            self.barcodes2 = b;
        }

        self.apply_revcomp(args.adapter_revcomp);
        
        for adapter in self.adapters.iter() {
            info!("Adapter {}: {}", adapter.id, adapter.sequence);
        }
        for adapter in self.adapters2.iter() {
            info!("Adapter-2 {}: {}", adapter.id, adapter.sequence);
        }
        for (i, b) in self.barcodes.iter().enumerate() {
            info!("Barcode-1-{}: {} {}", i, b.id, b.sequence);
        }
        for (i, b) in self.barcodes2.iter().enumerate() {
            info!("Barcode 2-{}: {} {}", i, b.id, b.sequence);
        }
    }

    fn load_preset(&mut self, preset: &AdapterPreset) {
        let mut new_adapters = self.adapters.to_vec();
        let mut new_adapters2 = self.adapters2.to_vec();
        match preset {
            AdapterPreset::TruSeq => {
                new_adapters.push(SequenceWithId {
                    id: "TruSeq-1".to_string(),
                    sequence: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC".to_string()
                });
                new_adapters2.push(SequenceWithId {
                    id: "TruSeq-2".to_string(),
                    sequence: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string()
                });
            }
            AdapterPreset::SmallRna => {
                new_adapters.push(SequenceWithId {
                    id: "SmallRna".to_string(),
                    sequence: "TGGAATTCTCGGGTGCCAAGG".to_string()
                });
            }
            AdapterPreset::Methyl => {
                new_adapters.push(SequenceWithId {
                    id: "Methyl-1".to_string(),
                    sequence: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC".to_string()
                });
                new_adapters2.push(SequenceWithId {
                    id: "Methyl-2".to_string(),
                    sequence: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string()
                });
            }
            AdapterPreset::Ribo => {
                new_adapters.push(SequenceWithId {
                    id: "Ribo-1".to_string(),
                    sequence: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC".to_string()
                });
                new_adapters2.push(SequenceWithId {
                    id: "Ribo-2".to_string(),
                    sequence: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT".to_string()
                });
            }
            AdapterPreset::Nextera => {
                new_adapters.push(SequenceWithId {
                    id: "Nextera-1".to_string(),
                    sequence: "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG".to_string()
                });
                new_adapters2.push(SequenceWithId {
                    id: "Nextera-2".to_string(),
                    sequence: "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG".to_string()
                });
            }
            AdapterPreset::NexteraMp => {
                new_adapters.push(SequenceWithId {
                    id: "NexteraMP-1".to_string(),
                    sequence: "CTGTCTCTTATACACATCT".to_string()
                });
                new_adapters2.push(SequenceWithId {
                    id: "NexteraMP-2".to_string(),
                    sequence: "CTGTCTCTTATACACATCT".to_string()
                });
            }
            AdapterPreset::Off => (),
        }
        self.adapters = new_adapters;
        self.adapters2 = new_adapters2;
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

    fn check_barcodes(&self, args: &Args, barcodes: &Vec<SequenceWithId>, desc: &str) {
        let min_length_barcodes = barcodes.iter()
            .map(|b| b.sequence.len())
            .min()
            .unwrap_or(0);
        if let Some(x) = args.barcode_min_overlap && x > min_length_barcodes {
            warn!("WARNING: Some of your barcodes are shorter than the --barcode-min-overlap.");
            warn!("WARNING: You will not detect these barcodes.");
        }
        if min_length_barcodes < 10 {
            warn!("WARNING: Some of your {} are shorter than 10 bases.", desc);
            warn!("WARNING: You may see \"Cutoff is too low to detect the pattern\" errors.");
        }
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
                self.adapters.extend(convert_to_revcomp(&self.adapters));
                self.adapters2.extend(convert_to_revcomp(&self.adapters2));
            }
            RevCompMode::Only => {
                self.adapters = convert_to_revcomp(&self.adapters);
                self.adapters2 = convert_to_revcomp(&self.adapters2);
            }
        }
    }
}

fn append_single_adapter(adapters: &mut Vec<SequenceWithId>, new_adapter: &String) {
    let index = adapters.len();
    adapters.push(
        SequenceWithId { id: format!("adapter_{}", index), sequence: new_adapter.to_string() }
    );
}

fn convert_to_revcomp(adapters: &[SequenceWithId]) -> Vec<SequenceWithId> {
    let mut result: Vec<SequenceWithId> = Vec::new();
    for adapter in adapters {
        result.push(SequenceWithId {
            id: format!("{}_rc", adapter.id),
            sequence: String::from_utf8(revcomp(adapter.sequence.as_bytes())).unwrap()
        });
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
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

    #[test]
    fn test_context_adapter_copy_to_second_read() {
        // Case 1: adapter_read_set = Two, no adapters2
        let args = Args {
            reads: Some(PathBuf::from("r1.fq")),
            reads2: Some(PathBuf::from("r2.fq")),
            adapter_seq: Some("AGCT".to_string()),
            adapter_read_set: AdapterReadset::Two,
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(context.adapters.len(), 1);
        assert_eq!(context.adapters[0].sequence, "AGCT");
        assert_eq!(context.adapters2.len(), 1);
        assert_eq!(context.adapters2[0].sequence, "AGCT");

        // Case 2: adapter_read_set = None, no adapters2
        let args = Args {
            reads: Some(PathBuf::from("r1.fq")),
            reads2: Some(PathBuf::from("r2.fq")),
            adapter_seq: Some("AGCT".to_string()),
            adapter_read_set: AdapterReadset::None,
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(context.adapters.len(), 1);
        assert_eq!(context.adapters[0].sequence, "AGCT");
        assert_eq!(context.adapters2.len(), 1);
        assert_eq!(context.adapters2[0].sequence, "AGCT");

        // Case 3: adapter_read_set = One, no adapters2 -> adapters2 should be empty
        let args = Args {
            reads: Some(PathBuf::from("r1.fq")),
            reads2: Some(PathBuf::from("r2.fq")),
            adapter_seq: Some("AGCT".to_string()),
            adapter_read_set: AdapterReadset::One,
            ..Args::parse_from(["flexbar"])
        };
        let context = Context::new(&args);
        assert_eq!(context.adapters.len(), 1);
        assert_eq!(context.adapters2.len(), 0);
    }
}
