use crate::config::Config;
use crate::types::{AdapterPreset, AdapterRemoval, RevCompMode};
use bio::io::fasta;
use bio::alphabets::dna::revcomp;
use log::{info};
use std::fs::File;
use std::io::BufReader;
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
}

impl Context {
    pub fn new() -> Self {
        Self {
            adapters: Vec::new(),
            adapters2: Vec::new(),
            barcodes: Vec::new(),
            barcodes2: Vec::new(),
        }
    }

    pub fn load(&mut self, config: &Config) {
        let args = &config.args;

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

        if config.adapter_removal == AdapterRemoval::Normal2 {
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
    use tempfile::NamedTempFile;

    #[test]
    fn test_context_load_preset() {
        let mut args = Args::default();
        args.adapter_preset = Some(AdapterPreset::TruSeq);
        let config = Config::new(args);
        let mut context = Context::new();
        context.load(&config);
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
        let config = Config::new(args);
        let mut context = Context::new();
        context.load(&config);
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
        let config = Config::new(args);
        let mut context = Context::new();
        context.load(&config);
        assert_eq!(context.barcodes.len(), 1);
        assert_eq!(context.barcodes[0].id, "bc1");
        assert_eq!(context.barcodes[0].sequence, "GGGG");
    }

    #[test]
    fn test_context_revcomp() {
        let mut args = Args::default();
        args.adapter_seq = Some("ACGT".to_string());
        args.adapter_revcomp = RevCompMode::On;
        let config = Config::new(args);
        let mut context = Context::new();
        context.load(&config);
        assert_eq!(context.adapters.len(), 2);
        assert_eq!(context.adapters[0], "ACGT");
        assert_eq!(context.adapters[1], "ACGT"); // RC of ACGT is ACGT
        
        let mut args = Args::default();
        args.adapter_seq = Some("AAAA".to_string());
        args.adapter_revcomp = RevCompMode::Only;
        let config = Config::new(args);
        let mut context = Context::new();
        context.load(&config);
        assert_eq!(context.adapters.len(), 1);
        assert_eq!(context.adapters[0], "TTTT");
    }
}
