use std::io::BufRead;
use rust_htslib::bam::{self, Read};
use crate::logging::log_fatal;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::records::unified_read::UnifiedRead;
use rust_htslib::bam::record::Record;
use log::warn;
use std::sync::atomic::{AtomicBool, Ordering};

static BAM_MAPPING_WARNING_SHOWN: AtomicBool = AtomicBool::new(false);

pub struct BamBackend<READER: BufRead + Send> {
    reader: bam::Reader,
    source: String,
    _phantom: std::marker::PhantomData<READER>,
}

impl<READER: BufRead + Send> BamBackend<READER> {
    pub fn new(_reader: READER, source: &str) -> Self {
        // rust-htslib's bam::Reader usually takes a path or uses stdin.
        // If we are given a BufRead, it might be tricky to use with htslib which wants a path for BAM.
        // However, if source is "-", we can use stdin.
        let reader = if source == "-" {
            bam::Reader::from_stdin().unwrap_or_else(|e| {
                log_fatal(&format!("Failed to open BAM from stdin: {:?}", e));
            })
        } else {
            bam::Reader::from_path(source).unwrap_or_else(|e| {
                log_fatal(&format!("Failed to open BAM file {}: {:?}", source, e));
            })
        };

        Self {
            reader,
            source: source.to_string(),
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<READER: BufRead + Send> GenericReaderBackend<READER> for BamBackend<READER> {
    fn read_record(&mut self) -> Option<UnifiedRead> {
        let mut record = Record::new();
        match self.reader.read(&mut record) {
            Some(Ok(())) => {
                if !record.is_unmapped() && !BAM_MAPPING_WARNING_SHOWN.load(Ordering::Relaxed) {
                    warn!("BAM record found that is not unmapped (0x4 flag not set). Mapping information may be inconsistent in the output.");
                    BAM_MAPPING_WARNING_SHOWN.store(true, Ordering::Relaxed);
                }

                Some(UnifiedRead::from_bam(record))
            },
            None => None,
            Some(Err(e)) => {
                log_fatal(&format!("Failed to read from BAM {}: {:?}", self.source, e));
            }
        }
    }
}
