use std::io::BufRead;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use crate::logging::log_fatal;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::types::SingeRecordDto;

pub struct FastqBackend<READER: BufRead + Send> {
    reader: fastq::Reader<READER>,
    source: String,
}

impl<Reader: BufRead + Send> FastqBackend<Reader> {
    pub fn new(reader: Reader, source: &str) -> Self {
        Self {
            reader: fastq::Reader::from_bufread(reader),
            source: source.to_string(),
        }
    }
}

impl<Reader: BufRead + Send> GenericReaderBackend<Reader> for FastqBackend<Reader> {
    fn read_record(&mut self) -> Option<SingeRecordDto> {
        let mut record = SingeRecordDto::new();
        match self.reader.read(&mut record) {
            Ok(()) => {
                if record.is_empty() {
                    None
                } else {
                    Some(record)
                }
            },
            Err(e) => {
                log_fatal(&format!("Failed to read from {}: {:?}", self.source, e));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    #[test]
    fn test_fastq_reader() {
        let data = b"@id1\nACGT\n+\n!!!!\n";
        let cursor = Cursor::new(data);
        let mut reader = FastqBackend::new(cursor,"test");
        
        let r1 = reader.read_record().unwrap();
        assert_eq!(r1.id(), "id1");
        assert_eq!(r1.seq(), b"ACGT");
        assert_eq!(r1.qual(), b"!!!!");
        let r2 = reader.read_record();
        assert!(r2.is_none());
    }
}
