use std::io::BufRead;
use bio::io::fastq;
use bio::io::fastq::{FastqRead};
use crate::types::{RecordDto, SingeRecordDto};
use crate::logging::log_fatal;
use crate::readers::generic::{GenericReader, GenericReaderBackend, GenericPairedBackend};

pub struct FastqBackend<R: BufRead + Send> {
    reader: fastq::Reader<R>,
}

impl<R: BufRead + Send> GenericReaderBackend<R> for FastqBackend<R> {
    type Error = bio::io::fastq::Error;
    fn new(mut readers: Vec<R>) -> Self {
        if readers.len() > 1 {
            log_fatal("FastqBackend only supports one input file.");
        }
        Self {
            reader: fastq::Reader::from_bufread(readers.remove(0)),
        }
    }

    fn read_record(&mut self) -> Result<Option<RecordDto>, Self::Error> {
        let mut record = SingeRecordDto::new();
        self.reader.read(&mut record)?;
        if record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(RecordDto::Single(record)))
        }
    }
}

pub type SingleReaderFastq<R> = GenericReader<R, FastqBackend<R>>;
pub type PairedFastqBackend<R> = GenericPairedBackend<R, FastqBackend<R>>;
pub type PairedReaderFastq<R> = GenericReader<R, PairedFastqBackend<R>>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::types::CompressionType;

    #[test]
    fn test_fastq_reader() {
        let data = b"@id1\nACGT\n+\n!!!!\n";
        let cursor = Cursor::new(data);
        let reader = SingleReaderFastq::new(
            vec![cursor], "test.fq".to_string(), CompressionType::Uncompressed, 0, None, true
        );
        
        let records: Vec<_> = reader.collect();
        assert_eq!(records.len(), 1);
        let ref r = records[0];
        assert_eq!(r.id(), "id1");
        assert_eq!(r.seq(), b"ACGT");
        assert_eq!(r.qual().unwrap(), b"!!!!");
    }
}
