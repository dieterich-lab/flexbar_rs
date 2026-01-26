use std::io::BufRead;
use bio::io::fastq;
use bio::io::fastq::{FastqRead};
use crate::types::RecordDto;
use crate::readers::generic::{GenericReader, GenericReaderBackend};

pub struct FastqBackend<R: BufRead + Send> {
    reader: fastq::Reader<R>,
}

impl<R: BufRead + Send> GenericReaderBackend<R> for FastqBackend<R> {
    type Error = bio::io::fastq::Error;
    fn new(reader: R) -> Self {
        Self {
            reader: fastq::Reader::from_bufread(reader),
        }
    }

    fn read_record(&mut self) -> Result<Option<RecordDto>, Self::Error> {
        let mut record = RecordDto::new();
        self.reader.read(&mut record)?;
        if record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(record))
        }
    }
}

pub type SingleReaderFastq<R> = GenericReader<R, FastqBackend<R>>;

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
            cursor, "test.fq".to_string(), CompressionType::Uncompressed, 0, None
        );
        
        let records: Vec<_> = reader.collect();
        assert_eq!(records.len(), 1);
        let ref r = records[0];
        assert_eq!(r.id(), "id1");
        assert_eq!(r.seq(), b"ACGT");
        assert_eq!(r.qual(), b"!!!!");
    }
}
