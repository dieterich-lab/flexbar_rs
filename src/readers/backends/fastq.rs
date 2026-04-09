use std::io::BufRead;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::types::SingeRecordDto;

pub struct FastqBackend<READER: BufRead + Send> {
    reader: fastq::Reader<READER>,
}

impl<Reader: BufRead + Send> FastqBackend<Reader> {
    pub fn new(reader: Reader) -> Self {
        Self {
            reader: fastq::Reader::from_bufread(reader),
        }
    }
}

impl<Reader: BufRead + Send> GenericReaderBackend<Reader> for FastqBackend<Reader> {
    type Error = bio::io::fastq::Error;

    fn read_record(&mut self) -> Result<Option<SingeRecordDto>, Self::Error> {
        let mut record = SingeRecordDto::new();
        self.reader.read(&mut record)?;
        if record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(record))
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
        let mut reader = FastqBackend::new(cursor);
        
        let r1 = reader.read_record().unwrap().unwrap();
        assert_eq!(r1.id(), "id1");
        assert_eq!(r1.seq(), b"ACGT");
        assert_eq!(r1.qual(), b"!!!!");
        let r2 = reader.read_record().unwrap();
        assert!(r2.is_none());
    }
}
