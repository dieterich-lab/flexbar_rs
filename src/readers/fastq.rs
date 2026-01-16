use std::io::{self, BufRead};
use std::path::Path;
use bio::io::fastq;
use bio::io::fastq::{FastqRead, Record};
use crate::readers::generic::{GenericReader, GenericReaderBackend};
use crate::utils::get_reader;
use crate::types::CompressionType;

pub struct FastqBackend<R: BufRead + Send> {
    reader: fastq::Reader<R>,
}

impl<R: BufRead + Send> GenericReaderBackend<R, Record> for FastqBackend<R> {
    type Error = bio::io::fastq::Error;
    fn new(reader: R) -> Self {
        Self {
            reader: fastq::Reader::from_bufread(reader),
        }
    }

    fn read_record(&mut self) -> Result<Option<Record>, Self::Error> {
        let mut record = Record::new();
        self.reader.read(&mut record)?;
        if record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(record))
        }
    }

    fn get_id(record: &Record) -> &str {
        record.id()
    }
}

pub type SingleReaderFastq<R> = GenericReader<R, Record, FastqBackend<R>>;

impl SingleReaderFastq<io::BufReader<Box<dyn io::Read + Send>>> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let (reader, source, compression) = get_reader(path)?;
        Ok(Self::new(reader, source, compression))
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
        let reader = SingleReaderFastq::new(cursor, "test.fq".to_string(), CompressionType::Uncompressed);
        
        let records: Vec<_> = reader.collect();
        assert_eq!(records.len(), 1);
        let ref r = records[0];
        assert_eq!(r.id(), "id1");
        assert_eq!(r.seq(), b"ACGT");
        assert_eq!(r.qual(), b"!!!!");
    }
}
