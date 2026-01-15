use std::io::{self, BufRead};
use std::path::Path;
use bio::io::fasta;
use bio::io::fasta::{FastaRead, Record};
use crate::readers::generic::{GenericReader, GenericReaderBackend};
use crate::utils::get_reader;
use crate::types::CompressionType;

pub struct FastaBackend<R: BufRead + Send> {
    reader: fasta::Reader<R>,
}

impl<R: BufRead + Send> GenericReaderBackend<R, Record> for FastaBackend<R> {
    type Error = io::Error;
    fn create(reader: R) -> Self {
        Self {
            reader: fasta::Reader::from_bufread(reader),
        }
    }

    fn read_record(&mut self, record: &mut Record) -> Result<(), Self::Error> {
        self.reader.read(record)
    }

    fn is_empty(record: &Record) -> bool {
        record.is_empty()
    }

    fn get_id(record: &Record) -> &str {
        record.id()
    }
}

pub type SingleReaderFasta<R> = GenericReader<R, Record, FastaBackend<R>>;

impl SingleReaderFasta<io::BufReader<Box<dyn io::Read + Send>>> {
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
    fn test_fasta_reader() {
        let data = b">id1 desc1\nACGT\n>id2 desc2\nTGCA\n";
        let cursor = Cursor::new(data);
        let reader = SingleReaderFasta::new(cursor, "test.fa".to_string(), CompressionType::Uncompressed);
        
        let records: Vec<_> = reader.collect();
        assert_eq!(records.len(), 2);
        let ref r = records[0];
        assert_eq!(r.id(), "id1");
        assert_eq!(r.seq(), b"ACGT");
    }
}
