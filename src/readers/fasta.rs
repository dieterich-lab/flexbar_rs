use std::io::{self, BufRead};
use std::path::Path;
use bio::io::fasta;
use bio::io::fasta::FastaRead;
use bio::io::fastq;
use crate::readers::generic::{GenericReader, GenericReaderBackend};
use crate::utils::get_reader;

pub struct FastaBackend<R: BufRead + Send> {
    reader: fasta::Reader<R>,
    fasta_record: fasta::Record,
}

impl<R: BufRead + Send> GenericReaderBackend<R, fastq::Record> for FastaBackend<R> {
    type Error = io::Error;
    fn new(reader: R) -> Self {
        Self {
            reader: fasta::Reader::from_bufread(reader),
            fasta_record: fasta::Record::new(),
        }
    }

    fn read_record(&mut self) -> Result<Option<fastq::Record>, Self::Error> {
        self.reader.read(&mut self.fasta_record)?;
        if self.fasta_record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(fastq::Record::with_attrs(
                self.fasta_record.id(),
                self.fasta_record.desc(),
                self.fasta_record.seq(),
                &vec![b'I'; self.fasta_record.seq().len()],
            )))
        }
    }

    fn get_id(record: &fastq::Record) -> &str {
        record.id()
    }
}

pub type SingleReaderFasta<R> = GenericReader<R, fastq::Record, FastaBackend<R>>;

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
    use crate::types::CompressionType;

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
        assert_eq!(r.qual(), b"IIII");
    }
}
