use std::io::{self, BufRead};
use std::path::Path;
use bio::io::fasta;
use bio::io::fasta::FastaRead;
use crate::args::Args;
use crate::types::RecordDto;
use crate::readers::generic::{GenericReader, GenericReaderBackend};
use crate::utils::get_raw_reader;

pub struct FastaBackend<R: BufRead + Send> {
    reader: fasta::Reader<R>,
}

impl<R: BufRead + Send> GenericReaderBackend<R> for FastaBackend<R> {
    type Error = io::Error;
    fn new(reader: R) -> Self {
        Self {
            reader: fasta::Reader::from_bufread(reader),
        }
    }

    fn read_record(&mut self) -> Result<Option<RecordDto>, Self::Error> {
        let mut fasta_record = fasta::Record::default();
        self.reader.read(&mut fasta_record)?;
        if fasta_record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(RecordDto::with_attrs(
                fasta_record.id(),
                fasta_record.desc(),
                fasta_record.seq(),
                &vec![b'I'; fasta_record.seq().len()],
            )))
        }
    }
}

pub type SingleReaderFasta<R> = GenericReader<R, FastaBackend<R>>;

impl SingleReaderFasta<io::BufReader<Box<dyn io::Read + Send>>> {
    pub fn create_form_args<P: AsRef<Path>>(path: P, args: &Args) -> io::Result<Self> {
        let (reader, source, compression) = get_raw_reader(path)?;
        Ok(Self::new(reader, source, compression, args.skip, args.count))
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
        let reader = SingleReaderFasta::new(cursor, "test.fa".to_string(), CompressionType::Uncompressed, 0 ,None);
        
        let records: Vec<_> = reader.collect();
        assert_eq!(records.len(), 2);
        let ref r = records[0];
        assert_eq!(r.id(), "id1");
        assert_eq!(r.seq(), b"ACGT");
        assert_eq!(r.qual(), b"IIII");
    }
}
