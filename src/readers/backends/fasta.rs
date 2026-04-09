use std::io::{self, BufRead};
use bio::io::fasta;
use bio::io::fasta::FastaRead;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::types::SingeRecordDto;
pub struct FastaBackend<READER: BufRead + Send> {
    reader: fasta::Reader<READER>,
}

impl<Reader: BufRead + Send> FastaBackend<Reader> {
    pub fn new(reader: Reader) -> Self {
        Self {
            reader: fasta::Reader::from_bufread(reader),
        }
    }
}

impl<Reader: BufRead + Send> GenericReaderBackend<Reader> for FastaBackend<Reader> {
    type Error = io::Error;

    fn read_record(&mut self) -> Result<Option<SingeRecordDto>, Self::Error> {
        let mut fasta_record = fasta::Record::default();
        self.reader.read(&mut fasta_record)?;
        if fasta_record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(SingeRecordDto::with_attrs(
                fasta_record.id(),
                fasta_record.desc(),
                fasta_record.seq(),
                &vec![b'I'; fasta_record.seq().len()],
            )))
        }
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
        let mut reader = FastaBackend::new(cursor);
        
        let r1 = reader.read_record().unwrap().unwrap();
        assert_eq!(r1.id(), "id1");
        assert_eq!(r1.seq(), b"ACGT");
        assert_eq!(r1.qual(), b"IIII");
        let r2 = reader.read_record().unwrap().unwrap();
        assert_eq!(r2.id(), "id2");
        let r3 = reader.read_record().unwrap();
        assert!(r3.is_none());
    }
}
