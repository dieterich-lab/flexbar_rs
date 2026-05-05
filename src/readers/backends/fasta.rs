use std::io::{BufRead};
use bio::io::fasta;
use bio::io::fasta::FastaRead;
use crate::logging::log_fatal;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::records::unified::UnifiedRead;
use crate::records::header::SequenceHeader;
pub struct FastaBackend<READER: BufRead + Send> {
    reader: fasta::Reader<READER>,
    source: String,
}

impl<Reader: BufRead + Send> FastaBackend<Reader> {
    pub fn new(reader: Reader, source: &str) -> Self {
        Self {
            reader: fasta::Reader::from_bufread(reader),
            source: source.to_string(),
        }
    }
}

impl<Reader: BufRead + Send> GenericReaderBackend<Reader> for FastaBackend<Reader> {

    fn read_record(&mut self) -> Option<UnifiedRead> {
        let mut fasta_record = fasta::Record::default();
        match self.reader.read(&mut fasta_record) {
            Ok(()) => {
                if fasta_record.is_empty() {
                    None
                } else {
                    Some(UnifiedRead::from_fasta(fasta_record))
                }
            }
            Err(e) => {
                log_fatal(&format!("Failed to read from {}: {:?}", self.source, e));
            }
        }
    }

    fn get_header(&self) -> SequenceHeader {
        SequenceHeader::Default
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::records::generic_read::GenericRawRead;

    #[test]
    fn test_fasta_reader() {
        let data = b">id1 desc1\nACGT\n>id2 desc2\nTGCA\n";
        let cursor = Cursor::new(data);
        let mut reader = FastaBackend::new(cursor, "test.fa");
        
        let r1 = reader.read_record().unwrap();
        assert_eq!(r1.id(), "id1");
        assert_eq!(r1.seq(), b"ACGT");
        assert_eq!(r1.qual(), None);
        let r2 = reader.read_record().unwrap();
        assert_eq!(r2.id(), "id2");
        let r3 = reader.read_record();
        assert!(r3.is_none());
    }
}
