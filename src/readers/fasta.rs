use std::io::{self, BufRead};
use bio::io::fasta;
use bio::io::fasta::FastaRead;
use crate::types::{RecordDto, SingeRecordDto};
use crate::logging::log_fatal;
use crate::readers::generic::{GenericReader, GenericReaderBackend, GenericPairedBackend};

pub struct FastaBackend<R: BufRead + Send> {
    reader: fasta::Reader<R>,
}

impl<R: BufRead + Send> GenericReaderBackend<R> for FastaBackend<R> {
    type Error = io::Error;
    fn new(mut readers: Vec<R>) -> Self {
        if readers.len() > 1 {
            log_fatal("FastaBackend only supports one input file.");
        }
        Self {
            reader: fasta::Reader::from_bufread(readers.remove(0)),
        }
    }

    fn read_record(&mut self) -> Result<Option<RecordDto>, Self::Error> {
        let mut fasta_record = fasta::Record::default();
        self.reader.read(&mut fasta_record)?;
        if fasta_record.is_empty() {
            Ok(None)
        } else {
            Ok(Some(RecordDto::Single(SingeRecordDto::with_attrs(
                fasta_record.id(),
                fasta_record.desc(),
                fasta_record.seq(),
                &vec![b'I'; fasta_record.seq().len()],
            ))))
        }
    }
}

pub type SingleReaderFasta<R> = GenericReader<R, FastaBackend<R>>;
pub type PairedFastaBackend<R> = GenericPairedBackend<R, FastaBackend<R>>;
pub type PairedReaderFasta<R> = GenericReader<R, PairedFastaBackend<R>>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::types::CompressionType;

    #[test]
    fn test_fasta_reader() {
        let data = b">id1 desc1\nACGT\n>id2 desc2\nTGCA\n";
        let cursor = Cursor::new(data);
        let reader = SingleReaderFasta::new(
            vec![cursor], "test.fa".to_string(), CompressionType::Uncompressed, 0 ,None, true
        );
        
        let records: Vec<_> = reader.collect();
        assert_eq!(records.len(), 2);
        let ref r = records[0];
        assert_eq!(r.id(), "id1");
        assert_eq!(r.seq(), b"ACGT");
        assert_eq!(r.qual().unwrap(), b"IIII");
    }
}
