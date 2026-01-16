use std::io;
use bio::io::fasta;
use bio::io::fastq;
use crate::writers::generic::{GenericWriter, GenericWriterBackend};

type FastaBackend<'a> = fasta::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a, fastq::Record> for FastaBackend<'a> {
    fn new(writer: Box<dyn io::Write + Send + 'a>) -> Self {
        fasta::Writer::new(writer)
    }

    fn write_record(&mut self, record: &fastq::Record) -> io::Result<()> {
        let fasta_record = fasta::Record::with_attrs(
            record.id(),
            record.desc(),
            record.seq(),
        );
        self.write_record(&fasta_record)
    }

    fn get_extension() -> &'static str {
        "fasta"
    }
}

pub type SingleWriterFasta<'a> = GenericWriter<'a, fastq::Record, FastaBackend<'a>>;

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;
    use crate::writers::Writer;

    #[test]
    fn test_fasta_writer_cursor() {
        let mut buffer = Vec::new();
        {
            let mut writer = SingleWriterFasta::from_writer(Box::new(&mut buffer));
            let record = Record::with_attrs("id1", Some("desc1"), b"ACGT", b"!!!!");
            writer.write(&record).expect("Failed to write record");
        }
        
        let result = String::from_utf8(buffer).expect("Not valid UTF-8");
        assert!(result.contains(">id1 desc1"));
        assert!(result.contains("ACGT"));
        assert!(!result.contains("!!!!"));
    }
}


