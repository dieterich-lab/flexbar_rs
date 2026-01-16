use std::io;
use bio::io::fastq;
use crate::writers::generic::{GenericWriter, GenericWriterBackend};

type FastqBackend<'a> = fastq::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a, fastq::Record> for FastqBackend<'a> {
    fn new(writer: Box<dyn io::Write + Send + 'a>) -> Self {
        fastq::Writer::new(writer)
    }

    fn write_record(&mut self, record: &fastq::Record) -> io::Result<()> {
        self.write_record(record)
    }

    fn get_extension() -> &'static str {
        "fastq"
    }
}

pub type SingleWriterFastq<'a> = GenericWriter<'a, fastq::Record, FastqBackend<'a>>;


#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;
    use crate::writers::Writer;

    #[test]
    fn test_fastq_writer_cursor() {
        let mut buffer = Vec::new();
        {
            let mut writer = SingleWriterFastq::from_writer(Box::new(&mut buffer));
            let record = Record::with_attrs("id1", Some("desc1"), b"ACGT", b"!!!!");
            writer.write(&record).expect("Failed to write record");
        }
        
        let result = String::from_utf8(buffer).expect("Not valid UTF-8");
        assert!(result.contains("@id1 desc1"));
        assert!(result.contains("ACGT"));
        assert!(result.contains("!!!!"));
    }
}
