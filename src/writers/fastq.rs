use std::io;
use bio::io::fastq;
use crate::types::{RecordDto};
use crate::writers::generic::{GenericWriter, GenericWriterBackend};

type FastqBackend<'a> = fastq::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastqBackend<'a> {
    fn new(writer: Box<dyn io::Write + Send + 'a>) -> Self {
        fastq::Writer::new(writer)
    }

    fn write_record(&mut self, record: &RecordDto) -> io::Result<()> {
        match record {
            RecordDto::Single(rec) => self.write_record(rec),
            RecordDto::Paired(_) => panic!("Fastq write can handle only single records"),
        }
    }

    fn get_extension() -> &'static str {
        "fastq"
    }
}

pub type SingleWriterFastq<'a> = GenericWriter<'a, FastqBackend<'a>>;


#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::SingeRecordDto;
    use crate::writers::Writer;

    #[test]
    fn test_fastq_writer_cursor() {
        let buffers = std::sync::Arc::new(std::sync::Mutex::new(std::collections::HashMap::new()));
        {
            let mut writer: SingleWriterFastq = crate::writers::generic::get_test_writer(buffers.clone());
            let record = RecordDto::Single(SingeRecordDto::with_attrs(
                "id1", Some("desc1"), b"ACGT", b"!!!!"
            ));
            writer.write(&record, None).expect("Failed to write record");
        }
        
        let bufs = buffers.lock().unwrap();
        let result = String::from_utf8(bufs.get("/\\/default/\\/").unwrap().clone()).expect("Not valid UTF-8");
        assert!(result.contains("@id1 desc1"));
        assert!(result.contains("ACGT"));
        assert!(result.contains("!!!!"));
    }
}
