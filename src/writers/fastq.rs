use std::io;
use bio::io::fastq;
use crate::types::SingeRecordDto;
use crate::writers::generic::{GenericWriter, GenericWriterBackend};

type FastqBackend<'a> = fastq::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastqBackend<'a> {
    fn new(writer: Box<dyn io::Write + Send + 'a>) -> Self {
        fastq::Writer::new(writer)
    }

    fn write_record(&mut self, record: &SingeRecordDto) -> io::Result<()> {
        self.write_record(record)
    }

    fn get_extension() -> &'static str {
        "fastq"
    }
}

pub type SingleWriterFastq<'a> = GenericWriter<'a, FastqBackend<'a>>;


#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::{Arc, Mutex};
    use super::*;
    use crate::types::{RecordDto, SingeRecordDto};
    use crate::writers::Writer;

    #[test]
    fn test_fastq_writer_cursor() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
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
