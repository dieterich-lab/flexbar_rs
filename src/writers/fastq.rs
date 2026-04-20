use std::io;
use bio::io::fastq;
use crate::types::SingeRecordDto;
use crate::writers::generic::{GenericWriterBackend};

pub type FastqBackend<'a> = fastq::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastqBackend<'a> {
    fn write_read(&mut self, record: &SingeRecordDto) -> io::Result<()> {
        self.write_record(record)
    }
}

