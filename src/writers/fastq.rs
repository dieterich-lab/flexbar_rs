use std::io;
use bio::io::fastq;
use crate::records::unified::{UnifiedRead};
use crate::writers::generic::{GenericWriterBackend};

pub type FastqBackend<'a> = fastq::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastqBackend<'a> {
    fn write_read(&mut self, record: &UnifiedRead) -> io::Result<()> {
        match record {
            UnifiedRead::Fastq(r) => self.write_record(r),
            _ => panic!("expected fastq read"),
        }
    }
}

