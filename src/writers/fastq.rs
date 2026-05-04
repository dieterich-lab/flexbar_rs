use std::io;
use bio::io::fastq;
use crate::records::unified::{UnifiedRead, GenericRawRead};
use crate::writers::generic::{GenericWriterBackend};

pub type FastqBackend<'a> = fastq::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastqBackend<'a> {
    fn write_read(&mut self, record: &UnifiedRead) -> io::Result<()> {
        self.write_record(
            &fastq::Record::with_attrs(
                record.id(),
                record.desc(),
                record.seq(),
                record.qual().unwrap(),
            )
        )
    }
}

