use std::io;
use bio::io::fasta;
use crate::types::SingeRecordDto;
use crate::writers::generic::{GenericWriterBackend};

pub type FastaBackend<'a> = fasta::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastaBackend<'a> {
    fn new(writer: Box<dyn io::Write + Send + 'a>) -> Self {
        fasta::Writer::new(writer)
    }

    fn write_read(&mut self, record: &SingeRecordDto) -> io::Result<()> {
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



