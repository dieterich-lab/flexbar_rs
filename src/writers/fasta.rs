use std::io;
use bio::io::fasta;
use crate::records::unified::{UnifiedRead, GenericRawRead};
use crate::writers::generic::{GenericWriterBackend};

pub type FastaBackend<'a> = fasta::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastaBackend<'a> {
    fn write_read(&mut self, read: &UnifiedRead) -> io::Result<()> {
        let fasta_record = fasta::Record::with_attrs(
            read.id(),
            read.desc(),
            read.seq(),
        );
        self.write_record(&fasta_record)
    }
}



