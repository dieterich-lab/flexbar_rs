use std::io;
use bio::io::fasta;
use crate::records::unified::{UnifiedRead};
use crate::writers::generic::{GenericWriterBackend};

pub type FastaBackend<'a> = fasta::Writer<Box<dyn io::Write + Send + 'a>>;

impl<'a> GenericWriterBackend<'a> for FastaBackend<'a> {
    fn write_read(&mut self, read: &UnifiedRead) -> io::Result<()> {
        match read {
            UnifiedRead::Fasta(r) => self.write_record(r),
            _ => panic!("expected fasta read"),
        }
    }
}



