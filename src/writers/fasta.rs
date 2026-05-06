use std::io;
use bio::io::fasta;
use crate::records::unified::{UnifiedRead};
use crate::writers::generic::{GenericWriterBackend};
use crate::records::generic_read::GenericRawRead;

use rust_htslib::bam::record::Aux;
use crate::types::BamTag;

pub struct FastaBackend<'a> {
    writer: fasta::Writer<Box<dyn io::Write + Send + 'a>>,
    bam_desc_tag: BamTag,
}

impl<'a> FastaBackend<'a> {
    pub fn new(writer: Box<dyn io::Write + Send + 'a>, bam_desc_tag: BamTag) -> Self {
        Self {
            writer: fasta::Writer::new(writer),
            bam_desc_tag,
        }
    }
}

impl<'a> GenericWriterBackend<'a> for FastaBackend<'a> {
    fn write_read(&mut self, read: &UnifiedRead) -> io::Result<()> {
        match read {
            UnifiedRead::Fasta(r) => self.writer.write_record(r),
            UnifiedRead::Fastq(r) => {
                let fasta_rec = fasta::Record::with_attrs(r.id(), r.desc(), r.seq());
                self.writer.write_record(&fasta_rec)
            }
            UnifiedRead::Bam(r) => {
                let id = r.id();
                let seq = r.seq();
                let mut desc = None;
                if let Ok(Aux::String(d)) = r.record.aux(&self.bam_desc_tag.bytes) {
                    desc = Some(d);
                }
                let fasta_rec = fasta::Record::with_attrs(id, desc, seq);
                self.writer.write_record(&fasta_rec)
            }
        }
    }
}



