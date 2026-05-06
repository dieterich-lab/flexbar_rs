use std::io;
use bio::io::fastq;
use crate::records::unified::{UnifiedRead};
use crate::writers::generic::{GenericWriterBackend};
use crate::records::generic_read::GenericRawRead;

use rust_htslib::bam::record::Aux;
use crate::types::BamTag;

pub struct FastqBackend<'a> {
    writer: fastq::Writer<Box<dyn io::Write + Send + 'a>>,
    bam_desc_tag: BamTag,
}

impl<'a> FastqBackend<'a> {
    pub fn new(writer: Box<dyn io::Write + Send + 'a>, bam_desc_tag: BamTag) -> Self {
        Self {
            writer: fastq::Writer::new(writer),
            bam_desc_tag,
        }
    }
}

impl<'a> GenericWriterBackend<'a> for FastqBackend<'a> {
    fn write_read(&mut self, record: &UnifiedRead) -> io::Result<()> {
        match record {
            UnifiedRead::Fastq(r) => self.writer.write_record(r),
            UnifiedRead::Fasta(r) => {
                let qual = vec![b'I'; r.seq().len()];
                let fastq_rec = fastq::Record::with_attrs(r.id(), r.desc(), r.seq(), &qual);
                self.writer.write_record(&fastq_rec)
            }
            UnifiedRead::Bam(r) => {
                let id = r.id();
                let seq = r.seq();
                let qual = r.qual().map(|q| q.to_vec()).unwrap_or_else(|| vec![b'I'; seq.len()]);
                let desc = match r.record.aux(&self.bam_desc_tag.bytes) {
                    Ok(Aux::String(d)) => Some(d),
                    _ => None,
                };
                let fastq_rec = fastq::Record::with_attrs(id, desc, seq, &qual);
                self.writer.write_record(&fastq_rec)
            }
        }
    }
}

