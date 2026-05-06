use std::io;
use rust_htslib::bam::{self, Header, Writer};
use crate::records::unified::{UnifiedRead};
use crate::writers::generic::{GenericWriterBackend};
use crate::records::header::SequenceHeader;

use rust_htslib::bam::record::{Aux, Record};
use crate::types::BamTag;

pub struct BamWriterBackend {
    writer: Writer,
    bam_desc_tag: BamTag,
}

impl BamWriterBackend {
    pub fn new(path: &str, header: &SequenceHeader, bam_desc_tag: BamTag) -> io::Result<Self> {
        let hts_header = match header {
            SequenceHeader::Bam(hv) => Header::from_template(hv),
            SequenceHeader::Default => Header::new(),
        };
        
        // If path is "-", we use stdout
        let writer = if path == "STDOUT" || path == "-" {
            Writer::from_stdout(&hts_header, bam::Format::Bam).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
        } else {
            Writer::from_path(path, &hts_header, bam::Format::Bam).map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
        };

        Ok(Self { writer, bam_desc_tag })
    }

    fn create_bam_record(&mut self, id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) -> Record {
        let mut bam_rec = Record::new();
        
        bam_rec.set(id.as_bytes(), None, seq, qual);
        bam_rec.set_unmapped();
        if let Some(x) = desc {
            let tag: [u8; 2] = self.bam_desc_tag.bytes;
            bam_rec.push_aux(&tag, Aux::String(x)).map_err(|e| io::Error::new(io::ErrorKind::Other, e)).unwrap();
        }
        bam_rec
    }
}

impl<'a> GenericWriterBackend<'a> for BamWriterBackend {
    fn write_read(&mut self, record: &UnifiedRead) -> io::Result<()> {
        let bam_record = match record {
            UnifiedRead::Bam(wrapper) => {
                &wrapper.record
            },
            UnifiedRead::Fasta(r) => {
                let qual = vec![b'I'; r.seq().len()];
                &self.create_bam_record(r.id(), r.desc(), r.seq(), qual.as_slice())
            },
            UnifiedRead::Fastq(r) => {
                &self.create_bam_record(r.id(), r.desc(), r.seq(), r.qual())
            }
        };
        self.writer.write(bam_record).map_err(|e| io::Error::new(io::ErrorKind::Other, e))
    }
}
