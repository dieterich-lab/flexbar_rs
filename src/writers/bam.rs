use std::io;
use rust_htslib::bam::{self, Header, Writer};
use crate::records::unified::AnnotatedRecord;
use crate::records::unified_read::RawRead;
use crate::writers::generic::{GenericWriterBackend};
use crate::records::generic_read::GenericRawRead;
use crate::records::header::SequenceHeader;
use crate::writers::utils::get_qual;

use rust_htslib::bam::record::{Aux, Record};

pub struct BamWriterBackend {
    writer: Writer,
}

impl BamWriterBackend {
    pub fn new(path: &str, header: &SequenceHeader) -> io::Result<Self> {
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

        Ok(Self { writer })
    }
}

impl<'a> GenericWriterBackend<'a> for BamWriterBackend {
    fn write_read(&mut self, record: &AnnotatedRecord, read_index: usize) -> io::Result<()> {
        let read = record.record.reads[read_index].as_ref().expect("Read missing");
        let mut bam_record = match &read.raw {
            RawRead::Bam(wrapper) => wrapper.record.clone(),
            _ => {
                let qual: Vec<u8> = get_qual(record, read_index)
                    .into_iter()
                    .map(|q| q.saturating_sub(33))
                    .collect();
                let mut bam_rec = Record::new();
                bam_rec.set(read.id().as_bytes(), None, read.seq(), &qual);
                bam_rec.set_unmapped();
                bam_rec
            }
        };

        if !record.umi_tags.is_empty() {
            let umi_str = record.umi_tags.join("_");
            bam_record.push_aux(b"MI", Aux::String(&umi_str)).unwrap();
        }

        let mut co_parts: Vec<String> = Vec::new();
        if let Ok(Aux::String(old)) = bam_record.aux(b"CO") {
            co_parts.push(old.to_string());
        }
        for tag in  record.removal_tags[read_index].iter() {
            co_parts.push(tag.to_string());
        }
        if record.paired_read_dropped {
            co_parts.push("single".to_string());
        }
        if !co_parts.is_empty() {
            let co_str = co_parts.join("_");
            bam_record.push_aux(b"CO", Aux::String(&co_str)).unwrap();
        }
        if let Some(barcode_id) = &record.barcode_id {
            bam_record.push_aux(b"BC", Aux::String(barcode_id)).unwrap();
        }
        self.writer.write(&bam_record).map_err(|e| io::Error::new(io::ErrorKind::Other, e))
    }
}
