use std::io;
use rust_htslib::bam::{self, Header, Writer};
use crate::records::unified::{UnifiedRead};
use crate::writers::generic::{GenericWriterBackend};
use crate::records::header::SequenceHeader;

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
    fn write_read(&mut self, record: &UnifiedRead) -> io::Result<()> {
        match record {
            UnifiedRead::Bam(wrapper) => {
                self.writer.write(&wrapper.record).map_err(|e| io::Error::new(io::ErrorKind::Other, e))
            },
            _ => {
                // We might want to convert other reads to BAM records if needed, 
                // but for now, we follow the pattern of other backends.
                panic!("Expected BAM read for BAM writer");
            }
        }
    }
}
