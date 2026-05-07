use std::io;
use bio::io::fastq;
use crate::records::unified::AnnotatedRecord;
use crate::writers::generic::GenericWriterBackend;
use crate::records::generic_read::GenericRawRead;
use crate::writers::utils::{get_desc, get_qual};

pub struct FastqBackend<'a> {
    writer: fastq::Writer<Box<dyn io::Write + Send + 'a>>,
    barcode_mark: bool,
}

impl<'a> FastqBackend<'a> {
    pub fn new(writer: Box<dyn io::Write + Send + 'a>, barcode_mark: bool) -> Self {
        Self {
            writer: fastq::Writer::new(writer),
            barcode_mark,
        }
    }
}

impl<'a> GenericWriterBackend<'a> for FastqBackend<'a> {
    fn write_read(&mut self, record: &AnnotatedRecord, read_index: usize) -> io::Result<()> {
        let read = record.record.reads[read_index].as_ref().expect("Read missing");
        let id = read.id();
        let seq = read.seq();
        let qual_vec = get_qual(record, read_index);

        let mut final_desc = get_desc(record, read_index);

        if self.barcode_mark {
            if let Some(barcode_id) = &record.barcode_id {
                let bc_mark = format!("BC:{}", barcode_id);
                final_desc = match final_desc {
                    Some(d) => Some(format!("{}_{}", d, bc_mark)),
                    None => Some(bc_mark),
                };
            }
        }

        let fastq_rec = fastq::Record::with_attrs(id, final_desc.as_deref(), seq, &qual_vec);
        self.writer.write_record(&fastq_rec)
    }
}

