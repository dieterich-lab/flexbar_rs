use std::io;
use bio::io::fasta;
use crate::records::unified::AnnotatedRecord;
use crate::writers::generic::GenericWriterBackend;
use crate::records::generic_read::GenericRawRead;
use crate::writers::utils::get_desc;

pub struct FastaBackend<'a> {
    writer: fasta::Writer<Box<dyn io::Write + Send + 'a>>,
    barcode_mark: bool,
}

impl<'a> FastaBackend<'a> {
    pub fn new(writer: Box<dyn io::Write + Send + 'a>, barcode_mark: bool) -> Self {
        Self {
            writer: fasta::Writer::new(writer),
            barcode_mark,
        }
    }
}

impl<'a> GenericWriterBackend<'a> for FastaBackend<'a> {
    fn write_read(&mut self, record: &AnnotatedRecord, read_index: usize) -> io::Result<()> {
        let read = record.record.reads[read_index].as_ref().expect("Read missing");
        let id = read.id();
        let seq = read.seq();

        let mut desc = get_desc(record, read_index);
        
        if self.barcode_mark {
            if let Some(barcode_id) = &record.barcode_id {
                let bc_mark = format!("BC:{}", barcode_id);
                desc = match desc {
                    Some(d) => Some(format!("{}_{}", d, bc_mark)),
                    None => Some(bc_mark),
                };
            }
        }
        
        let rec = fasta::Record::with_attrs(id, desc.as_deref(), seq);
        self.writer.write_record(&rec)
    }
}



