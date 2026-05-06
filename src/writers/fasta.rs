use std::io;
use bio::io::fasta;
use crate::records::unified::AnnotatedRecord;
use crate::writers::generic::GenericWriterBackend;
use crate::records::generic_read::GenericRawRead;
use crate::writers::utils::get_desc;

pub struct FastaBackend<'a> {
    writer: fasta::Writer<Box<dyn io::Write + Send + 'a>>,
}

impl<'a> FastaBackend<'a> {
    pub fn new(writer: Box<dyn io::Write + Send + 'a>) -> Self {
        Self {
            writer: fasta::Writer::new(writer),
        }
    }
}

impl<'a> GenericWriterBackend<'a> for FastaBackend<'a> {
    fn write_read(&mut self, record: &AnnotatedRecord, read_index: usize) -> io::Result<()> {
        let read = record.record.reads[read_index].as_ref().expect("Read missing");
        let id = read.id();
        let seq = read.seq();

        let desc = get_desc(record, read_index);
        let rec = fasta::Record::with_attrs(id, desc.as_deref(), seq);
        self.writer.write_record(&rec)
    }
}



