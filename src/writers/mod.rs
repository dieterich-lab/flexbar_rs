use crate::types::{RecordDto};

pub mod fasta;
pub mod fastq;
pub mod generic;

pub trait Writer: Send {
    fn write(&mut self, record: &RecordDto, tag: Option<&str>) -> std::io::Result<()>;
     fn statistics(&self) -> crate::statistics::WriterStatistics;
}
