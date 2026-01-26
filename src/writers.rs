use crate::types::SingeRecordDto;

pub mod fasta;
pub mod fastq;
pub mod generic;

pub trait Writer: Send {
    fn write(&mut self, record: &SingeRecordDto, tag: Option<&str>) -> std::io::Result<()>;
     fn statistics(&self) -> crate::statistics::WriterStatistics;
}
