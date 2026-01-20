use crate::types::RecordDto;

pub mod fasta;
pub mod fastq;
pub mod generic;

pub trait Writer {
    fn write(&mut self, record: &RecordDto) -> std::io::Result<()>;
    fn statistics(&self) -> crate::statistics::WriterStatistics;
}
