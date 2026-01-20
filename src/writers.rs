pub mod fasta;
pub mod fastq;
pub mod generic;

pub trait Writer<R> {
    fn write(&mut self, record: &R) -> std::io::Result<()>;
    fn statistics(&self) -> crate::statistics::WriterStatistics;
}
