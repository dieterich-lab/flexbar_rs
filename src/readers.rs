use crate::types::RecordDto;

pub mod fasta;
pub mod fastq;
mod generic;

pub trait SequenceReader: Iterator<Item = RecordDto> {
    fn compression(&self) -> crate::types::CompressionType;
    fn statistics(&self) -> crate::statistics::ReaderStatistics;
}