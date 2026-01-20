pub mod fasta;
pub mod fastq;
mod generic;

pub trait SequenceReader: Iterator<Item = bio::io::fastq::Record> {
    fn compression(&self) -> crate::types::CompressionType;
}