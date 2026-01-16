pub mod fasta;
pub mod fastq;
mod generic;

pub trait SequenceReader: Iterator<Item = bio::io::fastq::Record> {
    fn source(&self) -> &str;
    fn compression(&self) -> crate::types::CompressionType;
}