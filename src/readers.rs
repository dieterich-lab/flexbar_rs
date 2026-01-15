pub mod fasta;
pub mod fastq;
mod generic;

pub trait SequenceReader {
    fn source(&self) -> &str;
    fn compression(&self) -> crate::types::CompressionType;
}