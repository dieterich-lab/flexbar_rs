use std::io::BufRead;
use crate::readers::backends::fasta::FastaBackend;
use crate::readers::backends::fastq::FastqBackend;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::types::SingeRecordDto;

pub mod fasta;
pub mod fastq;
pub mod generic;

pub enum ReaderBackend<Reader: BufRead + Send> {
    Fasta(FastaBackend<Reader>),
    Fastq(FastqBackend<Reader>),
}

impl<Reader: BufRead + Send> GenericReaderBackend<Reader> for ReaderBackend<Reader> {
    fn read_record(&mut self) -> Option<SingeRecordDto> {
        match self { 
            Self::Fasta(x) =>x.read_record(),
            Self::Fastq(x) => x.read_record()
        }
    }
}