use std::fmt::{Display, Formatter};
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

pub enum ReaderBackendError<Reader: BufRead + Send> {
    Fasta(<FastaBackend<Reader> as GenericReaderBackend<Reader>>::Error),
    Fastq(<FastqBackend<Reader> as GenericReaderBackend<Reader>>::Error),
}

impl<Reader: BufRead + Send> Display for ReaderBackendError<Reader> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fasta(x) => x.fmt(f),
            Self::Fastq(x) => x.fmt(f),
        }
    }
}
impl<Reader: BufRead + Send> GenericReaderBackend<Reader> for ReaderBackend<Reader> {
    type Error = ReaderBackendError<Reader>;

    fn read_record(&mut self) -> Result<Option<SingeRecordDto>, Self::Error> {
        match self { 
            Self::Fasta(x) => 
                match x.read_record() {
                    Ok(x) => Ok(x),
                    Err(e) => Err(ReaderBackendError::Fasta(e)),
                },
            Self::Fastq(x) => 
                match x.read_record() {
                    Ok(x) => Ok(x),
                    Err(e) => Err(ReaderBackendError::Fastq(e)),
                }
        }
    }
}