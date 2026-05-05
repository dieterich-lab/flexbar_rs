use std::io::BufRead;
use crate::readers::backends::fasta::FastaBackend;
use crate::readers::backends::fastq::FastqBackend;
use crate::readers::backends::bam::BamBackend;
use crate::readers::backends::generic::GenericReaderBackend;
use crate::records::unified::UnifiedRead;

pub mod fasta;
pub mod fastq;
pub mod bam;
pub mod generic;

pub enum ReaderBackend<Reader: BufRead + Send> {
    Fasta(FastaBackend<Reader>),
    Fastq(FastqBackend<Reader>),
    Bam(BamBackend<Reader>),
}

impl<Reader: BufRead + Send> GenericReaderBackend<Reader> for ReaderBackend<Reader> {
    fn read_record(&mut self) -> Option<UnifiedRead> {
        match self { 
            Self::Fasta(x) =>x.read_record(),
            Self::Fastq(x) => x.read_record(),
            Self::Bam(x) => x.read_record(),
        }
    }
}