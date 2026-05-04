use serde::{Deserialize, Serialize};
use crate::records::generic_read::GenericRawRead;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum UnifiedRead {
    Fasta(bio::io::fasta::Record),
    Fastq(bio::io::fastq::Record),
}

impl UnifiedRead {
    pub fn from_fasta(r: bio::io::fasta::Record) -> Self {
        UnifiedRead::Fasta(r)
    }

    pub fn from_fastq(r: bio::io::fastq::Record) -> Self {
        UnifiedRead::Fastq(r)
    }

    pub fn new_fasta(id: &str, desc: Option<&str>, seq: &[u8]) -> Self {
        UnifiedRead::Fasta(
            bio::io::fasta::Record::with_attrs(id, desc, seq)
        )
    }

    pub fn new_fastq(id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) -> Self {
        UnifiedRead::Fastq(
            bio::io::fastq::Record::with_attrs(id, desc, seq, qual)
        )
    }
}

impl GenericRawRead for UnifiedRead {
    fn seq(&self) -> &[u8] {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::seq(r),
            UnifiedRead::Fastq(r) => GenericRawRead::seq(r),
        }
    }

    fn qual(&self) -> Option<&[u8]> {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::qual(r),
            UnifiedRead::Fastq(r) => GenericRawRead::qual(r),
        }
    }

    fn id(&self) -> &str {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::id(r),
            UnifiedRead::Fastq(r) => GenericRawRead::id(r),
        }
    }

    fn desc(&self) -> Option<&str> {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::desc(r),
            UnifiedRead::Fastq(r) => GenericRawRead::desc(r),
        }
    }

    fn desc_str(&self) -> Option<&str> {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::desc_str(r),
            UnifiedRead::Fastq(r) => GenericRawRead::desc_str(r),
        }
    }
}
