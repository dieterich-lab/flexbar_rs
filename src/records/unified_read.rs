use std::ops::Range;
use crate::records::generic_read::GenericRawRead;

pub enum UnifiedRead {
    Fasta(bio::io::fasta::Record),
    Fastq(bio::io::fastq::Record),
}

impl UnifiedRead {
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

    fn desc(&self) -> &[u8] {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::desc(r),
            UnifiedRead::Fastq(r) => GenericRawRead::desc(r),
        }
    }

    fn transform(self, range: Range<usize>, extra_desc_parts: Vec<String>) -> Self {
        match self {
            UnifiedRead::Fasta(r) => UnifiedRead::Fasta(r.transform(range, extra_desc_parts)),
            UnifiedRead::Fastq(r) => UnifiedRead::Fastq(r.transform(range, extra_desc_parts)),
        }
    }
}