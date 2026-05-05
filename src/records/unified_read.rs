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

use std::ops::Range;

impl GenericRawRead for UnifiedRead {
    #[inline(always)]
    fn seq(&self) -> &[u8] {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::seq(r),
            UnifiedRead::Fastq(r) => GenericRawRead::seq(r),
        }
    }

    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::qual(r),
            UnifiedRead::Fastq(r) => GenericRawRead::qual(r),
        }
    }

    #[inline(always)]
    fn id(&self) -> &str {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::id(r),
            UnifiedRead::Fastq(r) => GenericRawRead::id(r),
        }
    }

    #[inline(always)]
    fn set_seq(&self, seq: &[u8]) -> Self {
        match self {
            UnifiedRead::Fasta(r) => UnifiedRead::Fasta(r.set_seq(seq)),
            UnifiedRead::Fastq(r) => UnifiedRead::Fastq(r.set_seq(seq)),
        }
    }

    #[inline(always)]
    fn finalize(
        &self,
        range: &Range<usize>,
        removal_tags: &Vec<String>,
        umi_tags: &Vec<String>,
        paired_read_dropped: bool,
    ) -> Self {
        match self {
            UnifiedRead::Fasta(r) => UnifiedRead::Fasta(r.finalize(range, removal_tags, umi_tags, paired_read_dropped)),
            UnifiedRead::Fastq(r) => UnifiedRead::Fastq(r.finalize(range, removal_tags, umi_tags, paired_read_dropped)),
        }
    }
}
