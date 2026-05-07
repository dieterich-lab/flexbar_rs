use crate::records::generic_read::GenericRawRead;
use crate::records::bam_read::BamRecordWrapper;

#[derive(Clone, Debug)]
pub enum UnifiedRead {
    Fasta(bio::io::fasta::Record),
    Fastq(bio::io::fastq::Record),
    Bam(BamRecordWrapper),
}

impl UnifiedRead {
    pub fn from_fasta(r: bio::io::fasta::Record) -> Self {
        UnifiedRead::Fasta(r)
    }

    pub fn from_fastq(r: bio::io::fastq::Record) -> Self {
        UnifiedRead::Fastq(r)
    }

    pub fn from_bam(r: rust_htslib::bam::record::Record) -> Self {
        UnifiedRead::Bam(BamRecordWrapper::new(r))
    }

    #[cfg(test)]
    pub fn new_fasta(id: &str, desc: Option<&str>, seq: &[u8]) -> Self {
        UnifiedRead::Fasta(
            bio::io::fasta::Record::with_attrs(id, desc, seq)
        )
    }

    #[cfg(test)]
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
            UnifiedRead::Bam(r) => GenericRawRead::seq(r),
        }
    }

    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::qual(r),
            UnifiedRead::Fastq(r) => GenericRawRead::qual(r),
            UnifiedRead::Bam(r) => GenericRawRead::qual(r),
        }
    }

    #[inline(always)]
    fn id(&self) -> &str {
        match self {
            UnifiedRead::Fasta(r) => GenericRawRead::id(r),
            UnifiedRead::Fastq(r) => GenericRawRead::id(r),
            UnifiedRead::Bam(r) => GenericRawRead::id(r),
        }
    }

    #[inline(always)]
    fn set_seq(&self, seq: &[u8]) -> Self {
        match self {
            UnifiedRead::Fasta(r) => UnifiedRead::Fasta(r.set_seq(seq)),
            UnifiedRead::Fastq(r) => UnifiedRead::Fastq(r.set_seq(seq)),
            UnifiedRead::Bam(r) => UnifiedRead::Bam(r.set_seq(seq)),
        }
    }

    #[inline(always)]
    fn trim(
        &self,
        range: &Range<usize>,
    ) -> Self {
        match self {
            UnifiedRead::Fasta(r) => UnifiedRead::Fasta(r.trim(range)),
            UnifiedRead::Fastq(r) => UnifiedRead::Fastq(r.trim(range)),
            UnifiedRead::Bam(r) => UnifiedRead::Bam(r.trim(range)),
        }
    }
}
