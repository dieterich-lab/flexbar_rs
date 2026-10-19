use crate::records::generic_read::GenericRawRead;
use crate::records::bam_read::BamRecordWrapper;
use std::ops::Range;

#[derive(Clone, Debug)]
pub enum RawRead {
    Fasta(bio::io::fasta::Record),
    Fastq(bio::io::fastq::Record),
    Bam(BamRecordWrapper),
}

#[derive(Clone, Debug)]
pub struct UnifiedRead {
    pub raw: RawRead,
    qual: Option<Vec<u8>>,
}

impl UnifiedRead {
    pub fn from_fasta(r: bio::io::fasta::Record) -> Self {
        UnifiedRead {
            raw: RawRead::Fasta(r),
            qual: None,
        }
    }

    pub fn from_fastq(r: bio::io::fastq::Record) -> Self {
        let qual = Some(r.qual().to_vec());
        UnifiedRead {
            raw: RawRead::Fastq(r),
            qual,
        }
    }

    pub fn from_bam(r: rust_htslib::bam::record::Record) -> Self {
        let qual = Some(r.qual().iter().map(|&q| q + 33).collect());
        UnifiedRead {
            raw: RawRead::Bam(BamRecordWrapper::new(r)),
            qual,
        }
    }

    #[cfg(test)]
    pub fn new_fasta(id: &str, desc: Option<&str>, seq: &[u8]) -> Self {
        Self::from_fasta(
            bio::io::fasta::Record::with_attrs(id, desc, seq)
        )
    }

    #[cfg(test)]
    pub fn new_fastq(id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) -> Self {
        Self::from_fastq(
            bio::io::fastq::Record::with_attrs(id, desc, seq, qual)
        )
    }
}

impl GenericRawRead for UnifiedRead {
    #[inline(always)]
    fn seq(&self) -> &[u8] {
        match &self.raw {
            RawRead::Fasta(r) => GenericRawRead::seq(r),
            RawRead::Fastq(r) => GenericRawRead::seq(r),
            RawRead::Bam(r) => GenericRawRead::seq(r),
        }
    }

    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        self.qual.as_deref()
    }

    #[inline(always)]
    fn id(&self) -> &str {
        match &self.raw {
            RawRead::Fasta(r) => GenericRawRead::id(r),
            RawRead::Fastq(r) => GenericRawRead::id(r),
            RawRead::Bam(r) => GenericRawRead::id(r),
        }
    }

    #[inline(always)]
    fn set_seq(&self, seq: &[u8]) -> Self {
        match &self.raw {
            RawRead::Fasta(r) => UnifiedRead {
                raw: RawRead::Fasta(r.set_seq(seq)),
                qual: self.qual.clone(),
            },
            RawRead::Fastq(r) => UnifiedRead {
                raw: RawRead::Fastq(r.set_seq(seq)),
                qual: self.qual.clone(),
            },
            RawRead::Bam(r) => UnifiedRead {
                raw: RawRead::Bam(r.set_seq(seq)),
                qual: self.qual.clone(),
            },
        }
    }

    #[inline(always)]
    fn trim(
        &self,
        range: &Range<usize>,
    ) -> Self {
        let qual = self.qual.as_ref().map(|q| q[range.clone()].to_vec());
        match &self.raw {
            RawRead::Fasta(r) => UnifiedRead {
                raw: RawRead::Fasta(r.trim(range)),
                qual,
            },
            RawRead::Fastq(r) => UnifiedRead {
                raw: RawRead::Fastq(r.trim(range)),
                qual,
            },
            RawRead::Bam(r) => UnifiedRead {
                raw: RawRead::Bam(r.trim(range)),
                qual,
            },
        }
    }
}
