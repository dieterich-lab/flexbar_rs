use serde::{Deserialize, Serialize};
use crate::records::generic_read::GenericRawRead;
use std::ops::Range;

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

    fn transform(self, range: Range<usize>, extra_desc_parts: Vec<String>, replace_with_dummy: bool) -> Self {
        match self {
            UnifiedRead::Fasta(r) => UnifiedRead::Fasta(r.transform(range, extra_desc_parts, replace_with_dummy)),
            UnifiedRead::Fastq(r) => UnifiedRead::Fastq(r.transform(range, extra_desc_parts, replace_with_dummy)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unified_fasta_transform() {
        let read = UnifiedRead::new_fasta("id1", Some("desc1"), b"ATGCATGC");
        let transformed = read.transform(2..6, vec!["extra".to_string()], false);
        
        match transformed {
            UnifiedRead::Fasta(r) => {
                assert_eq!(r.id(), "id1");
                assert_eq!(r.desc(), Some("desc1_extra"));
                assert_eq!(r.seq(), b"GCAT");
            },
            _ => panic!("Expected Fasta record"),
        }
    }

    #[test]
    fn test_unified_fastq_transform() {
        let read = UnifiedRead::new_fastq("id1", Some("desc1"), b"ATGCATGC", b"IIIIIIII");
        let transformed = read.transform(2..6, vec!["extra".to_string()], false);
        
        match transformed {
            UnifiedRead::Fastq(r) => {
                assert_eq!(r.id(), "id1");
                assert_eq!(r.desc(), Some("desc1_extra"));
                assert_eq!(r.seq(), b"GCAT");
                assert_eq!(r.qual(), b"IIII");
            },
            _ => panic!("Expected Fastq record"),
        }
    }
}