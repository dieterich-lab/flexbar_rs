use std::ops::Range;
use bio::io::fastq::Record;
use crate::records::generic_read::GenericRawRead;


impl GenericRawRead for Record {
    fn seq(&self) -> &[u8] {
        self.seq()
    }

    fn qual(&self) -> Option<&[u8]> {
        Some(self.qual())
    }

    fn id(&self) -> &str {
        self.id()
    }

    fn desc(&self) -> Option<&str> {
        self.desc()
    }

    fn desc_str(&self) -> Option<&str> {
        self.desc()
    }

    fn transform(self, range: Range<usize>, extra_desc_parts: Vec<String>, replace_with_dummy: bool) -> Self {
        let id = self.id().to_string();
        let mut desc_parts: Vec<String> = match self.desc() {
            Some(desc) => vec![desc.to_string()],
            _ => vec![],
        };
        desc_parts.extend(extra_desc_parts);
        let joined_desc = desc_parts.join("_");
        let desc = match joined_desc.len() > 0 {
            true => Some(joined_desc.as_ref()),
            false => None,
        };
        let (seq, qual) = if replace_with_dummy {
            (b"N".to_vec(), b"I".to_vec())
        } else {
            (self.seq()[range.clone()].to_vec(), self.qual()[range].to_vec())
        };
        Record::with_attrs(&id, desc, &seq, &qual)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;

    #[test]
    fn test_fastq_transform() {
        let record = Record::with_attrs("id1", Some("desc1"), b"ATGCATGC", b"IIIIIIII");
        let transformed = record.transform(2..6, vec!["extra".to_string()], false);
        
        assert_eq!(transformed.id(), "id1");
        assert_eq!(transformed.desc(), Some("desc1_extra"));
        assert_eq!(transformed.seq(), b"GCAT");
        assert_eq!(transformed.qual(), b"IIII");
    }

    #[test]
    fn test_fastq_transform_dummy() {
        let record = Record::with_attrs("id1", Some("desc1"), b"ATGCATGC", b"IIIIIIII");
        let transformed = record.transform(2..6, vec![], true);
        
        assert_eq!(transformed.id(), "id1");
        assert_eq!(transformed.desc(), Some("desc1"));
        assert_eq!(transformed.seq(), b"N");
        assert_eq!(transformed.qual(), b"I");
    }

    #[test]
    fn test_fastq_transform_no_desc() {
        let record = Record::with_attrs("id1", None, b"ATGCATGC", b"IIIIIIII");
        let transformed = record.transform(0..4, vec!["extra".to_string()], false);
        
        assert_eq!(transformed.id(), "id1");
        assert_eq!(transformed.desc(), Some("extra"));
        assert_eq!(transformed.seq(), b"ATGC");
        assert_eq!(transformed.qual(), b"IIII");
    }
}
