use rust_htslib::bam::record::Record;
use crate::records::generic_read::GenericRawRead;
use std::ops::Range;
use std::cell::OnceCell;

#[derive(Clone, Debug)]
pub struct BamRecordWrapper {
    pub record: Record,
    pub seq_cache: OnceCell<Vec<u8>>,
}

impl BamRecordWrapper {
    pub fn new(record: Record) -> Self {
        Self {
            record,
            seq_cache: OnceCell::new(),
        }
    }
}

impl GenericRawRead for BamRecordWrapper {
    #[inline(always)]
    fn seq(&self) -> &[u8] {
        self.seq_cache.get_or_init(|| self.record.seq().as_bytes())
    }

    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        Some(self.record.qual())
    }

    #[inline(always)]
    fn id(&self) -> &str {
        std::str::from_utf8(self.record.qname()).unwrap_or("unknown")
    }

    #[inline(always)]
    fn set_seq(&self, seq: &[u8]) -> Self {
        let mut new_record = self.record.clone();
        let cigar = self.record.cigar().take();
        new_record.set(self.record.qname(), Some(&cigar), seq, self.record.qual());
        Self::new(new_record)
    }

    fn finalize(
        &self,
        range: &Range<usize>,
        _removal_tags: &Vec<String>,
        _umi_tags: &Vec<String>,
        _paired_read_dropped: bool,
    ) -> Self {
        let mut new_record = Record::new();
        let seq = self.record.seq().as_bytes();
        let qual = self.record.qual();
        
        let new_seq = &seq[range.clone()];
        let new_qual = &qual[range.clone()];
        
        new_record.set(self.record.qname(), None, new_seq, new_qual);
        Self::new(new_record)
    }
}
