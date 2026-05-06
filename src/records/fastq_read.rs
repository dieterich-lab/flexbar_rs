use bio::io::fastq::Record;
use crate::records::generic_read::GenericRawRead;


use std::ops::Range;

impl GenericRawRead for Record {
    #[inline(always)]
    fn seq(&self) -> &[u8] {
        self.seq()
    }

    #[inline(always)]
    fn qual(&self) -> Option<&[u8]> {
        Some(self.qual())
    }

    #[inline(always)]
    fn id(&self) -> &str {
        self.id()
    }

    #[inline(always)]
    fn set_seq(&self, seq: &[u8]) -> Self {
        Record::with_attrs(self.id(), self.desc(), seq, self.qual())
    }

    fn trim(&self, range: &Range<usize>) -> Self {
        Record::with_attrs(self.id(), self.desc(), &self.seq()[range.clone()], &self.qual()[range.clone()])
    }
}

