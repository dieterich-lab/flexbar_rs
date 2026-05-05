
use std::ops::Range;

pub trait GenericRawRead {
    fn seq(&self) -> &[u8];
    fn qual(&self) -> Option<&[u8]>;
    fn id(&self) -> &str;
    fn set_seq(&self, seq: &[u8]) -> Self;

    fn finalize(
        &self,
        range: &Range<usize>,
        removal_tags: &Vec<String>,
        umi_tags: &Vec<String>,
        paired_read_dropped: bool,
        bam_desc_tag: &str,
    ) -> Self;
}
