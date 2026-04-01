use std::ops::Range;
use crate::context::SequenceWithId;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignerMode {
    Local,
    LeftOverlap,
    RightOverlap,
    LeftMin,
    RightMin,
}

#[derive(Clone)]
pub struct ReferenceWithRange {
    pub reference: SequenceWithId,
    pub range: Range<usize>,
    pub aligner_range: Range<usize>,
}

impl ReferenceWithRange {
    pub fn get_aligner_sequence(&self) -> &[u8] {
        &self.reference.sequence.as_bytes()[self.aligner_range.clone()]
    }

    pub fn get_aligner_start_offset(&self) -> usize {
        self.aligner_range.start - self.range.start
    }

    pub fn get_aligner_end_offset(&self) -> usize {
        self.range.end - self.aligner_range.end
    }
}

pub trait GenericRawResult {
    fn validate(
        &self, x: &[u8], y: &[u8], allowed_errors: usize, check_length: usize, x_offset: usize, y_range: Range<usize>
    ) -> bool;
    fn get_umi_tags(&self, x: &[u8], y: &[u8], x_offset: usize, y_range: Range<usize>) -> Vec<String>;
    fn format(&self, x: &[u8], y: &[u8], x_offset: usize, y_range: Range<usize>) -> String;
}
