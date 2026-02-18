use crate::multi_aligner::result::AlignmentResult;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignerMode {
    Local,
    LeftOverlap,
    RightOverlap,
}

pub trait GenericRawResult {
    fn validate(&self, mode: AlignerMode, y: &[u8], error_rate: f32, min_overlap: usize) -> bool;
    fn get_umi_tags(&self, x: &[u8], y: &[u8]) -> Vec<String>;
    fn format(&self, x: &[u8], y: &[u8]) -> String;
}

pub trait GenericAligner<'a> {
    fn align<'s>(&mut self, sequence: &'s [u8], reference_index: usize)
        -> AlignmentResult<'a, 's>;
}


