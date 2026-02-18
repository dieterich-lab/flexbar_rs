use crate::multi_aligner::result::MultiAlignerResult;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignerMode {
    Local,
    LeftOverlap,
    RightOverlap,
}

pub trait GenericRawResult {
    fn validate(&self, mode: AlignerMode, y: &[u8], error_rate: f32, min_overlap: Option<usize>)
        -> bool;
    fn get_umi_tags(&self, x: &[u8], y: &[u8]) -> Vec<String>;
    fn format(&self, x: &[u8], y: &[u8]) -> String;
}

pub trait GenericAligner {
    fn align(&mut self, sequence: &[u8], reference_index: usize)
        -> MultiAlignerResult;
}


