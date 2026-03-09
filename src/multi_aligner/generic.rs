use std::ops::Range;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignerMode {
    Local,
    LeftOverlap,
    RightOverlap,
    LeftMin,
    RightMin,
}

pub trait GenericRawResult {
    fn validate(
        &self, x: &[u8], y: &[u8], error_rate: f32, min_overlap: Option<usize>, x_offset: usize, y_range: Range<usize>
    ) -> bool;
    fn get_umi_tags(&self, x: &[u8], y: &[u8], x_offset: usize, y_range: Range<usize>) -> Vec<String>;
    fn format(&self, x: &[u8], y: &[u8], x_offset: usize, y_range: Range<usize>) -> String;
}
