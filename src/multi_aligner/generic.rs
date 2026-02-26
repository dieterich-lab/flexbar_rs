#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignerMode {
    Local,
    LeftOverlap,
    RightOverlap,
}

pub trait GenericRawResult {
    fn validate(
        &self, x: &[u8], y: &[u8], error_rate: f32, min_overlap: Option<usize>, x_offset: usize
    ) -> bool;
    fn get_umi_tags(&self, x: &[u8], y: &[u8], x_offset: usize) -> Vec<String>;
    fn format(&self, x: &[u8], y: &[u8], x_offset: usize) -> String;
}
