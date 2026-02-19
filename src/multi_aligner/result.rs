use bio::alignment::Alignment;
use crate::multi_aligner::generic::{AlignerMode, GenericRawResult};


pub enum RawAlignment {
    BioRust(Alignment),
}

pub struct MultiAlignerResult {
    pub raw: RawAlignment,
    pub mode: AlignerMode,
    pub offset: usize,
    pub reference_index: usize,
}

impl<'r> MultiAlignerResult {
    pub fn score(&self) -> i32 {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.score,
        }
    }

    pub fn validate(&self, error_rate: f32, y: &[u8], min_overlap: Option<usize>) -> bool {
        match &self.raw {
            RawAlignment::BioRust(raw) =>
                raw.validate(self.mode, y, error_rate, min_overlap),
        }
    }

    pub fn get_umi_tags(&self, x: &[u8], y: &[u8]) -> Vec<String> {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.get_umi_tags(x, y),
        }
    }
    pub fn format(&self, x: &[u8], y: &[u8]) -> String {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.format(x, y),
        }
    }

    pub fn get_start_index(&self) -> usize {
        match &self.raw {
            RawAlignment::BioRust(raw) => self.offset + raw.xstart,
        }
    }
    pub fn get_end_index(&self) -> usize {
        match &self.raw {
            RawAlignment::BioRust(raw) => self.offset + raw.xend,
        }
    }
}