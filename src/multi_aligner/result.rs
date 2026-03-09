use bio::alignment::Alignment;
use crate::multi_aligner::generic::{GenericRawResult};


pub enum RawAlignment {
    BioRust(Alignment),
}

pub struct MultiAlignerResult {
    pub raw: RawAlignment,
    pub sequence_offset: usize,
    pub reference_index: usize,
}

impl<'r> MultiAlignerResult {
    pub fn score(&self) -> i32 {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.score,
        }
    }

    pub fn validate(
        &self, error_rate: f32, x: &[u8], y: &[u8], min_overlap: Option<usize>, x_offset: usize
    ) -> bool {
        match &self.raw {
            RawAlignment::BioRust(raw) =>
                raw.validate(x, y, error_rate, min_overlap, x_offset),
        }
    }

    pub fn get_umi_tags(&self, x: &[u8], y: &[u8]) -> Vec<String> {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.get_umi_tags(x, y, self.sequence_offset),
        }
    }
    pub fn format(&self, x: &[u8], y: &[u8]) -> String {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.format(x, y, self.sequence_offset),
        }
    }

    pub fn get_start_index(&self) -> usize {
        match &self.raw {
            RawAlignment::BioRust(raw) => self.sequence_offset + raw.xstart,
        }
    }
    pub fn get_end_index(&self) -> usize {
        match &self.raw {
            RawAlignment::BioRust(raw) => self.sequence_offset + raw.xend,
        }
    }
}