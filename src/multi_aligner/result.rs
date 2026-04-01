use std::ops::Range;
use bio::alignment::Alignment;
use crate::multi_aligner::generic::{GenericRawResult};


#[derive(Clone)]
pub enum RawAlignment {
    BioRust(Alignment),
}

#[derive(Clone)]
pub struct MultiAlignerResult {
    pub raw: RawAlignment,
    pub sequence_offset: usize,
    pub reference_index: usize,
    pub reference_range: Range<usize>,
    pub stripped_leading_n_count: usize,
    pub stripped_trailing_n_count: usize,
    pub may_need_shrink: bool,
}

impl<'r> MultiAlignerResult {
    pub fn score(&self) -> i32 {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.score,
        }
    }

    pub fn validate(
        &self, allowed_errors: usize, check_length: usize, y: &[u8], min_overlap: Option<usize>, _x_offset: usize
    ) -> bool {
        if let Some(requested_min_overlap) = min_overlap {
            let corrected_min_overlap = requested_min_overlap.saturating_sub(
                self.stripped_leading_n_count + self.stripped_trailing_n_count
            );
            if corrected_min_overlap > check_length { return false }
        }
        match &self.raw {
            RawAlignment::BioRust(raw) =>
                raw.validate(
                    &y[self.reference_range.clone()], allowed_errors, check_length
                ),
        }
    }

    pub fn get_umi_tags(&self, x: &[u8], y: &[u8]) -> Vec<String> {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.get_umi_tags(x, y, self.sequence_offset, self.reference_range.clone()),
        }
    }
    pub fn format(&self, x: &[u8], y: &[u8]) -> String {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.format(x, y, self.sequence_offset, self.reference_range.clone()),
        }
    }

    pub fn get_start_index(&self) -> usize {
        match &self.raw {
            RawAlignment::BioRust(raw)
                => self.sequence_offset + raw.xstart - self.stripped_leading_n_count,
        }
    }
    pub fn get_end_index(&self) -> usize {
        match &self.raw {
            RawAlignment::BioRust(raw)
                => self.sequence_offset + raw.xend + self.stripped_trailing_n_count,
        }
    }
}