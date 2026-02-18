use bio::alignment::Alignment;
use crate::context::SequenceWithId;
use crate::multi_aligner::generic::{AlignerMode, GenericRawResult};


pub enum RawAlignment {
    BioRust(Alignment),
}

pub struct AlignmentResult<'r, 's> {
    pub raw: RawAlignment,
    pub mode: AlignerMode,
    pub offset: usize,
    pub reference: &'r SequenceWithId,
    pub sequence: &'s [u8],
}

impl<'r, 's> AlignmentResult<'r, 's> {
    pub fn new(
        sequence: &'s [u8],
        reference: &'r SequenceWithId,
        mode: AlignerMode,
        offset: usize,
        raw: RawAlignment
    ) -> Self {
        Self {
            raw,
            mode,
            offset,
            reference,
            sequence,
        }
    }

    pub fn score(&self) -> i32 {
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.score,
        }
    }
    pub fn reference_id(&self) -> &str { &self.reference.id }
    pub fn offset(&self) -> usize {
        self.offset
    }
    pub fn mode(&self) -> AlignerMode { self.mode }

    pub fn validate(&self, error_rate: f32, min_overlap: usize) -> bool {
        let y = self.reference.sequence.as_bytes();
        match &self.raw {
            RawAlignment::BioRust(raw) =>
                raw.validate(self.mode, y, error_rate, min_overlap),
        }
    }

    pub fn get_umi_tags(&self) -> Vec<String> {
        let x = self.sequence;
        let y = self.reference.sequence.as_bytes();
        match &self.raw {
            RawAlignment::BioRust(raw) => raw.get_umi_tags(x, y),
        }
    }
    pub fn format(&self) -> String {
        let x = self.sequence;
        let y = self.reference.sequence.as_bytes();
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