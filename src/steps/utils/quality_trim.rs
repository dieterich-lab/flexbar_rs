use log::trace;
use crate::args::Args;
use crate::types::{FileFormat, QualTrimType, QualityType, SingeRecordDto};
use crate::context::Context;
use crate::statistics::WorkerStatistics;
use crate::steps::generic::PerReadAlignmentResults;

pub struct QualityTrim {
    pub trim_type: QualTrimType,
    pub threshold: usize,
    pub win_size: usize,
    pub trace: bool,
}

impl QualityTrim {
    pub fn new(args: &Args, context: &Context) -> Self {
        let format = context.format;
        let threshold = if format == FileFormat::Fastq {
            args.qtrim_threshold
                + match args.qtrim_format {
                    QualityType::Sanger => 33,
                    QualityType::Solexa => 59,
                    QualityType::Illumina => 64,
                }
        } else {
            args.qtrim_threshold
        };

        Self {
            trace: args.trace,
            trim_type: args.qtrim,
            threshold,
            win_size: args.qtrim_win_size,
        }
    }
    
    pub fn trim_read(
        &self,
        read: &SingeRecordDto,
        read_index: Option<usize>,
        per_read: &mut PerReadAlignmentResults,
        _stats: &mut WorkerStatistics
    ) {
        let current_qual = &read.qual()[per_read.read_range.clone()];
        let trim_pos = self.get_trim_pos(current_qual);
        if trim_pos < current_qual.len() {
            per_read.bases_quality_trim += current_qual.len() - trim_pos;
            per_read.read_range.end = per_read.read_range.start + trim_pos;
            if self.trace {
                let read_index_str = match read_index {
                    Some(i) => format!(" (R{})", i+1),
                    None => "".to_string(),
                };
                trace!(
                    "Quality trim on {}{}: Shortened to {} bases",
                    read.id(), read_index_str, trim_pos
                );
            }
        }
    }

    pub fn get_trim_pos(&self, qual: &[u8]) -> usize {
        match self.trim_type {
            QualTrimType::Tail => self.get_tail_trim_pos(qual),
            QualTrimType::Win => self.get_win_trim_pos(qual),
            QualTrimType::Bwa => self.get_bwa_trim_pos(qual),
            QualTrimType::Off => qual.len(),
        }
    }

    fn get_tail_trim_pos(&self, qual: &[u8]) -> usize {
        let mut trim_pos = qual.len();
        while trim_pos > 0 && qual[trim_pos - 1] < self.threshold as u8 {
            trim_pos -= 1;
        }
        trim_pos
    }

    fn get_win_trim_pos(&self, qual: &[u8]) -> usize {
        if qual.len() < self.win_size {
            return self.get_tail_trim_pos(qual);
        }

        let mut current_sum: usize = qual.iter().take(self.win_size).map(|&q| q as usize).sum();
        let threshold_sum = self.threshold * self.win_size;

        if current_sum < threshold_sum {
            return 0;
        }

        for i in 0..=(qual.len() - self.win_size) {
            if current_sum < threshold_sum {
                return i;
            }
            if i + self.win_size < qual.len() {
                current_sum -= qual[i] as usize;
                current_sum += qual[i + self.win_size] as usize;
            }
        }
        qual.len()
    }

    fn get_bwa_trim_pos(&self, qual: &[u8]) -> usize {
        let mut max_score = 0i64;
        let mut current_score = 0i64;
        let mut trim_pos = qual.len();

        for (i, &q) in qual.iter().enumerate().rev() {
            current_score += (self.threshold as i64) - (q as i64);
            if current_score < 0 {
                current_score = 0;
            }
            if current_score > max_score {
                max_score = current_score;
                trim_pos = i;
            }
        }
        trim_pos
    }
}
