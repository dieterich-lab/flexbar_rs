use crate::config::Config;
use crate::types::{QualTrimType, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};

pub struct SingleQualityTrim {
    qtrim: QualTrimType,
    threshold: usize,
    win_size: usize,
}

impl SingleQualityTrim {
    pub fn new(config: &Config) -> Self {
        Self {
            qtrim: config.args.qtrim,
            threshold: config.corrected_qtrim_threshold,
            win_size: config.args.qtrim_win_size,
        }
    }

    pub fn get_trim_pos(&self, qual: &[u8]) -> usize {
        match self.qtrim {
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

impl GenericStep<SingeRecordDto> for SingleQualityTrim {
    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        trace!("Quality trim processing record: {}", rec.record.id());
        let current_qual = &rec.record.qual()[rec.alignment_results.seq_start..rec.alignment_results.seq_end];
        let trim_pos = self.get_trim_pos(current_qual);

        if trim_pos < current_qual.len() {
            let old_len = rec.alignment_results.seq_end - rec.alignment_results.seq_start;
            rec.alignment_results.seq_end = rec.alignment_results.seq_start + trim_pos;
            trace!(
                "Quality trim on {}: Shortened to {} bases - new sequence:\n{:?}",
                rec.record.id(), trim_pos, &rec.record.seq()[
                    rec.alignment_results.seq_start..rec.alignment_results.seq_end
                ]
            );
        } else {
            trace!("Quality trim on {}: NOP", rec.record.id());
        }
        rec
    }
}
