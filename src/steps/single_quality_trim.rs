use crate::config::Config;
use crate::types::{QualTrimType, RecordDto};
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

    pub fn trim(&self, record: &RecordDto) -> RecordDto {
        let trim_pos = match self.qtrim {
            QualTrimType::Tail => self.get_tail_trim_pos(record.qual()),
            QualTrimType::Win => self.get_win_trim_pos(record.qual()),
            QualTrimType::Bwa => self.get_bwa_trim_pos(record.qual()),
            QualTrimType::Off => unreachable!(),
        };

        if trim_pos < record.seq().len() {
            let mut new_seq = record.seq().to_vec();
            let mut new_qual = record.qual().to_vec();
            new_seq.truncate(trim_pos);
            new_qual.truncate(trim_pos);
            RecordDto::with_attrs(record.id(), record.desc(), &new_seq, &new_qual)
        } else {
            record.clone()
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

impl GenericStep<RecordDto> for SingleQualityTrim {
    fn apply(&self, mut rec: RecordInTransformation<RecordDto>) -> RecordInTransformation<RecordDto> {
        if let Some(record) = rec.record {
            rec.record = Some(self.trim(&record));
        }
        rec
    }
}
