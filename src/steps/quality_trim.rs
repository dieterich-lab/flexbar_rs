use crate::args::Args;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;
use log::trace;
use crate::types::{FileFormat, QualTrimType, QualityType};

pub struct QualityTrim {
    pub trim_type: QualTrimType,
    pub threshold: usize,
    pub win_size: usize,
    pub trace: bool,
}

impl GenericStep for QualityTrim {
    type State = ();

    fn new(args: &Args, context: &Context) -> Self {
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

    fn create_state(&self) -> Self::State { () }

    fn apply(&self, mut rec: RecordInTransformation, _state: &mut Self::State) -> RecordInTransformation {
        for i in 0..2 {
            if rec.record.reads[i].is_none() || rec.is_dropped_read(i) { continue }
            self.trim_read(&mut rec, i);
        }
        rec
    }
}

impl QualityTrim {
    fn trim_read(&self, rec: &mut RecordInTransformation, read_index: usize) {
        let current_qual = rec.get_qual(read_index)
            .unwrap_or_else(|| panic!("Quality trim can't be done on reads without quality information."));
        let trim_pos = self.get_trim_pos(current_qual);
        let current_length = current_qual.len();
        if trim_pos < current_length {
            let per_read = rec.alignment_results.per_read
                .get_mut(read_index).unwrap().as_mut().unwrap();
            per_read.bases_quality_trim += current_length - trim_pos;
            per_read.read_range.end = per_read.read_range.start + trim_pos;
            if self.trace {
                trace!(
                    "Quality trim on {}: Shortened to {} bases",
                    rec.read_id(read_index), trim_pos
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;
    use crate::types::{FileFormat, RunType};

    #[test]
    fn test_qtrim_threshold_calculation_single() {
        let args = Args::parse_from(&[
            "flexbar_rs",
            "--reads", "test.fq",
            "--qtrim", "TAIL",
            "--qtrim-format", "sanger",
            "--qtrim-threshold", "20"
        ]);
        
        // Fastq format should add quality offset
        let mut context_fastq = Context::empty();
        context_fastq.format = FileFormat::Fastq;
        context_fastq.run_type = RunType::Single;
        let setup_fastq = QualityTrim::new(&args, &context_fastq);
        assert_eq!(setup_fastq.threshold, 20 + 33);

        // Fasta format should NOT add quality offset
        let mut context_fasta = Context::empty();
        context_fasta.format = FileFormat::Fasta;
        context_fasta.run_type = RunType::Single;
        let setup_fasta = QualityTrim::new(&args, &context_fasta);
        assert_eq!(setup_fasta.threshold, 20);

        // Test other quality formats
        let args_solexa = Args::parse_from(&["flexbar_rs", "--qtrim-format", "solexa", "--qtrim-threshold", "20"]);
        let mut context_solexa = Context::empty();
        context_solexa.format = FileFormat::Fastq;
        context_solexa.run_type = RunType::Single;
        let setup_solexa = QualityTrim::new(&args_solexa, &context_solexa);
        assert_eq!(setup_solexa.threshold, 20 + 59);

        let args_illumina = Args::parse_from(&["flexbar_rs", "--qtrim-format", "i1.3", "--qtrim-threshold", "20"]);
        let mut context_illumina = Context::empty();
        context_illumina.format = FileFormat::Fastq;
        context_illumina.run_type = RunType::Single;
        let setup_illumina = QualityTrim::new(&args_illumina, &context_illumina);
        assert_eq!(setup_illumina.threshold, 20 + 64);
    }
}
