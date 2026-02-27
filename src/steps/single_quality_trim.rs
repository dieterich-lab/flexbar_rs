use crate::args::Args;
use crate::types::{FileFormat, QualTrimType, QualityType, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;

pub struct SingleQualityTrim {}

#[derive(Clone)]
pub struct SingleQualityTrimSetup {
    pub qtrim: QualTrimType,
    pub threshold: usize,
    pub win_size: usize,
    pub trace: bool,
}

impl SingleQualityTrim {
    pub fn get_trim_pos(&self, qual: &[u8], setup: &SingleQualityTrimSetup) -> usize {
        match setup.qtrim {
            QualTrimType::Tail => self.get_tail_trim_pos(qual, setup),
            QualTrimType::Win => self.get_win_trim_pos(qual, setup),
            QualTrimType::Bwa => self.get_bwa_trim_pos(qual, setup),
            QualTrimType::Off => qual.len(),
        }
    }

    fn get_tail_trim_pos(&self, qual: &[u8], setup: &SingleQualityTrimSetup) -> usize {
        let mut trim_pos = qual.len();
        while trim_pos > 0 && qual[trim_pos - 1] < setup.threshold as u8 {
            trim_pos -= 1;
        }
        trim_pos
    }

    fn get_win_trim_pos(&self, qual: &[u8], setup: &SingleQualityTrimSetup) -> usize {
        if qual.len() < setup.win_size {
            return self.get_tail_trim_pos(qual, setup);
        }

        let mut current_sum: usize = qual.iter().take(setup.win_size).map(|&q| q as usize).sum();
        let threshold_sum = setup.threshold * setup.win_size;

        if current_sum < threshold_sum {
            return 0;
        }

        for i in 0..=(qual.len() - setup.win_size) {
            if current_sum < threshold_sum {
                return i;
            }
            if i + setup.win_size < qual.len() {
                current_sum -= qual[i] as usize;
                current_sum += qual[i + setup.win_size] as usize;
            }
        }
        qual.len()
    }

    fn get_bwa_trim_pos(&self, qual: &[u8], setup: &SingleQualityTrimSetup) -> usize {
        let mut max_score = 0i64;
        let mut current_score = 0i64;
        let mut trim_pos = qual.len();

        for (i, &q) in qual.iter().enumerate().rev() {
            current_score += (setup.threshold as i64) - (q as i64);
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

impl<'s> GenericStep<'s, SingeRecordDto> for SingleQualityTrim {
    type Setup<'a> = SingleQualityTrimSetup;

    fn create_setup(args: &Args, context: &Context) -> SingleQualityTrimSetup {
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

        SingleQualityTrimSetup {
            qtrim: args.qtrim,
            threshold,
            win_size: args.qtrim_win_size,
            trace: args.trace,
        }
    }

    fn new(_setup: &SingleQualityTrimSetup) -> Self {
        Self {}
    }

    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>, setup: &SingleQualityTrimSetup) -> RecordInTransformation<SingeRecordDto> {
        let current_qual = &rec.record.qual()[rec.alignment_results.seq_start..rec.alignment_results.seq_end];
        let trim_pos = self.get_trim_pos(current_qual, setup);
        if trim_pos < current_qual.len() {
            rec.alignment_results.seq_end = rec.alignment_results.seq_start + trim_pos;
        }
        if setup.trace {
            if trim_pos < current_qual.len() {
                trace!(
                    "Quality trim on {}: Shortened to {} bases - new sequence:\n{}",
                    rec.record.id(), trim_pos, rec.get_sequence_as_str()
                );
            } else {
                trace!("Quality trim on {}: NOP", rec.record.id());
            }
        }
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;

    #[test]
    fn test_qtrim_threshold_calculation() {
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
        let setup_fastq = SingleQualityTrim::create_setup(&args, &context_fastq);
        assert_eq!(setup_fastq.threshold, 20 + 33);

        // Fasta format should NOT add quality offset
        let mut context_fasta = Context::empty();
        context_fasta.format = FileFormat::Fasta;
        let setup_fasta = SingleQualityTrim::create_setup(&args, &context_fasta);
        assert_eq!(setup_fasta.threshold, 20);

        // Test other quality formats
        let args_solexa = Args::parse_from(&["flexbar_rs", "--qtrim-format", "solexa", "--qtrim-threshold", "20"]);
        let mut context_solexa = Context::empty();
        context_solexa.format = FileFormat::Fastq;
        let setup_solexa = SingleQualityTrim::create_setup(&args_solexa, &context_solexa);
        assert_eq!(setup_solexa.threshold, 20 + 59);

        let args_illumina = Args::parse_from(&["flexbar_rs", "--qtrim-format", "i1.3", "--qtrim-threshold", "20"]);
        let mut context_illumina = Context::empty();
        context_illumina.format = FileFormat::Fastq;
        let setup_illumina = SingleQualityTrim::create_setup(&args_illumina, &context_illumina);
        assert_eq!(setup_illumina.threshold, 20 + 64);
    }
}
