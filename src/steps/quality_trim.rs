use crate::args::Args;
use crate::types::{RecordDto, RunType};
use crate::steps::generic::{GenericStep, RecordInTransformation, CombinedAligmentResults};
use crate::context::Context;
use crate::steps::utils::quality_trim::QualityTrim;

pub struct QualityTrimStep {
    pub inner: QualityTrim,
}

pub struct QualityTrimState {}

impl GenericStep<RecordDto> for QualityTrimStep {
    type State = QualityTrimState;
    type AlignmentResult = CombinedAligmentResults;

    fn new(args: &Args, context: &Context) -> Self {
        match context.run_type {
            RunType::SingleBarcoded | RunType::PairedBarcoded => {
                panic!("Barcoded run types not implemented yet");
            }
            _ => {}
        }
        Self {
            inner: QualityTrim::new(args, context),
        }
    }

    fn create_state(&self) -> Self::State {
        QualityTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<RecordDto, CombinedAligmentResults>, _state: &mut Self::State) -> RecordInTransformation<RecordDto, CombinedAligmentResults> {
        match &rec.record {
            RecordDto::Single(s) => {
                self.inner.trim_read(s, None, rec.alignment_results.per_read[0].as_mut().unwrap(), &mut rec.stats);
            }
            RecordDto::Paired(p) => {
                for i in 0..2 {
                    self.inner.trim_read(&p.reads[i], Some(i), rec.alignment_results.per_read[i].as_mut().unwrap(), &mut rec.stats);
                }
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
    use crate::types::FileFormat;

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
        let setup_fastq = QualityTrimStep::new(&args, &context_fastq);
        assert_eq!(setup_fastq.inner.threshold, 20 + 33);

        // Fasta format should NOT add quality offset
        let mut context_fasta = Context::empty();
        context_fasta.format = FileFormat::Fasta;
        context_fasta.run_type = RunType::Single;
        let setup_fasta = QualityTrimStep::new(&args, &context_fasta);
        assert_eq!(setup_fasta.inner.threshold, 20);

        // Test other quality formats
        let args_solexa = Args::parse_from(&["flexbar_rs", "--qtrim-format", "solexa", "--qtrim-threshold", "20"]);
        let mut context_solexa = Context::empty();
        context_solexa.format = FileFormat::Fastq;
        context_solexa.run_type = RunType::Single;
        let setup_solexa = QualityTrimStep::new(&args_solexa, &context_solexa);
        assert_eq!(setup_solexa.inner.threshold, 20 + 59);

        let args_illumina = Args::parse_from(&["flexbar_rs", "--qtrim-format", "i1.3", "--qtrim-threshold", "20"]);
        let mut context_illumina = Context::empty();
        context_illumina.format = FileFormat::Fastq;
        context_illumina.run_type = RunType::Single;
        let setup_illumina = QualityTrimStep::new(&args_illumina, &context_illumina);
        assert_eq!(setup_illumina.inner.threshold, 20 + 64);
    }
}
