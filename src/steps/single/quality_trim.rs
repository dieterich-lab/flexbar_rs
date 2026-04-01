use crate::args::Args;
use crate::types::SingeRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::single::SingleAlignmentResults;
use crate::context::Context;
use crate::steps::utils::quality_trim::QualityTrim;

pub struct SingleQualityTrim {
    pub inner: QualityTrim,
}

pub struct SingleQualityTrimState {}

impl GenericStep<SingeRecordDto> for SingleQualityTrim {
    type State = SingleQualityTrimState;
    type AlignmentResult = SingleAlignmentResults;

    fn new(args: &Args, context: &Context) -> Self {
        Self {
            inner: QualityTrim::new(args, context),
        }
    }

    fn create_state(&self) -> Self::State {
        SingleQualityTrimState {}
    }

    fn apply(
        &self,
        mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, 
        _state: &mut Self::State
    ) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        self.inner.trim_read(&rec.record, None, &mut rec.alignment_results.read, &mut rec.stats);
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
        let setup_fastq = SingleQualityTrim::new(&args, &context_fastq);
        assert_eq!(setup_fastq.inner.threshold, 20 + 33);

        // Fasta format should NOT add quality offset
        let mut context_fasta = Context::empty();
        context_fasta.format = FileFormat::Fasta;
        let setup_fasta = SingleQualityTrim::new(&args, &context_fasta);
        assert_eq!(setup_fasta.inner.threshold, 20);

        // Test other quality formats
        let args_solexa = Args::parse_from(&["flexbar_rs", "--qtrim-format", "solexa", "--qtrim-threshold", "20"]);
        let mut context_solexa = Context::empty();
        context_solexa.format = FileFormat::Fastq;
        let setup_solexa = SingleQualityTrim::new(&args_solexa, &context_solexa);
        assert_eq!(setup_solexa.inner.threshold, 20 + 59);

        let args_illumina = Args::parse_from(&["flexbar_rs", "--qtrim-format", "i1.3", "--qtrim-threshold", "20"]);
        let mut context_illumina = Context::empty();
        context_illumina.format = FileFormat::Fastq;
        let setup_illumina = SingleQualityTrim::new(&args_illumina, &context_illumina);
        assert_eq!(setup_illumina.inner.threshold, 20 + 64);
    }
}
