use crate::args::Args;
use crate::types::SingeRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::single::SingleAlignmentResults;
use crate::context::Context;
use crate::steps::utils::pre_trim::apply_pre_trim;

pub struct SinglePreTrim {
    pub pre_trim_left: usize,
    pub pre_trim_right: usize,
}

pub struct SinglePreTrimState {}

impl GenericStep<SingeRecordDto> for SinglePreTrim {
    type State = SinglePreTrimState;
    type AlignmentResult = SingleAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            pre_trim_left: args.pre_trim_left,
            pre_trim_right: args.pre_trim_right,
        }
    }

    fn create_state(&self) -> Self::State {
        SinglePreTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        apply_pre_trim(&mut rec.alignment_results.read, self.pre_trim_left, self.pre_trim_right);
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use crate::steps::generic::PerReadAlignmentResults;
    use crate::statistics::WorkerStatistics;
    use clap::Parser;

    #[test]
    fn test_single_pre_trim_apply() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--pre-trim-left", "5", "--pre-trim-right", "3"]);
        let context = Context::empty();
        let step = SinglePreTrim::new(&args, &context);
        let mut state = step.create_state();

        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: SingleAlignmentResults {
                read: PerReadAlignmentResults {
                    read_range: 0..15,
                    dropped: false,
                    ..Default::default()
                },
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };

        let res = step.apply(rec_in, &mut state);
        assert_eq!(res.alignment_results.read.read_range, 5..12);
        assert_eq!(res.alignment_results.read.bases_pre_trim, 8);
    }
}
