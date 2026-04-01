use crate::args::Args;
use crate::types::SingeRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::single::SingleAlignmentResults;
use crate::context::Context;
use crate::steps::utils::post_trim::apply_post_trim;

pub struct SinglePostTrim {
    pub post_trim_length: usize,
}

pub struct SinglePostTrimState {}

impl GenericStep<SingeRecordDto> for SinglePostTrim {
    type State = SinglePostTrimState;
    type AlignmentResult = SingleAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            post_trim_length: args.post_trim_length,
        }
    }

    fn create_state(&self) -> Self::State {
        SinglePostTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        apply_post_trim(&mut rec.alignment_results.read, self.post_trim_length);
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
    fn test_single_post_trim_apply() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--post-trim-length", "10"]);
        let context = Context::empty();
        let step = SinglePostTrim::new(&args, &context);
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
        assert_eq!(res.alignment_results.read.read_range, 0..10);
        assert_eq!(res.alignment_results.read.bases_post_trim, 5);
    }
}
