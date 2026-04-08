use crate::args::Args;
use crate::types::RecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;
use crate::steps::generic::PerReadAlignmentResults;

fn apply_post_trim(per_read: &mut PerReadAlignmentResults, post_trim_length: usize) {
    if per_read.dropped || post_trim_length == 0 {
        return;
    }
    let current_length = per_read.read_range.end - per_read.read_range.start;
    if current_length > post_trim_length {
        per_read.bases_post_trim += current_length - post_trim_length;
        per_read.read_range.end = per_read.read_range.start + post_trim_length;
    }
}

pub struct PostTrim {
    pub post_trim_length: usize,
}

pub struct PostTrimState {}

impl GenericStep for PostTrim {
    type State = PostTrimState;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            post_trim_length: args.post_trim_length,
        }
    }

    fn create_state(&self) -> Self::State {
        PostTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation, _state: &mut Self::State) -> RecordInTransformation {
        match rec.record {
            RecordDto::Single(_) => {
                apply_post_trim(rec.alignment_results.per_read[0].as_mut().unwrap(), self.post_trim_length);
            }
            RecordDto::Paired(_) => {
                for i in 0..2 {
                    apply_post_trim(rec.alignment_results.per_read[i].as_mut().unwrap(), self.post_trim_length);
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
    use crate::types::{RunType, SingeRecordDto, PairedRecordDto};

    #[test]
    fn test_post_trim_apply_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--post-trim-length", "10"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let step = PostTrim::new(&args, &context);
        let mut state = step.create_state();

        let record = RecordDto::Single(SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII"));
        let res = step.apply(record.get_record_in_transformation(), &mut state);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..10);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().bases_post_trim, 5);
    }

    #[test]
    fn test_post_trim_apply_paired() {
        let mut args = Args::parse_from(["flexbar", "--post-trim-length", "10"]);
        args.post_trim_length = 10;
        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        let step = PostTrim::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("id", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));

        let res = step.apply(rec.get_record_in_transformation(), &mut state);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..10);
        assert_eq!(res.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..10);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().bases_post_trim, 5);
        assert_eq!(res.alignment_results.per_read[1].as_ref().unwrap().bases_post_trim, 5);
    }
}
