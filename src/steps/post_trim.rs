use crate::args::Args;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;

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
        for i in 0..2 {
            if rec.record.reads[i].is_none() || rec.is_dropped_read(i) { continue }
            let read_length = rec.get_read_length(i);
            if read_length > self.post_trim_length {
                let start = rec.alignment_results.per_read[i].as_ref().unwrap().read_range.start;
                rec.alignment_results.per_read[i].as_mut().unwrap().read_range.end = start + self.post_trim_length;
                rec.alignment_results.per_read[i].as_mut().unwrap().bases_post_trim = read_length - self.post_trim_length;
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
    use crate::types::{RecordDto, RunType, SingeRecordDto};

    #[test]
    fn test_post_trim_apply_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--post-trim-length", "10"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let step = PostTrim::new(&args, &context);
        let mut state = step.create_state();

        let record = RecordDto::single(SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII"));
        let res = step.apply(RecordInTransformation::from_dto(record), &mut state);
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
        let rec = RecordDto::paired(r1, r2);

        let res = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range, 0..10);
        assert_eq!(res.alignment_results.per_read[1].as_ref().unwrap().read_range, 0..10);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().bases_post_trim, 5);
        assert_eq!(res.alignment_results.per_read[1].as_ref().unwrap().bases_post_trim, 5);
    }
}
