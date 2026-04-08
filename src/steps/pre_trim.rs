use crate::args::Args;
use crate::types::{RecordDto, RunType};
use crate::steps::generic::{GenericStep, RecordInTransformation, CombinedAligmentResults};
use crate::context::Context;
use crate::steps::utils::pre_trim::apply_pre_trim;

pub struct PreTrim {
    pub pre_trim_left: usize,
    pub pre_trim_right: usize,
}

pub struct PreTrimState {}

impl GenericStep<RecordDto> for PreTrim {
    type State = PreTrimState;
    type AlignmentResult = CombinedAligmentResults;

    fn new(args: &Args, context: &Context) -> Self {
        match context.run_type {
            RunType::SingleBarcoded | RunType::PairedBarcoded => {
                panic!("Barcoded run types not implemented yet");
            }
            _ => {}
        }
        Self {
            pre_trim_left: args.pre_trim_left,
            pre_trim_right: args.pre_trim_right,
        }
    }

    fn create_state(&self) -> Self::State {
        PreTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<RecordDto, CombinedAligmentResults>, _state: &mut Self::State) -> RecordInTransformation<RecordDto, CombinedAligmentResults> {
        match rec.record {
            RecordDto::Single(_) => {
                apply_pre_trim(rec.alignment_results.per_read[0].as_mut().unwrap(), self.pre_trim_left, self.pre_trim_right);
            }
            RecordDto::Paired(_) => {
                for i in 0..2 {
                    apply_pre_trim(rec.alignment_results.per_read[i].as_mut().unwrap(), self.pre_trim_left, self.pre_trim_right);
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
    use crate::types::SingeRecordDto;
    use crate::types::PairedRecordDto;
    use crate::steps::generic::Transformable;
    use clap::Parser;

    #[test]
    fn test_pre_trim_apply_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--pre-trim-left", "5", "--pre-trim-right", "3"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let step = PreTrim::new(&args, &context);
        let mut state = step.create_state();

        let record = RecordDto::Single(SingeRecordDto::with_attrs("read1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII"));
        let res = step.apply(record.get_record_in_transformation(), &mut state);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range, 5..12);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().bases_pre_trim, 8);
    }

    #[test]
    fn test_pre_trim_apply_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test_1.fq", "--reads2", "test_2.fq", "--pre-trim-left", "5", "--pre-trim-right", "3"]);
        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        let step = PreTrim::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("read1/1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("read1/2", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let record = RecordDto::Paired(PairedRecordDto::new(r1, r2));
        
        let res = step.apply(record.get_record_in_transformation(), &mut state);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().read_range, 5..12);
        assert_eq!(res.alignment_results.per_read[1].as_ref().unwrap().read_range, 5..12);
        assert_eq!(res.alignment_results.per_read[0].as_ref().unwrap().bases_pre_trim, 8);
        assert_eq!(res.alignment_results.per_read[1].as_ref().unwrap().bases_pre_trim, 8);
    }
}
