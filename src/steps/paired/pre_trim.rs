use crate::args::Args;
use crate::types::PairedRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::context::Context;
use crate::steps::utils::pre_trim::apply_pre_trim;

pub struct PairedPreTrim {
    pub pre_trim_left: usize,
    pub pre_trim_right: usize,
}

pub struct PairedPreTrimState {}

impl GenericStep<PairedRecordDto> for PairedPreTrim {
    type State = PairedPreTrimState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            pre_trim_left: args.pre_trim_left,
            pre_trim_right: args.pre_trim_right,
        }
    }

    fn create_state(&self) -> Self::State {
        PairedPreTrimState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        for i in 0..2 {
            apply_pre_trim(&mut rec.alignment_results.per_read[i], self.pre_trim_left, self.pre_trim_right);
        }
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use crate::steps::generic::PerReadAlignmentResults;
    use crate::statistics::WorkerStatistics;
    use crate::types::SingeRecordDto;
    use clap::Parser;

    #[test]
    fn test_paired_pre_trim_apply() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test_1.fq", "--reads2", "test_2.fq", "--pre-trim-left", "5", "--pre-trim-right", "3"]);
        let context = Context::empty();
        let step = PairedPreTrim::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("read1/1", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("read1/2", None, b"0123456789ABCDE", b"IIIIIIIIIIIIIII");
        let record = PairedRecordDto { records: [r1, r2] };
        
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: PairedAlignmentResults {
                barcode_tag: None,
                per_read: [
                    PerReadAlignmentResults {
                        read_range: 0..15,
                        dropped: false,
                        ..Default::default()
                    },
                    PerReadAlignmentResults {
                        read_range: 0..15,
                        dropped: false,
                        ..Default::default()
                    },
                ],
            },
            stats: WorkerStatistics::default(),
        };

        let res = step.apply(rec_in, &mut state);
        assert_eq!(res.alignment_results.per_read[0].read_range, 5..12);
        assert_eq!(res.alignment_results.per_read[1].read_range, 5..12);
        assert_eq!(res.alignment_results.per_read[0].bases_pre_trim, 8);
        assert_eq!(res.alignment_results.per_read[1].bases_pre_trim, 8);
    }
}
