use crate::args::Args;
use crate::types::SingeRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::single::SingleAlignmentResults;
use crate::context::Context;
use crate::steps::utils::uncalled_filter::count_uncalled_bases;

pub struct SingleUncalledFilter {
    pub max_uncalled: usize,
}

pub struct SingleUncalledFilterState {}

impl GenericStep<SingeRecordDto> for SingleUncalledFilter {
    type State = SingleUncalledFilterState;
    type AlignmentResult = SingleAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            max_uncalled: args.max_uncalled,
        }
    }

    fn create_state(&self) -> Self::State {
        SingleUncalledFilterState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<SingeRecordDto, SingleAlignmentResults> {
        if !rec.alignment_results.read.dropped && self.max_uncalled > 0 {
            let n_count = count_uncalled_bases(&rec.record, &rec.alignment_results.read);
            if n_count > self.max_uncalled {
                rec.alignment_results.read.dropped = true;
                rec.stats.total_discarded_uncalled += 1;
            }
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
    use clap::Parser;

    #[test]
    fn test_single_uncalled_filter() {
        let args = Args::parse_from(&["flexbar_rs", "--max-uncalled", "2"]);
        let context = Context::empty();
        let step = SingleUncalledFilter::new(&args, &context);
        let mut state = step.create_state();

        let record = SingeRecordDto::with_attrs("read1", None, b"ACNNNGT", b"IIIIIII");
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: SingleAlignmentResults {
                read: PerReadAlignmentResults {
                    read_range: 0..7,
                    dropped: false,
                    ..Default::default()
                },
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };

        let res = step.apply(rec_in, &mut state);
        assert!(res.alignment_results.read.dropped);
        assert_eq!(res.stats.total_discarded_uncalled, 1);
        
        let record2 = SingeRecordDto::with_attrs("read2", None, b"ACNNGT", b"IIIIII");
        let rec_in2 = RecordInTransformation {
            record: record2,
            alignment_results: SingleAlignmentResults {
                read: PerReadAlignmentResults {
                    read_range: 0..6,
                    dropped: false,
                    ..Default::default()
                },
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };
        let res2 = step.apply(rec_in2, &mut state);
        assert!(!res2.alignment_results.read.dropped);
        assert_eq!(res2.stats.total_discarded_uncalled, 0);
    }
}
