use crate::args::Args;
use crate::types::PairedRecordDto;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::context::Context;
use crate::steps::utils::uncalled_filter::count_uncalled_bases;

pub struct PairedUncalledFilter {
    pub max_uncalled: usize,
}

pub struct PairedUncalledFilterState {}

impl GenericStep<PairedRecordDto> for PairedUncalledFilter {
    type State = PairedUncalledFilterState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        if args.max_uncalled == 0 {
            panic!("max_uncalled must be greater than 0");
        }
        Self {
            max_uncalled: args.max_uncalled,
        }
    }

    fn create_state(&self) -> Self::State {
        PairedUncalledFilterState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        let mut drop = false;
        for i in 0..2 {
            if !rec.alignment_results.per_read[i].dropped {
                let n_count = count_uncalled_bases(&rec.record.records[i], &rec.alignment_results.per_read[i]);
                if n_count > self.max_uncalled {
                    drop = true;
                    break;
                }
            }
        }

        if drop {
            for i in 0..2 {
                if !rec.alignment_results.per_read[i].dropped {
                    rec.alignment_results.per_read[i].dropped = true;
                    rec.stats.total_discarded_uncalled += 1;
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
    use crate::steps::generic::PerReadAlignmentResults;
    use crate::statistics::WorkerStatistics;
    use clap::Parser;

    #[test]
    fn test_paired_uncalled_filter() {
        let args = Args::parse_from(&["flexbar_rs", "--max-uncalled", "2"]);
        let context = Context::empty();
        let step = PairedUncalledFilter::new(&args, &context);
        let mut state = step.create_state();

        let r1 = SingeRecordDto::with_attrs("read1/1", None, b"ACNNNGT", b"IIIIIII");
        let r2 = SingeRecordDto::with_attrs("read1/2", None, b"ACGT", b"IIII");
        let record = PairedRecordDto::new(r1, r2);
        
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: PairedAlignmentResults {
                per_read: [
                    PerReadAlignmentResults {
                        read_range: 0..7,
                        dropped: false,
                        ..Default::default()
                    },
                    PerReadAlignmentResults {
                        read_range: 0..4,
                        dropped: false,
                        ..Default::default()
                    },
                ],
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };

        let res = step.apply(rec_in, &mut state);
        assert!(res.alignment_results.per_read[0].dropped);
        assert!(res.alignment_results.per_read[1].dropped);
        assert_eq!(res.stats.total_discarded_uncalled, 2);
        
        let r1_2 = SingeRecordDto::with_attrs("read2/1", None, b"ACNNGT", b"IIIIII");
        let r2_2 = SingeRecordDto::with_attrs("read2/2", None, b"ACNNGT", b"IIIIII");
        let record2 = PairedRecordDto::new(r1_2, r2_2);
        let rec_in2 = RecordInTransformation {
            record: record2,
            alignment_results: PairedAlignmentResults {
                per_read: [
                    PerReadAlignmentResults {
                        read_range: 0..6,
                        dropped: false,
                        ..Default::default()
                    },
                    PerReadAlignmentResults {
                        read_range: 0..6,
                        dropped: false,
                        ..Default::default()
                    },
                ],
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };
        let res2 = step.apply(rec_in2, &mut state);
        assert!(!res2.alignment_results.per_read[0].dropped);
        assert!(!res2.alignment_results.per_read[1].dropped);
        assert_eq!(res2.stats.total_discarded_uncalled, 0);
    }
}
