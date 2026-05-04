use crate::args::Args;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;
use crate::steps::generic::PerReadAlignmentResults;
use crate::records::unified::{UnifiedRead, GenericRawRead};

fn count_uncalled_bases(read: &UnifiedRead, results: &PerReadAlignmentResults) -> usize {
    if results.dropped {
        return 0;
    }
    let seq = read.seq();
    let range = results.read_range.clone();
    seq[range]
        .iter()
        .filter(|&&b| b == b'N' || b == b'n')
        .count()
}

pub struct UncalledFilter {
    pub max_uncalled: usize,
}

impl GenericStep for UncalledFilter {
    type State = ();

    fn new(args: &Args, _context: &Context) -> Self {
        if args.max_uncalled == 0 {
            panic!("max_uncalled must be greater than 0");
        }
        Self {
            max_uncalled: args.max_uncalled,
        }
    }

    fn create_state(&self) -> Self::State {
        ()
    }

    fn apply(&self, mut rec: RecordInTransformation, _state: &mut Self::State) -> RecordInTransformation {
        if rec.record.reads[1].is_none() {
            let per_read = rec.alignment_results.per_read[0].as_ref().unwrap();
            if !per_read.dropped {
                if let Some(s) = &rec.record.reads[0] {
                    let n_count = count_uncalled_bases(s, per_read);
                    if n_count > self.max_uncalled {
                        let per_read_mut = rec.alignment_results.per_read[0].as_mut().unwrap();
                        per_read_mut.dropped = true;
                        rec.stats.total_discarded_uncalled += 1;
                    }
                }
            }
            return rec;
        }

        let mut drop = false;
        for i in 0..2 {
            if !rec.alignment_results.per_read[i].as_ref().unwrap().dropped {
                if let Some(r) = &rec.record.reads[i] {
                    let n_count = count_uncalled_bases(r, rec.alignment_results.per_read[i].as_ref().unwrap());
                    if n_count > self.max_uncalled {
                        drop = true;
                        break;
                    }
                }
            }
        }
        if drop {
            for i in 0..2 {
                if !rec.alignment_results.per_read[i].as_ref().unwrap().dropped {
                    rec.alignment_results.per_read[i].as_mut().unwrap().dropped = true;
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
    use crate::types::RunType;
    use clap::Parser;
    use crate::records::unified::{UnifiedRecord, UnifiedRead};

    #[test]
    fn test_uncalled_filter_single() {
        let args = Args::parse_from(&["flexbar_rs", "--max-uncalled", "2"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let step = UncalledFilter::new(&args, &context);
        let mut state = step.create_state();

        let record = UnifiedRecord::single(UnifiedRead::new_fastq("read1", None, b"ACNNNGT", b"IIIIIII"));
        let res = step.apply(RecordInTransformation::from_dto(record), &mut state);
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.stats.total_discarded_uncalled, 1);
        
        let record2 = UnifiedRecord::single(UnifiedRead::new_fastq("read2", None, b"ACNNGT", b"IIIIII"));
        let res2 = step.apply(RecordInTransformation::from_dto(record2), &mut state);
        assert!(!res2.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res2.stats.total_discarded_uncalled, 0);
    }

    #[test]
    fn test_uncalled_filter_paired() {
        let args = Args::parse_from(&["flexbar_rs", "--max-uncalled", "2"]);
        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        let step = UncalledFilter::new(&args, &context);
        let mut state = step.create_state();

        let r1 = UnifiedRead::new_fastq("read1/1", None, b"ACNNNGT", b"IIIIIII");
        let r2 = UnifiedRead::new_fastq("read1/2", None, b"ACGT", b"IIII");
        let record = UnifiedRecord::paired(r1, r2);
        let res = step.apply(RecordInTransformation::from_dto(record), &mut state);
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert!(res.alignment_results.per_read[1].as_ref().unwrap().dropped);
        assert_eq!(res.stats.total_discarded_uncalled, 2);
        
        let r1_2 = UnifiedRead::new_fastq("read2/1", None, b"ACNNGT", b"IIIIII");
        let r2_2 = UnifiedRead::new_fastq("read2/2", None, b"ACNNGT", b"IIIIII");
        let record2 = UnifiedRecord::paired(r1_2, r2_2);
        let res2 = step.apply(RecordInTransformation::from_dto(record2), &mut state);
        assert!(!res2.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert!(!res2.alignment_results.per_read[1].as_ref().unwrap().dropped);
        assert_eq!(res2.stats.total_discarded_uncalled, 0);
    }
}
