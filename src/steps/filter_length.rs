use crate::args::Args;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;
use crate::records::generic_read::GenericRawRead;
use crate::records::unified::UnifiedRead;

pub struct FilterLength {
    pub min_read_length: usize,
    pub single_reads: bool,
    pub single_reads_paired: bool,
}

pub struct FilterLengthState {}

impl GenericStep for FilterLength {
    type State = FilterLengthState;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            min_read_length: args.min_read_length,
            single_reads: args.single_reads,
            single_reads_paired: args.single_reads_paired,
        }
    }

    fn create_state(&self) -> Self::State {
        FilterLengthState {}
    }

    fn apply(&self, mut rec: RecordInTransformation, _state: &mut Self::State) -> RecordInTransformation {
        if self.min_read_length == 0 {
            return rec;
        }
        let mut good_count = 0;
        let mut bad_count = 0;
        let mut is_short = [false; 2];
    
        for i in 0..2 {
            if rec.record.reads[i].is_none() || rec.is_dropped_read(i) { continue }
            if rec.get_read_length(i) < self.min_read_length {
                bad_count += 1;
                is_short[i] = true;
            } else {
                good_count += 1;
            }
        }

        for i in 0..2 {
            if is_short[i] {
                if good_count > 0 && self.single_reads_paired {  // Replace read with N
                    trace!("Length filter replaced {} with 'N'", rec.read_id(i));
                    let start = rec.alignment_results.per_read[i].as_ref().unwrap().read_range.start;
                    rec.alignment_results.per_read[i].as_mut().unwrap().read_range.end = start + 1;
                    if let Some(old) = &rec.record.reads[i] {
                        rec.record.reads[i] = Some(get_place_holder_read(old, start));
                    }
                    rec.stats.total_single_reads_replaced += 1;

                } else {
                    rec.drop_read(i);
                    trace!("Length filter dropped {}: ", rec.read_id(i));
                    rec.stats.total_discarded_length += 1;
                    
                }
            } else {
                if rec.record.reads[i].is_some() && !rec.is_dropped_read(i) && bad_count > 0 {
                    if !self.single_reads_paired && !self.single_reads {
                        trace!("Length filter dropped {} (other read too short): ", rec.read_id(i));
                        rec.drop_read(i);
                        rec.stats.total_discarded_length += 1;
                    }
                }
            }
        }
        rec
    }
}

fn get_place_holder_read(read: &UnifiedRead, start: usize) -> UnifiedRead {
    let seq = read.seq();
    let len = seq.len();
    let mut new_sequence: Vec<u8> = Vec::with_capacity(len);
    new_sequence.extend_from_slice(&seq[..start]);
    new_sequence.push(b'N');
    new_sequence.extend_from_slice(&seq[start + 1..]);
    read.set_seq(&new_sequence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;
    use crate::records::unified::{UnifiedRecord, UnifiedRead};
    use crate::types::RunType;

    #[test]
    fn test_length_filter_keep_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let filter = FilterLength::new(&args, &context);
        let mut state = filter.create_state();
        
        let record = UnifiedRecord::single(UnifiedRead::new_fastq("read1", None, b"0123456789", b"IIIIIIIIII"));
        let res = filter.apply(RecordInTransformation::from_dto(record), &mut state);
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.stats.total_discarded_length, 0);
    }

    #[test]
    fn test_length_filter_drop_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let filter = FilterLength::new(&args, &context);
        let mut state = filter.create_state();
        
        let record = UnifiedRecord::single(UnifiedRead::new_fastq("read1", None, b"012345678", b"IIIIIIIII"));
        let res = filter.apply(RecordInTransformation::from_dto(record), &mut state);
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert_eq!(res.stats.total_discarded_length, 1);
    }

    #[test]
    fn test_length_filter_disabled_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "0"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let filter = FilterLength::new(&args, &context);
        let mut state = filter.create_state();
        
        let record = UnifiedRecord::single(UnifiedRead::new_fastq("read1", None, b"A", b"I"));
        let res = filter.apply(RecordInTransformation::from_dto(record), &mut state);
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
    }

    #[test]
    fn test_length_filter_paired_keep() {
        let mut args = Args::parse_from(["flexbar", "--min-read-length", "10"]);
        args.min_read_length = 10;
        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        let step = FilterLength::new(&args, &context);
        let mut state = step.create_state();

        let r1 = UnifiedRead::new_fastq("id", None, b"0123456789", b"IIIIIIIIII");
        let r2 = UnifiedRead::new_fastq("id", None, b"0123456789", b"IIIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);

        let res = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        assert!(!res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert!(!res.alignment_results.per_read[1].as_ref().unwrap().dropped);
    }

    #[test]
    fn test_length_filter_paired_drop_all() {
        let mut args = Args::parse_from(["flexbar", "--min-read-length", "10"]);
        args.min_read_length = 10;
        let mut context = Context::empty();
        context.run_type = RunType::Paired;
        let step = FilterLength::new(&args, &context);
        let mut state = step.create_state();

        let r1 = UnifiedRead::new_fastq("id", None, b"012345678", b"IIIIIIIII");
        let r2 = UnifiedRead::new_fastq("id", None, b"012345678", b"IIIIIIIII");
        let rec = UnifiedRecord::paired(r1, r2);

        let res = step.apply(RecordInTransformation::from_dto(rec), &mut state);
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert!(res.alignment_results.per_read[1].as_ref().unwrap().dropped);
        assert_eq!(res.stats.total_discarded_length, 2);
    }
}
