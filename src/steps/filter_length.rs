use bio::io::fastq::Record;
use crate::args::Args;
use crate::types::{RecordDto, RunType, SingeRecordDto};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation, CombinedAligmentResults};
use crate::context::Context;

pub struct FilterLength {
    pub min_read_length: usize,
    pub single_reads: bool,
    pub single_reads_paired: bool,
}

pub struct FilterLengthState {}

impl GenericStep<RecordDto> for FilterLength {
    type State = FilterLengthState;
    type AlignmentResult = CombinedAligmentResults;

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

    fn apply(&self, mut rec: RecordInTransformation<RecordDto, CombinedAligmentResults>, _state: &mut Self::State) -> RecordInTransformation<RecordDto, CombinedAligmentResults> {
        if self.min_read_length == 0 {
            return rec;
        }

        match &mut rec.record {
            RecordDto::Single(s) => {
                let per_read = rec.alignment_results.per_read[0].as_mut().unwrap();
                if !per_read.dropped {
                    let current_len = per_read.read_range.end - per_read.read_range.start;
                    if current_len < self.min_read_length {
                        trace!("Length filter dropped record: {} (length: {})", s.id(), current_len);
                        per_read.dropped = true;
                        rec.stats.total_discarded_length += 1;
                    }
                }
            }
            RecordDto::Paired(p) => {
                let mut short = [false; 2];
                let mut good_reads: u32 = 0;
                for i in 0..2 {
                    let per_read = rec.alignment_results.per_read[i].as_ref().unwrap();
                    if !per_read.dropped {
                        if per_read.read_range.end - per_read.read_range.start < self.min_read_length {
                            short[i] = true;
                        } else {
                            good_reads += 1;
                        }
                    }
                }
                if good_reads == 2 {
                    // All fine ...
                } else if good_reads == 0 {
                    for i in 0..2 {
                        if short[i] {
                            rec.stats.total_discarded_length += 1;
                            rec.alignment_results.per_read[i].as_mut().unwrap().dropped = true;
                        }
                    }
                    trace!("Length filter dropped paired record: {}", p.id());
                } else if self.single_reads_paired {
                    for i in 0..2 {
                        if short[i] {
                            let start = rec.alignment_results.per_read[i].as_ref().unwrap().read_range.start;
                            rec.alignment_results.per_read[i].as_mut().unwrap().read_range.end = start + 1;
                            p.reads[i] = get_place_holder_record(&p.reads[i], start);
                            rec.stats.total_single_reads_replaced += 1;
                            trace!("Length filter replaced R{} with 'N': {}", i + 1, p.id());
                        }
                    }
                } else if self.single_reads {
                    for i in 0..2 {
                        if short[i] {
                            rec.alignment_results.per_read[i].as_mut().unwrap().dropped = true;
                            rec.stats.total_single_reads_written += 1;
                            trace!("Length filter dropped R{}, keeping R{}: {}",
                                i + 1, 2 - i, p.id()
                            );
                        }
                    }
                } else {
                    for i in 0..2 {
                        if !rec.alignment_results.per_read[i].as_ref().unwrap().dropped {
                            rec.alignment_results.per_read[i].as_mut().unwrap().dropped = true;
                            rec.stats.total_single_reads_dropped += 1;
                            rec.stats.total_discarded_length += 1;
                        }
                    }
                    trace!("Length filter dropped paired record: {}", p.id());
                }
            }
        }
        rec
    }
}

fn get_place_holder_record(record: &Record, start: usize) -> Record {
    let seq = record.seq();
    let len = seq.len();
    let mut new_sequence: Vec<u8> = Vec::with_capacity(len);
    new_sequence.extend_from_slice(&seq[..start]);
    new_sequence.push(b'N');
    new_sequence.extend_from_slice(&seq[start + 1..]);
    Record::with_attrs(record.id(), record.desc(), &*new_sequence, record.qual())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;
    use crate::types::PairedRecordDto;
    use crate::steps::generic::Transformable;

    #[test]
    fn test_length_filter_keep_single() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let mut context = Context::empty();
        context.run_type = RunType::Single;
        let filter = FilterLength::new(&args, &context);
        let mut state = filter.create_state();
        
        let record = RecordDto::Single(SingeRecordDto::with_attrs("read1", None, b"0123456789", b"IIIIIIIIII"));
        let res = filter.apply(record.get_record_in_transformation(), &mut state);
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
        
        let record = RecordDto::Single(SingeRecordDto::with_attrs("read1", None, b"012345678", b"IIIIIIIII"));
        let res = filter.apply(record.get_record_in_transformation(), &mut state);
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
        
        let record = RecordDto::Single(SingeRecordDto::with_attrs("read1", None, b"A", b"I"));
        let res = filter.apply(record.get_record_in_transformation(), &mut state);
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

        let r1 = SingeRecordDto::with_attrs("id", None, b"0123456789", b"IIIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"0123456789", b"IIIIIIIIII");
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));

        let res = step.apply(rec.get_record_in_transformation(), &mut state);
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

        let r1 = SingeRecordDto::with_attrs("id", None, b"012345678", b"IIIIIIIII");
        let r2 = SingeRecordDto::with_attrs("id", None, b"012345678", b"IIIIIIIII");
        let rec = RecordDto::Paired(PairedRecordDto::new(r1, r2));

        let res = step.apply(rec.get_record_in_transformation(), &mut state);
        assert!(res.alignment_results.per_read[0].as_ref().unwrap().dropped);
        assert!(res.alignment_results.per_read[1].as_ref().unwrap().dropped);
        assert_eq!(res.stats.total_discarded_length, 2);
    }
}
