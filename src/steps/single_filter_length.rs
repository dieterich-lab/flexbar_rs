use crate::args::Args;
use crate::types::SingeRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::context::Context;

pub struct SingleFilterLength {}

#[derive(Clone)]
pub struct SingleFilterLengthSetup {
    pub min_read_length: usize,
}

impl GenericStep<SingeRecordDto> for SingleFilterLength {
    type Setup = SingleFilterLengthSetup;

    fn create_setup(args: &Args, _context: &Context) -> SingleFilterLengthSetup {
        SingleFilterLengthSetup {
            min_read_length: args.min_read_length,
        }
    }

    fn new(_setup: &SingleFilterLengthSetup) -> Self {
        Self {}
    }

    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>, setup: &SingleFilterLengthSetup) -> RecordInTransformation<SingeRecordDto> {
        if !rec.alignment_results.record_dropped && setup.min_read_length > 0 {
            let current_len = rec.alignment_results.seq_end - rec.alignment_results.seq_start;
            if current_len < setup.min_read_length {
                trace!("Length filter dropped record: {} (length: {})", rec.record.id(), current_len);
                rec.alignment_results.record_dropped = true;
                rec.stats.total_discarded_length += 1;
            }
        }
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use crate::steps::generic::AlignmentResults;
    use crate::statistics::WorkerStatistics;
    use clap::Parser;

    #[test]
    fn test_length_filter_keep() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let context = Context::empty();
        let setup = SingleFilterLength::create_setup(&args, &context);
        let mut filter = SingleFilterLength::new(&SingleFilterLength::create_setup(&args, &context));
        
        let record = SingeRecordDto::with_attrs("read1", None, b"0123456789", b"IIIIIIIIII");
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: AlignmentResults {
                seq_start: 0,
                seq_end: 10,
                record_dropped: false,
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };
        let res = filter.apply(rec_in, &setup);
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_length, 0);
    }

    #[test]
    fn test_length_filter_drop() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let context = Context::empty();
        let setup = SingleFilterLength::create_setup(&args, &context);
        let mut filter = SingleFilterLength::new(&SingleFilterLength::create_setup(&args, &context));
        
        let record = SingeRecordDto::with_attrs("read1", None, b"012345678", b"IIIIIIIII");
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: AlignmentResults {
                seq_start: 0,
                seq_end: 9,
                record_dropped: false,
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };
        let res = filter.apply(rec_in, &setup);
        assert!(res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_length, 1);
    }

    #[test]
    fn test_length_filter_disabled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "0"]);
        let context = Context::empty();
        let setup = SingleFilterLength::create_setup(&args, &context);
        let mut filter = SingleFilterLength::new(&SingleFilterLength::create_setup(&args, &context));
        
        let record = SingeRecordDto::with_attrs("read1", None, b"A", b"I");
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: AlignmentResults {
                seq_start: 0,
                seq_end: 1,
                record_dropped: false,
                ..Default::default()
            },
            stats: WorkerStatistics::default(),
        };
        let res = filter.apply(rec_in, &setup);
        assert!(!res.alignment_results.record_dropped);
    }
}
