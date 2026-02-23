use crate::args::Args;
use crate::types::SingeRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};

use std::any::Any;

pub struct SingleFilterLength {
    min_read_length: usize,
}

impl SingleFilterLength {
    pub fn new(args: &Args) -> Self {
        Self {
            min_read_length: args.min_read_length,
        }
    }
}

pub struct SingleFilterLengthSetup;

impl GenericStep<SingeRecordDto> for SingleFilterLength {
    type Setup = Box<dyn Any + Send + Sync>;
    fn create_setup(&self) -> Box<dyn Any + Send + Sync> {
        Box::new(SingleFilterLengthSetup)
    }

    fn clone_setup(&self, _setup: &Box<dyn Any + Send + Sync>) -> Box<dyn Any + Send + Sync> {
        Box::new(SingleFilterLengthSetup)
    }

    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto>, _setup: &mut Box<dyn Any + Send + Sync>) -> RecordInTransformation<SingeRecordDto> {
        if !rec.alignment_results.record_dropped && self.min_read_length > 0 {
            let current_len = rec.alignment_results.seq_end - rec.alignment_results.seq_start;
            if current_len < self.min_read_length {
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
        let filter = SingleFilterLength::new(&args);
        
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
        let mut setup = filter.create_setup();
        let res = filter.apply(rec_in, &mut setup);
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_length, 0);
    }

    #[test]
    fn test_length_filter_drop() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "10"]);
        let filter = SingleFilterLength::new(&args);
        
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
        let mut setup = filter.create_setup();
        let res = filter.apply(rec_in, &mut setup);
        assert!(res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_length, 1);
    }

    #[test]
    fn test_length_filter_disabled() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--min-read-length", "0"]);
        let filter = SingleFilterLength::new(&args);
        
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
        let mut setup = filter.create_setup();
        let res = filter.apply(rec_in, &mut setup);
        assert!(!res.alignment_results.record_dropped);
    }
}
