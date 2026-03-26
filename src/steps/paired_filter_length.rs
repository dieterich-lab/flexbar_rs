use bio::io::fastq::Record;
use crate::args::Args;
use crate::types::PairedRecordDto;
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::context::Context;

pub struct PairedFilterLength {
    pub min_read_length: usize,
    pub single_reads: bool,
    pub single_reads_paired: bool,
}

pub struct PairedFilterLengthState {}

impl GenericStep<PairedRecordDto> for PairedFilterLength {
    type State = PairedFilterLengthState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, _context: &Context) -> Self {
        Self {
            min_read_length: args.min_read_length,
            single_reads: args.single_reads,
            single_reads_paired: args.single_reads_paired,
        }
    }

    fn create_state(&self) -> Self::State {
        PairedFilterLengthState {}
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, _state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        if self.min_read_length == 0 {
            return rec;
        }
        let mut short = [false; 2];
        let mut good_reads: u32 = 0;
        for i in 0..2 {
            let per_read = &rec.alignment_results.per_read[i];
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
                    rec.alignment_results.per_read[i].dropped = true;
                }
            }
            trace!("Length filter dropped paired record: {}", rec.record.id());
        } else if self.single_reads_paired {
            for i in 0..2 {
                if short[i] {
                    let start = rec.alignment_results.per_read[i].read_range.start;
                    rec.alignment_results.per_read[i].read_range.end = start + 1;
                    rec.record.records[i] = get_place_holder_record(&rec.record.records[i], start);
                    rec.stats.total_single_reads_replaced += 1;
                    trace!("Length filter replaced R{} with 'N': {}", i + 1, rec.record.id());
                }
            }
        } else if self.single_reads {
            for i in 0..2 {
                if short[i] {
                    rec.alignment_results.per_read[i].dropped = true;
                    rec.stats.total_single_reads_written += 1;
                    trace!("Length filter dropped R{}, keeping R{}: {}",
                        i + 1, 2 - i, rec.record.id()
                    );
                }
            }
        } else {
            for i in 0..2 {
                if !rec.alignment_results.per_read[i].dropped {
                    rec.alignment_results.per_read[i].dropped = true;
                    rec.stats.total_single_reads_dropped += 1;
                    rec.stats.total_discarded_length += 1;
                }
            }
            trace!("Length filter dropped paired record: {}", rec.record.id());
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
