use crate::steps::generic::FinalRecord;
use crate::types::{RecordDto};
use crate::workers::single::{SingleWorker, SingleWorkerState};
use crate::workers::paired::{PairedWorker, PairedWorkerState};
use crate::args::Args;
use crate::context::Context;

pub mod single;
pub mod paired;

pub enum Worker {
    Single(SingleWorker),
    Paired(PairedWorker)
}

pub enum WorkerState {
    Single(SingleWorkerState),
    Paired(PairedWorkerState)
}

impl Worker {
    pub fn new(args: &Args, context: &Context) -> Self {
        if args.reads2.is_some() {
            Worker::Paired(PairedWorker::new(args, context))
        } else {
            Worker::Single(SingleWorker::new(args, context))
        }
    }

    pub fn create_state(&self) -> WorkerState {
        match self {
            Worker::Single(w) => WorkerState::Single(w.create_state()),
            Worker::Paired(w) => WorkerState::Paired(w.create_state()),
        }
    }

    pub fn process_record(&self, record: &RecordDto, state: &mut WorkerState) -> FinalRecord<RecordDto> {
        match (self, record, state) {
            (
                Worker::Single(w), 
                RecordDto::Single(r), 
                WorkerState::Single(s)
            ) => {
                let raw     = w.process_record(r, s);
                FinalRecord {
                    record: match raw.record {
                        Some(r) => Some(RecordDto::Single(r)),
                        None => None,
                    },
                    output_tag: raw.output_tag,
                    stats: raw.stats,
                }
            },
            (
                Worker::Paired(w), 
                RecordDto::Paired(r), 
                WorkerState::Paired(s)
            ) => {
                let raw     = w.process_record(r, s);
                FinalRecord {
                    record: match raw.record {
                        Some(r) => Some(RecordDto::Paired(r)),
                        None => None,
                    },
                    output_tag: raw.output_tag,
                    stats: raw.stats,
                }
            },
            (Worker::Single(_), RecordDto::Paired(_), _) => {
                panic!("Single worker cannot process paired records");
            }
            (Worker::Paired(_), RecordDto::Single(_), _) => {
                panic!("Paired worker cannot process single records");
            }
            _ => panic!("Inconsistent worker and record types"),
        }
    }
}