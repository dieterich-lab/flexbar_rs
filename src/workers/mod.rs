use crate::steps::generic::FinalRecord;
use crate::types::{RecordDto};
use crate::workers::single::{SingleWorker, SingleWorkerState};
use crate::args::Args;
use crate::context::Context;

pub mod single;

pub enum Worker {
    Single(SingleWorker)
}

pub enum WorkerState {
    Single(SingleWorkerState)
}

impl Worker {
    pub fn new(args: &Args, context: &Context) -> Self {
        Worker::Single(SingleWorker::new(args, context))
    }

    pub fn create_state(&self) -> WorkerState {
        match self {
            Worker::Single(w) => WorkerState::Single(w.create_state()),
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
            (Worker::Single(_), RecordDto::Paired(_), _) => {
                panic!("Single worker cannot process paired records");
            }
        }
    }
}