use crate::steps::generic::FinalRecord;
use crate::types::{RecordDto};
use crate::workers::single::{SingleWorker, SingleWorkerSetup};

pub mod single;

pub enum WorkerSetup<'a> {
    Single(SingleWorkerSetup<'a>)
}

pub enum Worker {
    Single(SingleWorker)
}

pub struct WorkerFactory {
    cache: Option<Worker>,
}

impl WorkerFactory {
    pub fn new() -> Self {
        Self { cache: None }
    }
    
    pub fn get_or_create(&mut self, setup: &WorkerSetup) -> &mut Worker {
        match &self.cache {
            Some(_) => (),
            None => {
                self.cache = Some(
                    match setup {
                        WorkerSetup::Single(setup)
                            => Worker::Single(SingleWorker::new(setup)),
                    }
                )
            }
        }
        self.cache.as_mut().unwrap()
    }
}

impl Worker {
    pub fn process_record(&mut self, record: &RecordDto, setup: &WorkerSetup) -> FinalRecord<RecordDto> {
        match (self, record, setup) {
            (
                Worker::Single(w), 
                RecordDto::Single(r), 
                WorkerSetup::Single(s)
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
            }
        }
    }
}