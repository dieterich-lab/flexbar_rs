use std::io::BufRead;
use crate::types::SingeRecordDto;

pub trait GenericReaderBackend<Reader: BufRead + Send> {
    fn read_record(&mut self) -> Option<SingeRecordDto>;
}