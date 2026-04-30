use std::io::BufRead;
use crate::records::unified::SingeRecordDto;

pub trait GenericReaderBackend<Reader: BufRead + Send> {
    fn read_record(&mut self) -> Option<SingeRecordDto>;
}