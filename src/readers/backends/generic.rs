use std::io::BufRead;
use crate::types::SingeRecordDto;

pub trait GenericReaderBackend<Reader: BufRead + Send> {
    type Error: std::fmt::Display;
    fn read_record(&mut self) -> Result<Option<SingeRecordDto>, Self::Error>;
}