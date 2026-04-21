use std::io;
use crate::types::SingeRecordDto;
pub trait GenericWriterBackend<'a>: Send {
    fn write_read(&mut self, record: &SingeRecordDto) -> io::Result<()>;
}


pub const DEFAULT_TAG: &str = "/\\/default/\\/";
