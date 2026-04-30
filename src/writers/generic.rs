use std::io;
use crate::records::unified::SingeRecordDto;
pub trait GenericWriterBackend<'a>: Send {
    fn write_read(&mut self, record: &SingeRecordDto) -> io::Result<()>;
}


pub const DEFAULT_TAG: &str = "/\\/default/\\/";
