use std::io;
use crate::records::unified::AnnotatedRecord;
pub trait GenericWriterBackend<'a>: Send {
    fn write_read(&mut self, record: &AnnotatedRecord, read_index: usize) -> io::Result<()>;
}


pub const DEFAULT_TAG: &str = "/\\/default/\\/";
