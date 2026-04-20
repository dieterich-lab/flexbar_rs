use std::io;
use std::io::Write;
use crate::types::{SingeRecordDto};

pub trait GenericWriterBackend<'a>: Send {
    fn write_read(&mut self, record: &SingeRecordDto) -> io::Result<()>;
}


pub const DEFAULT_TAG: &str = "/\\/default/\\/";

// (bar_tag, read_index, extension) -> (writer, description)
pub type WriterFactory<'a> = Box<dyn
    Fn(Option<&str>, Option<usize>, &'static str)
        -> io::Result<(Box<dyn Write + Send + 'a>, String)
> + Send + 'a>;

