use std::io;
use crate::records::unified::{UnifiedRead};
pub trait GenericWriterBackend<'a>: Send {
    fn write_read(&mut self, record: &UnifiedRead) -> io::Result<()>;
}


pub const DEFAULT_TAG: &str = "/\\/default/\\/";
