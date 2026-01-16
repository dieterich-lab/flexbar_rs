use std::io;
use std::io::Write;
use crate::types::CompressionType;
use crate::utils::{get_writer_from_path, get_writer_path};
use crate::writers::Writer;

pub trait GenericWriterBackend<'a, R> {
    fn new(writer: Box<dyn Write + Send + 'a>) -> Self;
    fn write_record(&mut self, record: &R) -> io::Result<()>;
    fn get_extension() -> &'static str;
}


pub struct GenericWriter<'a, R, B: GenericWriterBackend<'a, R>> {
    backend: B,
    _phantom: std::marker::PhantomData<&'a (R, B)>,
}

impl<'a, R, B: GenericWriterBackend<'a, R>> GenericWriter<'a, R, B> {
    pub fn new(prefix: &str, compression: CompressionType) -> io::Result<Self> {
        let path = get_writer_path(prefix, compression, B::get_extension());
        let writer = get_writer_from_path(path.as_deref(), compression)?;
        Ok(Self::from_writer(writer))
    }

    pub fn from_writer(writer: Box<dyn Write + Send + 'a>) -> Self {
        Self {
            backend: B::new(writer),
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<'a, R, B: GenericWriterBackend<'a, R>> Writer<R> for GenericWriter<'a, R, B> {
    fn write(&mut self, record: &R) -> io::Result<()>
    {
        self.backend.write_record(record)
    }
}
