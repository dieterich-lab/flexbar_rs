use std::io;
use std::io::Write;
use crate::types::{CompressionType, RecordDto};
use crate::utils::{get_writer_from_path, get_writer_path};
use crate::writers::Writer;

pub trait GenericWriterBackend<'a>: Send {
    fn new(writer: Box<dyn Write + Send + 'a>) -> Self;
    fn write_record(&mut self, record: &RecordDto) -> io::Result<()>;
    fn get_extension() -> &'static str;
}


pub struct GenericWriter<'a, B: GenericWriterBackend<'a>> {
    backend: B,
    statistics: crate::statistics::WriterStatistics,
    _phantom: std::marker::PhantomData<&'a ()>,
}

impl<'a, B: GenericWriterBackend<'a>> GenericWriter<'a, B> {
    pub fn new(prefix: &str, compression: CompressionType) -> io::Result<Self> {
        let path = get_writer_path(prefix, compression, B::get_extension());
        let writer = get_writer_from_path(path.as_deref(), compression)?;
        Ok(Self::from_writer(writer))
    }

    pub fn from_writer(writer: Box<dyn Write + Send + 'a>) -> Self {
        Self {
            backend: B::new(writer),
            statistics: crate::statistics::WriterStatistics::default(),
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<'a, B: GenericWriterBackend<'a>> Writer for GenericWriter<'a, B> {
    fn write(&mut self, record: &RecordDto) -> io::Result<()>
    {
        self.backend.write_record(record)?;
        self.statistics.total_records_written += 1;
        Ok(())
    }

    fn statistics(&self) -> crate::statistics::WriterStatistics {
        self.statistics
    }
}
