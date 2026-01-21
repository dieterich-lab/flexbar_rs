use std::collections::HashMap;
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


const DEFAULT_TAG: &str = "/\\/default/\\/";

pub struct GenericWriter<'a, B: GenericWriterBackend<'a>> {
    prefix: Option<String>,
    compression: CompressionType,
    backends_by_tag: HashMap<String, B>,
    statistics: crate::statistics::WriterStatistics,
    _phantom: std::marker::PhantomData<&'a ()>,
}

impl<'a, B: GenericWriterBackend<'a>> GenericWriter<'a, B> {
    pub fn new(prefix: &str, compression: CompressionType) -> io::Result<Self> {
        let path = get_writer_path(prefix, compression, B::get_extension(), None);
        let writer = get_writer_from_path(path.as_deref(), compression)?;
        let mut tagged_backends = HashMap::new();
        tagged_backends.insert(DEFAULT_TAG.to_string(), B::new(writer));
        Ok(Self {
            prefix: Some(prefix.to_string()),
            compression,
            backends_by_tag: tagged_backends,
            statistics: crate::statistics::WriterStatistics::default(),
            _phantom: std::marker::PhantomData,
        })
    }

    fn get_backend(&mut self, tag: Option<&str>) -> io::Result<&mut B> {
        let t = tag.unwrap_or(DEFAULT_TAG);
        if !self.backends_by_tag.contains_key(t) {
            if let Some(ref prefix) = self.prefix {
                let path = get_writer_path(prefix, self.compression, B::get_extension(), Some(t));
                let writer = get_writer_from_path(path.as_deref(), self.compression)?;
                self.backends_by_tag.insert(t.to_string(), B::new(writer));
            } else {
                return Err(io::Error::new(io::ErrorKind::Other, "Cannot create tagged writer without prefix"));
            }
        }
        Ok(self.backends_by_tag.get_mut(t).unwrap())
    }
}

impl<'a, B: GenericWriterBackend<'a>> Writer for GenericWriter<'a, B> {
    fn write(&mut self, record: &RecordDto, tag: Option<&str>) -> io::Result<()>
    {
        let backend = self.get_backend(tag)?;
        backend.write_record(record)?;
        self.statistics.total_records_written += 1;
        Ok(())
    }

    fn statistics(&self) -> crate::statistics::WriterStatistics {
        self.statistics
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn test_generic_writer_tagging_unsafe() {
        let temp = tempfile::tempdir().unwrap();
        let prefix = temp.path().join("out").to_str().unwrap().to_string();
        
        struct TestBackend<'a> {
            writer: Box<dyn Write + Send + 'a>,
        }
        impl<'a> GenericWriterBackend<'a> for TestBackend<'a> {
            fn new(writer: Box<dyn Write + Send + 'a>) -> Self { Self { writer } }
            fn write_record(&mut self, record: &RecordDto) -> io::Result<()> {
                self.writer.write_all(record.id().as_bytes())
            }
            fn get_extension() -> &'static str { "txt" }
        }

        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = GenericWriter::new(&prefix, CompressionType::Uncompressed).unwrap();
            let rec1 = RecordDto::with_attrs("R1", None, b"A", b"I");
            
            // This tag contains characters that should be removed in the filename
            writer.write(&rec1, Some("tag/with\\slashes")).unwrap();
        }

        let tagged_path = temp.path().join("out_tagwithslashes.txt");
        assert!(tagged_path.exists());
        
        let mut s = String::new();
        std::fs::File::open(tagged_path).unwrap().read_to_string(&mut s).unwrap();
        assert_eq!(s, "R1");
    }

    #[test]
    fn test_generic_writer_tagging() {
        let temp = tempfile::tempdir().unwrap();
        let prefix = temp.path().join("out").to_str().unwrap().to_string();
        
        struct TestBackend<'a> {
            writer: Box<dyn Write + Send + 'a>,
        }
        impl<'a> GenericWriterBackend<'a> for TestBackend<'a> {
            fn new(writer: Box<dyn Write + Send + 'a>) -> Self { Self { writer } }
            fn write_record(&mut self, record: &RecordDto) -> io::Result<()> {
                self.writer.write_all(record.id().as_bytes())
            }
            fn get_extension() -> &'static str { "txt" }
        }

        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = GenericWriter::new(&prefix, CompressionType::Uncompressed).unwrap();
            let rec1 = RecordDto::with_attrs("R1", None, b"A", b"I");
            let rec2 = RecordDto::with_attrs("R2", None, b"A", b"I");
            
            writer.write(&rec1, None).unwrap();
            writer.write(&rec2, Some("tagA")).unwrap();
        }

        let default_path = temp.path().join("out.txt");
        let tagged_path = temp.path().join("out_tagA.txt");

        assert!(default_path.exists());
        assert!(tagged_path.exists());
        
        let mut s1 = String::new();
        std::fs::File::open(default_path).unwrap().read_to_string(&mut s1).unwrap();
        assert_eq!(s1, "R1");

        let mut s2 = String::new();
        std::fs::File::open(tagged_path).unwrap().read_to_string(&mut s2).unwrap();
        assert_eq!(s2, "R2");
    }
}
