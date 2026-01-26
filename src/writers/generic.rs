use std::collections::HashMap;
use std::io;
use std::io::Write;
use log::trace;
use crate::types::{CompressionType, SingeRecordDto};
use crate::utils::{get_writer_from_path, get_writer_path};
use crate::writers::Writer;

pub trait GenericWriterBackend<'a>: Send {
    fn new(writer: Box<dyn Write + Send + 'a>) -> Self;
    fn write_record(&mut self, record: &SingeRecordDto) -> io::Result<()>;
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
        Ok(Self {
            prefix: Some(prefix.to_string()),
            compression,
            backends_by_tag: HashMap::new(),
            statistics: crate::statistics::WriterStatistics::default(),
            _phantom: std::marker::PhantomData,
        })
    }

    pub fn from_writer(writer: Box<dyn Write + Send + 'a>) -> Self {
        let mut tagged_backends = HashMap::new();
        tagged_backends.insert(DEFAULT_TAG.to_string(), B::new(writer));
        Self {
            prefix: None,
            compression: CompressionType::Uncompressed,
            backends_by_tag: tagged_backends,
            statistics: crate::statistics::WriterStatistics::default(),
            _phantom: std::marker::PhantomData,
        }
    }

    fn get_backend(&mut self, tag: Option<&str>) -> io::Result<&mut B> {
        let t = tag.unwrap_or(DEFAULT_TAG);
        if !self.backends_by_tag.contains_key(t) {
            if let Some(ref prefix) = self.prefix {
                let actual_tag = if t == DEFAULT_TAG { None } else { Some(t) };
                let path = get_writer_path(prefix, self.compression, B::get_extension(), actual_tag);
                let writer = get_writer_from_path(path.as_deref(), self.compression)?;
                self.backends_by_tag.insert(t.to_string(), B::new(writer));
            } else {
                return Err(io::Error::new(io::ErrorKind::Other, format!("Cannot create tagged writer for tag '{}' without prefix", t)));
            }
        }
        Ok(self.backends_by_tag.get_mut(t).unwrap())
    }
}

impl<'a, B: GenericWriterBackend<'a>> Writer for GenericWriter<'a, B> {
    fn write(&mut self, record: &SingeRecordDto, tag: Option<&str>) -> io::Result<()>
    {
        let backend = self.get_backend(tag)?;
        backend.write_record(record)?;
        self.statistics.total_records_written += 1;
        self.statistics.total_bases_written += record.seq().len();
        if let Some(t) = tag {
            if t.starts_with("barcode_") && t != DEFAULT_TAG {
                let bc_id = &t["barcode_".len()..];
                if !bc_id.is_empty() {
                    *self.statistics.barcode_counts.entry(bc_id.to_string()).or_insert(0) += 1;
                }
            }
        }
        trace!("Writer processed record: {} with tag: {:?}", record.id(), tag);
        Ok(())
    }

    fn statistics(&self) -> crate::statistics::WriterStatistics {
        self.statistics.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn test_generic_writer_statistics() {
        let mut buffer = Vec::new();
        {
            struct TestBackend<'a> {
                _writer: Box<dyn Write + Send + 'a>,
            }
            impl<'a> GenericWriterBackend<'a> for TestBackend<'a> {
                fn new(writer: Box<dyn Write + Send + 'a>) -> Self { Self { _writer: writer } }
                fn write_record(&mut self, _record: &SingeRecordDto) -> io::Result<()> { Ok(()) }
                fn get_extension() -> &'static str { "txt" }
            }

            let mut writer: GenericWriter<'_, TestBackend<'_>> = GenericWriter::from_writer(Box::new(&mut buffer));
            let record = SingeRecordDto::with_attrs("id1", None, b"ACGT", b"!!!!");
            writer.write(&record, None).unwrap();
            
            let stats = writer.statistics();
            assert_eq!(stats.total_records_written, 1);
            assert_eq!(stats.total_bases_written, 4);
        }
    }

    #[test]
    fn test_generic_writer_barcode_statistics() {
        let temp = tempfile::tempdir().unwrap();
        let prefix = temp.path().join("out").to_str().unwrap().to_string();
        {
            struct TestBackend<'a> {
                _writer: Box<dyn Write + Send + 'a>,
            }
            impl<'a> GenericWriterBackend<'a> for TestBackend<'a> {
                fn new(writer: Box<dyn Write + Send + 'a>) -> Self { Self { _writer: writer } }
                fn write_record(&mut self, _record: &SingeRecordDto) -> io::Result<()> { Ok(()) }
                fn get_extension() -> &'static str { "txt" }
            }

            let mut writer: GenericWriter<'_, TestBackend<'_>> = GenericWriter::new(&prefix, CompressionType::Uncompressed).unwrap();
            let record = SingeRecordDto::with_attrs("id1", None, b"ACGT", b"!!!!");
            
            // Write with barcode tag
            writer.write(&record, Some("barcode_bc1")).unwrap();
            
            let stats = writer.statistics();
            assert_eq!(stats.total_records_written, 1);
            assert_eq!(stats.barcode_counts.get("bc1"), Some(&1));
            
            // Write another with same barcode
            writer.write(&record, Some("barcode_bc1")).unwrap();
            assert_eq!(writer.statistics().barcode_counts.get("bc1"), Some(&2));
            
            // Write another with different barcode
            writer.write(&record, Some("barcode_bc2")).unwrap();
            assert_eq!(writer.statistics().barcode_counts.get("bc2"), Some(&1));

            // Write without barcode tag
            writer.write(&record, None).unwrap();
            assert_eq!(writer.statistics().total_records_written, 4);
            assert_eq!(writer.statistics().barcode_counts.len(), 2);
        }
    }

    #[test]
    fn test_generic_writer_tagging_unsafe() {
        let temp = tempfile::tempdir().unwrap();
        let prefix = temp.path().join("out").to_str().unwrap().to_string();
        
        struct TestBackend<'a> {
            writer: Box<dyn Write + Send + 'a>,
        }
        impl<'a> GenericWriterBackend<'a> for TestBackend<'a> {
            fn new(writer: Box<dyn Write + Send + 'a>) -> Self { Self { writer } }
            fn write_record(&mut self, record: &SingeRecordDto) -> io::Result<()> {
                self.writer.write_all(record.id().as_bytes())
            }
            fn get_extension() -> &'static str { "txt" }
        }

        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = GenericWriter::new(&prefix, CompressionType::Uncompressed).unwrap();
            let rec1 = SingeRecordDto::with_attrs("R1", None, b"A", b"I");
            
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
    fn test_generic_writer_on_demand_default() {
        let temp = tempfile::tempdir().unwrap();
        let prefix = temp.path().join("out").to_str().unwrap().to_string();
        
        struct TestBackend<'a> {
            _writer: Box<dyn Write + Send + 'a>,
        }
        impl<'a> GenericWriterBackend<'a> for TestBackend<'a> {
            fn new(writer: Box<dyn Write + Send + 'a>) -> Self { Self { _writer: writer } }
            fn write_record(&mut self, _record: &SingeRecordDto) -> io::Result<()> { Ok(()) }
            fn get_extension() -> &'static str { "txt" }
        }

        let default_path = temp.path().join("out.txt");
        let tagged_path = temp.path().join("out_tagA.txt");

        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = GenericWriter::new(&prefix, CompressionType::Uncompressed).unwrap();
            let rec1 = SingeRecordDto::with_attrs("R1", None, b"A", b"I");
            
            // Write only with a tag
            writer.write(&rec1, Some("tagA")).unwrap();
            
            assert!(!default_path.exists(), "Default output file should not exist yet");
            assert!(tagged_path.exists(), "Tagged output file should exist");
        }
        
        assert!(!default_path.exists(), "Default output file should still not exist after drop");
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
            fn write_record(&mut self, record: &SingeRecordDto) -> io::Result<()> {
                self.writer.write_all(record.id().as_bytes())
            }
            fn get_extension() -> &'static str { "txt" }
        }

        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = GenericWriter::new(&prefix, CompressionType::Uncompressed).unwrap();
            let rec1 = SingeRecordDto::with_attrs("R1", None, b"A", b"I");
            let rec2 = SingeRecordDto::with_attrs("R2", None, b"A", b"I");
            
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
