use std::collections::HashMap;
use std::io;
use std::io::Write;
use log::trace;
use crate::types::RecordDto;
use crate::writers::Writer;

#[cfg(test)]
use std::sync::{Arc, Mutex};

pub trait GenericWriterBackend<'a>: Send {
    fn new(writer: Box<dyn Write + Send + 'a>) -> Self;
    fn write_record(&mut self, record: &RecordDto) -> io::Result<()>;
    fn get_extension() -> &'static str;
}


const DEFAULT_TAG: &str = "/\\/default/\\/";

// (bar_tag, read_index, extension) -> (writer, description)
pub type WriterFactory<'a> = Box<dyn
    Fn(Option<&str>, Option<usize>, &'static str)
        -> io::Result<(Box<dyn Write + Send + 'a>, String)
> + Send + 'a>;

pub struct GenericWriter<'a, B: GenericWriterBackend<'a>> {
    factory: WriterFactory<'a>,
    backends_by_tag: HashMap<String, B>,
    statistics: crate::statistics::WriterStatistics,
    _phantom: std::marker::PhantomData<&'a ()>,
}

impl<'a, B: GenericWriterBackend<'a>> GenericWriter<'a, B> {
    pub fn new(factory: WriterFactory<'a>) -> Self {
        Self {
            factory,
            backends_by_tag: HashMap::new(),
            statistics: crate::statistics::WriterStatistics::default(),
            _phantom: std::marker::PhantomData,
        }
    }

    fn get_backend(&mut self, tag: Option<&str>) -> io::Result<&mut B> {
        let t = tag.unwrap_or(DEFAULT_TAG);
        if !self.backends_by_tag.contains_key(t) {
            let actual_tag = if t == DEFAULT_TAG { None } else { Some(t) };
            let (writer, _description) = (self.factory)(actual_tag, None, B::get_extension())?;
            self.backends_by_tag.insert(t.to_string(), B::new(writer));
        }
        Ok(self.backends_by_tag.get_mut(t).unwrap())
    }
}

#[cfg(test)]
pub fn get_test_writer<'a, B: GenericWriterBackend<'a>>(
    buffers: Arc<Mutex<HashMap<String, Vec<u8>>>>,
) -> GenericWriter<'a, B> {
    let factory = Box::new(move |tag: Option<&str>, _read_index: Option<usize>, _extension: &'static str| {
        let t = tag.unwrap_or(DEFAULT_TAG).to_string();
        let mut bufs = buffers.lock().unwrap();
        bufs.entry(t.clone()).or_insert_with(Vec::new);
        
        // We need a writer that writes to the Vec<u8> inside the HashMap.
        // Since we want to keep the HashMap alive and accessible from the test,
        // and the factory is called multiple times, we have a problem with lifetimes if we just return &mut Vec<u8>.
        // A better way for tests is to use a custom writer that has a reference to the shared storage.
        
        struct SharedWriter {
            storage: Arc<Mutex<HashMap<String, Vec<u8>>>>,
            tag: String,
        }
        impl Write for SharedWriter {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                let mut storage = self.storage.lock().unwrap();
                storage.get_mut(&self.tag).unwrap().write(buf)
            }
            fn flush(&mut self) -> io::Result<()> {
                let mut storage = self.storage.lock().unwrap();
                storage.get_mut(&self.tag).unwrap().flush()
            }
        }

        let writer: Box<dyn Write + Send + 'a> = Box::new(SharedWriter {
            storage: buffers.clone(),
            tag: t.clone(),
        });
        
        Ok((writer, format!("memory:{}", t)))
    });

    GenericWriter::new(factory)
}

impl<'a, B: GenericWriterBackend<'a>> Writer for GenericWriter<'a, B> {
    fn write(&mut self, record: &RecordDto, tag: Option<&str>) -> io::Result<()>
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
    use crate::types::SingeRecordDto;

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

    #[test]
    fn test_generic_writer_statistics() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = get_test_writer(buffers.clone());
            let record = RecordDto::Single(
                SingeRecordDto::with_attrs("id1", None, b"ACGT", b"!!!!")
            );
            writer.write(&record, None).unwrap();
            
            let stats = writer.statistics();
            assert_eq!(stats.total_records_written, 1);
            assert_eq!(stats.total_bases_written, 4);
        }
    }

    #[test]
    fn test_generic_writer_barcode_statistics() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = get_test_writer(buffers.clone());
            let record = RecordDto::Single(SingeRecordDto::with_attrs(
                "id1", None, b"ACGT", b"!!!!"
            ));
            
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
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        
        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = get_test_writer(buffers.clone());
            let rec1 = RecordDto::Single(SingeRecordDto::with_attrs(
                "R1", None, b"A", b"I"
            ));
            
            // This tag contains characters that should be removed in the filename
            writer.write(&rec1, Some("tag/with\\slashes")).unwrap();
        }

        let bufs = buffers.lock().unwrap();
        let s = String::from_utf8(bufs.get("tag/with\\slashes").unwrap().clone()).unwrap();
        assert_eq!(s, "R1");
    }

    #[test]
    fn test_generic_writer_on_demand_default() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        
        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = get_test_writer(buffers.clone());
            let rec1 = RecordDto::Single(SingeRecordDto::with_attrs(
                "R1", None, b"A", b"I"
            ));
            
            // Write only with a tag
            writer.write(&rec1, Some("tagA")).unwrap();
            
            let bufs = buffers.lock().unwrap();
            assert!(!bufs.contains_key(DEFAULT_TAG), "Default output file should not exist yet");
            assert!(bufs.contains_key("tagA"), "Tagged output file should exist");
        }
        
        let bufs = buffers.lock().unwrap();
        assert!(!bufs.contains_key(DEFAULT_TAG), "Default output file should still not exist after drop");
    }

    #[test]
    fn test_generic_writer_tagging() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        
        {
            let mut writer: GenericWriter<'_, TestBackend<'_>> = get_test_writer(buffers.clone());
            let rec1 = RecordDto::Single(SingeRecordDto::with_attrs(
                "R1", None, b"A", b"I"
            ));
            let rec2 = RecordDto::Single(SingeRecordDto::with_attrs(
                "R2", None, b"A", b"I"
            ));
            
            writer.write(&rec1, None).unwrap();
            writer.write(&rec2, Some("tagA")).unwrap();
        }

        let bufs = buffers.lock().unwrap();
        assert!(bufs.contains_key(DEFAULT_TAG));
        assert!(bufs.contains_key("tagA"));
        
        let s1 = String::from_utf8(bufs.get(DEFAULT_TAG).unwrap().clone()).unwrap();
        assert_eq!(s1, "R1");

        let s2 = String::from_utf8(bufs.get("tagA").unwrap().clone()).unwrap();
        assert_eq!(s2, "R2");
    }
}
