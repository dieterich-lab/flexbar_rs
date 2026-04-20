use std::collections::HashMap;
use std::io;
use std::io::Write;
use log::trace;
use crate::types::{RecordDto, SingeRecordDto, FileFormat};
use crate::writers::generic::{WriterFactory, DEFAULT_TAG, GenericWriterBackend};
use crate::writers::fasta::FastaBackend;
use crate::writers::fastq::FastqBackend;

pub enum WriterBackend<'a> {
    Fasta(FastaBackend<'a>),
    Fastq(FastqBackend<'a>),
}

impl<'a> WriterBackend<'a> {
    fn write_read(&mut self, record: &SingeRecordDto) -> io::Result<()> {
        match self {
            WriterBackend::Fasta(b) => b.write_read(record),
            WriterBackend::Fastq(b) => b.write_read(record),
        }
    }

    fn get_extension(format: FileFormat) -> &'static str {
        match format {
            FileFormat::Fasta => "fasta",
            FileFormat::Fastq => "fastq",
        }
    }

    fn new(format: FileFormat, writer: Box<dyn Write + Send + 'a>) -> Self {
        match format {
            FileFormat::Fasta => WriterBackend::Fasta(FastaBackend::new(writer)),
            FileFormat::Fastq => WriterBackend::Fastq(FastqBackend::new(writer)),
        }
    }
}

pub struct SequenceWriter<'a> {
    format: FileFormat,
    factory: WriterFactory<'a>,
    backends_by_tag: HashMap<(String, Option<usize>), WriterBackend<'a>>,
    statistics: crate::statistics::WriterStatistics,
}

impl<'a> SequenceWriter<'a> {
    pub fn new(format: FileFormat, factory: WriterFactory<'a>) -> Self {
        Self {
            format,
            factory,
            backends_by_tag: HashMap::new(),
            statistics: crate::statistics::WriterStatistics::default(),
        }
    }

    fn get_backend(&mut self, tag: Option<&str>, read_index: Option<usize>) -> io::Result<&mut WriterBackend<'a>> {
        let t = tag.unwrap_or(DEFAULT_TAG).to_string();
        let key = (t, read_index);
        if !self.backends_by_tag.contains_key(&key) {
            let actual_tag = if key.0 == DEFAULT_TAG { None } else { Some(key.0.as_str()) };
            let extension = WriterBackend::get_extension(self.format);
            let (writer, _description) = (self.factory)(actual_tag, key.1, extension)?;
            self.backends_by_tag.insert(key.clone(), WriterBackend::new(self.format, writer));
        }
        Ok(self.backends_by_tag.get_mut(&key).unwrap())
    }

    pub fn write(&mut self, record: &RecordDto, tag: Option<&str>) -> io::Result<()> {
        let paired = record.is_paired();
        for i in 0..2 {
            if let Some(Some(read)) = record.reads.get(i) {
                if read.seq().is_empty() { continue; }
                let read_index = if paired { Some(i) } else { None };
                let backend = self.get_backend(tag, read_index)?;
                backend.write_read(&read)?;
                self.statistics.total_records_written += 1;
                self.statistics.total_bases_written += read.seq().len();
            }
        }

        if let Some(t) = tag {
            if t.starts_with("barcode_") && t != DEFAULT_TAG {
                let bc_id = &t["barcode_".len()..];
                if !bc_id.is_empty() {
                    *self.statistics.barcode_counts.entry(bc_id.to_string()).or_insert(0) += 1;
                }
            }
        }
        trace!("Writer processed record: {} with tag: {:?}", record.record_id(), tag);
        Ok(())
    }

    pub fn statistics(&self) -> crate::statistics::WriterStatistics {
        self.statistics.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};
    use crate::types::SingeRecordDto;

    pub fn get_test_writer<'a>(
        format: FileFormat,
        buffers: Arc<Mutex<HashMap<String, Vec<u8>>>>,
    ) -> SequenceWriter<'a> {
        let factory = Box::new(move |tag: Option<&str>, read_index: Option<usize>, _extension: &'static str| {
            let t = tag.unwrap_or(DEFAULT_TAG).to_string();
            let buffer_key = match read_index {
                Some(idx) => format!("{}_R{}", t, idx + 1),
                None => t.clone(),
            };
            let mut bufs = buffers.lock().unwrap();
            bufs.entry(buffer_key.clone()).or_insert_with(Vec::new);
            
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
                tag: buffer_key.clone(),
            });
            
            Ok((writer, format!("memory:{}", buffer_key)))
        });

        SequenceWriter::new(format, factory)
    }

    #[test]
    fn test_sequence_writer_statistics() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fasta, buffers.clone());
            let record = RecordDto::single(
                SingeRecordDto::with_attrs("id1", None, b"ACGT", b"!!!!")
            );
            writer.write(&record, None).unwrap();
            
            let stats = writer.statistics();
            assert_eq!(stats.total_records_written, 1);
            assert_eq!(stats.total_bases_written, 4);
        }
    }

    #[test]
    fn test_sequence_writer_barcode_statistics() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fasta, buffers.clone());
            let record = RecordDto::single(SingeRecordDto::with_attrs(
                "id1", None, b"ACGT", b"!!!!"
            ));
            
            writer.write(&record, Some("barcode_bc1")).unwrap();
            
            let stats = writer.statistics();
            assert_eq!(stats.total_records_written, 1);
            assert_eq!(stats.barcode_counts.get("bc1"), Some(&1));
            
            writer.write(&record, Some("barcode_bc1")).unwrap();
            assert_eq!(writer.statistics().barcode_counts.get("bc1"), Some(&2));
            
            writer.write(&record, Some("barcode_bc2")).unwrap();
            assert_eq!(writer.statistics().barcode_counts.get("bc2"), Some(&1));

            writer.write(&record, None).unwrap();
            assert_eq!(writer.statistics().total_records_written, 4);
            assert_eq!(writer.statistics().barcode_counts.len(), 2);
        }
    }

    #[test]
    fn test_fasta_writer_cursor() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fasta, buffers.clone());
            let record = RecordDto::single(SingeRecordDto::with_attrs(
                "id1", Some("desc1"), b"ACGT", b"!!!!"
            ));
            writer.write(&record, None).expect("Failed to write record");
        }
        
        let bufs = buffers.lock().unwrap();
        let result = String::from_utf8(bufs.get(DEFAULT_TAG).unwrap().clone()).expect("Not valid UTF-8");
        assert!(result.contains(">id1 desc1"));
        assert!(result.contains("ACGT"));
        assert!(!result.contains("!!!!"));
    }

    #[test]
    fn test_fastq_writer_cursor() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fastq, buffers.clone());
            let record = RecordDto::single(SingeRecordDto::with_attrs(
                "id1", Some("desc1"), b"ACGT", b"!!!!"
            ));
            writer.write(&record, None).expect("Failed to write record");
        }
        
        let bufs = buffers.lock().unwrap();
        let result = String::from_utf8(bufs.get(DEFAULT_TAG).unwrap().clone()).expect("Not valid UTF-8");
        assert!(result.contains("@id1 desc1"));
        assert!(result.contains("ACGT"));
        assert!(result.contains("!!!!"));
    }
}
