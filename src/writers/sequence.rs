use std::collections::HashMap;
use std::io;
use std::io::Write;
use log::trace;
use crate::records::unified::{AnnotatedRecord, UnifiedRead, GenericRawRead};
use crate::types::FileFormat;
use crate::utils::io::{OutputType, WriterFactory};
use crate::writers::generic::{GenericWriterBackend, DEFAULT_TAG};
use crate::writers::fasta::FastaBackend;
use crate::writers::fastq::FastqBackend;
use crate::writers::bam::BamWriterBackend;
use crate::records::header::SequenceHeader;

pub enum WriterBackend<'a> {
    Fasta(FastaBackend<'a>),
    Fastq(FastqBackend<'a>),
    Bam(BamWriterBackend),
}

impl<'a> WriterBackend<'a> {
    fn write_read(&mut self, record: &AnnotatedRecord, read_index: usize) -> io::Result<()> {
        match self {
            WriterBackend::Fasta(b) => b.write_read(record, read_index),
            WriterBackend::Fastq(b) => b.write_read(record, read_index),
            WriterBackend::Bam(b) => b.write_read(record, read_index),
        }
    }

    fn new(format: FileFormat, writer: Box<dyn Write + Send + 'a>, description: &str, header: &SequenceHeader, barcode_mark: bool) -> io::Result<Self> {
        match format {
            FileFormat::Fasta => Ok(WriterBackend::Fasta(FastaBackend::new(writer, barcode_mark))),
            FileFormat::Fastq => Ok(WriterBackend::Fastq(FastqBackend::new(writer, barcode_mark))),
            FileFormat::Bam => Ok(WriterBackend::Bam(BamWriterBackend::new(description, header)?)),
        }
    }
}

pub struct SequenceWriter<'a> {
    format: FileFormat,
    factory: WriterFactory<'a>,
    header: SequenceHeader,
    backends_by_tag: HashMap<(String, Option<usize>), WriterBackend<'a>>,
    statistics: crate::statistics::WriterStatistics,
    track_length_dist: bool,
    backend_descriptions: HashMap<(String, Option<usize>), String>,
    barcode_no_split: bool,
    barcode_mark: bool,
}

impl<'a> SequenceWriter<'a> {
    pub fn new(format: FileFormat, factory: WriterFactory<'a>, track_length_dist: bool, header: SequenceHeader, barcode_no_split: bool, barcode_mark: bool) -> Self {
        Self {
            format,
            factory,
            header,
            backends_by_tag: HashMap::new(),
            statistics: crate::statistics::WriterStatistics::default(),
            track_length_dist,
            backend_descriptions: HashMap::new(),
            barcode_no_split,
            barcode_mark,
        }
    }

    fn get_backend(&mut self, output_tag: &Option<String>, read_index: Option<usize>) -> io::Result<&mut WriterBackend<'a>> {
        let key = match output_tag {
            None => (DEFAULT_TAG.to_string(), read_index),
            Some(tag) => (tag.clone(), read_index),
        };
        if !self.backends_by_tag.contains_key(&key) {
            let actual_tag = if key.0 == DEFAULT_TAG { None } else { Some(key.0.as_str()) };
            let (writer, description) = (self.factory)(actual_tag, key.1, self.format, OutputType::Reads)?;
            self.backends_by_tag.insert(key.clone(), WriterBackend::new(self.format, writer, &description, &self.header, self.barcode_mark)?);
            self.backend_descriptions.insert(key.clone(), description);
        }
        Ok(self.backends_by_tag.get_mut(&key).unwrap())
    }

    pub fn write(&mut self, annotated: &AnnotatedRecord) -> io::Result<()> {
        let record = &annotated.record;
        let paired = record.is_paired();
        let output_tag = if self.barcode_no_split {
            None
        } else {
            annotated.barcode_id.clone()
        };

        for i in 0..2 {
            if let Some(Some(read)) = record.reads.get(i) {
                if read.seq().is_empty() { continue; }
                let read_index = if paired { Some(i) } else { None };
                let backend = self.get_backend(&output_tag, read_index)?;
                backend.write_read(annotated, i)?;
                self.update_statistics(read, &annotated.barcode_id, read_index);
            }
        }

        if let Some(bc) = &annotated.barcode_id {
            *self.statistics.barcode_counts.entry(bc.to_string()).or_insert(0) += 1;
        }
        trace!("Writer processed record: {} with barcode: {:?}", record.record_id(), annotated.barcode_id);
        Ok(())
    }

    fn update_statistics(&mut self, read: &UnifiedRead, barcode_id: &Option<String>, read_index: Option<usize>) {
        self.statistics.total_records_written += 1;
        let seq_len = read.seq().len();
        self.statistics.total_bases_written += seq_len;
        if self.track_length_dist {
            let t_str = match barcode_id {
                Some(bc) => bc,
                None => DEFAULT_TAG,
            };
            let key = (t_str.to_string(), read_index);
            let dist = self.statistics.length_distribution.entry(key).or_insert_with(HashMap::new);
            *dist.entry(seq_len).or_insert(0) += 1;
        }
    }

    pub fn write_length_distribution(&self) -> io::Result<()> {
        if !self.track_length_dist {
            return Ok(());
        }
        for (key, dist) in &self.statistics.length_distribution {
            let actual_tag = if key.0 == DEFAULT_TAG { None } else { Some(key.0.as_str()) };
            let (writer, description) = (self.factory)(actual_tag, key.1, self.format, OutputType::LengthDistribution)?;
            
            if description == "STDOUT" {
                continue;
            }

            let mut file = io::BufWriter::new(writer);
            writeln!(file, "Readlength\tCount")?;
            
            let mut lengths: Vec<_> = dist.keys().cloned().collect();
            lengths.sort_unstable();
            
            for len in lengths {
                let count = dist.get(&len).unwrap();
                if *count > 0 {
                    writeln!(file, "{}\t{}", len, count)?;
                }
            }
        }
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
    use crate::records::unified::{UnifiedRead, UnifiedRecord};

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


    pub fn get_test_writer<'a>(
        format: FileFormat,
        buffers: Arc<Mutex<HashMap<String, Vec<u8>>>>,
        track_length_dist: bool,
    ) -> SequenceWriter<'a> {
        let factory = Box::new(move |tag: Option<&str>, read_index: Option<usize>, _format: FileFormat, output_type: OutputType| {
            let t = tag.unwrap_or(DEFAULT_TAG).to_string();
            let base_key = match read_index {
                Some(idx) => format!("{}_R{}", t, idx + 1),
                None => t.clone(),
            };
            
            let buffer_key = match output_type {
                OutputType::LengthDistribution => format!("{}.lengthdist", base_key),
                _ => base_key,
            };

            let mut bufs = buffers.lock().unwrap();
            bufs.entry(buffer_key.clone()).or_insert_with(Vec::new);

            let writer: Box<dyn Write + Send + 'a> = Box::new(SharedWriter {
                storage: buffers.clone(),
                tag: buffer_key.clone(),
            });
            
            Ok((writer, format!("memory:{}", buffer_key)))
        });

        SequenceWriter::new(format, factory, track_length_dist, SequenceHeader::Default, false, false)
    }

    fn wrap(rec: UnifiedRecord, barcode_id: Option<&str>) -> AnnotatedRecord {
        AnnotatedRecord {
            record: rec,
            barcode_id: barcode_id.map(|s| s.to_string()),
            removal_tags: [vec![], vec![]],
            umi_tags: vec![],
            paired_read_dropped: false,
        }
    }

    #[test]
    fn test_sequence_writer_statistics() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fastq, buffers.clone(), false);
            let record = UnifiedRecord::single(
                UnifiedRead::new_fastq("id1", None, b"ACGT", b"!!!!")
            );
            writer.write(&wrap(record, None)).unwrap();
            
            let stats = writer.statistics();
            assert_eq!(stats.total_records_written, 1);
            assert_eq!(stats.total_bases_written, 4);
        }
    }

    #[test]
    fn test_sequence_writer_barcode_statistics() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fastq, buffers.clone(), false);
            let record = UnifiedRecord::single(UnifiedRead::new_fastq(
                "id1", None, b"ACGT", b"!!!!"
            ));
            
            writer.write(&wrap(record.clone(), Some("bc1"))).unwrap();
            
            let stats = writer.statistics();
            assert_eq!(stats.total_records_written, 1);
            assert_eq!(stats.barcode_counts.get("bc1"), Some(&1));
            
            writer.write(&wrap(record.clone(), Some("bc1"))).unwrap();
            assert_eq!(writer.statistics().barcode_counts.get("bc1"), Some(&2));
            
            writer.write(&wrap(record.clone(), Some("bc2"))).unwrap();
            assert_eq!(writer.statistics().barcode_counts.get("bc2"), Some(&1));

            writer.write(&wrap(record.clone(), None)).unwrap();
            assert_eq!(writer.statistics().total_records_written, 4);
            assert_eq!(writer.statistics().barcode_counts.len(), 2);
        }
    }

    #[test]
    fn test_fasta_writer_cursor() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fasta, buffers.clone(), false);
            let record = UnifiedRecord::single(UnifiedRead::new_fasta(
                "id1", Some("desc1"), b"ACGT"
            ));
            writer.write(&wrap(record, None)).expect("Failed to write record");
        }
        
        let bufs = buffers.lock().unwrap();
        let result = String::from_utf8(bufs.get(DEFAULT_TAG).unwrap().clone()).expect("Not valid UTF-8");
        assert!(result.contains(">id1 desc1"));
        assert!(result.contains("ACGT"));
        assert!(!result.contains("!!!!"));
    }

    #[test]
    fn test_fasta_to_fastq_conversion() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fastq, buffers.clone(), false);
            let record = UnifiedRecord::single(UnifiedRead::new_fasta(
                "id1", Some("desc1"), b"ACGT"
            ));
            writer.write(&wrap(record, None)).unwrap();
        }
        
        let bufs = buffers.lock().unwrap();
        let result = String::from_utf8(bufs.get(DEFAULT_TAG).unwrap().clone()).expect("Not valid UTF-8");
        assert!(result.contains("@id1 desc1"));
        assert!(result.contains("ACGT"));
        assert!(result.contains("IIII"));
    }

    #[test]
    fn test_bam_to_fastq_extraction() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fastq, buffers.clone(), false);
            
            use rust_htslib::bam::record::{Record, Aux};
            let mut bam_rec = Record::new();
            bam_rec.set(b"read1", None, b"ACGT", b"!!!!");
            bam_rec.push_aux(b"CO", Aux::String("packed_desc")).unwrap();
            
            let record = UnifiedRecord::single(UnifiedRead::from_bam(bam_rec));
            writer.write(&wrap(record, None)).unwrap();
        }
        
        let bufs = buffers.lock().unwrap();
        let result = String::from_utf8(bufs.get(DEFAULT_TAG).unwrap().clone()).expect("Not valid UTF-8");
        assert!(result.contains("ACGT"));
        assert!(result.contains("!!!!"));
    }

    #[test]
    fn test_barcode_no_split() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            // Set barcode_no_split = true
            let mut writer = SequenceWriter::new(
                FileFormat::Fastq,
                get_test_factory(buffers.clone()),
                false,
                SequenceHeader::Default,
                true,
                false
            );
            
            let record = UnifiedRecord::single(UnifiedRead::new_fastq(
                "id1", None, b"ACGT", b"!!!!"
            ));
            
            writer.write(&wrap(record, Some("bc1"))).unwrap();
            drop(writer); // Flush everything
            
            let bufs = buffers.lock().unwrap();
            // Should be in DEFAULT_TAG, not in "barcode_bc1"
            assert!(bufs.contains_key(DEFAULT_TAG));
            assert!(!bufs.contains_key("barcode_bc1"));
            
            let result = String::from_utf8(bufs.get(DEFAULT_TAG).unwrap().clone()).unwrap();
            println!("DEBUG: result: {:?}", result);
            assert!(result.contains("id1"));
        }
    }

    #[test]
    fn test_barcode_mark() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            // Set barcode_mark = true
            let mut writer = SequenceWriter::new(
                FileFormat::Fastq,
                get_test_factory(buffers.clone()),
                false,
                SequenceHeader::Default,
                false,
                true
            );
            
            let mut record = wrap(UnifiedRecord::single(UnifiedRead::new_fastq(
                "id1", None, b"ACGT", b"!!!!"
            )), None);
            record.barcode_id = Some("bc1".to_string());
            
            writer.write(&record).unwrap();
            drop(writer); // Flush everything
            
            let bufs = buffers.lock().unwrap();
            let result = String::from_utf8(bufs.get("bc1").unwrap().clone()).unwrap();
            // Should contain BC:bc1
            assert!(result.contains("@id1 barcode_bc1"));
        }
    }

    fn get_test_factory<'a>(buffers: Arc<Mutex<HashMap<String, Vec<u8>>>>) -> WriterFactory<'a> {
        Box::new(move |tag: Option<&str>, read_index: Option<usize>, _format: FileFormat, output_type: OutputType| {
            let t = tag.unwrap_or(DEFAULT_TAG).to_string();
            let base_key = match read_index {
                Some(idx) => format!("{}_R{}", t, idx + 1),
                None => t.clone(),
            };
            
            let buffer_key = match output_type {
                OutputType::LengthDistribution => format!("{}.lengthdist", base_key),
                _ => base_key,
            };

            let mut bufs = buffers.lock().unwrap();
            bufs.entry(buffer_key.clone()).or_insert_with(Vec::new);

            let writer = Box::new(SharedWriter {
                storage: buffers.clone(),
                tag: buffer_key.clone(),
            });
            
            Ok((writer, format!("memory:{}", buffer_key)))
        })
    }
}
