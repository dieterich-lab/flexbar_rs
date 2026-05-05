use std::collections::HashMap;
use std::io;
use std::io::Write;
use log::trace;
use crate::records::unified::{UnifiedRecord, UnifiedRead, GenericRawRead};
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
    fn write_read(&mut self, record: &UnifiedRead) -> io::Result<()> {
        match self {
            WriterBackend::Fasta(b) => b.write_read(record),
            WriterBackend::Fastq(b) => b.write_read(record),
            WriterBackend::Bam(b) => b.write_read(record),
        }
    }

    fn new(format: FileFormat, writer: Box<dyn Write + Send + 'a>, description: &str, header: &SequenceHeader) -> io::Result<Self> {
        match format {
            FileFormat::Fasta => Ok(WriterBackend::Fasta(FastaBackend::new(writer))),
            FileFormat::Fastq => Ok(WriterBackend::Fastq(FastqBackend::new(writer))),
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
}

impl<'a> SequenceWriter<'a> {
    pub fn new(format: FileFormat, factory: WriterFactory<'a>, track_length_dist: bool, header: SequenceHeader) -> Self {
        Self {
            format,
            factory,
            header,
            backends_by_tag: HashMap::new(),
            statistics: crate::statistics::WriterStatistics::default(),
            track_length_dist,
            backend_descriptions: HashMap::new(),
        }
    }

    fn get_backend(&mut self, tag: Option<&str>, read_index: Option<usize>) -> io::Result<&mut WriterBackend<'a>> {
        let t = tag.unwrap_or(DEFAULT_TAG).to_string();
        let key = (t, read_index);
        if !self.backends_by_tag.contains_key(&key) {
            let actual_tag = if key.0 == DEFAULT_TAG { None } else { Some(key.0.as_str()) };
            let (writer, description) = (self.factory)(actual_tag, key.1, self.format, OutputType::Reads)?;
            self.backends_by_tag.insert(key.clone(), WriterBackend::new(self.format, writer, &description, &self.header)?);
            self.backend_descriptions.insert(key.clone(), description);
        }
        Ok(self.backends_by_tag.get_mut(&key).unwrap())
    }

    pub fn write(&mut self, record: &UnifiedRecord, tag: Option<&str>) -> io::Result<()> {
        let paired = record.is_paired();
        for i in 0..2 {
            if let Some(Some(read)) = record.reads.get(i) {
                if read.seq().is_empty() { continue; }
                let read_index = if paired { Some(i) } else { None };
                let backend = self.get_backend(tag, read_index)?;
                backend.write_read(&read)?;
                self.update_statistics(read, tag, read_index);
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

    fn update_statistics(&mut self, read: &UnifiedRead, tag: Option<&str>, read_index: Option<usize>) {
        self.statistics.total_records_written += 1;
        let seq_len = read.seq().len();
        self.statistics.total_bases_written += seq_len;
        if self.track_length_dist {
            let t_str = tag.unwrap_or(DEFAULT_TAG).to_string();
            let key = (t_str.clone(), read_index);
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
    use crate::records::unified::{UnifiedRead};

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

        SequenceWriter::new(format, factory, track_length_dist, SequenceHeader::Default)
    }

    #[test]
    fn test_sequence_writer_statistics() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        {
            let mut writer = get_test_writer(FileFormat::Fastq, buffers.clone(), false);
            let record = UnifiedRecord::single(
                UnifiedRead::new_fastq("id1", None, b"ACGT", b"!!!!")
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
            let mut writer = get_test_writer(FileFormat::Fastq, buffers.clone(), false);
            let record = UnifiedRecord::single(UnifiedRead::new_fastq(
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
            let mut writer = get_test_writer(FileFormat::Fasta, buffers.clone(), false);
            let record = UnifiedRecord::single(UnifiedRead::new_fasta(
                "id1", Some("desc1"), b"ACGT"
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
    fn test_sequence_writer_length_distribution() {
        let buffers = Arc::new(Mutex::new(HashMap::new()));
        let mut writer = get_test_writer(FileFormat::Fasta, buffers.clone(), true);
        
        let record1 = UnifiedRecord::single(UnifiedRead::new_fasta("id1", None, b"ACGT")); // len 4
        let record2 = UnifiedRecord::single(UnifiedRead::new_fasta("id2", None, b"ACGTA")); // len 5
        let record3 = UnifiedRecord::single(UnifiedRead::new_fasta("id3", None, b"ACGT")); // len 4
        
        writer.write(&record1, None).unwrap();
        writer.write(&record2, None).unwrap();
        writer.write(&record3, None).unwrap();
        
        writer.write_length_distribution().unwrap();
        
        let bufs = buffers.lock().unwrap();
        let dist_key = format!("{}.lengthdist", DEFAULT_TAG);
        let content = String::from_utf8(bufs.get(&dist_key).unwrap().clone()).expect("Not valid UTF-8");
        
        let expected = "Readlength\tCount\n4\t2\n5\t1\n";
        assert_eq!(content, expected);
    }
}
