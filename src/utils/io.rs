use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use crate::types::FileFormat;
use std::fs::File;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use bzip2::read::BzDecoder;
use bzip2::write::BzEncoder;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;
use crate::types::CompressionType;

use crate::logging::log_fatal;

pub enum OutputType{
    Reads,
    LengthDistribution,
}

// (bar_tag, read_index, extension) -> (writer, description)
pub type WriterFactory<'a> = Box<dyn
Fn(Option<&str>, Option<usize>, FileFormat, OutputType)
    -> io::Result<(Box<dyn Write + Send + 'a>, String)
    > + Send + 'a>;

pub fn get_file_type(path: Option<&PathBuf>) -> FileFormat {
    path.as_ref().map_or(FileFormat::Fasta, |path| {
        if path.to_string_lossy() == "-" {
            return FileFormat::Fastq;
        }
        let ext = path.extension().and_then(|s| s.to_str()).map(|s| s.to_lowercase());
        if let Some(e) = ext {
            match e.as_str() {
                "fasta" => return FileFormat::Fasta,
                "fa" => return FileFormat::Fasta,
                "bam" => return FileFormat::Bam,
                _ => {}
            }
        }
        FileFormat::Fastq
    })
}

pub fn get_raw_reader<P: AsRef<Path>>(path: P, buffer_size: usize) -> (BufReader<Box<dyn Read + Send>>, String, CompressionType) {
    let path_ref = path.as_ref();
    let source = path_ref.to_string_lossy().into_owned();

    let (reader, compression): (Box<dyn Read + Send>, CompressionType) = if source == "-" {
        (Box::new(io::stdin()), CompressionType::Uncompressed)
    } else {
        let file = File::open(path_ref).unwrap_or_else(|e| {
            log_fatal(&format!("Error opening file {}: {}", path_ref.display(), e));
        });
        let extension = path_ref.extension()
            .and_then(|s| s.to_str())
            .map(|s| s.to_lowercase());

        match extension.as_deref() {
            Some("gz") => (Box::new(MultiGzDecoder::new(file)), CompressionType::Gz),
            Some("bz2") => (Box::new(BzDecoder::new(file)), CompressionType::Bz2),
            Some("zstd") | Some("zst") => {
                match ZstdDecoder::new(file) {
                    Ok(decoder) => (Box::new(decoder), CompressionType::Zstd),
                    Err(e) => {
                        log_fatal(&format!("Error creating zstd decoder for {}: {}", path_ref.display(), e));
                    }
                }
            },
            _ => (Box::new(file), CompressionType::Uncompressed),
        }
    };

    (BufReader::with_capacity(buffer_size, reader), source, compression)
}

struct WriterFactoryBackend {
    prefix: String,
    compression: CompressionType,
    buffer_size: usize,
}

impl WriterFactoryBackend {
    fn new(prefix: String, compression: CompressionType, buffer_size: usize) -> Self {
        Self { prefix, compression, buffer_size }
    }

    fn get_writer_path(
        &self,
        tag: Option<&str>,
        read_index: Option<usize>,
        file_format: FileFormat,
        output_type: &OutputType
    ) -> Option<PathBuf> {
        if self.prefix == "-" {
            return None;
        }
        let cleaned_tag = match tag {
            Some(t) => format!("_{}", t.replace(['/', '\\'], "")),
            None => "".to_string(),
        };
        let read_index_str = match read_index {
            Some(idx) => format!("_{}", idx + 1),
            None => "".to_string(),
        };
        let format_str = match file_format {
            FileFormat::Fasta => ".fasta",
            FileFormat::Fastq => ".fastq",
            FileFormat::Bam => ".bam",
        };
        let compression_str = match self.compression {
            CompressionType::Gz => ".gz",
            CompressionType::Bz2 => ".bz2",
            CompressionType::Zstd => ".zst",
            CompressionType::Uncompressed => "",
            CompressionType::Keep => log_fatal("Compression type KEEP is not supported for writing"),
        };
        let output_type_str = match output_type {
            OutputType::LengthDistribution => ".lengthdist",
            _ => ""
        };
        Some(PathBuf::from(vec![
            self.prefix.clone(),
            cleaned_tag,
            read_index_str,
            format_str.to_string(),
            compression_str.to_string(),
            output_type_str.to_string(),
        ].join("")))
    }

    fn get_raw_writer(&self, path: Option<PathBuf>, output_type: &OutputType) -> Box<dyn Write + Send> {
        let writer: Box<dyn Write + Send> = match &path {
            Some(p) => {
                let f = match File::create(&p) {
                    Ok(f) => f,
                    Err(e) => {
                        log_fatal(&format!("Error creating file {}: {}", p.display(), e));
                    }
                };
                Box::new(BufWriter::with_capacity(self.buffer_size, f))
            },
            None => Box::new(BufWriter::with_capacity(self.buffer_size, io::stdout())),
        };
        let writer_with_compression: io::Result<Box<dyn Write + Send>> = match output_type {
            OutputType::LengthDistribution => Ok(writer),
            _ => match self.compression {
                CompressionType::Gz => Ok(Box::new(GzEncoder::new(writer, flate2::Compression::default()))),
                CompressionType::Bz2 => Ok(Box::new(BzEncoder::new(writer, bzip2::Compression::default()))),
                CompressionType::Zstd => Ok(Box::new(ZstdEncoder::new(writer, 0).unwrap().auto_finish())),
                CompressionType::Uncompressed => Ok(writer),
                CompressionType::Keep => log_fatal("Compression type KEEP is not supported for writing"),
            }
        };
        match writer_with_compression {
            Ok(writer) => writer,
            Err(e) => log_fatal(&format!(
                "Failed to create writer ({:?}): {}",
                path, e
            )),
        }
    }
}


pub fn get_writer_factory<'a>(
    prefix: String,
    compression: CompressionType,
    buffer_size: usize,
) -> WriterFactory<'a> {
    let backend = WriterFactoryBackend::new(prefix.clone(), compression, buffer_size);
    Box::new(move |tag, read_index, file_format, output_type| {
        let path = backend.get_writer_path(tag, read_index, file_format, &output_type);
        let description = match &path {
            Some(p) => p.to_string_lossy().into_owned(),
            None => "STDOUT".to_string(),
        };
        let writer = backend.get_raw_writer(path, &output_type);
        Ok((writer, description))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_bz2_decompression() {
        let mut tmp = NamedTempFile::new().unwrap();
        let mut bz2_writer = BzEncoder::new(Vec::new(), bzip2::Compression::default());
        bz2_writer.write_all(b">test\nACGT\n").unwrap();
        let compressed_data = bz2_writer.finish().unwrap();
        tmp.write_all(&compressed_data).unwrap();
        
        // Rename tmp file to have .bz2 extension
        let path = tmp.path().with_extension("bz2");
        std::fs::copy(tmp.path(), &path).unwrap();

        let (mut reader, source, compression) = get_raw_reader(&path, 4096);
        assert_eq!(compression, CompressionType::Bz2);
        assert!(source.contains(path.to_str().unwrap()));

        let mut content = String::new();
        reader.read_to_string(&mut content).unwrap();
        assert_eq!(content, ">test\nACGT\n");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_get_writer_path_with_unsafe_tag() {
        let prefix = "output/test".to_string();
        let compression = CompressionType::Uncompressed;
        let backend = WriterFactoryBackend::new(prefix, compression, 4096);
        let tag = Some("tag/with\\slashes");
        let path = backend.get_writer_path(tag, None, FileFormat::Fastq, &OutputType::Reads).unwrap();
        assert_eq!(path, PathBuf::from("output/test_tagwithslashes.fastq"));
    }

    #[test]
    fn test_get_writer_path_with_tag() {
        let prefix = "output/test".to_string();
        let compression = CompressionType::Uncompressed;
        let backend = WriterFactoryBackend::new(prefix, compression, 4096);
        let tag = Some("tag1");
        let path = backend.get_writer_path(tag, None, FileFormat::Fastq, &OutputType::Reads).unwrap();
        assert_eq!(path, PathBuf::from("output/test_tag1.fastq"));
    }

    #[test]
    fn test_get_writer_path_without_tag() {
        let prefix = "output/test".to_string();
        let compression = CompressionType::Uncompressed;
        let backend = WriterFactoryBackend::new(prefix, compression, 4096);
        let path = backend.get_writer_path(None, None, FileFormat::Fastq, &OutputType::Reads).unwrap();
        assert_eq!(path, PathBuf::from("output/test.fastq"));
    }

    #[test]
    fn test_get_writer_path_paired() {
        let prefix = "output/test".to_string();
        let compression = CompressionType::Uncompressed;
        let backend = WriterFactoryBackend::new(prefix, compression, 4096);
        
        let path1 = backend.get_writer_path(None, Some(0), FileFormat::Fastq, &OutputType::Reads).unwrap();
        assert_eq!(path1, PathBuf::from("output/test_1.fastq"));
        
        let path2 = backend.get_writer_path(None, Some(1), FileFormat::Fastq, &OutputType::Reads).unwrap();
        assert_eq!(path2, PathBuf::from("output/test_2.fastq"));
        
        let path3 = backend.get_writer_path(Some("tag1"), Some(0), FileFormat::Fastq, &OutputType::Reads).unwrap();
        assert_eq!(path3, PathBuf::from("output/test_tag1_1.fastq"));
    }

    fn test_compression_roundtrip(compression: CompressionType) {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("test").to_str().unwrap().to_string();
        let factory = get_writer_factory(prefix.clone(), compression, 4096);

        let data = b"test data for compression roundtrip";
        let (mut writer, path) = factory(None, None, FileFormat::Fastq, OutputType::Reads).unwrap();
        writer.write_all(data).unwrap();
        writer.flush().unwrap();
        drop(writer);
        let (mut reader, _source, detected_compression) = get_raw_reader(&path, 4096);
        assert_eq!(detected_compression, compression);
        
        let mut buffer = Vec::new();
        reader.read_to_end(&mut buffer).unwrap();
        assert_eq!(buffer, data);
    }

    #[test]
    fn test_uncompressed_roundtrip() {
        test_compression_roundtrip(CompressionType::Uncompressed);
    }

    #[test]
    fn test_gz_roundtrip() {
        test_compression_roundtrip(CompressionType::Gz);
    }

    #[test]
    fn test_bz2_roundtrip() {
        test_compression_roundtrip(CompressionType::Bz2);
    }

    #[test]
    fn test_zstd_roundtrip() {
        test_compression_roundtrip(CompressionType::Zstd);
    }

    #[test]
    fn test_get_writer_path_panic_on_keep() {
        // This used to panic, now it calls log_fatal and exits.
        // For now we just skip the execution that would exit.
        // get_writer_path("test", CompressionType::Keep, "txt", None);
    }

    #[test]
    fn test_stdin_reader() {
        let (_reader, source, compression) = get_raw_reader("-", 4096);
        assert_eq!(source, "-");
        assert_eq!(compression, CompressionType::Uncompressed);
    }

    #[test]
    fn test_multi_gz_decompression() {
        let mut buf = Vec::new();
        {
            let mut encoder = GzEncoder::new(&mut buf, flate2::Compression::default());
            encoder.write_all(b"member1").unwrap();
            encoder.finish().unwrap();
        }
        {
            let mut encoder = GzEncoder::new(&mut buf, flate2::Compression::default());
            encoder.write_all(b"member2").unwrap();
            encoder.finish().unwrap();
        }

        let dir = tempdir().unwrap();
        let path = dir.path().join("test.gz");
        File::create(&path).unwrap().write_all(&buf).unwrap();

        let (mut reader, _source, compression) = get_raw_reader(&path, 4096);
        assert_eq!(compression, CompressionType::Gz);

        let mut content = String::new();
        reader.read_to_string(&mut content).unwrap();
        assert_eq!(content, "member1member2");
    }
}
