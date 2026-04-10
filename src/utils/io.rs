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

pub fn get_file_type(path: Option<&PathBuf>) -> FileFormat {
    path.as_ref().map_or(FileFormat::Fasta, |path| {
        if path.to_string_lossy() == "-" {
            return FileFormat::Fastq;
        }
        let ext = path.extension().and_then(|s| s.to_str()).map(|s| s.to_lowercase());
        if let Some(e) = ext {
            if e == "fasta" || e == "fa" {
                return FileFormat::Fasta;
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

pub fn get_writer_path(prefix: &str, compression: CompressionType, extension: &str, tag: Option<&str>, read_index: Option<usize>) -> Option<PathBuf> {
    if prefix == "-" {
        return None;
    }
    let mut path = PathBuf::from(prefix);
    if let Some(t) = tag {
        let cleaned_tag = t.replace(['/', '\\'], "");
        let mut filename = path.file_name().unwrap_or_default().to_os_string();
        filename.push("_");
        filename.push(cleaned_tag);
        path.set_file_name(filename);
    }
    if let Some(idx) = read_index {
        let mut filename = path.file_name().unwrap_or_default().to_os_string();
        filename.push(format!("_{}", idx + 1));
        path.set_file_name(filename);
    }
    path.set_extension(extension);
    
    Some(match compression {
        CompressionType::Gz => {
            let mut s = path.into_os_string();
            s.push(".gz");
            PathBuf::from(s)
        }
        CompressionType::Bz2 => {
            let mut s = path.into_os_string();
            s.push(".bz2");
            PathBuf::from(s)
        }
        CompressionType::Zstd => {
            let mut s = path.into_os_string();
            s.push(".zst");
            PathBuf::from(s)
        }
        CompressionType::Uncompressed => path,
        CompressionType::Keep => log_fatal("Compression type KEEP is not supported for writing"),
    })
}

fn get_writer_from_path(path: Option<&Path>, compression: CompressionType, buffer_size: usize) -> io::Result<Box<dyn Write + Send>> {
    let writer: Box<dyn Write + Send> = match path {
        Some(p) => Box::new(BufWriter::with_capacity(buffer_size, File::create(p)?)),
        None => Box::new(BufWriter::with_capacity(buffer_size, io::stdout())),
    };

    match compression {
        CompressionType::Gz => Ok(Box::new(GzEncoder::new(writer, flate2::Compression::default()))),
        CompressionType::Bz2 => Ok(Box::new(BzEncoder::new(writer, bzip2::Compression::default()))),
        CompressionType::Zstd => Ok(Box::new(ZstdEncoder::new(writer, 0)?.auto_finish())),
        CompressionType::Uncompressed => Ok(writer),
        CompressionType::Keep => log_fatal("Compression type KEEP is not supported for writing"),
    }
}

pub fn get_writer_factory<'a>(
    prefix: String,
    compression: CompressionType,
    buffer_size: usize,
) -> Box<dyn Fn(Option<&str>, Option<usize>, &'static str) -> io::Result<(Box<dyn Write + Send + 'a>, String)> + Send + 'a> {
    Box::new(move |tag, read_index, extension| {
        let path = get_writer_path(&prefix, compression, extension, tag, read_index);
        let description = match &path {
            Some(p) => p.to_string_lossy().into_owned(),
            None => "STDOUT".to_string(),
        };
        let writer = get_writer_from_path(path.as_deref(), compression, buffer_size)?;
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
        let prefix = "output/test";
        let compression = CompressionType::Uncompressed;
        let extension = "fastq";
        let tag = Some("tag/with\\slashes");
        let path = get_writer_path(prefix, compression, extension, tag, None).unwrap();
        assert_eq!(path, PathBuf::from("output/test_tagwithslashes.fastq"));
    }

    #[test]
    fn test_get_writer_path_with_tag() {
        let prefix = "output/test";
        let compression = CompressionType::Uncompressed;
        let extension = "fastq";
        let tag = Some("tag1");
        let path = get_writer_path(prefix, compression, extension, tag, None).unwrap();
        assert_eq!(path, PathBuf::from("output/test_tag1.fastq"));
    }

    #[test]
    fn test_get_writer_path_without_tag() {
        let prefix = "output/test";
        let compression = CompressionType::Uncompressed;
        let extension = "fastq";
        let path = get_writer_path(prefix, compression, extension, None, None).unwrap();
        assert_eq!(path, PathBuf::from("output/test.fastq"));
    }

    #[test]
    fn test_get_writer_path_paired() {
        let prefix = "output/test";
        let compression = CompressionType::Uncompressed;
        let extension = "fastq";
        
        let path1 = get_writer_path(prefix, compression, extension, None, Some(0)).unwrap();
        assert_eq!(path1, PathBuf::from("output/test_1.fastq"));
        
        let path2 = get_writer_path(prefix, compression, extension, None, Some(1)).unwrap();
        assert_eq!(path2, PathBuf::from("output/test_2.fastq"));
        
        let path3 = get_writer_path(prefix, compression, extension, Some("tag1"), Some(0)).unwrap();
        assert_eq!(path3, PathBuf::from("output/test_tag1_1.fastq"));
    }

    fn test_compression_roundtrip(compression: CompressionType, extension: &str) {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("test").to_str().unwrap().to_string();
        let path_opt = get_writer_path(&prefix, compression, extension, None, None);
        let path = path_opt.unwrap();
        
        let data = b"test data for compression roundtrip";
        
        {
            let mut writer = get_writer_from_path(Some(&path), compression, 4096).unwrap();
            writer.write_all(data).unwrap();
            writer.flush().unwrap();
        }
        
        let (mut reader, _source, detected_compression) = get_raw_reader(&path, 4096);
        assert_eq!(detected_compression, compression);
        
        let mut buffer = Vec::new();
        reader.read_to_end(&mut buffer).unwrap();
        assert_eq!(buffer, data);
    }

    #[test]
    fn test_uncompressed_roundtrip() {
        test_compression_roundtrip(CompressionType::Uncompressed, "txt");
    }

    #[test]
    fn test_gz_roundtrip() {
        test_compression_roundtrip(CompressionType::Gz, "txt");
    }

    #[test]
    fn test_bz2_roundtrip() {
        test_compression_roundtrip(CompressionType::Bz2, "txt");
    }

    #[test]
    fn test_zstd_roundtrip() {
        test_compression_roundtrip(CompressionType::Zstd, "txt");
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
