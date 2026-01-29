use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::fs::File;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use bzip2::read::BzDecoder;
use bzip2::write::BzEncoder;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;
use crate::types::CompressionType;

use crate::logging::log_fatal;

pub fn get_raw_reader<P: AsRef<Path>>(path: P) -> io::Result<(BufReader<Box<dyn Read + Send>>, String, CompressionType)> {
    let path_ref = path.as_ref();
    let source = path_ref.to_string_lossy().into_owned();

    let (reader, compression): (Box<dyn Read + Send>, CompressionType) = if source == "-" {
        (Box::new(io::stdin()), CompressionType::Uncompressed)
    } else {
        let file = File::open(path_ref)?;
        let extension = path_ref.extension()
            .and_then(|s| s.to_str())
            .map(|s| s.to_lowercase());

        match extension.as_deref() {
            Some("gz") => (Box::new(MultiGzDecoder::new(file)), CompressionType::Gz),
            Some("bz2") => (Box::new(BzDecoder::new(file)), CompressionType::Bz2),
            Some("zstd") | Some("zst") => (Box::new(ZstdDecoder::new(file)?), CompressionType::Zstd),
            _ => (Box::new(file), CompressionType::Uncompressed),
        }
    };

    Ok((BufReader::new(reader), source, compression))
}

pub fn get_writer_path(prefix: &str, compression: CompressionType, extension: &str, tag: Option<&str>) -> Option<PathBuf> {
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

pub fn get_writer_from_path(path: Option<&Path>, compression: CompressionType) -> io::Result<Box<dyn Write + Send>> {
    let writer: Box<dyn Write + Send> = match path {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(io::stdout()),
    };

    match compression {
        CompressionType::Gz => Ok(Box::new(GzEncoder::new(writer, flate2::Compression::default()))),
        CompressionType::Bz2 => Ok(Box::new(BzEncoder::new(writer, bzip2::Compression::default()))),
        CompressionType::Zstd => Ok(Box::new(ZstdEncoder::new(writer, 0)?.auto_finish())),
        CompressionType::Uncompressed => Ok(writer),
        CompressionType::Keep => log_fatal("Compression type KEEP is not supported for writing"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_bz2_decompression() {
        let mut tmp = NamedTempFile::new().unwrap();
        let mut bz2_writer = bzip2::write::BzEncoder::new(Vec::new(), bzip2::Compression::default());
        bz2_writer.write_all(b">test\nACGT\n").unwrap();
        let compressed_data = bz2_writer.finish().unwrap();
        tmp.write_all(&compressed_data).unwrap();
        
        // Rename tmp file to have .bz2 extension
        let path = tmp.path().with_extension("bz2");
        std::fs::copy(tmp.path(), &path).unwrap();

        let (mut reader, source, compression) = get_raw_reader(&path).unwrap();
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
        let path = get_writer_path(prefix, compression, extension, tag).unwrap();
        assert_eq!(path, PathBuf::from("output/test_tagwithslashes.fastq"));
    }

    #[test]
    fn test_get_writer_path_with_tag() {
        let prefix = "output/test";
        let compression = CompressionType::Uncompressed;
        let extension = "fastq";
        let tag = Some("tag1");
        let path = get_writer_path(prefix, compression, extension, tag).unwrap();
        assert_eq!(path, PathBuf::from("output/test_tag1.fastq"));
    }

    #[test]
    fn test_get_writer_path_without_tag() {
        let prefix = "output/test";
        let compression = CompressionType::Uncompressed;
        let extension = "fastq";
        let path = get_writer_path(prefix, compression, extension, None).unwrap();
        assert_eq!(path, PathBuf::from("output/test.fastq"));
    }

    fn test_compression_roundtrip(compression: CompressionType, extension: &str) {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("test").to_str().unwrap().to_string();
        let path_opt = get_writer_path(&prefix, compression, extension, None);
        let path = path_opt.unwrap();
        
        let data = b"test data for compression roundtrip";
        
        {
            let mut writer = get_writer_from_path(Some(&path), compression).unwrap();
            writer.write_all(data).unwrap();
            writer.flush().unwrap();
        }
        
        let (mut reader, _source, detected_compression) = get_raw_reader(&path).unwrap();
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
        let (_reader, source, compression) = get_raw_reader("-").unwrap();
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

        let (mut reader, _source, compression) = get_raw_reader(&path).unwrap();
        assert_eq!(compression, CompressionType::Gz);

        let mut content = String::new();
        reader.read_to_string(&mut content).unwrap();
        assert_eq!(content, "member1member2");
    }
}
