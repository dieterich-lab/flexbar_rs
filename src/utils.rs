use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::fs::File;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use bzip2::read::BzDecoder;
use bzip2::write::BzEncoder;
use zstd::stream::read::Decoder as ZstdDecoder;
use zstd::stream::write::Encoder as ZstdEncoder;
use crate::types::CompressionType;

pub fn get_reader<P: AsRef<Path>>(path: P) -> io::Result<(BufReader<Box<dyn Read + Send>>, String, CompressionType)> {
    let path_ref = path.as_ref();
    let source = path_ref.to_string_lossy().into_owned();
    let file = File::open(path_ref)?;
    
    let extension = path_ref.extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_lowercase());

    let (reader, compression): (Box<dyn Read + Send>, CompressionType) = match extension.as_deref() {
        Some("gz") => (Box::new(GzDecoder::new(file)), CompressionType::Gz),
        Some("bz2") => (Box::new(BzDecoder::new(file)), CompressionType::Bz2),
        Some("zstd") | Some("zst") => (Box::new(ZstdDecoder::new(file)?), CompressionType::Zstd),
        _ => (Box::new(file), CompressionType::Uncompressed),
    };

    Ok((BufReader::new(reader), source, compression))
}

pub fn get_writer_path(prefix: &str, compression: CompressionType, extension: &str) -> PathBuf {
    let mut path = PathBuf::from(prefix);
    path.set_extension(extension);
    
    match compression {
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
        CompressionType::Keep => panic!("Compression type KEEP is not supported for writing"),
    }
}

pub fn get_writer(path: &Path, compression: CompressionType) -> io::Result<Box<dyn Write + Send>> {
    let file = File::create(path)?;
    match compression {
        CompressionType::Gz => Ok(Box::new(GzEncoder::new(file, flate2::Compression::default()))),
        CompressionType::Bz2 => Ok(Box::new(BzEncoder::new(file, bzip2::Compression::default()))),
        CompressionType::Zstd => Ok(Box::new(ZstdEncoder::new(file, 0)?.auto_finish())),
        CompressionType::Uncompressed => Ok(Box::new(file)),
        CompressionType::Keep => panic!("Compression type KEEP is not supported for writing"),
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

        let (mut reader, source, compression) = get_reader(&path).unwrap();
        assert_eq!(compression, CompressionType::Bz2);
        assert!(source.contains(path.to_str().unwrap()));

        let mut content = String::new();
        reader.read_to_string(&mut content).unwrap();
        assert_eq!(content, ">test\nACGT\n");

        std::fs::remove_file(path).unwrap();
    }

    fn test_compression_roundtrip(compression: CompressionType, extension: &str) {
        let dir = tempdir().unwrap();
        let prefix = dir.path().join("test").to_str().unwrap().to_string();
        let path = get_writer_path(&prefix, compression, extension);
        
        let data = b"test data for compression roundtrip";
        
        {
            let mut writer = get_writer(&path, compression).unwrap();
            writer.write_all(data).unwrap();
            writer.flush().unwrap();
        }
        
        let (mut reader, _source, detected_compression) = get_reader(&path).unwrap();
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
    #[should_panic(expected = "Compression type KEEP is not supported for writing")]
    fn test_get_writer_path_panic_on_keep() {
        get_writer_path("test", CompressionType::Keep, "txt");
    }

    #[test]
    #[should_panic(expected = "Compression type KEEP is not supported for writing")]
    fn test_get_writer_panic_on_keep() {
        get_writer(Path::new("test"), CompressionType::Keep).unwrap();
    }
}
