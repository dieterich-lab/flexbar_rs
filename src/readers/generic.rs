use std::io::{self, BufRead};
use crate::readers::SequenceReader;
use crate::types::CompressionType;

pub trait GenericReaderBackend<R: BufRead + Send, Rec> {
    type Error: std::fmt::Display;
    fn new(reader: R) -> Self;
    fn read_record(&mut self) -> Result<Option<Rec>, Self::Error>;
    fn get_id(record: &Rec) -> &str;
}

pub struct GenericReader<R: BufRead + Send, Rec, B: GenericReaderBackend<R, Rec>> {
    backend: B,
    source: String,
    compression: CompressionType,
    last_id: String,
    _phantom: std::marker::PhantomData<(R, Rec)>,
}

impl<R: BufRead + Send, Rec, B: GenericReaderBackend<R, Rec>> GenericReader<R, Rec, B> {
    pub fn new(reader: R, source: String, compression: CompressionType) -> Self {
        Self {
            backend: B::new(reader),
            source,
            compression,
            last_id: String::new(),
            _phantom: std::marker::PhantomData,
        }
    }

    fn truncate_id(id: &str) -> String {
        if id.len() > 64 {
            id[..64].to_string() + "..."
        } else {
            id.to_string()
        }
    }
}

impl<R: BufRead + Send, Rec, B: GenericReaderBackend<R, Rec>> Iterator for GenericReader<R, Rec, B> {
    type Item = Rec;

    fn next(&mut self) -> Option<Self::Item> {
        match self.backend.read_record() {
            Ok(Some(record)) => {
                self.last_id = Self::truncate_id(B::get_id(&record));
                Some(record)
            }
            Ok(None) => None,
            Err(e) => {
                let context = if self.last_id.is_empty() {
                    " (no records processed yet)".to_string()
                } else {
                    format!(" (last successful record ID: {})", self.last_id)
                };
                panic!("Error reading from {}: {}{}", self.source, e, context);
            }
        }
    }
}

impl<R: BufRead + Send, B: GenericReaderBackend<R, bio::io::fastq::Record>> SequenceReader for GenericReader<R, bio::io::fastq::Record, B> {
    fn source(&self) -> &str {
        &self.source
    }

    fn compression(&self) -> CompressionType {
        self.compression
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    struct TestBackend {
        count: usize,
    }
    impl GenericReaderBackend<BufReader<&[u8]>, String> for TestBackend {
        type Error = io::Error;
        fn new(_reader: BufReader<&[u8]>) -> Self { TestBackend { count: 0 } }
        fn read_record(&mut self) -> Result<Option<String>, Self::Error> {
            self.count += 1;
            if self.count == 1 {
                Ok(Some("rec1".to_string()))
            } else {
                Err(io::Error::new(io::ErrorKind::Other, "test error"))
            }
        }
        fn get_id(record: &String) -> &str { record }
    }

    #[test]
    #[should_panic(expected = "Error reading from test_source: test error (last successful record ID: rec1)")]
    fn test_generic_reader_panic_with_id() {
        let data = b"dummy";
        let reader = BufReader::new(&data[..]);
        let mut gen_reader: GenericReader<_, String, TestBackend> = GenericReader::new(reader, "test_source".to_string(), CompressionType::Uncompressed);
        
        gen_reader.next(); // first call returns "rec1"
        gen_reader.next(); // second call panics
    }

    #[test]
    #[should_panic(expected = "Error reading from test_source: test error (last successful record ID: 0123456789012345678901234567890123456789012345678901234567890123...)")]
    fn test_generic_reader_panic_with_truncated_id() {
        struct LongIdBackend { count: usize }
        impl GenericReaderBackend<BufReader<&[u8]>, String> for LongIdBackend {
            type Error = io::Error;
            fn new(_reader: BufReader<&[u8]>) -> Self { LongIdBackend { count: 0 } }
            fn read_record(&mut self) -> Result<Option<String>, Self::Error> {
                self.count += 1;
                if self.count == 1 {
                    Ok(Some("0123456789".repeat(10))) // 100 chars
                } else {
                    Err(io::Error::new(io::ErrorKind::Other, "test error"))
                }
            }
            fn get_id(record: &String) -> &str { record }
        }

        let data = b"dummy";
        let reader = BufReader::new(&data[..]);
        let mut gen_reader: GenericReader<_, String, LongIdBackend> = GenericReader::new(reader, "test_source".to_string(), CompressionType::Uncompressed);
        
        gen_reader.next();
        gen_reader.next();
    }
}
