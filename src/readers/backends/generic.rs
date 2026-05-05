use std::io::BufRead;
use crate::records::unified::UnifiedRead;
use crate::records::header::SequenceHeader;

pub trait GenericReaderBackend<Reader: BufRead + Send> {
    fn read_record(&mut self) -> Option<UnifiedRead>;
    fn get_header(&self) -> SequenceHeader;
}