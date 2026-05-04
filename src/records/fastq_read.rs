use bio::io::fastq::Record;
use crate::records::generic_read::GenericRawRead;


impl GenericRawRead for Record {
    fn seq(&self) -> &[u8] {
        self.seq()
    }

    fn qual(&self) -> Option<&[u8]> {
        Some(self.qual())
    }

    fn id(&self) -> &str {
        self.id()
    }

    fn desc(&self) -> Option<&str> {
        self.desc()
    }

    fn desc_str(&self) -> Option<&str> {
        self.desc()
    }
}

