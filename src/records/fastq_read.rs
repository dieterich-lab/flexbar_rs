use std::ops::Range;
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

    fn desc(&self) -> &[u8] {
        self.desc().unwrap_or("").as_bytes()
    }

    fn transform(self, range: Range<usize>, extra_desc_parts: Vec<String>) -> Self {
        let id = self.id().to_string();
        let mut desc_parts: Vec<String> = match self.desc() {
            Some(desc) => vec![desc.to_string()],
            _ => vec![],
        };
        desc_parts.extend(extra_desc_parts);
        let joined_desc = desc_parts.join("_");
        let desc = match joined_desc.len() > 0 {
            true => Some(joined_desc.as_ref()),
            false => None,
        };
        let seq = self.seq()[range.clone()].to_vec();
        let qual = self.qual()[range].to_vec();
        Record::with_attrs(&id, desc, &seq, &qual)
    }
}
