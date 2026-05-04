use bio::io::fastq::Record;
use crate::records::generic_read::GenericRawRead;


use std::ops::Range;

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

    fn finalize(
        &self,
        range: &Range<usize>,
        removal_tags: &Vec<String>,
        umi_tags: &Vec<String>,
        paired_read_dropped: bool,
    ) -> Self {
        let mut desc_parts: Vec<String> = Vec::new();
        if let Some(x) = self.desc() {
            desc_parts.push(x.to_string());
        }
        for tag in removal_tags {
            desc_parts.push(tag.clone());
        }
        for tag in umi_tags {
            desc_parts.push(tag.clone());
        }
        if paired_read_dropped {
            desc_parts.push("single".to_string());
        }
        let new_desc = if desc_parts.is_empty() {
            None
        } else {
            Some(desc_parts.join("_"))
        };

        Record::with_attrs(self.id(), new_desc.as_deref(), &self.seq()[range.clone()], &self.qual()[range.clone()])
    }
}

