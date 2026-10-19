use crate::records::unified::AnnotatedRecord;
use crate::records::unified_read::RawRead;
use crate::records::generic_read::GenericRawRead;

pub fn get_desc(record: &AnnotatedRecord, read_index: usize, barcode_mark: bool) -> Option<String> {
    let read = record.record.reads[read_index].as_ref().expect("Read missing");

    let mut parts: Vec<String> = Vec::new();
    
    match &read.raw {
        RawRead::Fasta(r) => match r.desc() {
            Some(d) => parts.push(d.to_string()),
            _ => {}
        },
        RawRead::Fastq(r) => match r.desc() {
            Some(d) => parts.push(d.to_string()),
            _ => {}
        },
        _ => {}
    }
    for t in &record.removal_tags[read_index] {
        parts.push(t.clone());
    }
    for t in &record.umi_tags {
        parts.push(t.clone());
    }
    if record.paired_read_dropped {
        parts.push("single".to_string());
    }
    if barcode_mark && let Some(id) = &record.barcode_id {
        parts.push(format!("barcode_{}", id));
    }

    if parts.is_empty() {
        None
    } else {
        Some(parts.join("_"))
    }
}

pub fn get_qual(record: &AnnotatedRecord, read_index: usize) -> Vec<u8> {
    let read = record.record.reads[read_index].as_ref().expect("Read missing");
    match read.qual() {
        Some(q) => q.to_vec(),
        None => vec![b'I'; read.seq().len()],
    }
}
