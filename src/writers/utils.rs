use crate::records::unified::{UnifiedRead, AnnotatedRecord};
use crate::records::generic_read::GenericRawRead;

pub fn get_desc(record: &AnnotatedRecord, read_index: usize) -> Option<String> {
    let read = record.record.reads[read_index].as_ref().expect("Read missing");

    let mut parts: Vec<String> = Vec::new();
    
    match read {
        UnifiedRead::Fasta(r) => match r.desc() {
            Some(d) => parts.push(d.to_string()),
            _ => {}
        },
        UnifiedRead::Fastq(r) => match r.desc() {
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

    if parts.is_empty() {
        None
    } else {
        Some(parts.join("_"))
    }
}

pub fn get_qual(record: &AnnotatedRecord, read_index: usize) -> Vec<u8> {
    let read = record.record.reads[read_index].as_ref().expect("Read missing");
    match read {
        UnifiedRead::Fastq(r) => r.qual().to_vec(),
        UnifiedRead::Fasta(r) => vec![b'I'; r.seq().len()],
        UnifiedRead::Bam(r) => match r.qual() {
            Some(q) => q.to_vec(),
            None => vec![b'I'; r.seq().len()],
        },
    }
}
