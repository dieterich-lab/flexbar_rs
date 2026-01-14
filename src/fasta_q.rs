use bio::io::fasta;
use bio::alphabets::dna::revcomp;
use std::path::Path;
use std::io;

use std::fmt::Debug;
use crate::types::RevCompMode;

pub struct Bar {
    pub id: String,
    pub seq: String,
}

pub struct Bars(pub Vec<Bar>);

pub fn load_barcodes<P: AsRef<Path> + Debug>(path: P) -> Bars {
    let path_ref = path.as_ref();
    let source = path.as_ref().to_str().unwrap_or("<unknown>") ;
    let reader = fasta::Reader::from_file(path_ref).unwrap_or_else(|e| {
        panic!("Failed to open file {:?}: {}", path_ref, e);
    });
    let bars = load_sequence_from_reader(reader, source.as_ref(), true, false);
    check_unique_ids(&bars, path_ref);
    bars
}

pub fn load_adapters<P: AsRef<Path> + Debug>(path: P, rev_comp_mode: RevCompMode) -> Bars {
    let path_ref = path.as_ref();
    let source = path.as_ref().to_str().unwrap_or("<unknown>") ;
    let reader = fasta::Reader::from_file(path_ref).unwrap_or_else(|e| {
        panic!("Failed to open file {:?}: {}", path_ref, e);
    });

    let record_normal = rev_comp_mode != RevCompMode::Only;
    let record_reverse_complement = rev_comp_mode != RevCompMode::Off;

    let bars = load_sequence_from_reader(reader, source.as_ref(), record_normal, record_reverse_complement);
    check_unique_ids(&bars, path_ref);
    bars
}

fn check_unique_ids(bars: &Bars, path: &Path) {
    let mut ids = std::collections::HashSet::new();
    for bar in &bars.0 {
        if !ids.insert(&bar.id) {
            panic!(
                "Duplicate ID '{}' found in file: {:?}",
                bar.id, path
            );
        }
    }
}

fn load_sequence_from_reader<R: io::BufRead>(
    reader: fasta::Reader<R>,
    source: &str,
    record_normal: bool,
    record_reverse_complement: bool,
) -> Bars  {
    let mut bars = Vec::new();

    for result in reader.records() {
        match result {
            Ok(record) => {
                let raw_id = record.id().to_string();
                let raw_seq_bytes = record.seq();
                if record_normal {
                    let id = raw_id.clone();
                    let seq_bytes = raw_seq_bytes;
                    bars = process_sequence(bars, id, seq_bytes, source);
                }
                if record_reverse_complement {
                    let id_rc = format!("{}_rc", raw_id);
                    let rc_seq_bytes = revcomp(raw_seq_bytes);
                    bars = process_sequence(bars, id_rc, &rc_seq_bytes, source);
                }
            }
            Err(e) => {
                let context = format_last_record(bars.last());
                panic!(
                    "Error reading record in file {:?}: {}{}",
                    source, e, context
                );
            }
        }
    }

    Bars(bars)
}

fn process_sequence(mut bars: Vec<Bar>, id: String, seq_bytes: &[u8], source: &str) -> Vec<Bar> {
    match String::from_utf8(seq_bytes.to_vec()) {
        Ok(seq) => {
            bars.push(Bar { id: id.clone(), seq })
        },
        Err(e) => {
            let context = format_last_record(bars.last());
            panic!(
                "Failed to process ID {} as UTF-8 in file {:?}: {}{}",
                id.clone(), source, e, context
            );
        }
    }
    bars
}

fn format_last_record(last_record: Option<&Bar>) -> String {
    if let Some(bar) = last_record {
        let truncated_id = truncate_string(&bar.id, 20);
        let truncated_seq = truncate_string(&bar.seq, 50);
        format!(
            "\nLast successfully processed record: ID='{}', Seq='{}'",
            truncated_id, truncated_seq
        )
    } else {
        "\nNo records were successfully processed before the error.".to_string()
    }
}

fn truncate_string(s: &str, limit: usize) -> String {
    if s.len() > limit {
        format!("{}...", &s[..limit])
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_load_sequence_from_reader() {
        let data = ">seq1\nACGT\n>seq2\nTGCA";
        let cursor = Cursor::new(data);
        let reader = fasta::Reader::new(cursor);

        let bars = load_sequence_from_reader(reader, "test.fasta", true, false);
        assert_eq!(bars.0.len(), 2);
        assert_eq!(bars.0[0].id, "seq1");
        assert_eq!(bars.0[0].seq, "ACGT");
        assert_eq!(bars.0[1].id, "seq2");
        assert_eq!(bars.0[1].seq, "TGCA");
    }

    #[test]
    fn test_load_sequence_from_reader_revcomp() {
        let data = ">seq1\nAAAAAAAA";
        let cursor = Cursor::new(data);
        let reader = fasta::Reader::new(cursor);

        // Test normal and revcomp (corresponds to RevCompMode::On)
        let bars = load_sequence_from_reader(reader,"test.fasta", true, true);
        assert_eq!(bars.0.len(), 2);
        assert_eq!(bars.0[0].id, "seq1");
        assert_eq!(bars.0[0].seq, "AAAAAAAA");
        assert_eq!(bars.0[1].id, "seq1_rc");
        assert_eq!(bars.0[1].seq, "TTTTTTTT");

        // Test only revcomp (corresponds to RevCompMode::Only)
        let data2 = ">seq1\nGCTAGCTA"; // Revcomp: TAGCTAGC
        let cursor2 = Cursor::new(data2);
        let reader2 = fasta::Reader::new(cursor2);
        let bars2 = load_sequence_from_reader(reader2, "test.fasta", false, true);
        assert_eq!(bars2.0.len(), 1);
        assert_eq!(bars2.0[0].id, "seq1_rc");
        assert_eq!(bars2.0[0].seq, "TAGCTAGC");

        // Test normal only (corresponds to RevCompMode::Off)
        let data3 = ">seq1\nATGCATGC"; // Revcomp: GCATGCAT
        let cursor3 = Cursor::new(data3);
        let reader3 = fasta::Reader::new(cursor3);
        let bars3 = load_sequence_from_reader(reader3,"test.fasta", true, false);
        assert_eq!(bars3.0.len(), 1);
        assert_eq!(bars3.0[0].id, "seq1");
        assert_eq!(bars3.0[0].seq, "ATGCATGC");
    }

    #[test]
    #[should_panic(expected = "Error reading record in file \"test.fasta\": stream did not contain valid UTF-8\nNo records were successfully processed before the error.")]
    fn test_load_sequence_from_reader_panic() {
        // We trigger a read error by providing invalid UTF-8 bytes that the FASTA reader tries to interpret.
        let invalid_utf8_data = b">seq1\n\xFF\xFF\n";
        let cursor = Cursor::new(invalid_utf8_data);
        let reader = fasta::Reader::new(cursor);

        load_sequence_from_reader(reader,"test.fasta", true, false);
    }

    #[test]
    #[should_panic(expected = "Duplicate ID 'seq1' found in file: \"test.fasta\"")]
    fn test_check_unique_ids_panic() {
        let bars = Bars(vec![
            Bar { id: "seq1".to_string(), seq: "ACGT".to_string() },
            Bar { id: "seq1".to_string(), seq: "TGCA".to_string() },
        ]);
        check_unique_ids(&bars, Path::new("test.fasta"));
    }

    #[test]
    fn test_truncate_string() {
        assert_eq!(truncate_string("hello", 10), "hello");
        assert_eq!(truncate_string("hello world", 5), "hello...");
        assert_eq!(truncate_string("1234567890", 10), "1234567890");
    }
}
