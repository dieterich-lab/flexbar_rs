use bio::io::fasta;
use std::path::Path;
use std::io;

use std::fmt::Debug;

pub struct Bar {
    pub id: String,
    pub seq: String,
}

pub fn load_sequence<P: AsRef<Path> + Debug>(path: P) -> Vec<Bar> {
    let path_ref = path.as_ref();
    let reader = fasta::Reader::from_file(path_ref).unwrap_or_else(|e| {
        panic!("Failed to open file {:?}: {}", path_ref, e);
    });
    let bars = load_sequence_from_reader(reader, path_ref);
    check_unique_ids(&bars, path_ref);
    bars
}

fn check_unique_ids(bars: &[Bar], path: &Path) {
    let mut ids = std::collections::HashSet::new();
    for bar in bars {
        if !ids.insert(&bar.id) {
            panic!(
                "Duplicate ID '{}' found in file: {:?}",
                bar.id, path
            );
        }
    }
}

fn load_sequence_from_reader<R: io::BufRead>(reader: fasta::Reader<R>, path: &Path) -> Vec<Bar> {
    let mut bars = Vec::new();

    for result in reader.records() {
        match result {
            Ok(record) => {
                let id = record.id().to_string();
                let seq = match String::from_utf8(record.seq().to_vec()) {
                    Ok(s) => s,
                    Err(e) => {
                        let context = format_last_record(bars.last());
                        panic!(
                            "Failed to parse sequence as UTF-8 in file {:?}: {}{}",
                            path, e, context
                        );
                    }
                };
                bars.push(Bar { id, seq });
            }
            Err(e) => {
                let context = format_last_record(bars.last());
                panic!(
                    "Error reading FASTA record in file {:?}: {}{}",
                    path, e, context
                );
            }
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

        let bars = load_sequence_from_reader(reader, Path::new("test.fasta"));
        assert_eq!(bars.len(), 2);
        assert_eq!(bars[0].id, "seq1");
        assert_eq!(bars[0].seq, "ACGT");
        assert_eq!(bars[1].id, "seq2");
        assert_eq!(bars[1].seq, "TGCA");
    }

    #[test]
    #[should_panic(expected = "Error reading FASTA record in file \"test.fasta\": stream did not contain valid UTF-8\nNo records were successfully processed before the error.")]
    fn test_load_sequence_from_reader_panic() {
        // We trigger a read error by providing invalid UTF-8 bytes that the FASTA reader tries to interpret.
        let invalid_utf8_data = b">seq1\n\xFF\xFF\n";
        let cursor = Cursor::new(invalid_utf8_data);
        let reader = fasta::Reader::new(cursor);

        load_sequence_from_reader(reader, Path::new("test.fasta"));
    }

    #[test]
    #[should_panic(expected = "Duplicate ID 'seq1' found in file: \"test.fasta\"")]
    fn test_check_unique_ids_panic() {
        let bars = vec![
            Bar { id: "seq1".to_string(), seq: "ACGT".to_string() },
            Bar { id: "seq1".to_string(), seq: "TGCA".to_string() },
        ];
        check_unique_ids(&bars, Path::new("test.fasta"));
    }

    #[test]
    fn test_truncate_string() {
        assert_eq!(truncate_string("hello", 10), "hello");
        assert_eq!(truncate_string("hello world", 5), "hello...");
        assert_eq!(truncate_string("1234567890", 10), "1234567890");
    }
}
