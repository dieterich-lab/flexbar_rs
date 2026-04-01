use bio::io::fastq::Record;
use crate::steps::generic::PerReadAlignmentResults;

pub fn count_uncalled_bases(record: &Record, results: &PerReadAlignmentResults) -> usize {
    if results.dropped {
        return 0;
    }
    
    let seq = record.seq();
    let range = results.read_range.clone();
    
    seq[range].iter()
        .filter(|&&b| b == b'N' || b == b'n')
        .count()
}
