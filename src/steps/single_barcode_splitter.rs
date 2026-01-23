use crate::config::Config;
use crate::types::{TrimEnd, RecordDto};
use log::trace;
use crate::statistics::WorkerStatistics;
use crate::steps::generic::{GenericStep, RecordInTransformation, AlignmentResults};
use crate::steps::utils::{extract_umi_tags, get_trim_range};
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;

pub struct SingleBarcodeSplitter {
    barcodes: Vec<(String, String)>,
    trim_end: TrimEnd,
    min_overlap: usize,
    error_rate: f32,
    gap_open: i32,
    gap_extend: i32,
    barcode_match: i32,
    barcode_mismatch: i32,
    barcode_keep: bool,
    barcode_unassigned: bool,
    removal_tags: bool,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    umi_tags: bool,
    max_barcode_len: usize,
}

impl SingleBarcodeSplitter {
    pub fn new(config: &Config, barcodes: &[(String, String)]) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix)
            = match config.args.barcode_trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (0, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, 0, 0, 0),
        };

        let max_barcode_len = barcodes.iter().map(|(_, b)| b.len()).max().unwrap_or(0);
        let barcode_min_overlap = config.args.barcode_min_overlap.unwrap_or(max_barcode_len);

        Self {
            barcodes: barcodes.to_vec(),
            trim_end: config.args.barcode_trim_end,
            min_overlap: barcode_min_overlap,
            error_rate: config.args.barcode_error_rate,
            gap_open: config.args.barcode_gap_open,
            gap_extend: config.args.barcode_gap_extend,
            barcode_match: config.args.barcode_match,
            barcode_mismatch: config.args.barcode_mismatch,
            barcode_keep: config.args.barcode_keep,
            barcode_unassigned: config.args.barcode_unassigned,
            removal_tags: config.args.removal_tags,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            umi_tags: config.args.umi_tags,
            max_barcode_len,
        }
    }

    fn get_scoring(&self) -> Scoring<impl Fn(u8, u8) -> i32> {
        let m = self.barcode_match;
        let mm = self.barcode_mismatch;
        Scoring {
            match_fn: move |a, b| if a == b || a == b'N' { m } else { mm },
            match_scores: None,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    pub fn split_barcodes(&self, record: &RecordDto, results: &mut AlignmentResults, stats: &mut WorkerStatistics) {
        let current_seq = &record.seq()[results.seq_start..results.seq_end];
        let mut aligner = Aligner::with_capacity_and_scoring(
            current_seq.len(),
            self.max_barcode_len,
            self.get_scoring(),
        );

        let mut best_alignment: Option<(Alignment, &str, &str)> = None;

        for (id, barcode) in &self.barcodes {
            let alignment = aligner.semiglobal(barcode.as_bytes(), current_seq);
            let overlap = alignment.yend - alignment.ystart;
            let barcode_overlap = alignment.xend - alignment.xstart;
            
            if overlap >= self.min_overlap && barcode_overlap >= self.min_overlap {
                let matches = (alignment.score - overlap as i32 * self.barcode_mismatch) / (self.barcode_match - self.barcode_mismatch);
                let errors = overlap as i32 - matches;
                let error_rate = errors as f32 / overlap as f32;

                if error_rate <= self.error_rate {
                    if best_alignment.as_ref().map_or(true, |(best, _, _)| alignment.score > best.score) {
                        best_alignment = Some((alignment, id, barcode));
                    }
                }
            }
        }

        if let Some((alignment, id, barcode)) = best_alignment {
            if self.umi_tags {
                results.umi_tags.extend(extract_umi_tags(&alignment, barcode.as_bytes(), current_seq));
            }

            if !self.barcode_keep {
                let (rel_start, rel_end) = get_trim_range(current_seq.len(), &alignment, self.trim_end);
                let old_start = results.seq_start;
                results.seq_start = old_start + rel_start;
                results.seq_end = old_start + rel_end;

                if self.removal_tags {
                    results.removal_tags.push("Flexbar_removal".to_string());
                }
            }
            
            *stats.barcode_counts.entry(id.to_string()).or_insert(0) += 1;
            results.barcode_tag = Some(format!("barcode_{}", id));
        }
    }
}

impl GenericStep<RecordDto> for SingleBarcodeSplitter {
    fn apply(&self, mut rec: RecordInTransformation<RecordDto>) -> RecordInTransformation<RecordDto> {
        if !rec.alignment_results.record_dropped {
            trace!("Barcode splitter processing record: {}", rec.record.id());
            self.split_barcodes(&rec.record, &mut rec.alignment_results, &mut rec.stats);
            
            if rec.alignment_results.barcode_tag.is_none() {
                if self.barcode_unassigned {
                    rec.alignment_results.barcode_tag = Some("barcode_unassigned_".to_string());
                } else {
                    rec.stats.total_discarded_barcode += 1;
                    rec.alignment_results.record_dropped = true;
                    return rec;
                }
            }

            if let Some(ref t) = rec.alignment_results.barcode_tag {
                trace!("Barcode splitter assigned tag: {} to record: {}", t, rec.record.id());
            }
        }
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use crate::steps::generic::AlignmentResults;
    use clap::Parser;

    #[test]
    fn test_barcode_unassigned_false() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4"]);
        let config = Config::new(args);
        let barcodes = vec![("bc1".to_string(), "GGGG".to_string())];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);
        
        let record = RecordDto::with_attrs("read1", None, b"AAAAACGT", b"IIIIIIII");
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: AlignmentResults {
                barcode_tag: None,
                umi_tags: Vec::new(),
                removal_tags: Vec::new(),
                seq_start: 0,
                seq_end: record.seq().len(),
                record_dropped: false,
            },
            stats: WorkerStatistics::default(),
        };
        let res = splitter.apply(rec_in);
        
        assert!(res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_barcode, 1);
    }

    #[test]
    fn test_barcode_unassigned_true() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4", "--barcode-unassigned"]);
        let config = Config::new(args);
        let barcodes = vec![("bc1".to_string(), "GGGG".to_string())];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);
        
        let record = RecordDto::with_attrs("read1", None, b"AAAAACGT", b"IIIIIIII");
        let rec_in = RecordInTransformation {
            record: record.clone(),
            alignment_results: AlignmentResults {
                barcode_tag: None,
                umi_tags: Vec::new(),
                removal_tags: Vec::new(),
                seq_start: 0,
                seq_end: record.seq().len(),
                record_dropped: false,
            },
            stats: WorkerStatistics::default(),
        };
        let res = splitter.apply(rec_in);
        
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_unassigned_".to_string()));
        assert_eq!(res.alignment_results.seq_end - res.alignment_results.seq_start, record.seq().len());
    }
}
