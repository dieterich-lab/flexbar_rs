use crate::config::Config;
use crate::types::{TrimEnd, SingeRecordDto};
use log::trace;
use crate::statistics::WorkerStatistics;
use crate::steps::generic::{GenericStep, RecordInTransformation, AlignmentResults};
use crate::steps::utils::{extract_umi_tags, get_trim_range};
use bio::alignment::pairwise::*;
use bio::alignment::{Alignment, AlignmentOperation};
use crate::context::SequenceWithId;

pub struct SingleBarcodeSplitter {
    barcodes: Vec<SequenceWithId>,
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
    trace: bool,
}

struct BarCodeAlignment<'a> {
    alignment: Alignment,
    barcode: &'a SequenceWithId,
}

impl SingleBarcodeSplitter {
    pub fn new(config: &Config, barcodes: &Vec<SequenceWithId>) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix)
            = match config.args.barcode_trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (MIN_SCORE, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, MIN_SCORE, 0, 0),
        };

        let max_barcode_len = barcodes.iter().map(|b| b.sequence.len())
            .max().unwrap_or(0);
        let barcode_min_overlap = config.args.barcode_min_overlap.unwrap_or(max_barcode_len);

        Self {
            barcodes: (*barcodes).clone(),
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
            trace: config.args.trace
        }
    }

    fn get_scoring(&self) -> Scoring<impl Fn(u8, u8) -> i32> {
        let m = self.barcode_match;
        let mm = self.barcode_mismatch;
        Scoring {
            match_fn: move |a, b|
                if a == b'N' || b == b'N' { 0 }
                else if a == b { m }
                else { mm },
            match_scores: None,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    pub fn split_barcodes(&self, record: &SingeRecordDto, results: &mut AlignmentResults, _stats: &mut WorkerStatistics) {
        let current_seq = &record.seq()[results.seq_start..results.seq_end];
        let mut aligner = Aligner::with_capacity_and_scoring(
            current_seq.len(),
            self.max_barcode_len,
            self.get_scoring(),
        );
        let mut best_alignment: Option<BarCodeAlignment> = None;
        for barcode in &self.barcodes {
            let alignment = aligner.semiglobal(barcode.sequence.as_bytes(), current_seq);
            let candidate = BarCodeAlignment { alignment, barcode: &barcode };
            self.update_best_alignment(&mut best_alignment, candidate, current_seq);
        }

        match best_alignment {
            Some(best) => {
                results.barcode_tag = Some(format!("barcode_{}", best.barcode.id));
                if self.umi_tags {
                    let new_umi_tags = extract_umi_tags(
                        &best.alignment, best.barcode.sequence.as_bytes(), current_seq
                    );
                    trace!("Barcode splitter on {}: Extracted UMI tags {:?}", record.id(), new_umi_tags);
                    results.umi_tags.extend(new_umi_tags);
                }
                if !self.barcode_keep {
                    let (rel_start, rel_end) = get_trim_range(
                        current_seq.len(), &best.alignment, self.trim_end
                    );
                    let old_start = results.seq_start;
                    results.seq_start = old_start + rel_start;
                    results.seq_end = old_start + rel_end;
                    if self.removal_tags {
                        results.removal_tags.push("Flexbar_removal".to_string());
                    }
                }
            },
            None => {}
        }
    }

    fn update_best_alignment<'a>(
        &self,
        best_alignment: &mut Option<BarCodeAlignment<'a>>,
        candidate: BarCodeAlignment<'a>,
        current_seq: &[u8],
    ) {
        let overlap = candidate.alignment.yend - candidate.alignment.ystart;
        let barcode_overlap = candidate.alignment.xend - candidate.alignment.xstart;

        if overlap < self.min_overlap || barcode_overlap < self.min_overlap { return; }
        let matches = self.count_matches(
            &candidate.alignment,
            candidate.barcode.sequence.as_ref(),
            current_seq,
        );

        let errors = overlap - matches;
        let error_rate = errors as f32 / overlap as f32;
        if error_rate > self.error_rate { return; }

        match best_alignment {
            Some(best) => {
                if candidate.alignment.score > best.alignment.score {
                    *best_alignment = Some(candidate);
                }
            },
            None => { *best_alignment = Some(candidate); }
        }
    }

    fn count_matches(&self, alignment: &Alignment, seq_a: &[u8], seq_b: &[u8]) -> usize {
        let mut matches = 0;
        let mut xpos = alignment.xstart;
        let mut ypos = alignment.ystart;

        for op in &alignment.operations {
            match op {
                AlignmentOperation::Match | AlignmentOperation::Subst => {
                    let a = seq_a[xpos];
                    let b = seq_b[ypos];
                    if a == b || a == b'N' || b == b'N' { matches += 1; }
                    xpos += 1;
                    ypos += 1;
                }
                AlignmentOperation::Ins => {
                    xpos += 1;
                }
                AlignmentOperation::Del => {
                    ypos += 1;
                }
                _ => {}
            }
        }
        matches
    }
}

impl GenericStep<SingeRecordDto> for SingleBarcodeSplitter {
    fn apply(&self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        self.split_barcodes(&rec.record, &mut rec.alignment_results, &mut rec.stats);

        match &rec.alignment_results.barcode_tag {
            Some(tag) => {
                if self.trace {
                    if self.barcode_keep {
                        trace!("Barcode splitter on {}: Assigned tag {}", rec.record.id(), tag);
                    } else {
                        trace!(
                            "Barcode splitter on {}: Assigned tag: {} and shortened sequence to:\n{}",
                            rec.record.id(), tag, rec.get_sequence_as_str()
                        );
                    }
                }
            }
            None => {
                if self.barcode_unassigned {
                    rec.alignment_results.barcode_tag = Some("barcode_unassigned_".to_string());
                    trace!("Barcode splitter on {}: no tag assigned.", rec.record.id());
                } else {
                    rec.stats.total_discarded_barcode += 1;
                    rec.alignment_results.record_dropped = true;
                    trace!("Barcode splitter on {}: no tag assigned and DROPPED.", rec.record.id());
                }
            },
        }
        rec
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;
    use crate::steps::generic::Transformable;

    #[test]
    fn test_barcode_unassigned_false() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4"]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: "GGGG".to_string()}];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"AAAAACGT", b"IIIIIIII");
        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert!(res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_barcode, 1);
    }

    #[test]
    fn test_barcode_unassigned_true() {
        let args = Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-min-overlap", "4", "--barcode-unassigned"]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId { id: "bc1".to_string(), sequence: "GGGG".to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);
        
        let record = SingeRecordDto::with_attrs("read1", None, b"AAAAACGT", b"IIIIIIII");
        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_unassigned_".to_string()));
        assert_eq!(res.alignment_results.seq_end - res.alignment_results.seq_start, record.seq().len());
    }

    #[test]
    fn test_barcode_trimming_right() {
        let seq = b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA";
        let qual = b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII";
        let record = SingeRecordDto::with_attrs("read1", None, seq, qual);

        // Barcode is TAGCTCTGTAGGCACCATCAA (21 bp)
        // It's at the end of the sequence.
        let barcode_seq = "TAGCTCTGTAGGCACCATCAA";
        
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--barcode-trim-end", "RIGHT",
            "--barcode-min-overlap", "10"
        ]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: barcode_seq.to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);

        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // Current sequence length is 71.
        // Barcode starts at 50.
        // Trimming RIGHT should keep 0..50.
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 50);
        assert_eq!(res.get_sequence_as_str(), "CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCC");
    }

    #[test]
    fn test_barcode_trimming_left() {
        let seq = b"TAGCTCTGTAGGCACCATCAACATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCC";
        let qual = b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII";
        let record = SingeRecordDto::with_attrs("read1", None, seq, qual);

        let barcode_seq = "TAGCTCTGTAGGCACCATCAA";
        
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--barcode-trim-end", "LEFT",
            "--barcode-min-overlap", "10"
        ]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: barcode_seq.to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);

        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // Barcode is at the start (0..21).
        // Trimming LEFT should keep 21..71.
        assert_eq!(res.alignment_results.seq_start, 21);
        assert_eq!(res.alignment_results.seq_end, 71);
        assert_eq!(res.get_sequence_as_str(), "CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCC");
    }

    #[test]
    fn test_barcode_trimming_any_right() {
        let seq = b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA";
        let qual = b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII";
        let record = SingeRecordDto::with_attrs("read1", None, seq, qual);

        let barcode_seq = "TAGCTCTGTAGGCACCATCAA";
        
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--barcode-trim-end", "ANY",
            "--barcode-min-overlap", "10"
        ]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: barcode_seq.to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);

        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // Barcode is at the end. Trimming ANY should pick the side that removes less.
        // Wait, current ANY logic:
        // if alignment.ystart < len - alignment.yend { (alignment.yend, len) } else { (0, alignment.ystart) }
        // ystart = 50. len = 71. yend = 71.
        // 50 < 71 - 71 => 50 < 0 (false) => (0, 50)
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 50);
    }

    #[test]
    fn test_barcode_trimming_any_left() {
        let seq = b"TAGCTCTGTAGGCACCATCAACATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCC";
        let qual = b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII";
        let record = SingeRecordDto::with_attrs("read1", None, seq, qual);

        let barcode_seq = "TAGCTCTGTAGGCACCATCAA";
        
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--barcode-trim-end", "ANY",
            "--barcode-min-overlap", "10"
        ]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: barcode_seq.to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);

        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // Barcode is at the start.
        // ystart = 0. len = 71. yend = 21.
        // 0 < 71 - 21 => 0 < 50 (true) => (21, 71)
        assert_eq!(res.alignment_results.seq_start, 21);
        assert_eq!(res.alignment_results.seq_end, 71);
    }

    #[test]
    fn test_barcode_trimming_ltail() {
        let seq = b"TAGCTCTGTAGGCACCATCAACATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCC";
        let qual = b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII";
        let record = SingeRecordDto::with_attrs("read1", None, seq, qual);

        let barcode_seq = "TAGCTCTGTAGGCACCATCAA";
        
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--barcode-trim-end", "LTAIL",
            "--barcode-min-overlap", "10"
        ]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: barcode_seq.to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);

        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // Barcode is at the start (0..21).
        // LTAIL should keep 0..71 (it shouldn't really trim anything if it's perfectly at the start,
        // but it removes the matched part of the barcode from the read).
        // Actually LTAIL with (alignment.ystart, len) will keep (0, 71).
        // Let's see what happens if barcode is slightly shifted.
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 71);
    }

    #[test]
    fn test_barcode_trimming_rtail() {
        let seq = b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA";
        let qual = b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII";
        let record = SingeRecordDto::with_attrs("read1", None, seq, qual);

        let barcode_seq = "TAGCTCTGTAGGCACCATCAA";
        
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--barcode-trim-end", "RTAIL",
            "--barcode-min-overlap", "10"
        ]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: barcode_seq.to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);

        let rec_in = record.get_record_in_transformation();
        let res = splitter.apply(rec_in);
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        // Barcode is at the end (50..71).
        // RTAIL with (0, alignment.yend) will keep (0, 71).
        assert_eq!(res.alignment_results.seq_start, 0);
        assert_eq!(res.alignment_results.seq_end, 71);
    }

    #[test]
    fn test_issue_record() {
        // Record from issue:
        // @VH00748:119:AAGJ22MM5:1:1101:28608:1038 1:N:0:TACGGTCT+GCAATGGA
        // CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA
        // +
        // IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII

        let seq = b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA";
        let qual = b"IIIIIIIIIIIII-II9IIIIIIIIII9I9III9IIIII9IIIII9IIII99IIIIIIII9IIIII9IIII";
        let record = SingeRecordDto::with_attrs("VH00748:119:AAGJ22MM5:1:1101:28608:1038 1:N:0:TACGGTCT+GCAATGGA", None, seq, qual);

        // Barcode from issue:
        let barcode_seq = "NNNNNNTAGCTCTGTAGGCACCATCAAT";
        
        // Flexbar default is LTAIL for barcodes. 
        // LTAIL means it should match at the beginning of the read.
        let args = Args::parse_from(&[
            "flexbar_rs", 
            "--reads", "test.fq", 
            "--barcode-trim-end", "RIGHT",
            "--barcode-min-overlap", "11",
            "--post-trim-length", "248",
            "--umi-tags"
        ]);
        let config = Config::new(args);
        let barcodes = vec![
            SequenceWithId{ id: "bc1".to_string(), sequence: barcode_seq.to_string()}
        ];
        let splitter = SingleBarcodeSplitter::new(&config, &barcodes);

        let rec_in = record.get_record_in_transformation();
        
        let res = splitter.apply(rec_in);
        
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
        assert!(!res.alignment_results.record_dropped);
    }
}
