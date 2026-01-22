use crate::config::Config;
use crate::types::{TrimEnd, RecordDto};
use crate::statistics::WorkerStatistics;
use crate::steps::generic::{GenericStep, RecordInTransformation};
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
    removal_tags: bool,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
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
            removal_tags: config.args.removal_tags,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            max_barcode_len,
        }
    }

    fn get_scoring(&self) -> Scoring<impl Fn(u8, u8) -> i32> {
        let m = self.barcode_match;
        let mm = self.barcode_mismatch;
        Scoring {
            match_fn: move |a, b| if a == b { m } else { mm },
            match_scores: None,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    pub fn split_barcodes(&self, record: &RecordDto, stats: &mut WorkerStatistics) -> (RecordDto, Option<String>) {
        let mut aligner = Aligner::with_capacity_and_scoring(
            record.seq().len(),
            self.max_barcode_len,
            self.get_scoring(),
        );

        let mut best_alignment: Option<(Alignment, &str)> = None;

        for (id, barcode) in &self.barcodes {
            let alignment = aligner.semiglobal(barcode.as_bytes(), record.seq());
            let overlap = alignment.yend - alignment.ystart;
            let barcode_overlap = alignment.xend - alignment.xstart;
            
            if overlap >= self.min_overlap && barcode_overlap >= self.min_overlap {
                let matches = (alignment.score - overlap as i32 * self.barcode_mismatch) / (self.barcode_match - self.barcode_mismatch);
                let errors = overlap as i32 - matches;
                let error_rate = errors as f32 / overlap as f32;

                if error_rate <= self.error_rate {
                    if best_alignment.as_ref().map_or(true, |(best, _)| alignment.score > best.score) {
                        best_alignment = Some((alignment, id));
                    }
                }
            }
        }

        if let Some((alignment, id)) = best_alignment {
            let mut new_seq = record.seq().to_vec();
            let mut new_qual = record.qual().to_vec();

            if !self.barcode_keep {
                match self.trim_end {
                    TrimEnd::Right | TrimEnd::RTail => {
                        new_seq.truncate(alignment.ystart);
                        new_qual.truncate(alignment.ystart);
                    }
                    TrimEnd::Left | TrimEnd::LTail => {
                        new_seq = new_seq[alignment.yend..].to_vec();
                        new_qual = new_qual[alignment.yend..].to_vec();
                    }
                    TrimEnd::Any => {
                        if alignment.ystart < record.seq().len() - alignment.yend {
                             new_seq = new_seq[alignment.yend..].to_vec();
                             new_qual = new_qual[alignment.yend..].to_vec();
                        } else {
                             new_seq.truncate(alignment.ystart);
                             new_qual.truncate(alignment.ystart);
                        }
                    }
                }
            }
            
            *stats.barcode_counts.entry(id.to_string()).or_insert(0) += 1;
            
            let mut id_str = record.id().to_string();
            if self.removal_tags && !self.barcode_keep {
                id_str.push_str("_Flexbar_removal");
            }
            
            (RecordDto::with_attrs(&id_str, record.desc(), &new_seq, &new_qual), Some(id.to_string()))
        } else {
            (record.clone(), None)
        }
    }
}

impl GenericStep<RecordDto> for SingleBarcodeSplitter {
    fn apply(&self, mut rec: RecordInTransformation<RecordDto>) -> RecordInTransformation<RecordDto> {
        if let Some(record) = rec.record {
            let (new_record, tag) = self.split_barcodes(&record, &mut rec.stats);
            rec.record = Some(new_record);
            rec.barcode_tag = tag;
        }
        rec
    }
}
