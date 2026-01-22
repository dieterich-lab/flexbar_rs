use crate::config::Config;
use crate::types::{TrimEnd, RecordDto};
use crate::statistics::WorkerStatistics;
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;

pub struct SingleAdapterRemover<'a> {
    adapters: &'a [String],
    adapter_trim_end: TrimEnd,
    adapter_min_overlap: usize,
    adapter_error_rate: f32,
    removal_tags: bool,
    adapter_gap_open: i32,
    adapter_gap_extend: i32,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    max_adapter_len: usize,
}

impl<'a> SingleAdapterRemover<'a> {
    pub fn new(config: &Config, adapters: &'a [String]) -> Self {
        let (xclip_prefix, xclip_suffix, yclip_prefix, yclip_suffix)
            = match config.args.adapter_trim_end {
            TrimEnd::Any => (0, 0, 0, 0),
            TrimEnd::Left => (0, MIN_SCORE, MIN_SCORE, 0),
            TrimEnd::Right => (MIN_SCORE, 0, 0, MIN_SCORE),
            TrimEnd::LTail => (0, MIN_SCORE, 0, 0),
            TrimEnd::RTail => (MIN_SCORE, 0, 0, 0),
        };

        let max_adapter_len = adapters.iter().map(|a| a.len()).max().unwrap_or(0);

        Self {
            adapters,
            adapter_trim_end: config.args.adapter_trim_end,
            adapter_min_overlap: config.args.adapter_min_overlap,
            adapter_error_rate: config.args.adapter_error_rate,
            removal_tags: config.args.removal_tags,
            adapter_gap_open: config.args.adapter_gap_open,
            adapter_gap_extend: config.args.adapter_gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            max_adapter_len,
        }
    }

    fn get_scoring(&self) -> Scoring<fn(u8, u8) -> i32> {
        Scoring {
            match_fn: |a, b| if a == b { 1 } else { -1 },
            match_scores: None,
            gap_open: self.adapter_gap_open,
            gap_extend: self.adapter_gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    pub fn remove_adapters(&self, record: &RecordDto, stats: &mut WorkerStatistics) -> RecordDto {
        if self.adapters.is_empty() {
            return record.clone();
        }

        let mut aligner = Aligner::with_capacity_and_scoring(
            record.seq().len(),
            self.max_adapter_len,
            self.get_scoring(),
        );

        let mut best_alignment: Option<(Alignment, usize)> = None;

        for (i, adapter) in self.adapters.iter().enumerate() {
            let alignment = aligner.semiglobal(adapter.as_bytes(), record.seq());
            let overlap = alignment.yend - alignment.ystart;
            let adapter_overlap = alignment.xend - alignment.xstart;
            
            if overlap >= self.adapter_min_overlap && adapter_overlap >= self.adapter_min_overlap {
                let matches = (alignment.score + overlap as i32) / 2;
                let errors = overlap as i32 - matches;
                let error_rate = errors as f32 / overlap as f32;

                if error_rate <= self.adapter_error_rate {
                    if best_alignment.as_ref().map_or(true, |(best, _)| alignment.score > best.score) {
                        best_alignment = Some((alignment, i));
                    }
                }
            }
        }

        if let Some((alignment, i)) = best_alignment {
            let mut new_seq = record.seq().to_vec();
            let mut new_qual = record.qual().to_vec();

            match self.adapter_trim_end {
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
            
            stats.total_adapters_removed += 1;
            
            let mut id_str = record.id().to_string();
            if self.removal_tags {
                id_str.push_str(&format!("_Flexbar_removal_{}", i));
            }
            
            RecordDto::with_attrs(&id_str, record.desc(), &new_seq, &new_qual)
        } else {
            record.clone()
        }
    }
}
