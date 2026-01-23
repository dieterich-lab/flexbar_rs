use crate::config::Config;
use crate::types::{TrimEnd, RecordDto};
use log::trace;
use crate::statistics::WorkerStatistics;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::{extract_umi_tags, get_trim_range};
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;

pub struct SingleAdapterRemover {
    adapters: Vec<String>,
    adapter_trim_end: TrimEnd,
    adapter_min_overlap: usize,
    adapter_error_rate: f32,
    removal_tags: bool,
    adapter_match: i32,
    adapter_mismatch: i32,
    adapter_gap_open: i32,
    adapter_gap_extend: i32,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    umi_tags: bool,
    max_adapter_len: usize,
}

impl SingleAdapterRemover {
    pub fn new(config: &Config, adapters: &[String]) -> Self {
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
            adapters: adapters.to_vec(),
            adapter_trim_end: config.args.adapter_trim_end,
            adapter_min_overlap: config.args.adapter_min_overlap,
            adapter_error_rate: config.args.adapter_error_rate,
            removal_tags: config.args.removal_tags,
            adapter_match: config.args.adapter_match,
            adapter_mismatch: config.args.adapter_mismatch,
            adapter_gap_open: config.args.adapter_gap_open,
            adapter_gap_extend: config.args.adapter_gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            umi_tags: config.args.umi_tags,
            max_adapter_len,
        }
    }

    fn get_scoring(&self) -> Scoring<impl Fn(u8, u8) -> i32> {
        let m = self.adapter_match;
        let mm = self.adapter_mismatch;
        Scoring {
            match_fn: move |a, b| if a == b || a == b'N' { m } else { mm },
            match_scores: None,
            gap_open: self.adapter_gap_open,
            gap_extend: self.adapter_gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        }
    }

    pub fn remove_adapters(&self, record: &RecordDto, results: &mut crate::steps::generic::AlignmentResults, stats: &mut WorkerStatistics) {
        let current_seq = &record.seq()[results.seq_start..results.seq_end];
        let mut aligner = Aligner::with_capacity_and_scoring(
            current_seq.len(),
            self.max_adapter_len,
            self.get_scoring(),
        );

        let mut best_alignment: Option<(Alignment, usize)> = None;

        for (i, adapter) in self.adapters.iter().enumerate() {
            let alignment = aligner.semiglobal(adapter.as_bytes(), current_seq);
            let overlap = alignment.yend - alignment.ystart;
            let adapter_overlap = alignment.xend - alignment.xstart;
            
            if overlap >= self.adapter_min_overlap && adapter_overlap >= self.adapter_min_overlap {
                let matches = (alignment.score - overlap as i32 * self.adapter_mismatch) / (self.adapter_match - self.adapter_mismatch);
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
            if self.umi_tags {
                results.umi_tags.extend(extract_umi_tags(&alignment, self.adapters[i].as_bytes(), current_seq));
            }

            let (rel_start, rel_end) = get_trim_range(current_seq.len(), &alignment, self.adapter_trim_end);
            let old_start = results.seq_start;
            results.seq_start = old_start + rel_start;
            results.seq_end = old_start + rel_end;
            
            stats.total_adapters_removed += 1;
            
            if self.removal_tags {
                results.removal_tags.push(format!("Flexbar_removal_{}", i));
            }
        }
    }
}

impl GenericStep<RecordDto> for SingleAdapterRemover {
    fn apply(&self, mut rec: RecordInTransformation<RecordDto>) -> RecordInTransformation<RecordDto> {
        if !rec.alignment_results.record_dropped {
            trace!("Adapter remover processing record: {}", rec.record.id());
            self.remove_adapters(&rec.record, &mut rec.alignment_results, &mut rec.stats);
        }
        rec
    }
}
