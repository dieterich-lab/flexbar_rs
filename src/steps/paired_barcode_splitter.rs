use crate::args::Args;
use crate::types::{PairedRecordDto, TrimEnd};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::steps::utils::assign_barcode_tag;
use crate::context::{Context, SequenceWithId};
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup, MultiAlignerState};
use crate::steps::utils::barcode::get_barcode_aligner_config;

use crate::multi_aligner::result::MultiAlignerResult;

pub struct PairedBarcodeSplitter {
    pub barcodes: [Vec<SequenceWithId>; 2],
    pub multi_aligners: [Option<MultiAligner>; 2],
    pub multi_aligner_setups: [Option<MultiAlignerSetup>; 2],
    pub trim_end: TrimEnd,
    pub barcode_keep: bool,
    pub barcode_unassigned: bool,
    pub removal_tags: bool,
    pub umi_tags: bool,
    pub trace: bool,
}

pub struct PairedBarcodeSplitterState {
    pub multi_aligner_states: [Option<MultiAlignerState>; 2],
}

impl PairedBarcodeSplitter {
    fn create_state_by_index(&self, read_idx: usize) -> Option<MultiAlignerState> {
        match (&self.multi_aligner_setups[read_idx], &self.multi_aligners[read_idx]) {
            (Some(setup), Some(multi_aligner))
                => Some(setup.create_state(multi_aligner)),
            _ => None,
        }
    }
}

impl GenericStep<PairedRecordDto> for PairedBarcodeSplitter {
    type State = PairedBarcodeSplitterState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, context: &Context) -> Self {
        let barcodes = [
            context.barcodes.clone(),
            context.barcodes2.clone(),
        ];
        let mut multi_aligners: [Option<MultiAligner>; 2] = [None, None];
        let mut multi_aligner_setups: [Option<MultiAlignerSetup>; 2] = [None, None];
        for read_idx in 0..2 {
            if !barcodes[read_idx].is_empty() {
                let setup = MultiAlignerSetup::new(
                    format!("Barcode Splitter R{}", read_idx + 1),
                    args.barcode_trim_end,
                    &barcodes[read_idx],
                    get_barcode_aligner_config(args),
                    args.barcode_kmer
                );
                multi_aligners[read_idx] = Some(setup.create_aligner());
                multi_aligner_setups[read_idx] = Some(setup);
            }
        }

        Self {
            barcodes,
            multi_aligners,
            multi_aligner_setups,
            trim_end: args.barcode_trim_end,
            barcode_keep: args.barcode_keep,
            barcode_unassigned: args.barcode_unassigned,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    fn create_state(&self) -> Self::State {
        Self::State {
            multi_aligner_states: [
                self.create_state_by_index(0),
                self.create_state_by_index(1),
            ]
        }
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        let mut best_overall: Option<(usize, MultiAlignerResult)> = None;

        for read_idx in 0..2 {
            if self.multi_aligner_setups[read_idx].is_none() { continue; }
            let setup = self.multi_aligner_setups[read_idx].as_ref().unwrap();
            let aligner = self.multi_aligners[read_idx].as_ref().unwrap();
            let mut state = state.multi_aligner_states[read_idx].as_mut().unwrap();
            let seq = rec.get_sequence(read_idx);
            if let Some(best) = aligner.align(&mut state, &setup, &seq) {
                match &best_overall {
                    Some((_, current_best)) if best.score() > current_best.score() => {
                        best_overall = Some((read_idx, best));
                    }
                    None => {
                        best_overall = Some((read_idx, best));
                    }
                    _ => {}
                }
            }
        }

        match best_overall {
            Some((read_idx, best)) => {
                let barcode = &self.barcodes[read_idx][best.reference_index];
                let tag = format!("barcode_{}", barcode.id);
                assign_barcode_tag(&mut rec, tag, self.trace);
                if self.trace {
                    let bar_seq = barcode.sequence.as_bytes();
                    let sequence: &[u8] = rec.get_sequence(read_idx);
                    trace!(
                        "Barcode splitter on {} selected the following alignment from Read {}:",
                        rec.record.id(),
                        read_idx + 1
                    );
                    trace!("{}", best.format(sequence, bar_seq));
                }
                self.handle_umi_tags(&mut rec, &best, barcode, read_idx);
                self.remove_barcode(&mut rec, &best, read_idx);
            }
            None => {
                self.handle_unassigned(&mut rec);
            }
        }

        rec
    }
}

impl PairedBarcodeSplitter {
    fn handle_umi_tags(
        &self,
        rec: &mut RecordInTransformation<PairedRecordDto, PairedAlignmentResults>,
        best: &MultiAlignerResult,
        barcode: &SequenceWithId,
        read_idx: usize,
    ) {
        if !self.umi_tags { return; }
        
        let sequence: &[u8] = rec.get_sequence(read_idx);
        let bar_seq = barcode.sequence.as_bytes();
        let new_umi_tags = best.get_umi_tags(sequence, bar_seq);
        
        if self.trace {
            trace!(
                "Barcode splitter on {}: Extracted UMI tags {:?}",
                rec.record.id(), new_umi_tags.clone()
            );
        }
        rec.alignment_results.per_read[read_idx].umi_tags.extend(new_umi_tags.clone());
    }

    fn remove_barcode(
        &self,
        rec: &mut RecordInTransformation<PairedRecordDto, PairedAlignmentResults>,
        best: &MultiAlignerResult,
        read_idx: usize,
    ) {
        if self.barcode_keep { return; }
        
        let old_len = rec.alignment_results.per_read[read_idx].read_range.end - rec.alignment_results.per_read[read_idx].read_range.start;
        
        rec.trim(read_idx, self.trim_end, best.get_start_index(), best.get_end_index());
        
        let new_len = rec.alignment_results.per_read[read_idx].read_range.end - rec.alignment_results.per_read[read_idx].read_range.start;
        rec.stats.bases_barcode_removal += old_len - new_len;

        if self.removal_tags {
            rec.alignment_results.per_read[read_idx].removal_tags.push("Flexbar_removal".to_string());
        }

        if self.trace {
            trace!(
                "Barcode splitter on {} shortened sequence of Read {} to length {}",
                rec.record.id(), read_idx + 1, new_len
            );
        }
    }

    fn handle_unassigned(&self, rec: &mut RecordInTransformation<PairedRecordDto, PairedAlignmentResults>) {
        if self.barcode_unassigned {
            assign_barcode_tag(rec, "barcode_unassigned_".to_string(), self.trace);
        } else {
            for read_idx in 0..2 {
                if !rec.alignment_results.per_read[read_idx].dropped {
                    rec.alignment_results.per_read[read_idx].dropped = true;
                    rec.stats.total_discarded_barcode += 1;
                }
            }
            if self.trace {
                trace!("Barcode splitter on {}: no tag assigned and DROPPED.", rec.record.id());
            }
        }
    }
}

