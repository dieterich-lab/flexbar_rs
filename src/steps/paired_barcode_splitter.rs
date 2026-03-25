use crate::args::Args;
use crate::types::{PairedRecordDto, TrimEnd};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils::paired::PairedAlignmentResults;
use crate::steps::utils::assign_barcode_tag;
use crate::context::{Context, SequenceWithId};
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup, MultiAlignerState};
use crate::steps::utils::barcode::get_barcode_aligner_config;

pub struct PairedBarcodeSplitter {
    pub barcodes: Vec<SequenceWithId>,
    pub barcodes2: Vec<SequenceWithId>,
    pub multi_aligner: MultiAligner,
    pub multi_aligner_setup: MultiAlignerSetup,
    pub multi_aligner2: Option<MultiAligner>,
    pub multi_aligner_setup2: Option<MultiAlignerSetup>,
    pub trim_end: TrimEnd,
    pub barcode_keep: bool,
    pub barcode_unassigned: bool,
    pub removal_tags: bool,
    pub umi_tags: bool,
    pub trace: bool,
}

pub struct PairedBarcodeSplitterState {
    pub multi_aligner_state: MultiAlignerState,
    pub multi_aligner_state2: Option<MultiAlignerState>,
}

impl GenericStep<PairedRecordDto> for PairedBarcodeSplitter {
    type State = PairedBarcodeSplitterState;
    type AlignmentResult = PairedAlignmentResults;

    fn new(args: &Args, context: &Context) -> Self {
        let config = get_barcode_aligner_config(args);
        let multi_aligner_setup = MultiAlignerSetup::new(
            "Barcode Splitter R1".to_string(),
            args.barcode_trim_end,
            &context.barcodes,
            config.clone(),
            args.barcode_kmer
        );
        let multi_aligner = multi_aligner_setup.create_aligner();

        let (multi_aligner2, multi_aligner_setup2) = if !context.barcodes2.is_empty() {
            let setup = MultiAlignerSetup::new(
                "Barcode Splitter R2".to_string(),
                args.barcode_trim_end,
                &context.barcodes2,
                config,
                args.barcode_kmer
            );
            let aligner = setup.create_aligner();
            (Some(aligner), Some(setup))
        } else {
            (None, None)
        };

        Self {
            barcodes: context.barcodes.clone(),
            barcodes2: context.barcodes2.clone(),
            multi_aligner,
            multi_aligner_setup,
            multi_aligner2,
            multi_aligner_setup2,
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
            multi_aligner_state: self.multi_aligner_setup.create_state(&self.multi_aligner),
            multi_aligner_state2: self.multi_aligner2.as_ref().zip(self.multi_aligner_setup2.as_ref())
                .map(|(a, s)| s.create_state(a)),
        }
    }

    fn apply(&self, mut rec: RecordInTransformation<PairedRecordDto, PairedAlignmentResults>, state: &mut Self::State) -> RecordInTransformation<PairedRecordDto, PairedAlignmentResults> {
        let seq1 = rec.get_sequence(0);
        let best1 = self.multi_aligner.align(&mut state.multi_aligner_state, &self.multi_aligner_setup, &seq1);

        let mut best2 = None;
        if let (Some(aligner2), Some(setup2), Some(state2)) = (
            &self.multi_aligner2,
            &self.multi_aligner_setup2,
            &mut state.multi_aligner_state2
        ) {
            let seq2 = rec.get_sequence(1);
            best2 = aligner2.align(state2, setup2, &seq2);
        }

        let best_overall = match (best1, best2) {
            (Some(b1), Some(b2)) => {
                if b1.score() >= b2.score() {
                    Some((0, b1))
                } else {
                    Some((1, b2))
                }
            }
            (Some(b1), None) => Some((0, b1)),
            (None, Some(b2)) => Some((1, b2)),
            (None, None) => None,
        };

        match best_overall {
            Some((read_idx, best)) => {
                let barcode = if read_idx == 0 {
                    &self.barcodes[best.reference_index]
                } else {
                    &self.barcodes2[best.reference_index]
                };

                let tag = format!("barcode_{}", barcode.id);
                assign_barcode_tag(&mut rec, tag, self.trace);

                let other_idx = 1 - read_idx;

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

                if self.umi_tags {
                    let sequence: &[u8] = rec.get_sequence(read_idx);
                    let bar_seq = barcode.sequence.as_bytes();
                    let new_umi_tags = best.get_umi_tags(sequence, bar_seq);
                    rec.alignment_results.per_read[read_idx].umi_tags.extend(new_umi_tags.clone());
                    rec.alignment_results.per_read[other_idx].umi_tags.extend(new_umi_tags);
                }

                if !self.barcode_keep {
                    let old_len = rec.alignment_results.per_read[read_idx].read_range.end - rec.alignment_results.per_read[read_idx].read_range.start;
                    rec.trim(read_idx, self.trim_end, best.get_start_index(), best.get_end_index());
                    let new_len = rec.alignment_results.per_read[read_idx].read_range.end - rec.alignment_results.per_read[read_idx].read_range.start;
                    
                    if self.removal_tags {
                        let removal_tag = format!("BR:{}:{}:{}", barcode.id, best.get_start_index(), best.get_end_index());
                        rec.alignment_results.per_read[read_idx].removal_tags.push(removal_tag.clone());
                        rec.alignment_results.per_read[other_idx].removal_tags.push(removal_tag);
                    }
                    rec.stats.bases_barcode_removal += old_len - new_len;
                }
            }
            None => {
                if self.barcode_unassigned {
                    assign_barcode_tag(&mut rec, "barcode_unassigned_".to_string(), self.trace);
                } else {
                    for read_idx in 0..2 {
                        if !rec.alignment_results.per_read[read_idx].dropped {
                            rec.alignment_results.per_read[read_idx].dropped = true;
                            rec.stats.total_discarded_barcode += 1;
                        }
                    }
                }
            }
        }

        rec
    }
}

