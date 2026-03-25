use crate::args::Args;
use crate::types::{PairedRecordDto, SingeRecordDto, TrimEnd};
use log::trace;
use crate::steps::generic::{GenericStep, RecordInTransformation, Transformable, PairedAlignmentResults, SingleAlignmentResults};
use crate::context::{Context, SequenceWithId};
use crate::multi_aligner::{MultiAligner, MultiAlignerSetup, MultiAlignerState};

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
        let config = crate::multi_aligner::AlignerConfig {
            match_score: args.barcode_match,
            mismatch_score: args.barcode_mismatch,
            gap_open: args.barcode_gap_open,
            gap_extend: args.barcode_gap_extend,
            error_rate: args.barcode_error_rate,
            min_overlap: args.barcode_min_overlap,
            cycles: if args.barcode_keep { 1 } else { args.barcode_cycles },
        };
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
        let mut rec1_in = rec.record.records[0].get_record_in_transformation();
        let mut rec2_in = rec.record.records[1].get_record_in_transformation();

        let seq1 = rec1_in.get_sequence();
        let best1 = self.multi_aligner.align(&mut state.multi_aligner_state, &self.multi_aligner_setup, &seq1);

        let mut best2 = None;
        if let (Some(aligner2), Some(setup2), Some(state2)) = (
            &self.multi_aligner2,
            &self.multi_aligner_setup2,
            &mut state.multi_aligner_state2
        ) {
            let seq2 = rec2_in.get_sequence();
            best2 = aligner2.align(state2, setup2, &seq2);
        }

        let best_overall = match (best1, best2) {
            (Some(b1), Some(b2)) => {
                if b1.score() >= b2.score() {
                    Some((1, b1))
                } else {
                    Some((2, b2))
                }
            }
            (Some(b1), None) => Some((1, b1)),
            (None, Some(b2)) => Some((2, b2)),
            (None, None) => None,
        };

        match best_overall {
            Some((read_idx, best)) => {
                let barcode = if read_idx == 1 {
                    &self.barcodes[best.reference_index]
                } else {
                    &self.barcodes2[best.reference_index]
                };

                let tag = format!("barcode_{}", barcode.id);
                Self::assign_barcode_tag(&mut rec1_in, tag.clone(), self.trace);
                Self::assign_barcode_tag(&mut rec2_in, tag, self.trace);

                let (best_rec_in, other_rec_in) = if read_idx == 1 {
                    (&mut rec1_in, &mut rec2_in)
                } else {
                    (&mut rec2_in, &mut rec1_in)
                };

                if self.trace {
                    let bar_seq = barcode.sequence.as_bytes();
                    let sequence: &[u8] = best_rec_in.get_sequence();
                    trace!(
                        "Barcode splitter on {} selected the following alignment from Read {}:",
                        rec.record.id(),
                        read_idx
                    );
                    trace!("{}", best.format(sequence, bar_seq));
                }

                if self.umi_tags {
                    let sequence: &[u8] = best_rec_in.get_sequence();
                    let bar_seq = barcode.sequence.as_bytes();
                    let new_umi_tags = best.get_umi_tags(sequence, bar_seq);
                    best_rec_in.alignment_results.read.umi_tags.extend(new_umi_tags.clone());
                    other_rec_in.alignment_results.read.umi_tags.extend(new_umi_tags);
                }

                if !self.barcode_keep {
                    let old_len = best_rec_in.alignment_results.read.read_range.end - best_rec_in.alignment_results.read.read_range.start;
                    best_rec_in.trim(self.trim_end, best.get_start_index(), best.get_end_index());
                    let new_len = best_rec_in.alignment_results.read.read_range.end - best_rec_in.alignment_results.read.read_range.start;
                    
                    if self.removal_tags {
                        let removal_tag = format!("BR:{}:{}:{}", barcode.id, best.get_start_index(), best.get_end_index());
                        best_rec_in.alignment_results.read.removal_tags.push(removal_tag.clone());
                        other_rec_in.alignment_results.read.removal_tags.push(removal_tag);
                    }
                    rec.stats.bases_barcode_removal += old_len - new_len;
                }
            }
            None => {
                if self.barcode_unassigned {
                    let tag = "barcode_unassigned_".to_string();
                    Self::assign_barcode_tag(&mut rec1_in, tag.clone(), self.trace);
                    Self::assign_barcode_tag(&mut rec2_in, tag, self.trace);
                } else {
                    rec.alignment_results.per_read[0].dropped = true;
                    rec.alignment_results.per_read[1].dropped = true;
                    rec.stats.total_discarded_barcode += 1;
                }
            }
        }

        rec.record.records[0] = rec1_in.record;
        rec.record.records[1] = rec2_in.record;

        rec.alignment_results.barcode_tag = rec1_in.alignment_results.barcode_tag;

        rec.alignment_results.per_read[0] = rec1_in.alignment_results.read;
        rec.alignment_results.per_read[1] = rec2_in.alignment_results.read;

        rec
    }
}

impl PairedBarcodeSplitter {
    fn assign_barcode_tag(rec: &mut RecordInTransformation<SingeRecordDto, SingleAlignmentResults>, tag: String, trace: bool) {
        if trace {
            trace!("Barcode splitter on {}: Assigned tag {}", rec.record.id(), tag);
        }
        rec.alignment_results.barcode_tag = Some(tag);
    }
}
