use crate::args::Args;
use crate::types::{TrimEnd, SingeRecordDto};
use crate::steps::generic::{GenericStep, RecordInTransformation};
use crate::steps::utils_sigalign::{extract_umi_tags_sig, get_local_aligner_from_parameters, get_semiglobal_aligner_from_parameters, trace_alignment_sig, usize_to_u32, AlignerAlgorithmParameters};
use crate::context::SequenceWithId;
use sigalign::{ReferenceBuilder, Aligner, algorithms::Local, algorithms::SemiGlobal, Reference};
use sigalign::results::{Alignment, QueryAlignment};
use log::trace;
use std::cell::RefCell;
use std::ops::DerefMut;
use thread_local::ThreadLocal;
use crate::logging::log_fatal;

struct BarCodeAlignment {
    barcode_index: usize,
    penalty: u32,
    start: usize,
    end: usize,
    alignment_for_trace: Option<Alignment>
}

struct BarcodeAligners {
    prefix: Option<Aligner<SemiGlobal>>,
    middle: Option<Aligner<Local>>,
    suffix: Option<Aligner<SemiGlobal>>,
}

pub struct SingleBarcodeSplitterV2 {
    barcodes: Vec<SequenceWithId>,
    max_barcode_len: usize,
    trim_end: TrimEnd,
    min_overlap: usize,
    error_rate: f32,
    barcode_keep: bool,
    barcode_unassigned: bool,
    removal_tags: bool,
    umi_tags: bool,
    trace: bool,
    reference: Reference,
    aligner_params: AlignerAlgorithmParameters,
    aligners: ThreadLocal<RefCell<BarcodeAligners>>,
}

impl GenericStep<SingeRecordDto> for SingleBarcodeSplitterV2 {
    fn apply(&mut self, mut rec: RecordInTransformation<SingeRecordDto>) -> RecordInTransformation<SingeRecordDto> {
        let r_seq = rec.get_sequence();
        let best = self.try_alignment(r_seq);
        self.apply_best_alignment(&mut rec, &best);
        rec
    }
}

impl SingleBarcodeSplitterV2 {
    pub fn new(args: &Args, barcodes: &Vec<SequenceWithId>) -> Self {
        let max_barcode_len = barcodes.iter().map(|b| b.sequence.len())
            .max().unwrap_or(0);
        let barcode_min_overlap = args.barcode_min_overlap.unwrap_or(max_barcode_len);
        let mut builder  = ReferenceBuilder::new()
            .ignore_base(b'N');
        for barcode in barcodes {
            builder = builder.add_target(&*barcode.id, barcode.sequence.as_bytes());
        }
        let reference = builder.build().unwrap();

        Self {
            barcodes: (*barcodes).clone(),
            max_barcode_len,
            reference,
            aligner_params: get_algorithm_parameters_from_args(args, max_barcode_len),
            aligners: ThreadLocal::new(),
            min_overlap: barcode_min_overlap,
            trim_end: args.barcode_trim_end,
            error_rate: args.barcode_error_rate,
            barcode_keep: args.barcode_keep,
            barcode_unassigned: args.barcode_unassigned,
            removal_tags: args.removal_tags,
            umi_tags: args.umi_tags,
            trace: args.trace,
        }
    }

    fn try_alignment(&self, r_seq: &[u8]) -> Option<BarCodeAlignment> {
        let mut aligners = self.get_aligners();
        if let Some(aligner) = aligners.middle.as_mut() {
            let result = aligner.align(r_seq, &self.reference);
            if let Some(x) = self.get_best_valid_alignment(
                &result, 0, false, false
            ) {
                return Some(x);
            }
        }

        if let Some(aligner) = aligners.prefix.as_mut() {
            let len = r_seq.len().min(self.max_barcode_len);
            let prefix = &r_seq[..len];
            let result = aligner.align(prefix, &self.reference);
            if let Some(x) = self.get_best_valid_alignment(
                &result, 0,false, false
            ) {
                return Some(x);
            }
        }

        if let Some(aligner) = aligners.suffix.as_mut() {
            let start_index = if r_seq.len() > self.max_barcode_len {
                r_seq.len() - self.max_barcode_len
            } else { 0 };
            let suffix = &r_seq[start_index..];
            let result = aligner.align(suffix, &self.reference);
            if let Some(x) = self.get_best_valid_alignment(
                &result, start_index, false, false
            ) {
                return Some(x);
            }
        }
        None
    }

    fn get_aligners(&self) -> impl DerefMut<Target=BarcodeAligners> + '_ {
        self.aligners
            .get_or(|| RefCell::new(self.create_new_aligners()))
            .borrow_mut()
    }

    fn create_new_aligners(&self) -> BarcodeAligners {
        let new_aligners = BarcodeAligners {
            prefix: match self.trim_end {
                TrimEnd::Left | TrimEnd::LTail | TrimEnd::Any => Some(
                    get_semiglobal_aligner_from_parameters(&self.aligner_params)
                ),
                _ => None
            },
            middle: match self.trim_end {
                TrimEnd::Any | TrimEnd::Right | TrimEnd::Left => Some(
                    get_local_aligner_from_parameters(&self.aligner_params)
                ),
                _ => None
            },
            suffix: match self.trim_end {
                TrimEnd::Right | TrimEnd::RTail | TrimEnd::Any => Some(
                get_semiglobal_aligner_from_parameters(&self.aligner_params)
                ),
                _ => None
            },
        };
        new_aligners
    }

    fn get_best_valid_alignment(
        &self, all_alignments: &QueryAlignment, offset: usize, is_prefix: bool, is_suffix: bool
    ) -> Option<BarCodeAlignment> {
        let mut best_alignment: Option<BarCodeAlignment> = None;
        for (barcode_index, barcode) in self.barcodes.iter().enumerate() {
            if barcode_index >= all_alignments.0.len() { continue; }
            for alignment in all_alignments.0[barcode_index].alignments.iter() {
                if !self.validate_alignment(alignment, &barcode, is_prefix, is_suffix) { continue; }
                match best_alignment {
                    Some(ref best) => {
                        if alignment.penalty >= best.penalty { continue }
                    },
                    None => ()
                }
                best_alignment = Some(BarCodeAlignment {
                    barcode_index,
                    penalty: alignment.penalty,
                    start: offset + alignment.position.query.0 as usize,
                    end: offset + alignment.position.query.1 as usize,
                    alignment_for_trace: if self.trace { Some(alignment.clone()) } else { None },
                });
            }
        }
        best_alignment
    }

    fn validate_alignment(
        &self, alignment: &Alignment, barcode: &SequenceWithId, is_prefix: bool, is_suffix: bool
    ) -> bool {
        let overlap = (alignment.position.query.1 - alignment.position.query.0) as usize;
        let barcode_overlap = (alignment.position.target.1 - alignment.position.target.0) as usize;
        if overlap < self.min_overlap || barcode_overlap < self.min_overlap {
            return false;
        }
        if is_prefix {
            if alignment.position.query.0 > 0 { return false }
            if alignment.penalty as f32 > overlap as f32 * self.error_rate { return false }
        }
        if is_suffix {
            let query_length = self.max_barcode_len as u32;
            if alignment.position.query.1 < query_length { return false }
            if alignment.penalty as f32 > overlap as f32 * self.error_rate { return false }
        }

        if !is_prefix && !is_suffix {
            let length = barcode.sequence.len();
            let total_penalty = (length - overlap) as u32 + alignment.penalty;
            if total_penalty as f32 > length as f32 * self.error_rate {
                return false;
            }
        }
        true
    }

    fn apply_best_alignment(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        best_alignment: &Option<BarCodeAlignment>
    ) {
        if self.trace {
            trace!("Record {}:", rec.record.id());
        }
        match best_alignment {
            Some(best) => {
                let barcode = &self.barcodes[best.barcode_index];
                if let Some(alignment) = best.alignment_for_trace.as_ref() {
                    trace!(
                            "Barcode splitter on {} selected the following alignment:", 
                            rec.record.id()
                    );
                    trace_alignment_sig(
                        alignment, rec.get_sequence(), barcode.sequence.as_bytes()
                    );
                }
                self.assign_barcode_tag(rec, format!("barcode_{}", barcode.id));
                self.handle_umi_tags(rec, &best);
                self.remove_barcode(rec, &best);
            },
            None => {
                if self.trace && log::log_enabled!(log::Level::Trace) {
                    trace!("No barcode match found.");
                }
                if self.barcode_unassigned {
                    self.assign_barcode_tag(rec, "barcode_unassigned_".to_string());
                } else {
                    self.drop_record(rec);
                }
            }
        }
    }

    fn assign_barcode_tag(&self, rec: &mut RecordInTransformation<SingeRecordDto>, tag: String) {
        if self.trace {
            trace!("Barcode splitter on {}: Assigned tag {}", rec.record.id(), tag);
        }
        rec.alignment_results.barcode_tag = Some(tag);
    }

    fn handle_umi_tags(
        &self,
        rec: &mut RecordInTransformation<SingeRecordDto>,
        best: &BarCodeAlignment
    ) {
        if !self.umi_tags { return; }

        let barcode = &self.barcodes[best.barcode_index];
        let new_umi_tags = extract_umi_tags_sig(
            best.alignment_for_trace.as_ref().unwrap(), barcode.sequence.as_bytes(), rec.get_sequence()
        );
        if self.trace {
            trace!(
                "Barcode splitter on {}: Extracted UMI tags {:?}",
                rec.record.id(), new_umi_tags.clone()
            );
        }
        rec.alignment_results.umi_tags.extend(new_umi_tags);
    }

    fn remove_barcode(&self, rec: &mut RecordInTransformation<SingeRecordDto>, best: &BarCodeAlignment) {
        if self.barcode_keep { return; }
        let rel_start = best.start - rec.alignment_results.seq_start;
        let rel_end = best.end - rec.alignment_results.seq_start;
        let old_start = rec.alignment_results.seq_start;
        match self.trim_end {
            TrimEnd::Right | TrimEnd::RTail => {
                rec.alignment_results.seq_end = old_start + rel_start;
            }
            TrimEnd::Left | TrimEnd::LTail => {
                rec.alignment_results.seq_start = old_start + rel_end;
            }
            TrimEnd::Any => {
                let length = rec.get_sequence().len();
                if rel_start < length - rel_end {
                    rec.alignment_results.seq_start = old_start + rel_start;
                } else {
                    rec.alignment_results.seq_end = old_start + rel_end;
                }
            }
        }
        if self.removal_tags {
            rec.alignment_results.removal_tags.push("Flexbar_removal".to_string());
        }
        if self.trace {
            trace!(
                "Barcode splitter on {} shortened sequence to:\n{}",
                rec.record.id(), rec.get_sequence_as_str()
            );
            if self.removal_tags {
                trace!("Barcode splitter on {}: Added removal tag", rec.record.id());
            }
        }
    }

    fn drop_record(&self, rec: &mut RecordInTransformation<SingeRecordDto>) {
        rec.stats.total_discarded_barcode += 1;
        rec.alignment_results.record_dropped = true;
        if self.trace {
            trace!("Barcode splitter on {}: no tag assigned and DROPPED.", rec.record.id());
        }
    }
}

fn get_algorithm_parameters_from_args(args: &Args, barcode_len: usize)
    -> AlignerAlgorithmParameters
{
    // Some parameters need to be carefully auto-tuned to avoid
    // "Cutoff is too low to detect the pattern" errors.
    let min_length = usize_to_u32(
        args.barcode_min_overlap.unwrap_or(barcode_len),
        "barcode_min_overlap"
    );
    let mismatch = 1; // max(1, 64u32 / min_length);
    let maximum_penalty_per_length = if args.barcode_error_rate < 0.0001f32 {
        0.0001f32  // sigalign needs a positive value here - we sort this out later!
    } else {
        args.barcode_error_rate
    };
    AlignerAlgorithmParameters {
        mismatch,
        gap_open: args.barcode_gap_open.unsigned_abs(),
        gap_extend: args.barcode_gap_extend.unsigned_abs(),
        min_length,
        max_penalty_per_length: maximum_penalty_per_length,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;
    use crate::steps::generic::Transformable;

    fn process_record(record: SingeRecordDto, args: Args) -> RecordInTransformation<SingeRecordDto>
    {
        let barcodes = vec![
            SequenceWithId {
                id: "bc1".to_string(),
                sequence: "TAGCTCTGTAGGCACCATCAAT".to_string()
            }
        ];
        let mut splitter = SingleBarcodeSplitterV2::new(&args, &barcodes);
        splitter.apply(
            record.get_record_in_transformation()
        )
    }

    fn get_record_not_matching_barcode() -> SingeRecordDto {
        SingeRecordDto::with_attrs("read1", None, b"AAAAACGT", b"IIIIIIII")
    }

    #[test]
    fn test_barcode_unassigned_false() {
        let res = process_record(
            get_record_not_matching_barcode(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq"])
        );
        
        assert!(res.alignment_results.record_dropped);
        assert_eq!(res.stats.total_discarded_barcode, 1);
    }

    #[test]
    fn test_barcode_unassigned_true() {
        let res = process_record(
            get_record_not_matching_barcode(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-match", "1", "--barcode-mismatch", "1", "--barcode-gap-open", "1", "--barcode-gap-extend", "1", "--barcode-unassigned"])
        );
        
        assert!(!res.alignment_results.record_dropped);
        assert_eq!(res.alignment_results.barcode_tag, Some("barcode_unassigned_".to_string()));
        assert_eq!(res.alignment_results.seq_end - res.alignment_results.seq_start, 8);
    }

    fn get_record_with_barcode_at_end() -> SingeRecordDto {
        SingeRecordDto::with_attrs(
            "read1", None,
            b"AAAAAATAGCTCTGTAGGCACCATCAAT",
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIII"
        )
    }

    #[test]
    fn test_barcode_trimming_right() {
        let res = process_record(
            get_record_with_barcode_at_end(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "RIGHT", "--trace"])
        );

        if !res.alignment_results.record_dropped {
            assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
            assert_eq!(res.get_sequence(), b"AAAAAA");
        }
    }

    fn get_record_with_barcode_at_start() -> SingeRecordDto {
        SingeRecordDto::with_attrs(
            "read1", None,
            b"TAGCTCTGTAGGCACCATCAATAAAAAA",
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIII"
        )
    }

    #[test]
    fn test_barcode_trimming_left() {
        let res = process_record(
            get_record_with_barcode_at_start(),
            Args::parse_from(&["flexbar_rs", "--reads", "test.fq", "--barcode-trim-end", "LEFT", "--barcode-match", "1", "--barcode-mismatch", "1", "--barcode-gap-open", "1", "--barcode-gap-extend", "1"])
        );
        
        if !res.alignment_results.record_dropped {
            assert_eq!(res.alignment_results.barcode_tag, Some("barcode_bc1".to_string()));
            assert_eq!(res.get_sequence(), b"AAAAAA");
        }
    }
}
