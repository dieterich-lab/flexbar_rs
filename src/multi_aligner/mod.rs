
mod generic;
mod bio_rust;
pub mod result;

use std::cmp::min;
use crate::types::TrimEnd;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::aligner::{AlignerState, BioRustPairwiseAligner, BioRustPairwiseRawAligner};
use crate::multi_aligner::bio_rust::aligner_banded::{BandedExtraConfig, BandedStateExtra, BioRustBandedPairwiseAligner, BioRustBandedRawAligner};
use crate::multi_aligner::generic::AlignerMode;
use crate::multi_aligner::result::MultiAlignerResult;
use log::{debug, info, warn};
use bio_rust::generic_aligner::BioRustAlignerParameters;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustAlignerState};


fn get_minimum_kmer_k() -> usize {
    // If k becomes too low, two problems may arise:
    // - The hash may become huge.
    // - The kmer approach may become inefficient.
    4
}

pub enum SpecificAligner {
    BioRust(BioRustPairwiseAligner),
    BioRustBanded(BioRustBandedPairwiseAligner),
}

impl SpecificAligner {
    pub fn setup_aligner(&self, state: &mut SpecificAlignerState, sequence_length: usize) {
        match (self, state) {
            (SpecificAligner::BioRust(a), SpecificAlignerState::BioRust(s))
                => a.setup_aligner(s, sequence_length),
            (SpecificAligner::BioRustBanded(a), SpecificAlignerState::BioRustBanded(s))
                => a.setup_aligner(s, sequence_length),
            _ => panic!("Unexpected combination of aligner and state"),
        }
    }
}

pub enum SpecificAlignerState {
    BioRust(BioRustAlignerState<BioRustPairwiseRawAligner, AlignerState>),
    BioRustBanded(BioRustAlignerState<BioRustBandedRawAligner, BandedStateExtra>),
}

pub struct MultiAligner {
    aligners: Vec<SpecificAligner>,
    shrink_aligner: Option<SpecificAligner>,
}

pub struct MultiAlignerState {
    aligner_states: Vec<SpecificAlignerState>,
    shrink_state: Option<SpecificAlignerState>,
}

#[derive(Clone)]
pub struct BioRustPairwiseAlignerSetup {
    config: BioRustAlignerConfig<()>,
}

#[derive(Clone)]
pub struct BioRustBandedPairwiseAlignerSetup {
    config: BioRustAlignerConfig<BandedExtraConfig>,
}

#[derive(Clone)]
pub enum SpecificAlignerSetup {
    BioRust(BioRustPairwiseAlignerSetup),
    BioRustBanded(BioRustBandedPairwiseAlignerSetup),
}

impl SpecificAlignerSetup {
    fn get_mode(&self) -> AlignerMode {
        match self {
            SpecificAlignerSetup::BioRust(s) => s.config.mode,
            SpecificAlignerSetup::BioRustBanded(s) => s.config.mode,
        }
    }

    fn get_match_score(&self) -> i32 {
        match self {
            SpecificAlignerSetup::BioRust(s) => s.config.match_score,
            SpecificAlignerSetup::BioRustBanded(s) => s.config.match_score,
        }
    }

    fn is_perfect_alignment(&self, result: &MultiAlignerResult) -> bool {
        match self.get_mode() {
            AlignerMode::Local => {
                    let match_score = self.get_match_score();
                    let ref_length = result.reference_range.end - result.reference_range.start;
                    let perfect_score = match_score * ref_length as i32;
                    result.score() >= perfect_score
            },
            _ => false,
        }
    }
}

pub type AlignerConfig = BioRustAlignerParameters;

#[derive(Clone)]
pub struct MultiAlignerSetup {
    references: Vec<SequenceWithId>,
    config: AlignerConfig,
    aligner_setups: Vec<SpecificAlignerSetup>,
    shrink_setup: Option<SpecificAlignerSetup>,
    trim_end: TrimEnd,
}

impl MultiAlignerSetup {
    pub fn new(
        name: String,
        trim_end: TrimEnd,
        references: &[SequenceWithId],
        config: BioRustAlignerParameters,
        k_override: Option<usize>,
    ) -> Self {
        let kmer_k = determine_kmer_k(
            &name, references, config.error_rate, config.min_overlap, k_override
        );
        if kmer_k > 0 {
            debug!("{}: Using banded aligner with k={} with trim end {:?}.", name, kmer_k, trim_end);
        } else {
            debug!("{}: Using non-banded aligner with trim end {:?}.", name, trim_end);
        }

        let required_modes = Self::get_aligner_mode_order(trim_end, &config);

        let mut aligner_setups:Vec<SpecificAlignerSetup> = Vec::new();
        for mode in required_modes {
            // We do not trust the banded aligner with semi-global alignments.
            if kmer_k > 0 && mode != AlignerMode::LeftOverlap && mode != AlignerMode::RightOverlap {
                aligner_setups.push(Self::create_banded_setup(
                    &name, &config, mode, kmer_k, references
                ));
            } else {
                aligner_setups.push(Self::create_standard_setup(
                    &name, &config, mode, references
                ));
            }
        }
        let shrink_setup = if config.iterate {
            Some(Self::create_standard_setup(&name, &config, AlignerMode::Local, references))
        } else {
            None
        };
        MultiAlignerSetup { references: references.to_vec(), config, aligner_setups, shrink_setup, trim_end }
    }

    fn get_aligner_mode_order(trim_end: TrimEnd, config: &BioRustAlignerParameters) -> Vec<AlignerMode> {
        match trim_end {
            TrimEnd::Any => if config.min_overlap.is_some() {
                // We need to do three alignments in this case:
                // 1. A local one to find alignments embedded in the read.
                // 2. A semi-global one to find alignments overlapping with the left end.
                // 3. Another one for the right end.
                // It seems that the semi-global alignment is slow (why? Have not found out yet).
                // As a workaround we try first a local alignment on the end with the
                // minimal overlap on just the end of the read and skip the semi-global
                // one if it fails.
                vec![
                    AlignerMode::Local,
                    AlignerMode::LeftMin, AlignerMode::LeftOverlap,
                    AlignerMode::RightMin, AlignerMode::RightOverlap
                ]
            } else {
                // If no min_overlap is specified, the reference must be completely present in
                // the read. So we only need a local alignment.
                vec![AlignerMode::Local]
            }
            TrimEnd::Left => if config.min_overlap.is_some() {
                vec![AlignerMode::Local, AlignerMode::LeftMin, AlignerMode::LeftOverlap]
            } else {
                vec![AlignerMode::Local]
            },
            TrimEnd::Right => if config.min_overlap.is_some() {
                vec![AlignerMode::Local, AlignerMode::RightMin, AlignerMode::RightOverlap]
            } else {
                vec![AlignerMode::Local]
            },
            TrimEnd::LTail => vec![AlignerMode::LeftMin, AlignerMode::LeftOverlap],
            TrimEnd::RTail => vec![AlignerMode::RightMin, AlignerMode::RightOverlap],
        }
    }

    fn create_banded_setup(
        name: &str,
        config: &BioRustAlignerParameters,
        mode: AlignerMode,
        kmer_k: usize,
        _references: &[SequenceWithId],
    ) -> SpecificAlignerSetup {
        debug!(
            "{}: Using banded aligner for mode {:?}, k={}.", name, mode, kmer_k
        );
        let aligner_config = get_aligner_config(
            config, mode, BandedExtraConfig {
                k: kmer_k,
                min_overlap: config.min_overlap
            }
        );
        SpecificAlignerSetup::BioRustBanded(
            BioRustBandedPairwiseAlignerSetup { config: aligner_config }
        )
    }

    fn create_standard_setup(
        name: &str,
        config: &BioRustAlignerParameters,
        mode: AlignerMode,
        _references: &[SequenceWithId],
    ) -> SpecificAlignerSetup {
        debug!("{}: Using standard aligner for mode {:?}.", name, mode);
        let aligner_config = get_aligner_config(
            config, mode, ()
        );
        SpecificAlignerSetup::BioRust(
            BioRustPairwiseAlignerSetup { config: aligner_config }
        )
    }

    pub fn create_aligner(&self) -> MultiAligner {
        let aligners: Vec<SpecificAligner> = self.aligner_setups.iter()
        .map(|setup| match setup {
                SpecificAlignerSetup::BioRust(x) 
                    => SpecificAligner::BioRust(
                        BioRustPairwiseAligner::new(
                            &x.config, &self.references, self.config.error_rate
                        )
                ),
                SpecificAlignerSetup::BioRustBanded(x) 
                    => SpecificAligner::BioRustBanded(
                    BioRustBandedPairwiseAligner::new(
                        &x.config, &self.references, self.config.error_rate
                    )
                ),
        })
        .collect();
        let shrink_aligner = self.shrink_setup.as_ref().map(|setup| match setup {
                SpecificAlignerSetup::BioRust(x) 
                    => SpecificAligner::BioRust(
                        BioRustPairwiseAligner::new(
                            &x.config, &self.references, self.config.error_rate
                        )
                ),
                SpecificAlignerSetup::BioRustBanded(x) 
                    => SpecificAligner::BioRustBanded(
                    BioRustBandedPairwiseAligner::new(
                        &x.config, &self.references, self.config.error_rate
                    )
                ),
        });
        MultiAligner { aligners, shrink_aligner }
    }

    pub fn create_state(&self, aligner: &MultiAligner) -> MultiAlignerState {
        let aligner_states = aligner.aligners.iter().map(|a| match a {
            SpecificAligner::BioRust(a) => SpecificAlignerState::BioRust(a.create_state(&self.references)),
            SpecificAligner::BioRustBanded(a) => SpecificAlignerState::BioRustBanded(a.create_state(&self.references)),
        }).collect();
        let shrink_state = aligner.shrink_aligner.as_ref().map(|a| match a {
            SpecificAligner::BioRust(a) => SpecificAlignerState::BioRust(a.create_state(&self.references)),
            SpecificAligner::BioRustBanded(a) => SpecificAlignerState::BioRustBanded(a.create_state(&self.references)),
        });
        MultiAlignerState { aligner_states, shrink_state }
    }
}

fn get_aligner_config<Extra: Clone>(
    config: &BioRustAlignerParameters, mode: AlignerMode, extra: Extra
) -> BioRustAlignerConfig<Extra> {
    BioRustAlignerConfig {
        mode,
        match_score: config.match_score,
        mismatch_score: config.mismatch_score,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        min_overlap: config.min_overlap,
        extra,
    }
}


macro_rules! define_align_mode_func {
    ($name:ident) => {
        #[cfg_attr(not(debug_assertions), inline(always))]
        fn $name(
            aligner: &SpecificAligner,
            state: &mut SpecificAlignerState,
            sequence: &[u8],
            index: usize,
        ) -> MultiAlignerResult {
            match (aligner, state) {
                (SpecificAligner::BioRust(aligner), SpecificAlignerState::BioRust(s)) => {
                    aligner.align(s, sequence, index)
                },
                (SpecificAligner::BioRustBanded(aligner), SpecificAlignerState::BioRustBanded(s)) => {
                    aligner.align(s, sequence, index)
                }
                _ => panic!("Unexpected combination of aligner and setup"),
            }
        }
    };
}

impl MultiAligner {
    pub fn align(&self, state: &mut MultiAlignerState, setup: &MultiAlignerSetup, sequence: &[u8]) -> Option<MultiAlignerResult> {
        if sequence.len() == 0 { return None }  // The aligner seems to go berserk on sequences of size 0
        for (aligner, state) in self.aligners.iter().zip(state.aligner_states.iter_mut()) {
            aligner.setup_aligner(state, sequence.len());
        }
        let best = self.find_best_match(state, setup, sequence);
        if let Some(best) = best {
            if !best.may_need_shrink || !setup.config.iterate {
                return Some(best);
            }
            // Prepare the shrink_state
            if let Some(shrink_state) = state.shrink_state.as_mut() {
                self.shrink_aligner.as_ref().unwrap().setup_aligner(shrink_state, sequence.len());
            }
            return self.shrink(state, setup, sequence, best);
        }
        best
    }

    fn find_best_match(&self, state: &mut MultiAlignerState, setup: &MultiAlignerSetup, sequence: &[u8]) -> Option<MultiAlignerResult> {
        let mut best: Option<MultiAlignerResult> = None;

        'outer: for index in 0..setup.references.len() {
            let mut skip_next = false;
            let mut may_need_shrink = true;
            for ((some_setup, some_aligner), some_state) in
                setup.aligner_setups.iter().zip(self.aligners.iter())
                .zip(state.aligner_states.iter_mut())
            {
                if skip_next { skip_next = false; continue }
                let mut result= Self::do_alignment(
                    some_setup, some_aligner, some_state, sequence, index
                );
                let mode = some_setup.get_mode();
                let is_precheck = match mode {
                    AlignerMode::LeftMin | AlignerMode::RightMin => true,
                    _ => false,
                };
                let reference = &setup.references[index].sequence.as_bytes();
                if result.validate(
                    setup.config.error_rate, sequence, reference, setup.config.min_overlap, result.sequence_offset
                ) {
                    if !is_precheck {
                        result.may_need_shrink = may_need_shrink;
                        if get_score_from_best(&best) < result.score() {
                            best = Some(result);
                        }
                    }
                } else {
                    if is_precheck {
                        skip_next = true;
                    }
                    if mode == AlignerMode::Local {
                        may_need_shrink = false;
                    }
                }
                if let Some(res) = &best && some_setup.is_perfect_alignment(res) {
                    break 'outer;
                }
            }
        }
        best
    }

    fn shrink(&self, state: &mut MultiAlignerState, setup: &MultiAlignerSetup, sequence: &[u8], initial_best: MultiAlignerResult) -> Option<MultiAlignerResult> {
        match setup.trim_end {
            TrimEnd::Right => self.shrink_right(state, setup, sequence, initial_best),
            TrimEnd::RTail => self.shrink_rtail(state, setup, sequence, initial_best),
            TrimEnd::Left => self.shrink_left(state, setup, sequence, initial_best),
            TrimEnd::LTail => self.shrink_ltail(state, setup, sequence, initial_best),
            TrimEnd::Any => {
                let res_right = self.shrink_right(state, setup, sequence, initial_best.clone());
                let res_left = self.shrink_left(state, setup, sequence, initial_best);
                match (res_right, res_left) {
                    (Some(r), Some(l)) => {
                        let remaining_right = r.get_start_index();
                        let remaining_left = sequence.len() - l.get_end_index();
                        if remaining_right >= remaining_left { Some(r) } else { Some(l) }
                    },
                    (Some(r), None) => Some(r),
                    (None, Some(l)) => Some(l),
                    (None, None) => None,
                }
            }
        }
    }

    fn shrink_right(&self, state: &mut MultiAlignerState, setup: &MultiAlignerSetup, sequence: &[u8], mut best: MultiAlignerResult) -> Option<MultiAlignerResult> {
        let aligner = self.shrink_aligner.as_ref().unwrap();
        let s_state = state.shrink_state.as_mut().unwrap();
        let ref_index = best.reference_index;
        let reference = &setup.references[ref_index].sequence.as_bytes();
        loop {
            let next_len = best.get_start_index() + best.stripped_leading_n_count;
            if next_len == 0 { break }
            let next_seq = &sequence[..next_len];
            aligner.setup_aligner(s_state, next_seq.len());
            let result = Self::align_local(aligner, s_state, next_seq, ref_index);
            if result.validate(setup.config.error_rate, next_seq, reference, setup.config.min_overlap, result.sequence_offset) {
                if result.score() < best.score() {
                    break;
                }
                if result.get_start_index() == best.get_start_index() {
                    break;
                }
                best = result;
            } else {
                break;
            }
        }
        Some(best)
    }

    fn shrink_rtail(&self, state: &mut MultiAlignerState, setup: &MultiAlignerSetup, sequence: &[u8], mut best: MultiAlignerResult) -> Option<MultiAlignerResult> {
        let aligner = self.shrink_aligner.as_ref().unwrap();
        let s_state = state.shrink_state.as_mut().unwrap();
        let ref_index = best.reference_index;
        let reference_bytes = &setup.references[ref_index].sequence.as_bytes();
        let ref_len = reference_bytes.len();
        loop {
            let start = best.get_start_index();
            let window_start = start.saturating_sub(ref_len);
            let window_end = (start + best.stripped_leading_n_count).min(sequence.len());
            if window_start >= window_end { break }
            let next_seq = &sequence[window_start..window_end];
            aligner.setup_aligner(s_state, next_seq.len());
            let mut result = Self::align_local(aligner, s_state, next_seq, ref_index);
            result.sequence_offset += window_start;
            if result.validate(setup.config.error_rate, sequence, reference_bytes, setup.config.min_overlap, result.sequence_offset) {
                if result.get_start_index() >= best.get_start_index() {
                    break;
                }
                best = result;
            } else {
                break;
            }
        }
        Some(best)
    }

    fn shrink_left(&self, state: &mut MultiAlignerState, setup: &MultiAlignerSetup, sequence: &[u8], mut best: MultiAlignerResult) -> Option<MultiAlignerResult> {
        let aligner = self.shrink_aligner.as_ref().unwrap();
        let s_state = state.shrink_state.as_mut().unwrap();
        let ref_index = best.reference_index;
        let reference = &setup.references[ref_index].sequence.as_bytes();
        loop {
            let next_start = best.get_end_index() - best.stripped_trailing_n_count;
            if next_start >= sequence.len() { break }
            let next_seq = &sequence[next_start..];
            aligner.setup_aligner(s_state, next_seq.len());
            let mut result = Self::align_local(aligner, s_state, next_seq, ref_index);
            result.sequence_offset += next_start;
            if result.validate(setup.config.error_rate, sequence, reference, setup.config.min_overlap, result.sequence_offset) {
                if result.get_end_index() == best.get_end_index() {
                    break;
                }
                best = result;
            } else {
                break;
            }
        }
        Some(best)
    }

    fn shrink_ltail(&self, state: &mut MultiAlignerState, setup: &MultiAlignerSetup, sequence: &[u8], mut best: MultiAlignerResult) -> Option<MultiAlignerResult> {
        let aligner = self.shrink_aligner.as_ref().unwrap();
        let s_state = state.shrink_state.as_mut().unwrap();
        let ref_index = best.reference_index;
        let reference_bytes = &setup.references[ref_index].sequence.as_bytes();
        let ref_len = reference_bytes.len();
        loop {
            let end = best.get_end_index();
            let window_start = end.saturating_sub(ref_len);
            let window_end = (end + ref_len).min(sequence.len());
            if window_start >= window_end { break }
            let next_seq = &sequence[window_start..window_end];
            aligner.setup_aligner(s_state, next_seq.len());
            let mut result = Self::align_local(aligner, s_state, next_seq, ref_index);
            result.sequence_offset += window_start;
            if result.validate(setup.config.error_rate, sequence, reference_bytes, setup.config.min_overlap, result.sequence_offset) {
                if result.get_end_index() <= best.get_end_index() {
                    break;
                }
                best = result;
            } else {
                break;
            }
        }
        Some(best)
    }

    #[inline(always)]
    fn do_alignment(setup: &SpecificAlignerSetup, aligner: &SpecificAligner, state: &mut SpecificAlignerState, sequence: &[u8], reference_index: usize) -> MultiAlignerResult {
        let mode = setup.get_mode();
        // The split by mode is technically completely unnecessary. Its only purpose
        // is to allow more meaningful profiling runs. Hopefully, the optimizer will
        // inline the three functions called below and remove the match because all
        // do exactly the same.
        match mode {
            AlignerMode::Local => Self::align_local(aligner, state, sequence, reference_index),
            AlignerMode::LeftOverlap => Self::align_left_overlap(aligner, state, sequence, reference_index),
            AlignerMode::RightOverlap => Self::align_right_overlap(aligner, state, sequence, reference_index),
            AlignerMode::LeftMin => Self::align_left_min(aligner, state, sequence, reference_index),
            AlignerMode::RightMin => Self::align_right_min(aligner, state, sequence, reference_index),
        }
    }

    define_align_mode_func!(align_local);
    define_align_mode_func!(align_left_overlap);
    define_align_mode_func!(align_right_overlap);
    define_align_mode_func!(align_left_min);
    define_align_mode_func!(align_right_min);
}

fn get_score_from_best(best: &Option<MultiAlignerResult>) -> i32 {
    match best {
        Some(b) => b.score(),
        None => i32::MIN,
    }
}

fn determine_kmer_k(
    name: &str,
    references: &[SequenceWithId],
    error_rate: f32,
    min_overlap: Option<usize>,
    k_override: Option<usize>
) -> usize {
    let mut highest_possible_k = usize::MAX;
    for r in references.iter() {
        let possible_k = determine_kmer_k_for_sequence(&r.sequence, error_rate, min_overlap);
        if possible_k < highest_possible_k {
            highest_possible_k = possible_k;
        }
    }
    match k_override {
        Some(0) => {
            info!("{}: Using non-banded aligner on user request.", name);
            0
        },
        Some(val) => {
            if val > highest_possible_k {
                warn!(
                        "Warning: The k-mer size ({}, <{}) for {} might be too big.",
                        val, highest_possible_k, name
                    );
                warn!("Warning: Valid alignments for the given error rate/overlap may be missed.");
            }
            val
        },
        None => {
            if highest_possible_k > get_minimum_kmer_k() {
                highest_possible_k
            } else {
                warn!(
                    "Warning: The k-mer size ({}) for {} is too small - using non-banded aligner.",
                    highest_possible_k, name
                );
                warn!("Warning: Performance might be improved using a lower error rate and/or a larger overlap");
                0
            }
        },
    }
}

fn determine_kmer_k_for_sequence(
    sequence: &str,
    error_rate: f32,
    min_overlap: Option<usize>,
) -> usize {
    let mut max_clean_length: usize = 0;
    let mut present_clean_length: usize = 0;
    let mut n_count = 0;
    for c in sequence.chars() {
        if c == 'N' {
            if present_clean_length > max_clean_length {
                max_clean_length = present_clean_length;
            }
            present_clean_length = 0;
            n_count += 1;
        } else {
            present_clean_length += 1;
        }
    }
    if present_clean_length > max_clean_length {
        max_clean_length = present_clean_length;
    }
    let allowed_errors = (max_clean_length as f32 * error_rate) as usize;
    let min_match_length = max_clean_length / (allowed_errors + 1);
    match min_overlap {
        Some(min_overlap) => {
            if min_overlap < n_count {
                // That is a pathological case, which will give a warning anyway.
                0
            } else {
                min(min_match_length, min_overlap - n_count)
            }
        },
        None => min_match_length,
    }
}

#[cfg(test)]
mod tests {
    use crate::context::SequenceWithId;
    use crate::multi_aligner::bio_rust::generic_aligner::BioRustAlignerParameters;
    use crate::types::TrimEnd;
    use super::{MultiAlignerSetup, SpecificAlignerSetup};
    use super::result::MultiAlignerResult;

    fn run_alignment(
        trim_end: TrimEnd,
        config: BioRustAlignerParameters,
        k_overwrite: Option<usize>,
        seq_x: &[u8],
        references: &[SequenceWithId]
    ) -> Option<MultiAlignerResult> {
        let setup = MultiAlignerSetup::new(
            "Test Aligner".to_string(),
            trim_end,
            references,
            config,
            k_overwrite,
        );
        let aligner = setup.create_aligner();
        let mut state = setup.create_state(&aligner);
        aligner.align(&mut state, &setup, seq_x)
    }

    fn get_as_references(raw: &[u8]) -> Vec<SequenceWithId> {
        let mut references = Vec::new();
        references.push(
            SequenceWithId {
                id: "ref".to_string(),
                sequence: String::from_utf8_lossy(raw).to_string(),
            }
        );
        references
    }

    fn get_default_config() -> BioRustAlignerParameters {
        BioRustAlignerParameters {
            error_rate: 0.0,
            min_overlap: Some(3),
            match_score: 1,
            mismatch_score: -1,
            gap_open: -1,
            gap_extend: -1,
            iterate: false,
        }
    }

    fn get_overlap5_config() -> BioRustAlignerParameters {
        let mut config = get_default_config();
        config.min_overlap = Some(5);
        config
    }

    fn get_benchmark_config() -> BioRustAlignerParameters {
        BioRustAlignerParameters {
            error_rate: 0.0,
            min_overlap: Some(11),
            match_score: 1,
            mismatch_score: -1,
            gap_open: -9,
            gap_extend: -1,
            iterate: false,
        }
    }

    fn get_with_tolerance_config() -> BioRustAlignerParameters {
        let mut config = get_default_config();
        config.error_rate = 0.2;
        config
    }

    fn is_banded_setup(setup: &MultiAlignerSetup) -> bool {
        match setup.aligner_setups.get(0) {
            Some(SpecificAlignerSetup::BioRustBanded(_)) => true,
            _ => false,
        }
    }

    #[test]
    fn test_aligner_trivial() {
        let references = get_as_references(b"ACGT");
        let result = run_alignment(
            TrimEnd::Right, get_default_config(), None,
            b"ACGT", &references
        );
        assert!(result.is_some());
        let alignment = result.unwrap();
        assert_eq!(alignment.score(), 4);
        assert_eq!(alignment.get_start_index(), 0);
        assert_eq!(alignment.get_end_index(), 4);
    }

    #[test]
    fn test_aligner_right_overhang() {
        let references = get_as_references(b"NNNTTTTCCC");
        let result = run_alignment(
            TrimEnd::Right, get_overlap5_config(), None,
            b"AAAAATTT", &references
        );
        let alignment = result.unwrap();
        // Matching prefix NNNTTT at position 2
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 2);
        assert_eq!(alignment.get_end_index(), 8);
    }

    #[test]
    fn test_aligner_right_overhang_ambiguous() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Right, get_default_config(), None,
            b"AAAAATTT", &references
        );
        let alignment = result.unwrap();
        // Matching prefix TTT at position 5
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 8);
    }

    #[test]
    fn test_aligner_right_overhang_with_error() {
        let references = get_as_references(b"TTTTTTT");
        let result = run_alignment(
            TrimEnd::Right, get_with_tolerance_config(), None,
            b"AAAAATTATT", &references
        );
        let alignment = result.unwrap();
        // Matching prefix TTxTT with an error starting at x[5]
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_embedded() {
        let references = get_as_references(b"TTTTT");
        let result = run_alignment(
            TrimEnd::Right, get_overlap5_config(), None,
            b"AAAAATTTTTAAA", &references
        );
        let alignment = result.unwrap();
        // Prefect match att position 5
        assert_eq!(alignment.score(), 5);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_embedded_with_error() {
        let references = get_as_references(b"TTTTT");
        let result = run_alignment(
            TrimEnd::Right, get_with_tolerance_config(), None,
            b"AAAAATTATTAAA", &references
        );
        let alignment = result.unwrap();
        // Match with error TTxTT at position 5
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_embedded_with_wildcard() {
        let references = get_as_references(b"TTNTT");
        let result = run_alignment(
            TrimEnd::Right, get_overlap5_config(), None,
            b"AAAAATTATTAAA", &references
        );
        let alignment = result.unwrap();
        // Matching TTNTT at position 5
        assert_eq!(alignment.score(), 5);
        assert_eq!(alignment.get_start_index(), 5);
        assert_eq!(alignment.get_end_index(), 10);
    }

    #[test]
    fn test_aligner_right_wrong_side_overhang() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Right, get_with_tolerance_config(), None,
            b"TTTAAAA", &references
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_aligner_right_partial_embedded() {
        // This a strange case: The multi_aligner will score an invalid barcode match
        // in the middle of the string higher than an exact prefix match at the end.
        // So the same read may be mapped to a barcode with TrimEnd::RTail but not
        // with TrimEnd::Right.
        let references = get_as_references(b"NNNTTTTTCCCC");
        let mut config = get_default_config();
        config.min_overlap = Some(6);
        let result = run_alignment(
            TrimEnd::Right, config, None,
            b"AAAAATTTTTAAATTT", &references
        );
        let alignment = result.unwrap();
        // Matching NNNTTTTTT (prefix) at position 2
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 10);
        assert_eq!(alignment.get_end_index(), 16);
    }

    #[test]
    fn test_aligner_left_overhang() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Left, get_default_config(), None,
            b"TTTAAAAA", &references
        );
        let alignment = result.unwrap();
        // Match TTT (suffix) starting x[0]
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 0);
        assert_eq!(alignment.get_end_index(), 3);
    }

    #[test]
    fn test_aligner_rtail_overhang() {
        // See also: test_aligner_right_partial_embedded()!!!
        let references = get_as_references(b"NNNTTTTTCCCC");
        let result = run_alignment(
            TrimEnd::RTail, get_overlap5_config(), None,
            b"AAAAATTTTTAAATTT", &references
        );
        let alignment = result.unwrap();
        // Match NNNTTT starting x[10]
        assert_eq!(alignment.score(), 3);
        assert_eq!(alignment.get_start_index(), 10);
        assert_eq!(alignment.get_end_index(), 16);
    }

    #[test]
    fn test_aligner_k_override_zero() {
        let references = get_as_references(b"ACGTACGTACGT"); // Length 12
        // Default estimate_k for error_rate 0.0 is max_no_n = 12
        // 12 > 8, so default k would be 8 (banded aligner used since 8 >= 4)
        
        let setup = MultiAlignerSetup::new(
            "Test Aligner".to_string(),
            TrimEnd::Right,
            &references,
            get_overlap5_config(),
            Some(0)  // Override k to 0 should force the old aligner
        );

        assert!(!is_banded_setup(&setup));
    }

    #[test]
    fn test_aligner_k_override_banded() {
        let references = get_as_references(b"ACGT");
        // Default estimate_k for length 4 is 4.
        
        let setup = MultiAlignerSetup::new(
            "Test Aligner".to_string(),
            TrimEnd::Right,
            &references,
            get_default_config(),
            Some(5), // Override k to 5 should use banded aligner even if k=4 was estimated
        );
        assert!(is_banded_setup(&setup));
    }

    #[test]
    fn test_determine_k() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, None, None
        );
        // We need to detect matches of length 22 (28-6)
        assert_eq!(k, 22);
    }

    #[test]
    fn test_determine_k_small_overlap() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, Some(11), None
        );
        // With an overlap of 11 and 6 'N's we must be able to detect
        // matches of length 5. So K can be at most 5.
        assert_eq!(k, 5);
    }

    #[test]
    fn test_determine_k_with_min_overlap_low() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTGC".to_string(), // length 9, 6 Ns, max_no_n = 3
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, Some(8), None
        );
        // With only 3 bases other than 'N', we need a k below 4 - so we disable
        // the banded aligner.
        assert_eq!(k, 0);
    }

    #[test]
    fn test_determine_k_with_min_overlap_high() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.0, None, Some(0)
        );
        // Normale 22 - but the user disabled the banded aligner.
        assert_eq!(k, 0);
    }

    #[test]
    fn test_determine_k_with_error_rate() {
        let references = vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: "NNNNNNTAGCTCTGTAGGCACCATCAAT".to_string(), // length 28, 6 Ns
        }];
        let k = super::determine_kmer_k(
            "Test", &references, 0.1, None, None
        );
        // With error rate 0.0 we would need to detect matches of length 22 (28-6).
        // With an error rate of 0.1, we may see 2 errors in 22 bases, splitting
        // the alignment worst case into 3 parts of about 7 bases.
        assert_eq!(k, 7);
    }

    #[test]
    fn test_format_alignment_wildcard() {
        let x = b"ACGTATCG";
        let y = b"ACGTNTCG";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);
        assert_eq!(formatted, concat!(
            "Alignment with score 8:\n",
            "  ACGTATCG\n",
            "  ||||.|||\n",
            "  ACGTNTCG\n"
        ));
    }

    #[test]
    fn test_format_alignment_embedded() {
        let x = b"AAAAACCCCTTTTT";
        let y = b"CCCC";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);
        // Expected alignment:
        // AAAAACCCCTTTTT
        //      ||||
        //      CCCC
        assert!(formatted.contains("AAAAACCCCTTTTT"));
        assert!(formatted.contains("||||"));
        assert!(formatted.contains("CCCC"));
    }

    #[test]
    fn test_format_alignment_prefix_suffix() {
        let x = b"AAACCC";
        let y = b"CCCGGG";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        // AAACCC---
        //    |||
        // ---CCCGGG

        assert!(formatted.contains("AAACCC"));
        assert!(formatted.contains("|||"));
        assert!(formatted.contains("CCCGGG"));
    }

    #[test]
    fn test_format_alignment_suffix_prefix() {
        let x = b"CCCGGG";
        let y = b"AAACCC";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        // ---CCCGGG
        //    |||
        // AAACCC---

        assert!(formatted.contains("CCCGGG"));
        assert!(formatted.contains("|||"));
        assert!(formatted.contains("AAACCC"));
    }

    #[test]
    fn test_format_alignment_partial_with_error() {
        let x = b"CCCCCGGG";
        let y = b"AAACCACC";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_with_tolerance_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert!(formatted.contains("CCCCCGGG"));
        assert!(formatted.contains("|\\|"));
        assert!(formatted.contains("AAACCACC"));
        assert_eq!(formatted, concat!(
        "Alignment with score 3:\n",
        "     CCCCCGGG\n",
        "  ,,,||\\||´´´\n",
        "  AAACCACC   \n"
        ));
    }

    #[test]
    fn test_format_alignment_partial_with_deletion() {
        let x = b"CAATCCTTTT";
        let y = b"CCCCAACCT";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_with_tolerance_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert_eq!(formatted, concat!(
        "Alignment with score 4:\n",
        "     CAATCCTTTT\n",
        "  ,,,|||´|||´´´\n",
        "  CCCCAA CCT   \n"
        ));
    }

    #[test]
    fn test_format_alignment_partial_with_insertion() {
        let x = b"CAACCTTTT";
        let y = b"CCCCAATCCT";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_with_tolerance_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert_eq!(formatted, concat!(
        "Alignment with score 4:\n",
        "     CAA CCTTTT\n",
        "  ,,,|||,|||´´´\n",
        "  CCCCAATCCT   \n"
        ));
    }

    #[test]
    fn test_umi_tag() {
        let x = b"CAATCCTTT";
        let y = b"NNNTTT";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, &references
        );
        let alignment = result.unwrap();
        assert_eq!(alignment.get_umi_tags(x, y), vec![String::from("TCC")])
    }

    #[test]
    fn test_offset_handling() {
        let x = b"CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA";
        let y =b"NNNNNNTAGCTCTGTAGGCACCATCAAT";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, &references
        );
        assert!(result.is_some());
        let alignment = result.unwrap();
        let umi_tags = alignment.get_umi_tags(x, y);
        
        assert_eq!(umi_tags.len(), 1);
        assert_eq!(umi_tags[0], "GTGCCC");

        // Check format
        let formatted = alignment.format(x, y);
        assert_eq!(formatted, concat!(
            "Alignment with score 21:\n",
            "  CATCAGTGCGCCCCGGGCGGGTCGCGCCGTGTGGGGGGGCCGGGGTGCCCTAGCTCTGTAGGCACCATCAA \n",
            "  ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´......|||||||||||||||||||||,\n",
            "                                              NNNNNNTAGCTCTGTAGGCACCATCAAT\n"
        ));
    }

    #[test]
    fn test_small_overlap_alignment() {
        let x = b"TATCGTTTTTTCACTGACCCGGTGAGGCGGGGTGTACGCATACTGTAGGCACCATCGTAGAGTAGCTCTGT";
        let y = b"NNNNNNTAGCTCTGTAGGCACCATCAAT";
        let references = get_as_references(y);

        let result = run_alignment(
            TrimEnd::Right, get_benchmark_config(), None,
            x, &references
        );
        
        assert!(result.is_some(), "Alignment should be found");
        let alignment = result.unwrap();
        assert_eq!(alignment.get_start_index(), 56);
        
        let formatted = alignment.format(x, y);
        println!("{}", formatted);

        let umi_tags = alignment.get_umi_tags(x, y);
        assert_eq!(umi_tags[0], "GTAGAG");
        
        // Verify it matches the end of x
        assert_eq!(alignment.get_end_index(), x.len());
    }

    #[test]
    fn test_min_overlap_alignment() {
        let x = b"TGCTCTCGCTTCTGGCGCCAAGCGCTGCACCGACTGTAAGCCACGTCAATCACCCGGAAGAGCACATAGCT";
        let y = b"NNNNNNTAGCTCTGTAGGCACCATCAAT";
        let references = get_as_references(y);

        let result = run_alignment(
            TrimEnd::Right, get_benchmark_config(), None,
            x, &references
        );

        assert!(result.is_some(), "Alignment should be found");
        let alignment = result.unwrap();
        assert_eq!(alignment.get_start_index(), 60);
        assert_eq!(alignment.get_end_index(), 71);

        let formatted = alignment.format(x, y);
        println!("{}", formatted);

        let umi_tags = alignment.get_umi_tags(x, y);
        assert_eq!(umi_tags[0], "AGCACA");


    }

    #[test]
    fn invalid_alignment_at_end() {
        let x = b"TATTCGCTTTGGGCGCGCGCCGCGGCTGGACGAGGCGCCGCCGGCATCATAGCTCTGTAGGACCATCAATA";
        let y = b"NNNNNNTAGCTCTGTAGGCACCATCAAT";
        let references = get_as_references(y);

        let result = run_alignment(
            TrimEnd::Right, get_benchmark_config(), None,
            x, &references
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_shrinking_right() {
        // Two occurrences of AGCTCTGTAGGCACCATCAA, we want the leftmost one (shortest remaining seq)
        // according to the plan which says "focus on remaining sequence, left of the reference"
        // and "repeat that until the shrunken sequence no longer contains the reference".
        let x = b"AGCTCTGTAGGCACCATCAA_AGCTCTGTAGGCACCATCAA_TRAILING";
        let y = b"AGCTCTGTAGGCACCATCAA";
        let references = get_as_references(y);
        let mut config = get_default_config();
        config.iterate = true;
        
        let result = run_alignment(
            TrimEnd::Right, config, None,
            x, &references
        ).unwrap();
        
        // The first match is at 0, the second at 21.
        // If it starts with the one at 21, it will shrink to the left and find the one at 0.
        assert_eq!(result.get_start_index(), 0);
    }

    #[test]
    fn test_shrinking_left() {
        // Two occurrences, we want the rightmost one (shortest remaining suffix)
        let x = b"LEADING_AGCTCTGTAGGCACCATCAA_AGCTCTGTAGGCACCATCAA";
        let y = b"AGCTCTGTAGGCACCATCAA";
        let references = get_as_references(y);
        let mut config = get_default_config();
        config.iterate = true;

        let result = run_alignment(
            TrimEnd::Left, config, None,
            x, &references
        ).unwrap();

        // The first match is at 8, the second at 29.
        // If it starts at 8, it will shrink to the right and find the one at 29.
        assert_eq!(result.get_end_index(), 29 + 20);
    }

    #[test]
    fn test_shrinking_any() {
        let x = b"AGCTCTGTAGGCACCATCAA_MIDDLE_AGCTCTGTAGGCACCATCAA";
        let y = b"AGCTCTGTAGGCACCATCAA";
        let references = get_as_references(y);
        let mut config = get_default_config();
        config.iterate = true;

        let result = run_alignment(
            TrimEnd::Any, config.clone(), None,
            x, &references
        ).unwrap();

        // RIGHT shrinking finds the one at 0 (remaining 0)
        // LEFT shrinking finds the one at 28 (remaining 48 - 48 = 0)
        // Instructed to pick the one with longer remaining sequence.
        // Here both give 0.
        assert!(result.get_start_index() == 0 || result.get_end_index() == 48);
    }

    #[test]
    fn test_shrinking_rtail() {
        let x = b"AGCTCTGTAGGCACCATCAA_MIDDLE_AGCTCTGTAGGCACCATCAA";
        let y = b"AGCTCTGTAGGCACCATCAA"; 
        let references = get_as_references(y);
        let mut config = get_default_config();
        config.iterate = true;
        config.min_overlap = Some(5);

        let setup = MultiAlignerSetup::new("test".to_string(), TrimEnd::RTail, &references, config.clone(), None);
        let aligner = setup.create_aligner();
        let mut state = setup.create_state(&aligner);
        
        // Setup main aligners
        for (a, s) in aligner.aligners.iter().zip(state.aligner_states.iter_mut()) {
            a.setup_aligner(s, x.len());
        }
        
        let initial_best = aligner.find_best_match(&mut state, &setup, x).unwrap();
        
        // Setup shrink aligner
        if let Some(shrink_state) = state.shrink_state.as_mut() {
            aligner.shrink_aligner.as_ref().unwrap().setup_aligner(shrink_state, x.len());
        }
        
        let shrunken = aligner.shrink(&mut state, &setup, x, initial_best.clone()).unwrap();
        
        // Shrunken result should at least be as good as initial.
        assert!(shrunken.score() >= initial_best.score());
    }

    #[test]
    fn test_shrinking_n_handling() {
        // Reference has leading Ns
        let x = b"MIDDLE_AGCTCTGTAGGCACCATCAA_AGCTCTGTAGGCACCATCAA";
        let y = b"NNNNNNAGCTCTGTAGGCACCATCAA"; // 6 Ns
        let references = get_as_references(y);
        let mut config = get_default_config();
        config.iterate = true;

        let result = run_alignment(
            TrimEnd::Right, config, None,
            x, &references
        ).unwrap();

        // First occurrence is at 7 (AGCTCT...). get_start_index() = 7 - 6 = 1.
        // Second is at 28. get_start_index() = 28 - 6 = 22.
        // Shrink right should find the one at 1.
        assert_eq!(result.get_start_index(), 1);
    }
}