
mod generic;
mod bio_rust;
pub mod result;

use std::cmp::min;
use crate::types::TrimEnd;
use crate::context::SequenceWithId;
use crate::multi_aligner::bio_rust::aligner::BioRustPairwiseAligner;
use crate::multi_aligner::bio_rust::aligner_banded::{BandedAlignerData, BandedExtraConfig, BioRustBandedPairwiseAligner, BioRustBandedRawAligner};
use crate::multi_aligner::generic::AlignerMode;
use crate::multi_aligner::result::MultiAlignerResult;
use log::{debug, info, warn};
use bio_rust::generic_aligner::BioRustAlignerParameters;
use crate::multi_aligner::bio_rust::generic_aligner::{BioRustAlignerConfig, BioRustRawAligner};


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

pub struct MultiAligner {
    aligners: Vec<SpecificAligner>,
}

#[derive(Clone)]
pub struct BioRustPairwiseAlignerSetup {
    config: BioRustAlignerConfig<()>,
}

#[derive(Clone)]
pub struct BioRustBandedPairwiseAlignerSetup {
    config: BioRustAlignerConfig<BandedExtraConfig>,
    aux_data: BandedAlignerData
}

#[derive(Clone)]
pub enum SpecificAlignerSetup {
    BioRust(BioRustPairwiseAlignerSetup),
    BioRustBanded(BioRustBandedPairwiseAlignerSetup),
}

pub type AlignerConfig = BioRustAlignerParameters;

#[derive(Clone)]
pub struct MultiAlignerSetup {
    references: Vec<SequenceWithId>,
    config: AlignerConfig,
    kmer_k: usize,
    aligner_setups: Vec<SpecificAlignerSetup>
}

impl MultiAlignerSetup {
    pub fn new(
        name: String,
        trim_end: TrimEnd,
        references: Vec<SequenceWithId>,
        config: BioRustAlignerParameters,
        k_override: Option<usize>,
    ) -> Self {
        let kmer_k = determine_kmer_k(
            &name, &references, config.error_rate, config.min_overlap, k_override
        );
        if kmer_k > 0 {
            debug!("{}: Using banded aligner with k={} with trim end {:?}.", name, kmer_k, trim_end);
        } else {
            debug!("{}: Using non-banded aligner with trim end {:?}.", name, trim_end);
        }

        let required_modes = match trim_end {
            TrimEnd::Any => vec![
                AlignerMode::Local, AlignerMode::LeftOverlap, AlignerMode::RightOverlap
            ],
            TrimEnd::Left => vec![AlignerMode::Local, AlignerMode::LeftOverlap],
            TrimEnd::Right => vec![AlignerMode::Local, AlignerMode::RightOverlap],
            TrimEnd::LTail => vec![AlignerMode::LeftOverlap],
            TrimEnd::RTail => vec![AlignerMode::RightOverlap],
        };

        let mut aligner_setups:Vec<SpecificAlignerSetup> = Vec::new();
        for mode in required_modes {
            if kmer_k > 0 {
                debug!(
                    "{}: Using banded aligner for mode {:?}, k={}.", name, mode, kmer_k
                );
                let aligner_config = get_aligner_config(
                    &config, mode, BandedExtraConfig { k: kmer_k }
                );
                let aux_data = BioRustBandedRawAligner::create_aux_data(
                    &aligner_config, &references
                );
                aligner_setups.push(
                    SpecificAlignerSetup::BioRustBanded(
                        BioRustBandedPairwiseAlignerSetup { config: aligner_config, aux_data }
                    )
                )
            } else {
                debug!("{}: Using standard aligner for mode {:?}.", name, mode);
                let aligner_config = get_aligner_config(&config, mode, ());
                aligner_setups.push(SpecificAlignerSetup::BioRust(
                    BioRustPairwiseAlignerSetup { config: aligner_config }
                ));
            }

        }
        MultiAlignerSetup { references, config, kmer_k, aligner_setups }
    }

    pub fn create_aligner(&self) -> MultiAligner {
        let aligners: Vec<SpecificAligner> = self.aligner_setups.iter()
        .map(|setup| match setup {
                SpecificAlignerSetup::BioRust(x) 
                    => SpecificAligner::BioRust(
                        BioRustPairwiseAligner::new(&x.config, self.references.clone())
                ),
                SpecificAlignerSetup::BioRustBanded(x) 
                    => SpecificAligner::BioRustBanded(
                    BioRustBandedPairwiseAligner::new(&x.config, self.references.clone())
                ),
        })
        .collect();
        MultiAligner { aligners }
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
        extra,
    }
}


impl MultiAligner {
    pub fn align(&mut self, setup: &MultiAlignerSetup, sequence: &[u8]) -> Option<MultiAlignerResult> {
        if sequence.len() == 0 { return None }  // The aligner seems to go berserk on sequences of size 0
        let mut best: Option<MultiAlignerResult> = None;
        for index in 0..setup.references.len() {
            for (some_setup, some_aligner) in 
                setup.aligner_setups.iter().zip(self.aligners.iter_mut()
            ) {
                let result = match (some_aligner, some_setup) {
                    (SpecificAligner::BioRust(aligner), SpecificAlignerSetup::BioRust(_)) => {
                        aligner.align(sequence, index, &())
                    },
                    (SpecificAligner::BioRustBanded(aligner), SpecificAlignerSetup::BioRustBanded(a_setup)) => {
                        aligner.align(sequence, index, &a_setup.aux_data)
                    }
                    _ => panic!("Unexpected combination of aligner and setup"),
                };
                let reference = &setup.references[index].sequence.as_bytes();
                if !result.validate(
                    setup.config.error_rate, sequence, reference, setup.config.min_overlap, result.offset
                ) {
                    continue;
                }
                let old_score = match best.as_ref() {
                    Some(b) => b.score(),
                    None => i32::MIN,
                };
                if old_score >= result.score() { continue; }
                best = Some(result);
                break;
            }
        }
        best
    }


}

fn determine_kmer_k(
    name: &str,
    references: &Vec<SequenceWithId>,
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
    use super::{MultiAlignerSetup};
    use super::result::MultiAlignerResult;

    fn run_alignment(
        trim_end: TrimEnd,
        config: BioRustAlignerParameters,
        k_overwrite: Option<usize>,
        seq_x: &[u8],
        references: Vec<SequenceWithId>
    ) -> Option<MultiAlignerResult> {
        let setup = MultiAlignerSetup::new(
            "Test Aligner".to_string(),
            trim_end,
            references,
            config,
            k_overwrite,
        );
        let mut aligner = setup.create_aligner();
        aligner.align(&setup, seq_x)
    }

    fn get_as_references(raw: &[u8]) -> Vec<SequenceWithId> {
       vec![SequenceWithId {
            id: "ref".to_string(),
            sequence: String::from_utf8_lossy(raw).to_string(),
       }]
    }

    fn get_default_config() -> BioRustAlignerParameters {
        BioRustAlignerParameters {
            error_rate: 0.0,
            min_overlap: Some(3),
            match_score: 1,
            mismatch_score: -1,
            gap_open: -3,
            gap_extend: -1,
        }
    }

    fn get_overlap5_config() -> BioRustAlignerParameters {
        let mut config = get_default_config();
        config.min_overlap = Some(5);
        config
    }

    fn get_with_tolerance_config() -> BioRustAlignerParameters {
        let mut config = get_default_config();
        config.error_rate = 0.2;
        config
    }

    #[test]
    fn test_aligner_trivial() {
        let references = get_as_references(b"ACGT");
        let result = run_alignment(
            TrimEnd::Right, get_default_config(), None,
            b"ACGT", references
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
            TrimEnd::Right, get_default_config(), None,
            b"AAAAATTT", references
        );
        let alignment = result.unwrap();
        // Matching prefix NNNTTT at position 2
        assert_eq!(alignment.score(), 6);
        assert_eq!(alignment.get_start_index(), 2);
        assert_eq!(alignment.get_end_index(), 8);
    }

    #[test]
    fn test_aligner_right_overhang_ambiguous() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Right, get_default_config(), None,
            b"AAAAATTT", references
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
            TrimEnd::Right, get_overlap5_config(), None,
            b"AAAAATTATT", references
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
            b"AAAAATTTTTAAA", references
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
            TrimEnd::Right, get_overlap5_config(), None,
            b"AAAAATTATTAAA", references
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
            b"AAAAATTATTAAA", references
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
            b"TTTAAAA", references
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
            b"AAAAATTTTTAAATTT", references
        );
        let alignment = result.unwrap();
        // Matching NNNTTTTTT (prefix) at position 2
        assert_eq!(alignment.score(), 6);
        assert_eq!(alignment.get_start_index(), 10);
        assert_eq!(alignment.get_end_index(), 16);
    }

    #[test]
    fn test_aligner_left_overhang() {
        let references = get_as_references(b"TTTT");
        let result = run_alignment(
            TrimEnd::Left, get_default_config(), None,
            b"TTTAAAAA", references
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
            b"AAAAATTTTTAAATTT", references
        );
        let alignment = result.unwrap();
        // Match NNNTTT starting x[10]
        assert_eq!(alignment.score(), 6);
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
            references,
            get_overlap5_config(),
            Some(0)  // Override k to 0 should force the old aligner
        );
        assert_eq!(setup.kmer_k, 0);
    }

    #[test]
    fn test_aligner_k_override_banded() {
        let references = get_as_references(b"ACGT");
        // Default estimate_k for length 4 is 4.
        
        let setup = MultiAlignerSetup::new(
            "Test Aligner".to_string(),
            TrimEnd::Right,
            references,
            get_default_config(),
            Some(5), // Override k to 5 should use banded aligner even if k=4 was estimated
        );
        assert_eq!(setup.kmer_k, 5);
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
            x, references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);
        assert!(formatted.contains("ACGTATCG"));
        assert!(formatted.contains("||||.|||"));
        assert!(formatted.contains("ACGTNTCG"));
    }

    #[test]
    fn test_format_alignment_embedded() {
        let x = b"AAAAACCCCTTTTT";
        let y = b"CCCC";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, references
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
            x, references
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
            x, references
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
            TrimEnd::Any, get_default_config(), None,
            x, references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert!(formatted.contains("CCCCCGGG"));
        assert!(formatted.contains("|\\|"));
        assert!(formatted.contains("AAACCACC"));
    }

    #[test]
    fn test_format_alignment_partial_with_deletion() {
        let x = b"CAATCCTTTT";
        let y = b"CCCCAACCT";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert!(formatted.contains("CAATCCTTT"));
        assert!(formatted.contains("CCCCAA CCT"));
    }

    #[test]
    fn test_format_alignment_partial_with_insertion() {
        let x = b"CAACCTTTT";
        let y = b"CCCCAATCCT";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, references
        );
        let alignment = result.unwrap();
        let formatted = alignment.format(x, y);

        assert!(formatted.contains("CAA CCTTTT"));
        assert!(formatted.contains("CCCCAATCCT"));
    }

    #[test]
    fn test_umi_tag() {
        let x = b"CAATCCTTT";
        let y = b"NNNTTT";
        let references = get_as_references(y);
        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, references
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
            x, references
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
            "  ´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´......||||||||||||||||||||| \n",
            "                                              NNNNNNTAGCTCTGTAGGCACCATCAAT\n"
        ));
    }

    #[test]
    fn invalid_alignment_at_end() {
        let x = b"TATTCGCTTTGGGCGCGCGCCGCGGCTGGACGAGGCGCCGCCGGCATCATAGCTCTGTAGGACCATCAATA";
        let y = b"NNNNNNTAGCTCTGTAGGCACCATCAAT";
        let references = get_as_references(y);

        let result = run_alignment(
            TrimEnd::Any, get_default_config(), None,
            x, references
        );
        assert!(result.is_none());
    }
}