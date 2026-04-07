use serde::{Deserialize, Serialize};
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum TrimEnd {
    Any,
    Left,
    Right,
    LTail,
    RTail,
}

impl FromStr for TrimEnd {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "ANY" => Ok(TrimEnd::Any),
            "LEFT" => Ok(TrimEnd::Left),
            "RIGHT" => Ok(TrimEnd::Right),
            "LTAIL" => Ok(TrimEnd::LTail),
            "RTAIL" => Ok(TrimEnd::RTail),
            _ => Err(format!("Unknown trim-end: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum FileFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum QualityType {
    Sanger,
    Solexa,
    Illumina,
}

impl FromStr for QualityType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "sanger" | "i1.8" => Ok(QualityType::Sanger),
            "solexa" => Ok(QualityType::Solexa),
            "i1.3" | "i1.5" => Ok(QualityType::Illumina),
            _ => Err(format!("Unknown quality format: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum QualTrimType {
    Off,
    Tail,
    Win,
    Bwa,
}

impl FromStr for QualTrimType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "OFF" => Ok(QualTrimType::Off),
            "TAIL" => Ok(QualTrimType::Tail),
            "WIN" => Ok(QualTrimType::Win),
            "BWA" => Ok(QualTrimType::Bwa),
            _ => Err(format!("Unknown quality trim mode: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum LogAlign {
    None,
    All,
    Mod,
    Tab,
}

impl FromStr for LogAlign {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "NONE" => Ok(LogAlign::None),
            "ALL" => Ok(LogAlign::All),
            "MOD" => Ok(LogAlign::Mod),
            "TAB" => Ok(LogAlign::Tab),
            _ => Err(format!("Unknown log alignment mode: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum CompressionType {
    Uncompressed,
    Gz,
    Bz2,
    Zstd,
    Keep,
}

impl FromStr for CompressionType {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "GZ" => Ok(CompressionType::Gz),
            "BZ2" => Ok(CompressionType::Bz2),
            "ZSTD" => Ok(CompressionType::Zstd),
            "KEEP" => Ok(CompressionType::Keep),
            _ => Ok(CompressionType::Uncompressed),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum RunType {
    Single,
    Paired,
    SingleBarcoded,
    PairedBarcoded,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum BarcodeDetect {
    Off,
    WithinRead,
    WithinReadRemoval,
    WithinRead2,
    WithinReadRemoval2,
    BarcodeRead,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum AdapterRemoval {
    Off,
    Normal,
    Normal2,
    One,
    Two,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum AdapterReadset {
    One,
    Two,
    None
}

impl FromStr for AdapterRemoval {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1" => Ok(AdapterRemoval::One),
            "2" => Ok(AdapterRemoval::Two),
            "Off" => Ok(AdapterRemoval::Off),
            "Normal" => Ok(AdapterRemoval::Normal),
            "Normal2" => Ok(AdapterRemoval::Normal2),
            _ => Err(format!("Unknown adapter removal mode: {}", s)),
        }
    }
}

impl FromStr for AdapterReadset {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1" => Ok(AdapterReadset::One),
            "2" => Ok(AdapterReadset::Two),
            "None" => Ok(AdapterReadset::None),
            _ => Err(format!("Unknown adapter read set: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum RevCompMode {
    Off,
    On,
    Only,
}

impl FromStr for RevCompMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "ON" => Ok(RevCompMode::On),
            "ONLY" => Ok(RevCompMode::Only),
            _ => Ok(RevCompMode::Off),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum PairOverlapMode {
    Off,
    On,
    Short,
    Only,
}

impl FromStr for PairOverlapMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "ON" => Ok(PairOverlapMode::On),
            "SHORT" => Ok(PairOverlapMode::Short),
            "ONLY" => Ok(PairOverlapMode::Only),
            _ => Ok(PairOverlapMode::Off),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum AdapterPreset {
    Off,
    TruSeq,
    SmallRna,
    Methyl,
    Ribo,
    Nextera,
    NexteraMp,
}

impl FromStr for AdapterPreset {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "TruSeq" => Ok(AdapterPreset::TruSeq),
            "SmallRNA" => Ok(AdapterPreset::SmallRna),
            "Methyl" => Ok(AdapterPreset::Methyl),
            "Ribo" => Ok(AdapterPreset::Ribo),
            "Nextera" => Ok(AdapterPreset::Nextera),
            "NexteraMP" => Ok(AdapterPreset::NexteraMp),
            _ => Err(format!("Unknown adapter preset: {}", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum AdapterTrimmed {
    Off,
    On,
    Only,
}

impl FromStr for AdapterTrimmed {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "ON" => Ok(AdapterTrimmed::On),
            "OFF" => Ok(AdapterTrimmed::Off),
            "ONLY" => Ok(AdapterTrimmed::Only),
            _ => Err(format!("Unknown adapter trimmed out mode: {}", s)),
        }
    }
}

pub type SingeRecordDto = bio::io::fastq::Record;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedRecordDto {
    pub reads: [SingeRecordDto; 2],
}

impl PairedRecordDto {
    pub fn new(r1: SingeRecordDto, r2: SingeRecordDto) -> Self {
        Self { reads: [r1, r2] }
    }

    pub fn id(&self) -> &str {
        self.reads[0].id()
    }

    pub fn validate(&self) -> bool {
        let id1 = self.reads[0].id();
        let id2 = self.reads[1].id();
        if id1 == id2 {
            return true;
        }
        if id1.ends_with("/1") && id2.ends_with("/2") {
            return id1[..id1.len() - 2] == id2[..id2.len() - 2]
        }
        false
    }
}

pub enum RecordDto {
    Single(SingeRecordDto),
    Paired(PairedRecordDto),
}

macro_rules! implement_by_record_index {
    ($name: ident, $ReturnType: ty) => {
        #[inline(always)]
        pub fn $name(&self, read_index: Option<usize>) -> $ReturnType {
            match (self, read_index) {
                (RecordDto::Single(record), None) => record.$name(),
                (RecordDto::Paired(paired), Some(0)) => paired.reads[0].$name(),
                (RecordDto::Paired(paired), Some(1)) => paired.reads[1].$name(),
                _ => panic!("Invalid read index {:?}", read_index),
            }
        }
    };
}

impl RecordDto {
    pub fn id(&self, read_index: Option<usize>) -> &str {
        match (self, read_index) {
            (RecordDto::Single(record), None) => record.id(),
            (RecordDto::Paired(paired), Some(0)) => paired.reads[0].id(),
            (RecordDto::Paired(paired), Some(1)) => paired.reads[1].id(),
            // We sometimes need an ID for error messages without knowing if it's a single or
            // paired record. So we allow None here too and just return the first ID.
            (RecordDto::Paired(paired), None) => paired.reads[0].id(),
            _ => panic!("Invalid read index {:?}", read_index),
        }
    }
    implement_by_record_index!(desc, Option<&str>);
    implement_by_record_index!(seq, &[u8]);
    implement_by_record_index!(qual, &[u8]);
}
