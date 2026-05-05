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
    Bam,
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

