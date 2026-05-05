use rust_htslib::bam::HeaderView;

#[derive(Clone, Debug)]
pub enum SequenceHeader {
    Default,
    Bam(HeaderView),
}

impl Default for SequenceHeader {
    fn default() -> Self {
        SequenceHeader::Default
    }
}
