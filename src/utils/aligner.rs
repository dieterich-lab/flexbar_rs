use bio::alignment::Alignment;
use bio::alignment::pairwise::*;
use log::trace;
use crate::types::TrimEnd;


struct Capacity {
    n: usize,
    m: usize,
}

struct CachedAligner {
    aligner: Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>,
    capacity: Capacity,
}

pub struct AlignerFactory {
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    xclip_prefix: i32,
    xclip_suffix: i32,
    yclip_prefix: i32,
    yclip_suffix: i32,
    cache: thread_local::ThreadLocal<std::cell::RefCell<Option<CachedAligner>>>,
}

impl AlignerFactory {
    pub fn new(
        trim_end: TrimEnd,
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Self {
        let xclip_prefix = match trim_end {
            TrimEnd::Right => MIN_SCORE,
            TrimEnd::RTail => MIN_SCORE,
            _ => 0,
        };
        let xclip_suffix = match trim_end {
            TrimEnd::Left => MIN_SCORE,
            TrimEnd::LTail => MIN_SCORE,
            _ => 0,
        };
        let yclip_prefix = 0;
        let yclip_suffix = 0;

        Self {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            xclip_prefix,
            xclip_suffix,
            yclip_prefix,
            yclip_suffix,
            cache: thread_local::ThreadLocal::new(),
        }
    }

    pub fn get(&self, n: usize, m: usize) -> std::cell::RefMut<'_, Aligner<Box<dyn Fn(u8, u8) -> i32 + Send + Sync>>> {
        let cache_cell = self.cache.get_or(|| std::cell::RefCell::new(None));
        let mut cache = cache_cell.borrow_mut();
        let new_cap = self.get_capacity_requirement(&cache, n, m);

        match new_cap {
            Some(cap) => {
                let cached_aligner = self.build_aligner(cap);
                *cache = Some(cached_aligner);
            },
            None => ()
        }

        std::cell::RefMut::map(cache, |opt| {
            &mut opt.as_mut().unwrap().aligner
        })
    }

    fn get_capacity_requirement(&self, cache: &Option<CachedAligner>, n: usize, m: usize)
                                -> Option<Capacity> {
        match cache {
            Some(c) => {
                if n <= c.capacity.n && m <= c.capacity.m {
                    None
                } else {
                    let new_n = if n > c.capacity.n { n + 16 } else { c.capacity.n };
                    let new_m = if m > c.capacity.m { m + 16 } else { c.capacity.m };
                    Some(Capacity { n: new_n, m: new_m }) 
                }
            },
            None => Some(Capacity { n, m })
        }
    }

    fn build_aligner(&self, capacity: Capacity) -> CachedAligner {
        let match_score = self.match_score;
        let mismatch_score = self.mismatch_score;

        let match_fn: Box<dyn Fn(u8, u8) -> i32 + Send + Sync> = Box::new(move |a, b| {
            if a == b || a == b'N' || b == b'N' {
                match_score
            } else {
                mismatch_score
            }
        });

        let scoring = Scoring {
            match_fn,
            match_scores: None,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            xclip_prefix: self.xclip_prefix,
            xclip_suffix: self.xclip_suffix,
            yclip_prefix: self.yclip_prefix,
            yclip_suffix: self.yclip_suffix,
        };

        let aligner = Aligner::with_capacity_and_scoring(capacity.m, capacity.n, scoring);
        CachedAligner { aligner, capacity }
    }
}

pub fn trace_alignment(alignment: &Alignment, x: &[u8], y: &[u8]) {
    trace!(
        "Alignment with score {}:\n{}",
        alignment.score, alignment.pretty(x, y, 80)
    );
}