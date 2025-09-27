//! A dense_hash_set for u64 keys.
//!
//! Compared to std::collections::HashSet<u64>, this uses a different layout: no metadata table, just plain data.
//! This is similar to Google's dense_hash_map, which predates the SwissTable design. By avoiding a metadata table,
//! we may need to do longer probe sequences (each probe is 8 bytes, not 1 byte), but on the other hand we only take
//! 1 cache miss per access, not 2.

use std::arch::x86_64::{_mm_prefetch, _MM_HINT_T0};
use std::hint::select_unpredictable;
use std::mem::transmute;
type S = wide::i64x4;

use rustc_hash::FxHashMap;
use wide::CmpEq;

fn mul_high(a: u64, b: usize) -> usize {
    ((a as u128) * (b as u128) >> 64) as usize
}

pub struct U64HashSet {
    buckets: usize,
    table: Box<[Bucket]>,
    len: usize,
    has_zero: bool,
    hits: usize,
    skips: usize,
    skips2: usize,
    last_b: usize,
    last_j: usize,
    last_key: u64,
}

const PADDING: usize = 1000;

const BUCKET_SIZE: usize = 8;

#[derive(Clone, Copy, PartialEq, Eq)]
#[repr(align(64))] // Cache line alignment
struct Bucket([u64; BUCKET_SIZE]);

impl U64HashSet {
    pub fn with_capacity(n: usize) -> Self {
        eprintln!("N        {n:>10}");
        let capacity = n * 15 / 10;
        eprintln!("CAPACITY {capacity:>10}");
        let buckets = capacity.div_ceil(BUCKET_SIZE);
        let table = vec![Bucket([0u64; BUCKET_SIZE]); buckets + PADDING].into_boxed_slice();
        Self {
            buckets,
            table,
            len: 0,
            has_zero: false,
            hits: 0,
            skips: 0,
            skips2: 0,
            last_b: 0,
            last_j: 0,
            last_key: 0,
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len + self.has_zero as usize
    }

    pub fn iter(&self) -> impl Iterator<Item = u64> {
        std::iter::repeat_n(0, self.has_zero as usize).chain(
            self.table
                .iter()
                .flat_map(|b| b.0.iter().copied())
                .filter(|x| *x != 0),
        )
    }

    pub fn test(&self) {
        for x in self.iter() {
            if !self.contains(x) {
                eprintln!("Did not find {x}!");
                let hash64 = x;
                let bucket_i = mul_high(hash64, self.buckets);

                let [h1, h2]: &[S; 2] = unsafe { transmute(&self.table[bucket_i]) };
                let c0 = h1.cmp_eq(S::ZERO).move_mask().count_ones() as usize;
                let c1 = h2.cmp_eq(S::ZERO).move_mask().count_ones() as usize;
                let elems = BUCKET_SIZE - c0 - c1;
                eprintln!("Intended bucket {bucket_i} of size {elems}");
                let flat = unsafe { self.table.align_to::<u64>().1 };
                let pos = flat.iter().position(|y| *y == x);
                eprintln!("Found at {pos:?}");
                if let Some(p) = pos {
                    let bucket = p / BUCKET_SIZE;
                    eprintln!("Actual bucket {bucket}");
                }

                panic!();
            }
        }
    }

    pub fn stats(&self) {
        let mut counts = [0; 9];

        eprintln!("Size    : {}", self.len);
        eprintln!("hits    : {}", self.hits);
        eprintln!("Skips   : {}", self.skips);
        eprintln!("Skips/el: {}", self.skips as f32 / self.hits as f32);
        eprintln!("Skips2/el {}", self.skips2 as f32 / self.hits as f32);
        // return;

        let mut sum = 0;
        let mut cnt = 0;
        for bucket in &self.table {
            let [h1, h2]: &[S; 2] = unsafe { transmute(&bucket.0) };
            let c0 = h1.cmp_eq(S::ZERO).move_mask().count_ones() as usize;
            let c1 = h2.cmp_eq(S::ZERO).move_mask().count_ones() as usize;
            let elems = BUCKET_SIZE - c0 - c1;
            counts[elems] += 1;
            cnt += 1;
            sum += elems;
        }
        for i in 0..=8 {
            eprintln!("{i}: {:>9}", counts[i]);
        }
        eprintln!("buckets {cnt}");
        eprintln!("slots   {}", cnt * BUCKET_SIZE);
        eprintln!("sum {sum}");
        eprintln!("avg {}", sum as f32 / cnt as f32);

        self.test();
    }

    #[inline(always)]
    fn bucket_idx(&self, key: u64) -> usize {
        mul_high(key, self.buckets)
    }

    #[inline(always)]
    pub fn prefetch(&self, key: u64) {
        let hash64 = key;
        let bucket_i = self.bucket_idx(hash64);
        // Safety: bucket_mask is correct because the number of buckets is a power of 2.
        unsafe {
            _mm_prefetch::<_MM_HINT_T0>(
                self.table.get_unchecked(bucket_i) as *const Bucket as *const i8
            )
        };
    }

    #[inline(always)]
    pub fn contains(&self, key: u64) -> bool {
        if key == 0 {
            return self.has_zero;
        }
        let hash64 = key;
        let mut bucket_i = self.bucket_idx(hash64);

        // type S = wide::u64x4;
        type S = wide::i64x4;
        let keys = S::splat(key as i64);

        loop {
            use std::mem::transmute;
            // Safety: bucket_mask is correct because the number of buckets is a power of 2.
            let bucket = unsafe { self.table.get_unchecked(bucket_i) };
            let [h1, h2]: &[S; 2] = unsafe { transmute(&bucket.0) };
            let mask = (h1.cmp_eq(keys) | h2.cmp_eq(keys)).move_mask() as u8;
            if mask > 0 {
                return true;
            }
            let has_zero = (h1.cmp_eq(S::ZERO) | h2.cmp_eq(S::ZERO)).move_mask() as u8;
            if has_zero > 0 {
                return false;
            }

            bucket_i += 1;
        }
    }

    #[inline(always)]
    pub fn insert(&mut self, key: u64) -> bool {
        if key == 0 {
            self.len += !self.has_zero as usize;
            let ret = !self.has_zero;
            self.has_zero = true;
            return ret;
        }
        let hash64 = key;
        let mut bucket_i = self.bucket_idx(hash64);
        let keys = S::splat(key as i64);

        let mut i = 0;
        loop {
            // Safety: bucket_mask is correct because the number of buckets is a power of 2.
            let bucket = &mut self.table[bucket_i];
            let [h1, h2]: &[S; 2] = unsafe { transmute(&bucket.0) };

            let mask = (h1.cmp_eq(keys) | h2.cmp_eq(keys)).move_mask() as u8;
            if mask > 0 {
                // Element already exists.
                return false;
            }

            let c0 = h1.cmp_eq(S::ZERO).move_mask().count_ones() as usize;
            let c1 = h2.cmp_eq(S::ZERO).move_mask().count_ones() as usize;
            let taken = BUCKET_SIZE - c0 - c1;

            if taken < BUCKET_SIZE {
                let element_i = taken;
                let element = &mut bucket.0[element_i];
                if *element == 0 {
                    self.hits += 1;
                    *element = key;
                    self.len += 1;
                    self.skips = i;
                    self.skips2 += i * i;
                    return true;
                }
                panic!();
            }

            bucket_i += 1;
            i += 1;
            self.skips += 1;
        }
    }

    #[inline(always)]
    pub fn insert_in_order(&mut self, key: u64) {
        self.len += 1;
        assert!(
            self.len <= self.buckets * BUCKET_SIZE,
            "Inserted too many keys!"
        );
        if key == 0 {
            assert!(self.last_key == 0, "Keys must be inserted in order");
            self.has_zero = true;
            return;
        }
        assert!(key > self.last_key, "Keys must be inserted in order");
        self.last_key = key;
        let hash64 = key;
        let bucket_i = self.bucket_idx(hash64);
        // same bucket?
        self.last_j = select_unpredictable(bucket_i > self.last_b, 0, self.last_j + 1);
        self.last_b = self.last_b.max(bucket_i);
        self.last_b += self.last_j >> 3;
        self.last_j &= 7;
        self.hits += 1;
        self.skips += self.last_b - bucket_i;
        self.skips2 += (self.last_b - bucket_i).pow(2);
        self.table[self.last_b].0[self.last_j] = key;
    }
}
