//! A dense_hash_set for u64 keys.
//!
//! Compared to std::collections::HashSet<u64>, this uses a different layout: no metadata table, just plain data.
//! This is similar to Google's dense_hash_map, which predates the SwissTable design. By avoiding a metadata table,
//! we may need to do longer probe sequences (each probe is 8 bytes, not 1 byte), but on the other hand we only take
//! 1 cache miss per access, not 2.

use std::hash::{BuildHasher, BuildHasherDefault};

use wide::CmpEq;

type Hasher = BuildHasherDefault<rustc_hash::FxHasher>;

pub struct U64HashSet {
    table: Box<[Bucket]>,
    len: usize,
    has_zero: bool,
}

const BUCKET_SIZE: usize = 8;

#[derive(Clone, Copy)]
#[repr(align(64))] // Cache line alignment
struct Bucket([u64; BUCKET_SIZE]);

impl U64HashSet {
    pub fn with_capacity(capacity: usize) -> Self {
        // TODO: integer overflow...
        let num_buckets = (capacity.next_power_of_two() * 2).div_ceil(BUCKET_SIZE);
        let table = vec![Bucket([0u64; BUCKET_SIZE]); num_buckets].into_boxed_slice();
        Self {
            table,
            len: 0,
            has_zero: false,
        }
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len + self.has_zero as usize
    }

    pub fn iter(&self) -> impl Iterator<Item = u64>{
        std::iter::repeat_n(0, self.has_zero as usize).chain(
            self.table
                .iter()
                .flat_map(|b| b.0.iter().copied())
                .filter(|x| *x != 0),
        )
    }

    #[inline(always)]
    pub fn prefetch(&self, key: u64) {
        let hash64 = Hasher::default().hash_one(key);
        let bucket_mask = self.table.len() - 1;
        let bucket_i = hash64 as usize;
        // Safety: bucket_mask is correct because the number of buckets is a power of 2.
        unsafe {
            std::intrinsics::prefetch_write_data::<_, 0>(
                self.table.get_unchecked(bucket_i & bucket_mask) as *const Bucket as *const u8,
            )
        };
    }

    #[inline(always)]
    pub fn contains(&self, key: u64) -> bool {
        if key == 0 {
            return self.has_zero;
        }
        let hash64 = Hasher::default().hash_one(key);
        let bucket_mask = self.table.len() - 1;
        let mut bucket_i = hash64 as usize;

        // type S = wide::u64x4;
        type S = wide::i64x4;
        let keys = S::splat(key as i64);

        loop {
            use std::mem::transmute;
            // Safety: bucket_mask is correct because the number of buckets is a power of 2.
            let bucket = unsafe { self.table.get_unchecked(bucket_i & bucket_mask) };
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
    pub fn insert(&mut self, key: u64) {
        if key == 0 {
            self.len += !self.has_zero as usize;
            self.has_zero = true;
            return;
        }
        let hash64 = Hasher::default().hash_one(key);
        let bucket_mask = self.table.len() - 1;
        let element_offset_in_bucket = (hash64 >> 61) as usize;
        let mut bucket_i = hash64 as usize;

        loop {
            // Safety: bucket_mask is correct because the number of buckets is a power of 2.
            let bucket = unsafe { self.table.get_unchecked_mut(bucket_i & bucket_mask) };
            for element_i in 0..BUCKET_SIZE {
                let element = &mut bucket.0[(element_i + element_offset_in_bucket) % BUCKET_SIZE];
                if *element == 0 {
                    *element = key;
                    self.len += 1;
                    return;
                }
                if *element == key {
                    return;
                }
            }
            bucket_i += 1;
        }
    }
}
