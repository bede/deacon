/// A vector of values that deduplicates every time it doubles in size.
///
/// This way, the space overhead of duplicate elements is at most a factor 2 compared to the unique elements.
pub struct DedupVec<T: Ord + Clone> {
    data: Vec<T>,
    limit: usize,
}

impl<T: Ord + Clone> Default for DedupVec<T> {
    fn default() -> Self {
        DedupVec {
            data: Vec::new(),
            // Only start dedupping once shards hit 1M elements, or ~8MB, for 8GB total memory.
            limit: 1 << 22,
        }
    }
}

impl<T: Ord + Clone> DedupVec<T> {
    /// Hand a full staging buffer to its shard, compacting the shard if it has grown
    #[inline]
    pub fn extend(&mut self, values: &mut Vec<T>) {
        self.data.extend_from_slice(values);
        values.clear();

        if self.data.len() >= self.limit {
            self.data.sort_unstable();
            self.data.dedup();
            self.limit = self.limit.max(2 * self.data.len());
        }
    }

    /// Sort and deduplicate every shard in parallel
    pub fn finish(mut self) -> Vec<T> {
        self.data.sort_unstable();
        self.data.dedup();
        self.data
    }
}
