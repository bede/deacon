/// A vector of values that deduplicates every time it doubles in size.
///
/// This way, the space overhead of duplicate elements is at most a factor 2 compared to the unique elements.
pub struct DeduppingVec<T: Ord + Clone> {
    data: Vec<T>,
    limit: usize,
}

impl<T: Ord + Clone> Default for DeduppingVec<T> {
    fn default() -> Self {
        DeduppingVec {
            data: Vec::new(),
            // Only start dedupping once shards hit 1M elements, or ~8MB, for 8GB total memory.
            limit: 1 << 20,
        }
    }
}

impl<T: Ord + Clone> DeduppingVec<T> {
    /// Hand a full staging buffer to its shard, compacting the shard if it has grown
    #[inline]
    pub fn extend(&mut self, values: impl IntoIterator<Item = T>) {
        self.data.extend(values);

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
