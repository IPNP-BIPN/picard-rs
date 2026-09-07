//! `java.util.HashMap`'s iteration order, for the reports that are written in it.
//!
//! A Java `HashMap` iterates its buckets in index order and each bucket in insertion order, so the
//! order it hands out is deterministic, reproducible, and nothing like insertion order overall. It
//! reaches an output whenever a tool loops over a map to write: `SamFileValidator` keeps the reads
//! still waiting for a mate in a `HashMap<String, PairEndInfo>` and reports the leftovers by
//! iterating it, so `ValidateSamFile`'s `MATE_NOT_FOUND` lines come out in bucket order. The names
//! are the same either way; the ORDER is the file, and a report in a different order is a different
//! file.
//!
//! What is reproduced here is the whole of what decides that order:
//!
//! * `String.hashCode`, `h = 31 * h + c` over UTF-16 code units, wrapping at 32 bits;
//! * `HashMap.hash`, which spreads the high bits down: `h ^ (h >>> 16)`;
//! * the table: 16 buckets to start, doubling when the size passes three quarters of the capacity,
//!   and NEVER shrinking, so a map that grew and emptied keeps the width it reached;
//! * `resize`'s split, which preserves each bucket's relative order: an entry stays at `j` or moves
//!   to `j + oldCapacity` depending on one bit, and the two lists keep the order they had.
//!
//! Not reproduced: treeification. A bucket of eight or more entries in a table of at least 64
//! becomes a red-black tree and iterates in tree order instead. Reaching it takes eight names
//! colliding in one bucket, which no corpus here comes close to, and a map that treeifies would
//! iterate differently -- so this is a limit of the emulation and not an approximation of it.

/// `java.lang.String.hashCode`, over UTF-16 code units.
pub fn string_hash_code(s: &str) -> i32 {
    let mut h: i32 = 0;
    for unit in s.encode_utf16() {
        h = h.wrapping_mul(31).wrapping_add(unit as i32);
    }
    h
}

/// `HashMap.hash`: the key's hash with its high half folded down.
fn spread(key: &str) -> u32 {
    let h = string_hash_code(key) as u32;
    h ^ (h >> 16)
}

/// A `java.util.HashMap<String, V>`, kept only for the order it iterates in.
pub struct JavaHashMap<V> {
    table: Vec<Vec<(String, V)>>,
    size: usize,
}

impl<V> Default for JavaHashMap<V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<V> JavaHashMap<V> {
    /// An empty map. Java allocates the table lazily at the first put, with 16 buckets.
    pub fn new() -> Self {
        JavaHashMap {
            table: Vec::new(),
            size: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    fn index(&self, key: &str) -> usize {
        (self.table.len() as u32 - 1) as usize & spread(key) as usize
    }

    /// `HashMap.put`: replace in place when the key is there, else append to the bucket's tail and
    /// grow the table if the size has passed the threshold.
    pub fn put(&mut self, key: &str, value: V) {
        if self.table.is_empty() {
            self.table = (0..16).map(|_| Vec::new()).collect();
        }
        let index = self.index(key);
        if let Some(slot) = self.table[index].iter_mut().find(|(k, _)| k == key) {
            slot.1 = value;
            return;
        }
        self.table[index].push((key.to_string(), value));
        self.size += 1;
        if self.size > self.table.len() * 3 / 4 {
            self.resize();
        }
    }

    /// `HashMap.remove`. The table keeps its width: Java never shrinks one.
    pub fn remove(&mut self, key: &str) -> Option<V> {
        if self.table.is_empty() {
            return None;
        }
        let index = self.index(key);
        let position = self.table[index].iter().position(|(k, _)| k == key)?;
        self.size -= 1;
        Some(self.table[index].remove(position).1)
    }

    pub fn contains_key(&self, key: &str) -> bool {
        !self.table.is_empty() && self.table[self.index(key)].iter().any(|(k, _)| k == key)
    }

    /// `resize`: double the table and split each bucket in two, preserving relative order.
    fn resize(&mut self) {
        let old_capacity = self.table.len();
        let mut grown: Vec<Vec<(String, V)>> = (0..old_capacity * 2).map(|_| Vec::new()).collect();
        for (j, bucket) in std::mem::take(&mut self.table).into_iter().enumerate() {
            for (key, value) in bucket {
                // `(e.hash & oldCap) == 0` keeps the entry where it was; otherwise it moves by
                // exactly the old capacity.
                let high = spread(&key) as usize & old_capacity != 0;
                grown[if high { j + old_capacity } else { j }].push((key, value));
            }
        }
        self.table = grown;
    }

    /// The map's entries in ITERATION order: buckets by index, each bucket in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &V)> {
        self.table
            .iter()
            .flat_map(|bucket| bucket.iter().map(|(k, v)| (k.as_str(), v)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The hashes are Java's. `"read0000".hashCode()` and friends were taken from a JVM.
    #[test]
    fn string_hash_code_matches_java() {
        assert_eq!(string_hash_code(""), 0);
        assert_eq!(string_hash_code("a"), 97);
        assert_eq!(string_hash_code("Aa"), 2112);
        assert_eq!(string_hash_code("BB"), 2112); // the classic collision
        assert_eq!(string_hash_code("hello"), 99162322);
    }

    /// Sixteen keys in a table of sixteen, iterated in bucket order rather than insertion order.
    #[test]
    fn iteration_is_bucket_order_not_insertion_order() {
        let mut map = JavaHashMap::new();
        for i in 0..8 {
            map.put(&format!("key{i}"), i);
        }
        let seen: Vec<&str> = map.iter().map(|(k, _)| k).collect();
        let inserted: Vec<String> = (0..8).map(|i| format!("key{i}")).collect();
        assert_eq!(seen.len(), inserted.len());
        assert_ne!(
            seen,
            inserted.iter().map(String::as_str).collect::<Vec<_>>(),
            "these keys do not land in insertion order, which is the point"
        );
    }

    #[test]
    fn a_removed_key_is_gone_and_the_table_keeps_its_width() {
        let mut map = JavaHashMap::new();
        for i in 0..20 {
            map.put(&format!("k{i}"), i);
        }
        assert_eq!(map.len(), 20);
        assert_eq!(map.remove("k7"), Some(7));
        assert_eq!(map.remove("k7"), None);
        assert_eq!(map.len(), 19);
        assert!(!map.contains_key("k7"));
        assert!(map.contains_key("k8"));
    }
}
