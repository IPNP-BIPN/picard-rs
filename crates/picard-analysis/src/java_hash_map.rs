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
//! Not reproduced: what a long bucket does. When an insertion makes one bucket long enough
//! (`TREEIFY_THRESHOLD`, which OpenJDK 17's `put` reaches on a bucket's ninth entry and its
//! `computeIfAbsent` on the eighth), `HashMap.treeifyBin` runs. In a table of fewer than 64 buckets
//! (`MIN_TREEIFY_CAPACITY`) it doubles the table early; from 64 up it turns the bucket into a tree
//! and moves the tree's root to the head of the bucket's list, which changes the order the bucket
//! iterates in. Reaching either takes that many names colliding in one bucket, which no corpus here
//! comes close to, and a map that did would iterate differently, so this is a limit of the
//! emulation and not an approximation of it.
//!
//! The JDK is GPL2, so nothing here transcribes its source: the behaviour above is named, not
//! quoted.
//!
//! # Why not the obvious Rust
//!
//! Rust's `std::collections::HashMap` seeds its hasher randomly per process, so its iteration order
//! is not even stable between two runs of the same binary, let alone equal to Java's. An
//! insertion-ordered or sorted map would be stable but would still be a different order. The only
//! way to get Java's order is to rebuild the table Java builds, which is what this file does.
//!
//! # Reading this file without reading Rust
//!
//! - `&str` is a borrowed, read-only string; Rust strings are UTF-8 bytes, not UTF-16 like Java's,
//!   which is why the hash re-encodes before it hashes.
//! - `wrapping_mul` and `wrapping_add` are Java's `int` `*` and `+`: they wrap on overflow. Plain
//!   `*` and `+` in Rust panic on overflow in a debug build, so Java's arithmetic has to be asked
//!   for by name.
//! - `x as u32` reinterprets the 32 bits of an `i32` as unsigned; `>>` on a `u32` shifts in zeros,
//!   which is Java's `>>>`.
//! - `JavaHashMap<V>` is generic over the value type, like Java's `HashMap<String, V>`.
//! - `Option<V>` is "a value or nothing", Java's nullable return made explicit; `?` after an
//!   `Option` returns `None` from the function at once when there is nothing.
//! - A body's last expression, written without `;`, is its return value.

/// `java.lang.String.hashCode`, over UTF-16 code units.
///
/// `encode_utf16` yields the same 16-bit `char` values a Java `String` holds, so a name outside
/// ASCII hashes as it does in Java. Hashing the UTF-8 bytes of the `&str` would agree only on
/// ASCII. Also called on its own by `fingerprint.rs` and `crosscheck_run.rs`.
pub fn string_hash_code(s: &str) -> i32 {
    let mut h: i32 = 0;
    for unit in s.encode_utf16() {
        h = h.wrapping_mul(31).wrapping_add(unit as i32);
    }
    h
}

/// `HashMap.hash`: the key's hash with its high half folded down.
///
/// Computed as unsigned so that `>>` is Java's `>>>`; the bits are the same as Java's `int`.
fn spread(key: &str) -> u32 {
    let h = string_hash_code(key) as u32;
    h ^ (h >> 16)
}

/// A `java.util.HashMap<String, V>`, kept only for the order it iterates in.
///
/// `table` is Java's bucket array: one `Vec` per bucket, holding that bucket's entries in the order
/// of Java's linked list. An empty `table` is Java's not-yet-allocated `table == null`. Lookups
/// walk a bucket comparing keys, as Java walks a chain. `size` is Java's `size`, the entry count
/// the growth threshold is tested against.
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

    /// The bucket a key falls in: `(n - 1) & hash`, which for a power-of-two width `n` keeps the
    /// hash's low bits. Only called once the table is allocated.
    fn index(&self, key: &str) -> usize {
        (self.table.len() as u32 - 1) as usize & spread(key) as usize
    }

    /// `HashMap.put`: replace in place when the key is there, else append to the bucket's tail and
    /// grow the table if the size has passed the threshold.
    ///
    /// The threshold is `capacity * 3 / 4`: Java's load factor is 0.75 and its threshold starts at
    /// 12 for 16 buckets and doubles with the table, so integer arithmetic gives the same number
    /// at every width. A replaced value keeps its entry's place in the bucket. Used where the
    /// reference calls `put`, `ValidateSamFile`'s map of reads waiting for a mate among them.
    pub fn put(&mut self, key: &str, value: V) {
        // The first insertion allocates the sixteen buckets, as Java's first `resize` does.
        // `(0..16).map(..).collect()` builds a list of sixteen empty buckets, like a Java stream.
        if self.table.is_empty() {
            self.table = (0..16).map(|_| Vec::new()).collect();
        }
        let index = self.index(key);
        // `if let Some(slot) = ...` runs the block only when the search found an entry;
        // `iter_mut` hands out writable references, so `slot.1 = value` (the pair's second
        // field) replaces the value inside the bucket.
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

    /// `HashMap.computeIfAbsent` for a key that is not there: the new node goes at the HEAD of its
    /// bucket, not the tail as `put` appends it, and the table is grown BEFORE the insertion, when
    /// the size already exceeds the threshold, not after it as `put` does. The two differ in the
    /// order a bucket iterates in, which is the order a report written from the map comes out in.
    ///
    /// A key already present is left as it is, as `computeIfAbsent` returns the existing value
    /// without calling the function. Used where the reference fills its map with
    /// `computeIfAbsent`: the per-sample details of AccumulateVariantCallingMetrics and the UMI
    /// and library maps of UmiAwareMarkDuplicatesWithMateCigar.
    pub fn insert_front_if_absent(&mut self, key: &str, value: V) {
        if self.table.is_empty() {
            self.table = (0..16).map(|_| Vec::new()).collect();
        } else if self.size > self.table.len() * 3 / 4 {
            self.resize();
        }
        let index = self.index(key);
        if self.table[index].iter().any(|(k, _)| k == key) {
            return;
        }
        // `insert(0, ..)` puts the entry first and shifts the rest along, keeping their order.
        self.table[index].insert(0, (key.to_string(), value));
        self.size += 1;
    }

    /// `HashMap.remove`. The table keeps its width: Java never shrinks one.
    ///
    /// The entry is unlinked from its bucket and the others keep their order, as Java's unlinking
    /// of a chain node does. `None` when the key is absent, Java's `null`.
    pub fn remove(&mut self, key: &str) -> Option<V> {
        if self.table.is_empty() {
            return None;
        }
        let index = self.index(key);
        // `position` is the index of the first matching entry, or `None`; the `?` returns `None`
        // from `remove` when the key is not in its bucket.
        let position = self.table[index].iter().position(|(k, _)| k == key)?;
        self.size -= 1;
        Some(self.table[index].remove(position).1)
    }

    pub fn contains_key(&self, key: &str) -> bool {
        !self.table.is_empty() && self.table[self.index(key)].iter().any(|(k, _)| k == key)
    }

    /// `resize`: double the table and split each bucket in two, preserving relative order.
    ///
    /// Each old bucket `j` is walked in order and each entry appended to `j` or `j + oldCap` of
    /// the new table, so the entries that land together keep the order they had, as they do in
    /// Java's split.
    fn resize(&mut self) {
        let old_capacity = self.table.len();
        let mut grown: Vec<Vec<(String, V)>> = (0..old_capacity * 2).map(|_| Vec::new()).collect();
        // `std::mem::take` moves the old table out and leaves an empty one in its place, so its
        // entries can be moved into `grown` rather than copied; `into_iter` consumes the buckets
        // and `enumerate` pairs each with its index `j`.
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
    ///
    /// The order of Java's `entrySet().iterator()`: the bucket array from index 0, each bucket
    /// head to tail. `impl Iterator` returns a lazy sequence, like a Java `Iterator`, without
    /// naming its type; `flat_map` is a stream's `flatMap`, joining the buckets end to end.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &V)> {
        self.table
            .iter()
            .flat_map(|bucket| bucket.iter().map(|(k, v)| (k.as_str(), v)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The hashes are Java's `String.hashCode`, each of which can be checked by hand from
    /// `h = 31 * h + c`.
    #[test]
    fn string_hash_code_matches_java() {
        assert_eq!(string_hash_code(""), 0);
        assert_eq!(string_hash_code("a"), 97);
        assert_eq!(string_hash_code("Aa"), 2112);
        assert_eq!(string_hash_code("BB"), 2112); // the classic collision
        assert_eq!(string_hash_code("hello"), 99162322);
    }

    /// Eight keys in a table of sixteen, iterated in bucket order rather than insertion order.
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
