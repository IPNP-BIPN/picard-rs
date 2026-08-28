//! `MakeVcfSampleNameMap`: a TSV from VCF path to sample name, one line per input.
//!
//! Reading the VCF headers is not ported. What is ported is what the tool does with the names once
//! it has them, which is the whole of the output: the map, its order, and the line it writes.
//!
//! Ported from `picard.vcf.MakeVcfSampleNameMap` in Picard 3.4.0.

/// `doWork`, on an input whose header does not name exactly one sample.
pub fn wrong_sample_count_message(path: &str, count: usize) -> String {
    format!(
        "Input: {path} was expected to contain a single sample but actually contained {count} samples."
    )
}

/// What one input contributes: the path as GIVEN and the one sample its header names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry {
    pub path: String,
    pub sample: String,
}

/// `String.hashCode`: `s[0]*31^(n-1) + ... + s[n-1]`, over the UTF-16 code units, wrapping.
pub fn java_string_hash(text: &str) -> i32 {
    let mut hash: i32 = 0;
    for unit in text.encode_utf16() {
        hash = hash.wrapping_mul(31).wrapping_add(i32::from(unit));
    }
    hash
}

/// `HashMap.hash`: the key's hash exclusive-ored with its own top half, which is what spreads the
/// high bits into the low ones a small table looks at.
pub fn hash_spread(hash: i32) -> u32 {
    let hash = hash as u32;
    hash ^ (hash >> 16)
}

/// `HashMap.tableSizeFor`: the least power of two at or above the given capacity.
pub fn table_size_for(capacity: usize) -> usize {
    let mut size = 1usize;
    while size < capacity {
        size <<= 1;
    }
    size
}

/// The order a `java.util.HashMap` built by `new HashMap<>(initialCapacity)` iterates its keys in.
///
/// The table is walked bucket by bucket and each bucket in insertion order, so the order is
/// neither the arguments' nor sorted: it is the spread hashes modulo the table's capacity. The
/// capacity starts at `tableSizeFor(initialCapacity)` and doubles whenever the size passes three
/// quarters of it, which is why two keys can be ordered one way and three another.
pub fn hash_map_order(keys: &[String], initial_capacity: usize) -> Vec<String> {
    let mut capacity = table_size_for(initial_capacity.max(1));
    let mut threshold = (capacity as f64 * 0.75) as usize;
    let mut table: Vec<Vec<String>> = vec![Vec::new(); capacity];
    let mut size = 0usize;
    for key in keys {
        let bucket = (capacity - 1) & hash_spread(java_string_hash(key)) as usize;
        if table[bucket].iter().any(|held| held == key) {
            continue;
        }
        table[bucket].push(key.clone());
        size += 1;
        if size > threshold {
            capacity *= 2;
            threshold = (capacity as f64 * 0.75) as usize;
            let mut resized: Vec<Vec<String>> = vec![Vec::new(); capacity];
            for held in table.into_iter().flatten() {
                let bucket = (capacity - 1) & hash_spread(java_string_hash(&held)) as usize;
                resized[bucket].push(held);
            }
            table = resized;
        }
    }
    table.into_iter().flatten().collect()
}

/// The map the tool builds, keyed by the path STRING, in the order it iterates.
///
/// The key is the argument as written and never a normalised or absolute path, so the same file
/// named `a.vcf` and `./a.vcf` is two entries while the same string twice is one. A sample name
/// that turns up twice is only warned about: both entries are kept.
///
/// The map is sized from the INPUT count, so three inputs start in a table of four and two start
/// in one of two and are resized to four by the second: the order is the table's and not the
/// arguments', and naming the same three inputs backwards changes nothing about it.
pub fn build(entries: &[Entry]) -> Vec<(String, String)> {
    let paths: Vec<String> = entries.iter().map(|entry| entry.path.clone()).collect();
    hash_map_order(&paths, entries.len())
        .into_iter()
        .map(|path| {
            let sample = entries
                .iter()
                .find(|entry| entry.path == path)
                .expect("the path was one of the entries")
                .sample
                .clone();
            (path, sample)
        })
        .collect()
}

/// One line of the file: the PATH first and the name second, which is the opposite way round from
/// the tool's own name.
pub fn line(path: &str, sample: &str) -> String {
    format!("{path}\t{sample}")
}

/// The whole file. `Files.write` writes a newline after EVERY line, so the file always ends on
/// one, and an empty map would leave a file of no bytes at all.
pub fn render(lines: &[(String, String)]) -> String {
    lines
        .iter()
        .map(|(path, sample)| line(path, sample) + "\n")
        .collect()
}
