//! The fingerprinting engine the five fingerprint binaries share: the haplotype map as
//! `HaplotypeMap` builds it, the `HaplotypeProbabilities` family, `Fingerprint`, the pileup
//! `FingerprintChecker.fingerprintSamFile` walks, the VCF loader, and `calculateMatchResults`.
//!
//! [`crate::haplotype_map`], [`crate::check_fingerprint`] and their siblings hold the arithmetic
//! of each tool in isolation. This is the machinery that turns files into those numbers, ported
//! line for line because three things in it leak into the answer:
//!
//! * **Hash order.** The haplotype blocks come out of a `HashMap` keyed by anchor name, a block's
//!   SNPs out of one keyed by SNP name, and every per-read-group fingerprint out of a `HashMap`
//!   keyed by `FingerprintIdDetails`, whose hash includes the file's URI. The rows of the
//!   crosscheck tools are in that last order, so it is reproduced (see [`JavaMap`]) and the URI is
//!   the one the reference saw (see [`reference_view`]).
//! * **Normalisation after every base.** `HaplotypeProbabilitiesUsingLogLikelihoods.setLogLikelihoods`
//!   renormalises the three log-likelihoods each time it is called, so the order bases are added
//!   in is part of the result's last bits.
//! * **One base per read name**, across every locus and every read group, in the order the locus
//!   iterator meets them.
//!
//! Ported from Picard 3.4.0 `picard.fingerprint` (`HaplotypeMap`, `HaplotypeBlock`, `Snp`,
//! `Fingerprint`, `FingerprintChecker`, `FingerprintIdDetails`, the `HaplotypeProbabilities*`
//! classes, `CappedHaplotypeProbabilities`, `MatchResults`, `LocusResult`, `DiploidGenotype`),
//! `picard.util.MathUtil` and `picard.util.AlleleSubsettingUtils`, and htsjdk 4.2.0
//! `SamLocusIterator`/`AbstractLocusIterator`.

use std::collections::BTreeMap;

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

use crate::java_hash_map::string_hash_code;
use crate::theoretical_sensitivity::JavaRandom;

// ---------------------------------------------------------------------------------------------
// Java collections and math.
// ---------------------------------------------------------------------------------------------

/// `HashMap.hash`: the high half folded down.
fn spread(h: i32) -> u32 {
    let h = h as u32;
    h ^ (h >> 16)
}

/// `HashMap.tableSizeFor`.
fn table_size_for(capacity: usize) -> usize {
    capacity.max(1).next_power_of_two()
}

/// A `java.util.HashMap` with arbitrary keys, kept for the order it iterates in.
///
/// Buckets in index order, each in insertion order; 16 buckets to start (or `tableSizeFor` of a
/// requested capacity), doubling past three quarters, with `resize`'s order-preserving split.
#[derive(Debug, Clone)]
pub struct JavaMap<K, V> {
    table: Vec<Vec<(i32, K, V)>>,
    initial: usize,
    size: usize,
}

impl<K: PartialEq + Clone, V> Default for JavaMap<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: PartialEq + Clone, V> JavaMap<K, V> {
    pub fn new() -> Self {
        Self::with_capacity(16)
    }

    /// `new HashMap<>(initialCapacity)`.
    pub fn with_capacity(initial: usize) -> Self {
        JavaMap {
            table: Vec::new(),
            initial: table_size_for(initial),
            size: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.size
    }

    pub fn is_empty(&self) -> bool {
        self.size == 0
    }

    fn index(&self, hash: i32) -> usize {
        (self.table.len() - 1) & spread(hash) as usize
    }

    /// `put`: replaces the value in place when the key is there (the key object is kept).
    pub fn put(&mut self, hash: i32, key: K, value: V) {
        if self.table.is_empty() {
            self.table = (0..self.initial).map(|_| Vec::new()).collect();
        }
        let index = self.index(hash);
        if let Some(slot) = self.table[index]
            .iter_mut()
            .find(|(h, k, _)| *h == hash && *k == key)
        {
            slot.2 = value;
            return;
        }
        self.table[index].push((hash, key, value));
        self.size += 1;
        if self.size > self.table.len() * 3 / 4 {
            let old = self.table.len();
            let mut grown: Vec<Vec<(i32, K, V)>> = (0..old * 2).map(|_| Vec::new()).collect();
            for (j, bucket) in std::mem::take(&mut self.table).into_iter().enumerate() {
                for entry in bucket {
                    let high = spread(entry.0) as usize & old != 0;
                    grown[if high { j + old } else { j }].push(entry);
                }
            }
            self.table = grown;
        }
    }

    pub fn get(&self, hash: i32, key: &K) -> Option<&V> {
        if self.table.is_empty() {
            return None;
        }
        self.table[self.index(hash)]
            .iter()
            .find(|(h, k, _)| *h == hash && k == key)
            .map(|(_, _, v)| v)
    }

    pub fn get_mut(&mut self, hash: i32, key: &K) -> Option<&mut V> {
        if self.table.is_empty() {
            return None;
        }
        let index = self.index(hash);
        self.table[index]
            .iter_mut()
            .find(|(h, k, _)| *h == hash && k == key)
            .map(|(_, _, v)| v)
    }

    pub fn contains_key(&self, hash: i32, key: &K) -> bool {
        self.get(hash, key).is_some()
    }

    pub fn remove(&mut self, hash: i32, key: &K) -> Option<V> {
        if self.table.is_empty() {
            return None;
        }
        let index = self.index(hash);
        let position = self.table[index]
            .iter()
            .position(|(h, k, _)| *h == hash && k == key)?;
        self.size -= 1;
        Some(self.table[index].remove(position).2)
    }

    /// Entries in iteration order.
    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.table
            .iter()
            .flat_map(|bucket| bucket.iter().map(|(_, k, v)| (k, v)))
    }

    pub fn into_entries(self) -> Vec<(i32, K, V)> {
        self.table.into_iter().flatten().collect()
    }
}

/// `java.util.concurrent.ConcurrentHashMap`, for the order `putAll` into one leaves.
///
/// Bins are appended at their tail like a `HashMap`'s, but a resize builds each half of a split
/// bin by PREPENDING every node before the bin's last run, so the order within a bin can reverse.
/// That is only observable for keys that share a bin, which is also the only time it reaches a
/// later `HashMap`'s order.
#[derive(Debug, Clone)]
pub struct JavaConcurrentMap<K, V> {
    table: Vec<Vec<(i32, K, V)>>,
    size_ctl: usize,
    size: usize,
}

impl<K: PartialEq + Clone, V> JavaConcurrentMap<K, V> {
    /// `new ConcurrentHashMap<>(initialCapacity)`.
    pub fn with_capacity(initial: usize) -> Self {
        let cap = table_size_for(initial + (initial >> 1) + 1);
        JavaConcurrentMap {
            table: Vec::new(),
            size_ctl: cap,
            size: 0,
        }
    }

    fn index(len: usize, hash: i32) -> usize {
        // `spread`: `(h ^ (h >>> 16)) & HASH_BITS`.
        (len - 1) & ((spread(hash) & 0x7fff_ffff) as usize)
    }

    /// `putAll`, which first presizes the table for the incoming map (`tryPresize`).
    pub fn put_all(&mut self, entries: Vec<(i32, K, V)>) {
        let size = entries.len();
        let c = table_size_for(size + (size >> 1) + 1);
        loop {
            if self.table.is_empty() {
                let n = self.size_ctl.max(c);
                self.table = (0..n).map(|_| Vec::new()).collect();
                self.size_ctl = n - (n >> 2);
            } else if c <= self.size_ctl || self.table.len() >= (1 << 30) {
                break;
            } else {
                self.transfer();
            }
        }
        for (h, k, v) in entries {
            self.put(h, k, v);
        }
    }

    /// `transfer`: double the table, each bin split with its prefix prepended.
    fn transfer(&mut self) {
        let n = self.table.len();
        let mut grown: Vec<Vec<(i32, K, V)>> = (0..n * 2).map(|_| Vec::new()).collect();
        for (i, bin) in std::mem::take(&mut self.table).into_iter().enumerate() {
            if bin.is_empty() {
                continue;
            }
            let bit = |h: i32| (spread(h) & 0x7fff_ffff) as usize & n;
            let mut last_run = 0;
            let mut run_bit = bit(bin[0].0);
            for (j, node) in bin.iter().enumerate().skip(1) {
                let b = bit(node.0);
                if b != run_bit {
                    run_bit = b;
                    last_run = j;
                }
            }
            let mut nodes = bin;
            let tail = nodes.split_off(last_run);
            let (mut low, mut high) = if run_bit == 0 {
                (tail, Vec::new())
            } else {
                (Vec::new(), tail)
            };
            for node in nodes {
                if bit(node.0) == 0 {
                    low.insert(0, node);
                } else {
                    high.insert(0, node);
                }
            }
            grown[i] = low;
            grown[i + n] = high;
        }
        self.table = grown;
        let len = self.table.len();
        self.size_ctl = len - (len >> 2);
    }

    pub fn put(&mut self, hash: i32, key: K, value: V) {
        if self.table.is_empty() {
            let n = self.size_ctl;
            self.table = (0..n).map(|_| Vec::new()).collect();
            self.size_ctl = n - (n >> 2);
        }
        let index = Self::index(self.table.len(), hash);
        if let Some(slot) = self.table[index]
            .iter_mut()
            .find(|(h, k, _)| *h == hash && *k == key)
        {
            slot.2 = value;
            return;
        }
        self.table[index].push((hash, key, value));
        self.size += 1;
        while self.size >= self.size_ctl && self.table.len() < (1 << 30) {
            self.transfer();
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = (&K, &V)> {
        self.table
            .iter()
            .flat_map(|bucket| bucket.iter().map(|(_, k, v)| (k, v)))
    }

    pub fn into_entries(self) -> Vec<(i32, K, V)> {
        self.table.into_iter().flatten().collect()
    }
}

/// `Integer.hashCode`.
pub fn integer_hash(value: i32) -> i32 {
    value
}

/// `Math.log10`.
pub fn log10(x: f64) -> f64 {
    jmath::math::log10(x)
}

/// `Math.pow`: a HotSpot intrinsic, approximated by the platform's correctly rounded-in-practice
/// `pow`. Every value it produces here is summed, logged and printed to at most six places.
pub fn pow(x: f64, y: f64) -> f64 {
    x.powf(y)
}

/// `MathUtil.MAX_PROB_BELOW_ONE`.
const MAX_PROB_BELOW_ONE: f64 = 0.9999999999999999;

/// `MathUtil.max`: a strict `>` scan, so the first of equal maxima wins.
pub fn java_max(nums: &[f64]) -> f64 {
    let mut max = nums[0];
    for &n in &nums[1..] {
        if n > max {
            max = n;
        }
    }
    max
}

/// `MathUtil.pNormalizeLogProbability`.
pub fn p_normalize_log_probability(values: &[f64; 3]) -> [f64; 3] {
    let max = java_max(values);
    let bump = 300.0 - max;
    let mut tmp = [0.0; 3];
    let mut total = 0.0;
    for i in 0..3 {
        tmp[i] = pow(10.0, values[i] + bump);
        total += tmp[i];
    }
    let max_p = MAX_PROB_BELOW_ONE;
    let min_p = (1.0 - MAX_PROB_BELOW_ONE) / 2.0;
    for value in &mut tmp {
        *value /= total;
        if *value > max_p {
            *value = max_p;
        } else if *value < min_p {
            *value = min_p;
        }
    }
    tmp
}

/// `MathUtil.pNormalizeVector`.
pub fn p_normalize_vector(values: &[f64]) -> Vec<f64> {
    let total: f64 = values.iter().fold(0.0, |a, b| a + b);
    let max_p = MAX_PROB_BELOW_ONE;
    let min_p = (1.0 - MAX_PROB_BELOW_ONE) / (values.len() as f64 - 1.0);
    values
        .iter()
        .map(|v| {
            let p = v / total;
            if p > max_p {
                max_p
            } else if p < min_p {
                min_p
            } else {
                p
            }
        })
        .collect()
}

fn multiply(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] * b[0], a[1] * b[1], a[2] * b[2]]
}

fn sum3(a: &[f64; 3]) -> f64 {
    let mut r = 0.0;
    for v in a {
        r += v;
    }
    r
}

/// `MathUtil.LOG_10_MATH.sum`.
fn log10_sum(values: &[f64; 3]) -> f64 {
    let scaling = java_max(values);
    let mut simple = 0.0;
    for v in values {
        simple += pow(10.0, v - scaling);
    }
    log10(simple) + scaling
}

/// `QualityUtil.getErrorProbabilityFromPhredScore`.
fn error_probability(quality: u8) -> f64 {
    htsjdk_bam::quality_util::error_probability_from_phred_score(i32::from(quality))
        .unwrap_or_else(|| 1.0 / 10f64.powf(f64::from(quality) / 10.0))
}

/// `StringUtil.toUpperCase(byte)`.
pub fn to_upper(b: u8) -> u8 {
    b.to_ascii_uppercase()
}

// ---------------------------------------------------------------------------------------------
// The haplotype map.
// ---------------------------------------------------------------------------------------------

/// `Snp`.
#[derive(Debug, Clone, PartialEq)]
pub struct Snp {
    pub name: String,
    pub chrom: String,
    pub pos: i32,
    pub allele1: u8,
    pub allele2: u8,
    pub maf: f64,
}

impl Snp {
    /// `getAlleleString`.
    pub fn allele_string(&self) -> String {
        format!(
            "{}{}",
            self.allele1 as char,
            (self.allele2 as char).to_ascii_lowercase()
        )
    }

    /// `compareTo`: contig name, then position.
    pub fn compare(&self, other: &Snp) -> std::cmp::Ordering {
        self.chrom.cmp(&other.chrom).then(self.pos.cmp(&other.pos))
    }

    /// `getGenotype(haplotype)`: the three genotypes `DiploidGenotype.fromBases` builds.
    pub fn genotype(&self, index: usize) -> Result<DiploidGenotype, String> {
        let (a, b) = match index {
            0 => (self.allele1, self.allele1),
            1 => (self.allele1, self.allele2),
            _ => (self.allele2, self.allele2),
        };
        DiploidGenotype::from_bases(a, b)
    }

    /// `toString`.
    pub fn display(&self) -> String {
        format!("{}:{}", self.chrom, self.pos)
    }
}

/// `DiploidGenotype`, named by its two bases in `ACGT` order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct DiploidGenotype {
    /// The enum's ordinal: AA AC AG AT CC CG CT GG GT TT.
    pub ordinal: u8,
}

const GENOTYPE_NAMES: [&str; 10] = ["AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"];

impl DiploidGenotype {
    /// `fromBases`: looked up by the SUM of the two upper-cased bases, so the order is irrelevant.
    pub fn from_bases(a: u8, b: u8) -> Result<DiploidGenotype, String> {
        let sum = u32::from(to_upper(a)) + u32::from(to_upper(b));
        for (ordinal, name) in GENOTYPE_NAMES.iter().enumerate() {
            let bytes = name.as_bytes();
            if u32::from(bytes[0]) + u32::from(bytes[1]) == sum {
                return Ok(DiploidGenotype {
                    ordinal: ordinal as u8,
                });
            }
        }
        Err(format!(
            "Unknown genotype string [{}{}], any pair of ACTG case insensitive is acceptable",
            a as char, b as char
        ))
    }

    pub fn name(self) -> &'static str {
        GENOTYPE_NAMES[self.ordinal as usize]
    }

    pub fn is_heterozygous(self) -> bool {
        let b = self.name().as_bytes();
        b[0] != b[1]
    }

    pub fn is_homozygous(self) -> bool {
        !self.is_heterozygous()
    }
}

/// `HaplotypeBlock`.
#[derive(Debug, Clone)]
pub struct Block {
    pub maf: f64,
    pub frequencies: [f64; 3],
    /// Indices into [`HaplotypeMap::snps`], in the iteration order of `snpsByName`.
    pub snps: Vec<usize>,
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub first_snp: usize,
}

impl Block {
    fn new(maf: f64) -> Block {
        let major = 1.0 - maf;
        Block {
            maf,
            frequencies: [major * major, major * maf * 2.0, maf * maf],
            snps: Vec::new(),
            chrom: String::new(),
            start: 0,
            end: 0,
            first_snp: 0,
        }
    }

    /// The key a `Fingerprint` (a `TreeMap`) orders and identifies blocks by.
    pub fn key(&self) -> BlockKey {
        BlockKey {
            chrom: self.chrom.clone(),
            start: self.start,
            end: self.end,
        }
    }

    /// `toString`.
    pub fn display(&self) -> String {
        format!("{}[{}-{}]", self.chrom, self.start, self.end)
    }
}

/// `HaplotypeBlock.compareTo`: contig name, start, end.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BlockKey {
    pub chrom: String,
    pub start: i32,
    pub end: i32,
}

/// `HaplotypeMap`.
#[derive(Debug, Clone)]
pub struct HaplotypeMap {
    /// The header's sequence dictionary.
    pub dictionary: Vec<(String, i32)>,
    pub snps: Vec<Snp>,
    /// `getHaplotypes()`, in the order `fromHaplotypes` added them.
    pub blocks: Vec<Block>,
    /// `haplotypesBySnpLocus` and `snpsByPosition`.
    by_locus: std::collections::HashMap<(String, i32), (usize, usize)>,
}

/// A refusal from building the map, as the class and message the reference throws.
pub type Thrown = (String, String);

impl HaplotypeMap {
    /// `new HaplotypeMap(file)` for a haplotype database (not a VCF).
    pub fn from_database(text: &str, path_for_message: &str) -> Result<HaplotypeMap, Thrown> {
        let mut lines = text.split('\n').map(|l| l.strip_suffix('\r').unwrap_or(l));
        let mut header = String::new();
        let mut first_body: Option<&str> = None;
        for line in lines.by_ref() {
            if line.starts_with('@') {
                header.push_str(line);
                header.push('\n');
            } else {
                first_body = Some(line);
                break;
            }
        }
        if header.is_empty() {
            return Err((
                "java.lang.IllegalStateException".into(),
                format!("Haplotype map file must contain header: {path_for_message}"),
            ));
        }
        let (sam_header, _) = htsjdk_bam::sam_file::read_sam(&header).map_err(|e| {
            (
                "htsjdk.samtools.SAMFormatException".to_string(),
                format!("{e:?}"),
            )
        })?;
        let dictionary: Vec<(String, i32)> = sam_header
            .sequences
            .iter()
            .map(|s| (s.name.clone(), s.length))
            .collect();
        let mut map = HaplotypeMap {
            dictionary,
            snps: Vec::new(),
            blocks: Vec::new(),
            by_locus: std::collections::HashMap::new(),
        };
        // `text.split('\n')` leaves an empty last element for a file ending in a newline, which
        // `readLine` never returns; it is blank and skipped either way.
        let mut anchor_to_block: crate::java_hash_map::JavaHashMap<usize> =
            crate::java_hash_map::JavaHashMap::new();
        let mut blocks: Vec<Block> = Vec::new();
        let mut held: Vec<(Snp, String)> = Vec::new();
        let Some(first) = first_body else {
            return Ok(map);
        };
        let body = std::iter::once(first).chain(lines);
        for line in body {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = java_split_tab(line);
            if fields.len() < 6 || fields.len() > 8 {
                return Err((
                    "picard.PicardException".into(),
                    format!(
                        "Invalid haplotype map record contains {} fields: {line}",
                        fields.len()
                    ),
                ));
            }
            let pos: i32 = fields[1].parse().map_err(|_| {
                (
                    "java.lang.NumberFormatException".to_string(),
                    format!("For input string: \"{}\"", fields[1]),
                )
            })?;
            let maf: f64 = fields[5].parse().map_err(|_| {
                (
                    "java.lang.NumberFormatException".to_string(),
                    format!("For input string: \"{}\"", fields[5]),
                )
            })?;
            let snp = Snp {
                name: fields[2].to_string(),
                chrom: fields[0].to_string(),
                pos,
                allele1: to_upper(fields[3].as_bytes()[0]),
                allele2: to_upper(fields[4].as_bytes()[0]),
                maf,
            };
            let anchor = fields.get(6).copied();
            match anchor {
                Some(a) if !a.trim().is_empty() && a != snp.name => {
                    held.push((snp, a.to_string()));
                }
                _ => {
                    let mut block = Block::new(maf);
                    let name = snp.name.clone();
                    map.snps.push(snp);
                    add_snp(&mut block, &map.snps, map.snps.len() - 1)?;
                    blocks.push(block);
                    anchor_to_block.put(&name, blocks.len() - 1);
                }
            }
        }
        for (snp, anchor) in held {
            let Some(&index) = anchor_to_block
                .iter()
                .find(|(k, _)| *k == anchor)
                .map(|(_, v)| v)
            else {
                return Err((
                    "picard.PicardException".into(),
                    format!("No haplotype found for anchor snp {anchor}"),
                ));
            };
            map.snps.push(snp);
            let s = map.snps.len() - 1;
            add_snp(&mut blocks[index], &map.snps, s)?;
        }
        let order: Vec<usize> = anchor_to_block.iter().map(|(_, v)| *v).collect();
        for index in order {
            map.add_haplotype(blocks[index].clone())?;
        }
        Ok(map)
    }

    /// `addHaplotype`.
    fn add_haplotype(&mut self, block: Block) -> Result<(), Thrown> {
        let b = self.blocks.len();
        for &s in &block.snps {
            let snp = &self.snps[s];
            let key = (snp.chrom.clone(), snp.pos);
            if self.by_locus.contains_key(&key) {
                return Err((
                    "java.lang.IllegalStateException".into(),
                    format!("Same snp name cannot be used twice{}", snp.display()),
                ));
            }
            self.by_locus.insert(key, (b, s));
        }
        self.blocks.push(block);
        Ok(())
    }

    /// `getHaplotype(chrom, pos)` and `getSnp(chrom, pos)`.
    pub fn at(&self, chrom: &str, pos: i32) -> Option<(usize, usize)> {
        self.by_locus.get(&(chrom.to_string(), pos)).copied()
    }

    /// `getSequenceIndex` in the map's own dictionary.
    pub fn sequence_index(&self, chrom: &str) -> i32 {
        self.dictionary
            .iter()
            .position(|(n, _)| n == chrom)
            .map_or(-1, |i| i as i32)
    }

    /// `getActiveDictionary`: the dictionary up to the last contig any SNP is on.
    pub fn active_dictionary(&self) -> Vec<(String, i32)> {
        let max = self
            .snps
            .iter()
            .filter(|s| self.by_locus.contains_key(&(s.chrom.clone(), s.pos)))
            .map(|s| self.sequence_index(&s.chrom))
            .max();
        match max {
            None => self.dictionary.clone(),
            Some(m) => self.dictionary[..(m + 1).max(0) as usize].to_vec(),
        }
    }

    /// `getIntervalList()` as the locus iterator uses it: every SNP's position, sorted by the
    /// dictionary order and uniqued. Only the positions matter: the intervals are single bases.
    pub fn loci(&self, sam_dictionary: &[String]) -> Vec<(i32, i32)> {
        let mut loci: Vec<(i32, i32)> = self
            .blocks
            .iter()
            .flat_map(|b| b.snps.iter())
            .filter_map(|&s| {
                let snp = &self.snps[s];
                sam_dictionary
                    .iter()
                    .position(|n| *n == snp.chrom)
                    .map(|i| (i as i32, snp.pos))
            })
            .collect();
        loci.sort();
        loci.dedup();
        loci
    }
}

/// `String.split("\\t")`: trailing empty strings removed.
fn java_split_tab(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    while fields.len() > 1 && fields.last() == Some(&"") {
        fields.pop();
    }
    if fields.len() == 1 && fields[0].is_empty() {
        return fields;
    }
    fields
}

/// `HaplotypeBlock.addSnp`.
fn add_snp(block: &mut Block, snps: &[Snp], index: usize) -> Result<(), Thrown> {
    let snp = &snps[index];
    if block.snps.is_empty() {
        block.chrom = snp.chrom.clone();
        block.start = snp.pos;
        block.end = snp.pos;
        block.first_snp = index;
    } else if block.chrom != snp.chrom {
        return Err((
            "picard.PicardException".into(),
            format!(
                "Snp chromosome {} does not agree with chromosome of existing snp(s): {}",
                snp.chrom, block.chrom
            ),
        ));
    } else {
        if snp.pos < block.start {
            block.start = snp.pos;
            block.first_snp = index;
        }
        if snp.pos > block.end {
            block.end = snp.pos;
        }
    }
    // `snpsByName.put`: a name already present keeps its place and takes the new SNP.
    let names: Vec<String> = block.snps.iter().map(|&s| snps[s].name.clone()).collect();
    let mut order: crate::java_hash_map::JavaHashMap<usize> =
        crate::java_hash_map::JavaHashMap::new();
    for (n, &s) in names.iter().zip(block.snps.iter()) {
        order.put(n, s);
    }
    order.put(&snp.name, index);
    block.snps = order.iter().map(|(_, v)| *v).collect();
    Ok(())
}

// ---------------------------------------------------------------------------------------------
// HaplotypeProbabilities.
// ---------------------------------------------------------------------------------------------

/// The log-likelihood state `HaplotypeProbabilitiesUsingLogLikelihoods` keeps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LogLikelihoods {
    pub ll: [f64; 3],
}

impl LogLikelihoods {
    fn zero() -> Self {
        LogLikelihoods { ll: [0.0; 3] }
    }

    /// `setLogLikelihoods`: shifted so the largest is 0, then normalised so they sum to one.
    fn set(&mut self, ll: [f64; 3]) {
        let max = java_max(&ll);
        let removed = [ll[0] + -max, ll[1] + -max, ll[2] + -max];
        // `getProbabilityFromLog` refuses an underflow; the maximum here is always zero.
        let mut sum = 0.0;
        for v in removed {
            sum += pow(10.0, v);
        }
        let shift = -log10(sum);
        self.ll = [removed[0] + shift, removed[1] + shift, removed[2] + shift];
    }
}

/// The concrete `HaplotypeProbabilities` classes.
#[derive(Debug, Clone, PartialEq)]
pub enum Kind {
    /// `HaplotypeProbabilitiesFromSequence`.
    Sequence,
    /// `HaplotypeProbabilitiesFromContaminatorSequence`.
    Contaminator {
        contamination: f64,
        map: [[f64; 3]; 3],
    },
    /// `HaplotypeProbabilitiesFromGenotypeLikelihoods`.
    GenotypeLikelihoods,
    /// `HaplotypeProbabilitiesFromGenotype`.
    Genotype { snp: usize, likelihoods: [f64; 3] },
    /// `CappedHaplotypeProbabilities`.
    Capped,
}

/// One block's evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct Probabilities {
    pub block: usize,
    pub kind: Kind,
    pub log: LogLikelihoods,
    pub obs1: i32,
    pub obs2: i32,
    pub other: i32,
}

impl Probabilities {
    pub fn sequence(block: usize) -> Probabilities {
        Probabilities {
            block,
            kind: Kind::Sequence,
            log: LogLikelihoods::zero(),
            obs1: 0,
            obs2: 0,
            other: 0,
        }
    }

    pub fn contaminator(block: usize, contamination: f64) -> Probabilities {
        Probabilities {
            kind: Kind::Contaminator {
                contamination,
                map: [[0.0; 3]; 3],
            },
            ..Probabilities::sequence(block)
        }
    }

    pub fn genotype_likelihoods(block: usize) -> Probabilities {
        Probabilities {
            kind: Kind::GenotypeLikelihoods,
            ..Probabilities::sequence(block)
        }
    }

    pub fn genotype(block: usize, snp: usize, likelihoods: [f64; 3]) -> Probabilities {
        Probabilities {
            kind: Kind::Genotype { snp, likelihoods },
            ..Probabilities::sequence(block)
        }
    }

    /// `new CappedHaplotypeProbabilities(probabilities, cap)`.
    pub fn capped(other: &Probabilities, map: &HaplotypeMap, cap: f64) -> Probabilities {
        let ll = other.log_likelihoods(map);
        let max = java_max(&ll);
        let shifted = [ll[0] + -max, ll[1] + -max, ll[2] + -max];
        let floored = [
            java_math_max(shifted[0], cap),
            java_math_max(shifted[1], cap),
            java_math_max(shifted[2], cap),
        ];
        let mut log = LogLikelihoods::zero();
        log.set(floored);
        Probabilities {
            block: other.block,
            kind: Kind::Capped,
            log,
            obs1: 0,
            obs2: 0,
            other: 0,
        }
    }

    /// `getPriorProbablities`.
    pub fn priors<'a>(&self, map: &'a HaplotypeMap) -> &'a [f64; 3] {
        &map.blocks[self.block].frequencies
    }

    /// `getLogLikelihoods`.
    pub fn log_likelihoods(&self, map: &HaplotypeMap) -> [f64; 3] {
        match &self.kind {
            Kind::Genotype { likelihoods, .. } => [
                log10(likelihoods[0]),
                log10(likelihoods[1]),
                log10(likelihoods[2]),
            ],
            Kind::Contaminator { map: m, .. } => {
                // `updateLikelihoods`, which the getter runs first.
                let priors = map.blocks[self.block].frequencies;
                let log_priors = [log10(priors[0]), log10(priors[1]), log10(priors[2])];
                let mut ll = [0.0; 3];
                for c in 0..3 {
                    let summed = [
                        log_priors[0] + m[c][0],
                        log_priors[1] + m[c][1],
                        log_priors[2] + m[c][2],
                    ];
                    ll[c] = log10_sum(&summed);
                }
                let mut log = LogLikelihoods::zero();
                log.set(ll);
                log.ll
            }
            _ => self.log.ll,
        }
    }

    /// `getLikelihoods`.
    pub fn likelihoods(&self, map: &HaplotypeMap) -> [f64; 3] {
        match &self.kind {
            Kind::Genotype { likelihoods, .. } => *likelihoods,
            _ => p_normalize_log_probability(&self.log_likelihoods(map)),
        }
    }

    /// `getPosteriorLikelihoods`: likelihoods times priors, not normalised.
    pub fn posterior_likelihoods(&self, map: &HaplotypeMap) -> [f64; 3] {
        multiply(&self.likelihoods(map), self.priors(map))
    }

    fn shifted_log_posterior(&self, map: &HaplotypeMap) -> [f64; 3] {
        let ll = self.log_likelihoods(map);
        let f = self.priors(map);
        [
            ll[0] + log10(f[0]),
            ll[1] + log10(f[1]),
            ll[2] + log10(f[2]),
        ]
    }

    /// `getPosteriorProbabilities`.
    pub fn posterior_probabilities(&self, map: &HaplotypeMap) -> [f64; 3] {
        match &self.kind {
            Kind::Genotype { .. } => {
                let v = p_normalize_vector(&self.posterior_likelihoods(map));
                [v[0], v[1], v[2]]
            }
            _ => p_normalize_log_probability(&self.shifted_log_posterior(map)),
        }
    }

    /// `hasEvidence`.
    pub fn has_evidence(&self, map: &HaplotypeMap) -> bool {
        match &self.kind {
            Kind::Genotype { .. } => true,
            Kind::Sequence | Kind::Contaminator { .. } => {
                self.log_likelihoods(map).iter().any(|v| *v != 0.0)
                    || self.obs1 > 0
                    || self.obs2 > 0
            }
            _ => self.log_likelihoods(map).iter().any(|v| *v != 0.0),
        }
    }

    pub fn total_obs(&self) -> i32 {
        match self.kind {
            Kind::Sequence | Kind::Contaminator { .. } => self.obs1 + self.obs2 + self.other,
            _ => 0,
        }
    }

    pub fn obs_allele1(&self) -> i32 {
        match self.kind {
            Kind::Sequence | Kind::Contaminator { .. } => self.obs1,
            _ => 0,
        }
    }

    pub fn obs_allele2(&self) -> i32 {
        match self.kind {
            Kind::Sequence | Kind::Contaminator { .. } => self.obs2,
            _ => 0,
        }
    }

    /// `getRepresentativeSnp`.
    pub fn representative_snp(&self, map: &HaplotypeMap) -> usize {
        match &self.kind {
            Kind::Genotype { snp, .. } => *snp,
            _ => map.blocks[self.block].first_snp,
        }
    }

    /// `getMostLikelyGenotype(snp)`, by `getMostLikelyIndex`.
    pub fn most_likely_index(&self, map: &HaplotypeMap) -> usize {
        let p = self.posterior_probabilities(map);
        if p[0] > p[1] && p[0] > p[2] {
            0
        } else if p[1] > p[2] {
            1
        } else {
            2
        }
    }

    /// `getLodMostProbableGenotype`.
    pub fn lod_most_probable_genotype(&self, map: &HaplotypeMap) -> f64 {
        match &self.kind {
            Kind::Genotype { .. } => {
                let probs = self.posterior_probabilities(map);
                let mut biggest = 0.0;
                let mut second = 0.0;
                for p in probs {
                    if p > biggest {
                        second = biggest;
                        biggest = p;
                        continue;
                    }
                    if p > second {
                        second = p;
                    }
                }
                log10(biggest) - log10(second)
            }
            _ => {
                let logs = self.shifted_log_posterior(map);
                let mut biggest = -f64::MAX;
                let mut second = biggest;
                for p in logs {
                    if p > biggest {
                        second = biggest;
                        biggest = p;
                        continue;
                    }
                    if p > second {
                        second = p;
                    }
                }
                biggest - second
            }
        }
    }

    /// `HaplotypeProbabilitiesFromSequence.addToProbs` and its contaminator override.
    pub fn add_to_probs(&mut self, snp: &Snp, base: u8, quality: u8) {
        match &mut self.kind {
            Kind::Contaminator { contamination, map } => {
                let alt = if base == snp.allele1 {
                    self.obs1 += 1;
                    false
                } else if base == snp.allele2 {
                    self.obs2 += 1;
                    true
                } else {
                    self.other += 1;
                    return;
                };
                let p_err = error_probability(quality);
                let c = *contamination;
                for (cont, row) in map.iter_mut().enumerate() {
                    for (main, cell) in row.iter_mut().enumerate() {
                        let theta = 0.5 * ((1.0 - c) * main as f64 + c * cont as f64);
                        let matching = if alt { theta } else { 1.0 - theta };
                        let opposing = if !alt { theta } else { 1.0 - theta };
                        *cell += log10(matching * (1.0 - p_err) + opposing * p_err);
                    }
                }
            }
            _ => {
                let mut ll = self.log.ll;
                let p_error = error_probability(quality);
                if base == snp.allele1 {
                    self.obs1 += 1;
                    for (g, value) in ll.iter_mut().enumerate() {
                        let p_alt = g as f64 / 2.0;
                        *value += log10((1.0 - p_alt) * (1.0 - p_error) + p_alt * p_error);
                    }
                } else if base == snp.allele2 {
                    self.obs2 += 1;
                    for (g, value) in ll.iter_mut().enumerate() {
                        let p_alt = 1.0 - g as f64 / 2.0;
                        *value += log10((1.0 - p_alt) * (1.0 - p_error) + p_alt * p_error);
                    }
                } else {
                    self.other += 1;
                }
                self.log.set(ll);
            }
        }
    }

    /// `HaplotypeProbabilitiesFromGenotypeLikelihoods.addToLogLikelihoods`.
    pub fn add_to_log_likelihoods(&mut self, snp: &Snp, alleles: (u8, u8), gls: [f64; 3]) {
        let ll = self.log.ll;
        if snp.allele1 == alleles.0 && snp.allele2 == alleles.1 {
            self.log
                .set([ll[0] + gls[0], ll[1] + gls[1], ll[2] + gls[2]]);
        } else if snp.allele2 == alleles.0 && snp.allele1 == alleles.1 {
            self.log
                .set([ll[0] + gls[2], ll[1] + gls[1], ll[2] + gls[0]]);
        }
    }

    /// `merge`, which a merged fingerprint applies to a deep copy of its first member.
    pub fn merge(&mut self, other: &Probabilities, map: &HaplotypeMap) -> Result<(), Thrown> {
        let same_class = std::mem::discriminant(&self.kind) == std::mem::discriminant(&other.kind);
        match (&mut self.kind, &other.kind) {
            (Kind::Genotype { likelihoods, .. }, Kind::Genotype { likelihoods: o, .. }) => {
                for g in 0..3 {
                    likelihoods[g] *= o[g];
                }
                return Ok(());
            }
            (Kind::Genotype { .. }, _) => {
                return Err((
                    "java.lang.IllegalArgumentException".into(),
                    "Can only merge HaplotypeProbabilities of same class.".into(),
                ))
            }
            _ => {}
        }
        if matches!(other.kind, Kind::Genotype { .. }) {
            return Err((
                "java.lang.IllegalArgumentException".into(),
                format!(
                    "Can only merge HaplotypeProbabilities of same class. Found {} and {}",
                    self.java_class(),
                    other.java_class()
                ),
            ));
        }
        let a = self.log_likelihoods(map);
        let b = other.log_likelihoods(map);
        self.log.set([a[0] + b[0], a[1] + b[1], a[2] + b[2]]);
        if matches!(self.kind, Kind::Sequence | Kind::Contaminator { .. }) {
            if !matches!(other.kind, Kind::Sequence | Kind::Contaminator { .. }) {
                return Err((
                    "java.lang.IllegalArgumentException".into(),
                    format!(
                        "Can only merge() HaplotypeProbabilities of same class: Tried to merge a {} with a {}.",
                        self.java_class().trim_start_matches("class "),
                        other.java_class().trim_start_matches("class ")
                    ),
                ));
            }
            self.obs1 += other.obs1;
            self.obs2 += other.obs2;
            self.other += other.other;
        }
        if let Kind::Contaminator { map: m, .. } = &mut self.kind {
            let Kind::Contaminator { map: o, .. } = &other.kind else {
                return Err((
                    "java.lang.IllegalArgumentException".into(),
                    "Can only merge HaplotypeProbabilities of same class.".into(),
                ));
            };
            for c in 0..3 {
                for g in 0..3 {
                    m[c][g] += o[c][g];
                }
            }
        }
        let _ = same_class;
        Ok(())
    }

    fn java_class(&self) -> &'static str {
        match self.kind {
            Kind::Sequence => "class picard.fingerprint.HaplotypeProbabilitiesFromSequence",
            Kind::Contaminator { .. } => {
                "class picard.fingerprint.HaplotypeProbabilitiesFromContaminatorSequence"
            }
            Kind::GenotypeLikelihoods => {
                "class picard.fingerprint.HaplotypeProbabilitiesFromGenotypeLikelihoods"
            }
            Kind::Genotype { .. } => "class picard.fingerprint.HaplotypeProbabilitiesFromGenotype",
            Kind::Capped => "class picard.fingerprint.CappedHaplotypeProbabilities",
        }
    }
}

/// `Math.max(double, double)`.
fn java_math_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a > b {
        a
    } else if b > a {
        b
    } else if a == 0.0 && b == 0.0 {
        if a.is_sign_positive() {
            a
        } else {
            b
        }
    } else {
        a
    }
}

/// `HaplotypeProbabilityOfNormalGivenTumor.getLikelihoods`.
fn tumor_likelihoods(likelihoods: [f64; 3], p_loh: f64) -> [f64; 3] {
    let t = [
        [1.0, 0.0, 0.0],
        [p_loh / 2.0, 1.0 - p_loh, p_loh / 2.0],
        [0.0, 0.0, 1.0],
    ];
    let mut out = [0.0; 3];
    for n in 0..3 {
        out[n] = 0.0;
        for g in 0..3 {
            out[n] += likelihoods[g] * t[n][g];
        }
    }
    out
}

// ---------------------------------------------------------------------------------------------
// Fingerprints and their identities.
// ---------------------------------------------------------------------------------------------

/// `Fingerprint`: a `TreeMap` from block to its evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct Fingerprint {
    pub sample: Option<String>,
    /// The source's URI string, as the reference would print it.
    pub source: Option<String>,
    pub info: Option<String>,
    pub map: BTreeMap<BlockKey, Probabilities>,
}

impl Fingerprint {
    pub fn new(sample: Option<String>, source: Option<String>, info: Option<String>) -> Self {
        Fingerprint {
            sample,
            source,
            info,
            map: BTreeMap::new(),
        }
    }

    pub fn add(&mut self, map: &HaplotypeMap, p: Probabilities) {
        self.map.insert(map.blocks[p.block].key(), p);
    }

    /// `merge`: a key set gathered in a `HashSet` and then visited, so a missing block is a deep
    /// copy of the other's and a shared one is merged into this one's.
    pub fn merge(&mut self, other: &Fingerprint, map: &HaplotypeMap) -> Result<(), Thrown> {
        for (key, theirs) in &other.map {
            match self.map.get_mut(key) {
                None => {
                    self.map.insert(key.clone(), theirs.clone());
                }
                Some(mine) => mine.merge(theirs, map)?,
            }
        }
        Ok(())
    }
}

/// `FingerprintIdDetails`.
#[derive(Debug, Clone, Default)]
pub struct IdDetails {
    pub platform_unit: Option<String>,
    pub run_barcode: Option<String>,
    pub run_lane: Option<i32>,
    pub molecular_barcode: Option<String>,
    pub library: Option<String>,
    pub file: Option<String>,
    pub sample: Option<String>,
    /// Not part of equality or the hash.
    pub group: Option<String>,
}

impl PartialEq for IdDetails {
    fn eq(&self, o: &Self) -> bool {
        self.platform_unit == o.platform_unit
            && self.run_barcode == o.run_barcode
            && self.run_lane == o.run_lane
            && self.molecular_barcode == o.molecular_barcode
            && self.library == o.library
            && self.file == o.file
            && self.sample == o.sample
    }
}

pub const MULTIPLE_VALUES: &str = "<MULTIPLE_VALUES>";

impl IdDetails {
    /// `new FingerprintIdDetails(platformUnit, file)`.
    pub fn from_platform_unit(pu: Option<&str>, file: &str) -> IdDetails {
        let mut id = IdDetails {
            run_barcode: Some("?".into()),
            run_lane: Some(-1),
            molecular_barcode: Some("?".into()),
            ..IdDetails::default()
        };
        if let Some(pu) = pu {
            let parts = java_split_regex_dot(pu);
            if parts.len() == 3 || parts.len() == 2 {
                id.run_barcode = Some(parts[0].clone());
                id.molecular_barcode = Some(if parts.len() == 3 {
                    parts[2].clone()
                } else {
                    String::new()
                });
                if let Ok(lane) = parts[1].parse::<i32>() {
                    if is_java_int(&parts[1]) {
                        id.run_lane = Some(lane);
                    }
                }
            }
        }
        id.platform_unit = pu.map(str::to_string);
        id.file = Some(file.to_string());
        id
    }

    /// `hashCode`.
    pub fn hash(&self) -> i32 {
        let h = |s: &Option<String>| s.as_deref().map_or(0, string_hash_code);
        let mut r = h(&self.platform_unit);
        r = r.wrapping_mul(31).wrapping_add(h(&self.run_barcode));
        r = r
            .wrapping_mul(31)
            .wrapping_add(self.run_lane.map_or(0, integer_hash));
        r = r.wrapping_mul(31).wrapping_add(h(&self.molecular_barcode));
        r = r.wrapping_mul(31).wrapping_add(h(&self.library));
        r = r.wrapping_mul(31).wrapping_add(h(&self.file));
        r = r.wrapping_mul(31).wrapping_add(h(&self.sample));
        r
    }

    /// `merge`.
    pub fn merge(&mut self, o: &IdDetails) {
        fn pick<T: Clone + PartialEq>(l: &Option<T>, r: &Option<T>, or: T) -> Option<T> {
            match (l, r) {
                (_, None) => l.clone(),
                (None, _) => r.clone(),
                (Some(a), Some(b)) => Some(if a == b { a.clone() } else { or }),
            }
        }
        self.platform_unit = pick(
            &self.platform_unit,
            &o.platform_unit,
            MULTIPLE_VALUES.into(),
        );
        self.run_barcode = pick(&self.run_barcode, &o.run_barcode, MULTIPLE_VALUES.into());
        self.run_lane = pick(&self.run_lane, &o.run_lane, i32::MIN);
        self.library = pick(&self.library, &o.library, MULTIPLE_VALUES.into());
        self.file = pick(&self.file, &o.file, MULTIPLE_VALUES.into());
        self.sample = pick(&self.sample, &o.sample, MULTIPLE_VALUES.into());
        self.molecular_barcode = pick(
            &self.molecular_barcode,
            &o.molecular_barcode,
            MULTIPLE_VALUES.into(),
        );
    }
}

/// `Integer.parseInt` accepts an optional sign and decimal digits, nothing else.
fn is_java_int(s: &str) -> bool {
    let digits = s.strip_prefix(['-', '+']).unwrap_or(s);
    !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())
}

/// `String.split("\\.")`: trailing empty strings removed.
fn java_split_regex_dot(s: &str) -> Vec<String> {
    if s.is_empty() {
        return vec![String::new()];
    }
    let mut parts: Vec<String> = s.split('.').map(str::to_string).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

/// `CrosscheckMetric.DataType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    File,
    Sample,
    Library,
    ReadGroup,
}

impl DataType {
    pub fn parse(s: &str) -> Option<DataType> {
        match s {
            "FILE" => Some(DataType::File),
            "SAMPLE" => Some(DataType::Sample),
            "LIBRARY" => Some(DataType::Library),
            "READGROUP" => Some(DataType::ReadGroup),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            DataType::File => "FILE",
            DataType::Sample => "SAMPLE",
            DataType::Library => "LIBRARY",
            DataType::ReadGroup => "READGROUP",
        }
    }
}

/// `Fingerprint.getFingerprintIdDetailsStringFunction`.
pub fn group_of(id: &IdDetails, by: DataType) -> String {
    let s = |v: &Option<String>| v.clone().unwrap_or_else(|| "null".into());
    let value = match by {
        DataType::ReadGroup => id.platform_unit.clone(),
        DataType::Library => Some(format!("{}::{}", s(&id.sample), s(&id.library))),
        DataType::File => Some(format!("{}::{}", s(&id.file), s(&id.sample))),
        DataType::Sample => id.sample.clone(),
    };
    value.unwrap_or_else(|| id.hash().to_string())
}

/// A map of fingerprints keyed by their identity, in `HashMap` order.
pub type FingerprintMap = JavaMap<IdDetails, Fingerprint>;

/// `Fingerprint.mergeFingerprintsBy`.
///
/// `groupingBy` collects the entries into lists in a `HashMap<String, List>` in the order the
/// source map iterates; `toMap` then builds the result in the order THAT map iterates. A group of
/// one keeps its own key and fingerprint object; a larger group gets a fresh key merged from the
/// members and a fresh fingerprint merged from them.
///
/// The members are merged through a `HashSet<Fingerprint>`, whose order is identity hashes and
/// not reproducible. Merging is commutative for two members, which is all the corpus groups.
pub fn merge_fingerprints_by(
    fingerprints: &FingerprintMap,
    by: DataType,
    map: &HaplotypeMap,
) -> Result<FingerprintMap, Thrown> {
    let entries: Vec<(IdDetails, Fingerprint)> = fingerprints
        .iter()
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect();
    merge_entries_by(&entries, by, map)
}

/// [`merge_fingerprints_by`] over entries already in their source map's iteration order.
pub fn merge_entries_by(
    fingerprints: &[(IdDetails, Fingerprint)],
    by: DataType,
    map: &HaplotypeMap,
) -> Result<FingerprintMap, Thrown> {
    let mut index: crate::java_hash_map::JavaHashMap<usize> =
        crate::java_hash_map::JavaHashMap::new();
    let mut lists: Vec<Vec<(IdDetails, Fingerprint)>> = Vec::new();
    for (id, fp) in fingerprints {
        let key = group_of(id, by);
        let slot = match index.get(&key) {
            Some(&i) => i,
            None => {
                lists.push(Vec::new());
                index.put_front_if_absent(&key, lists.len() - 1);
                lists.len() - 1
            }
        };
        lists[slot].push((id.clone(), fp.clone()));
    }
    let collected: Vec<(String, &Vec<(IdDetails, Fingerprint)>)> = index
        .iter()
        .map(|(k, &i)| (k.to_string(), &lists[i]))
        .collect();
    let mut out: FingerprintMap = JavaMap::new();
    for (key, list) in collected {
        let mut final_id = if list.len() == 1 {
            list[0].0.clone()
        } else {
            let mut id = IdDetails::default();
            for (member, _) in list {
                id.merge(member);
            }
            id
        };
        final_id.group = Some(key.to_string());
        let fingerprint = if list.len() == 1 {
            list[0].1.clone()
        } else {
            let first = &list[0].0;
            let mut merged =
                Fingerprint::new(first.sample.clone(), None, Some(group_of(first, by)));
            for (_, fp) in list.iter() {
                merged.merge(fp, map)?;
            }
            merged
        };
        let hash = final_id.hash();
        if out.contains_key(hash, &final_id) {
            return Err((
                "java.lang.IllegalStateException".into(),
                format!("Duplicate key {}", key),
            ));
        }
        out.put(hash, final_id, fingerprint);
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// Reading sequence data.
// ---------------------------------------------------------------------------------------------

/// `ValidationStringency`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stringency {
    Strict,
    Lenient,
    Silent,
}

/// What a fingerprinting walk needs besides the file.
pub struct SamOptions {
    pub minimum_base_quality: u8,
    pub minimum_mapping_quality: u8,
    pub allow_duplicate_reads: bool,
    pub locus_max_reads: i32,
    pub stringency: Stringency,
    pub default_sample: String,
}

impl Default for SamOptions {
    fn default() -> Self {
        SamOptions {
            minimum_base_quality: 20,
            minimum_mapping_quality: 10,
            allow_duplicate_reads: false,
            locus_max_reads: 0,
            stringency: Stringency::Strict,
            default_sample: "<UNKNOWN>".into(),
        }
    }
}

/// A line the reference logs through htsjdk's `Log`, which the harness reads its refusal from.
pub fn log(level: &str, class: &str, message: &str) {
    eprintln!("{level}\t2026-01-01 00:00:00\t{class}\t{message}");
}

/// The SAM flags the locus iterator reads.
const UNMAPPED: u16 = 0x4;
const SECONDARY: u16 = 0x100;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// `FingerprintChecker.fingerprintSamFile`.
///
/// `path` is the path as given (it is what `Path.toString` prints in messages) and `uri` the
/// `toUri().toString()` the reference would have built. `random` is the checker's static
/// `Random(42)`, which a run shares across every file it reads.
#[allow(clippy::too_many_arguments)]
pub fn fingerprint_sam(
    header: &SamHeader,
    records: &[BamRecord],
    path: &str,
    uri: &str,
    map: &HaplotypeMap,
    options: &SamOptions,
    random: &mut JavaRandom,
    make: &dyn Fn(usize) -> Probabilities,
) -> Result<FingerprintMap, Thrown> {
    let names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    check_dictionary(
        map,
        &header
            .sequences
            .iter()
            .map(|s| (s.name.clone(), s.length))
            .collect::<Vec<_>>(),
    )?;
    let sort = header.attributes.get("SO");
    match sort {
        None | Some("unsorted") | Some("coordinate") => {}
        Some(_) => {
            return Err((
                "htsjdk.samtools.SAMException".into(),
                "SamLocusIterator cannot operate on a SAM file that is not coordinate sorted."
                    .into(),
            ))
        }
    }

    // One fingerprint per read group, keyed by its identity.
    let mut by_group: Vec<(String, IdDetails)> = Vec::new();
    let mut fingerprints: FingerprintMap = JavaMap::new();
    for rg in &header.read_groups {
        let mut id = IdDetails::from_platform_unit(rg.attributes.get("PU"), uri);
        id.library = rg.attributes.get("LB").map(str::to_string);
        id.sample = rg.attributes.get("SM").map(str::to_string);
        let mut fp = Fingerprint::new(
            id.sample.clone(),
            Some(uri.to_string()),
            id.platform_unit.clone(),
        );
        for b in 0..map.blocks.len() {
            fp.add(map, make(b));
        }
        fingerprints.put(id.hash(), id.clone(), fp);
        by_group.push((rg.id.clone(), id));
    }
    let mut unknown: Option<IdDetails> = None;

    // The records the iterator would accumulate, by locus, in file order.
    let loci = map.loci(&names);
    let mut pileups: BTreeMap<(i32, i32), Vec<(usize, usize)>> = BTreeMap::new();
    for locus in &loci {
        pileups.insert(*locus, Vec::new());
    }
    for (r, rec) in records.iter().enumerate() {
        let filtered = if options.allow_duplicate_reads {
            rec.flags & SECONDARY != 0
        } else {
            rec.flags & (SECONDARY | SUPPLEMENTARY) != 0 || rec.flags & DUPLICATE != 0
        };
        if filtered {
            continue;
        }
        if rec.reference_index < 0 {
            // `finishedAlignedReads`: the first unplaced record ends the walk.
            break;
        }
        if rec.flags & UNMAPPED != 0 || rec.mapping_quality < options.minimum_mapping_quality {
            continue;
        }
        for block in alignment_blocks(&rec.cigar, rec.alignment_start) {
            for i in 0..block.length {
                let offset = (block.read_start + i - 1) as usize;
                let quality_ok = options.minimum_base_quality == 0
                    || rec.base_qualities.is_empty()
                    || rec.base_qualities[offset] >= options.minimum_base_quality;
                if !quality_ok {
                    continue;
                }
                if let Some(list) =
                    pileups.get_mut(&(rec.reference_index, block.reference_start + i))
                {
                    list.push((r, offset));
                }
            }
        }
    }

    let mut used: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut found = false;
    for ((contig, pos), list) in &pileups {
        let chrom = &names[*contig as usize];
        let Some((block, snp)) = map.at(chrom, *pos) else {
            continue;
        };
        let chosen: Vec<(usize, usize)> =
            if options.locus_max_reads == 0 || list.len() <= options.locus_max_reads as usize {
                list.clone()
            } else {
                // `MathUtil.randomSublist`.
                let mut needed = options.locus_max_reads as usize;
                let mut available = list.len();
                let mut short = Vec::new();
                for item in list {
                    if random.next_double() < needed as f64 / available as f64 {
                        short.push(*item);
                        needed -= 1;
                    }
                    if needed == 0 {
                        break;
                    }
                    available -= 1;
                }
                short
            };
        for (r, offset) in chosen {
            let rec = &records[r];
            let group = match rec.tags.get(Tag::new(b"RG")) {
                Some(TagValue::Str(id)) => by_group.iter().find(|(g, _)| g == id).map(|(_, d)| d),
                _ => None,
            };
            let details = match group {
                Some(d) => d.clone(),
                None => {
                    if unknown.is_none() {
                        let message = format!(
                            "Found read with no readgroup: {} in file: {path}",
                            rec.read_name
                        );
                        if options.stringency == Stringency::Strict {
                            log("ERROR", "FingerprintChecker", &message);
                            return Err(("picard.PicardException".into(), message));
                        }
                        let mut id = IdDetails::from_platform_unit(Some("<UNKNOWN>.0.ZZZ"), uri);
                        id.sample = Some(options.default_sample.clone());
                        id.library = Some("<UNKNOWN>".into());
                        let mut fp = Fingerprint::new(
                            id.sample.clone(),
                            Some(uri.to_string()),
                            id.platform_unit.clone(),
                        );
                        for b in 0..map.blocks.len() {
                            fp.add(map, make(b));
                        }
                        fingerprints.put(id.hash(), id.clone(), fp);
                        unknown = Some(id);
                    }
                    unknown.clone().expect("set above")
                }
            };
            found = true;
            if used.contains(&rec.read_name) {
                continue;
            }
            let key = map.blocks[block].key();
            let fp = fingerprints
                .get_mut(details.hash(), &details)
                .expect("every identity has a fingerprint");
            let probs = fp.map.get_mut(&key).expect("every block has evidence");
            let base = to_upper(rec.read_bases[offset]);
            let quality = rec.base_qualities.get(offset).copied().unwrap_or(0);
            probs.add_to_probs(&map.snps[snp], base, quality);
            used.insert(rec.read_name.clone());
        }
    }
    if !found && sort != Some("coordinate") {
        return Err((
            "picard.PicardException".into(),
            format!(
                "Couldn't even find one locus with reads to fingerprint in file {path}, which in \
                 addition isn't coordinate-sorted. Please sort the file and try again."
            ),
        ));
    }
    Ok(fingerprints)
}

/// `checkDictionaryGoodForFingerprinting`.
pub fn check_dictionary(map: &HaplotypeMap, other: &[(String, i32)]) -> Result<(), Thrown> {
    let active = map.active_dictionary();
    if other.len() < active.len() {
        return Err((
            "htsjdk.samtools.util.SequenceUtil$SequenceListsDifferException".into(),
            "Dictionary on fingerprinted file smaller than that on Haplotype Database!".into(),
        ));
    }
    for (i, (name, length)) in active.iter().enumerate() {
        if other[i].0 != *name || other[i].1 != *length {
            return Err((
                "picard.PicardException".into(),
                "Dictionary on fingerprinted file does not match dictionary in Haplotype Database."
                    .into(),
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------------------------
// Reading genotypes.
// ---------------------------------------------------------------------------------------------

/// One sample's GT allele indices (`None` for a no-call allele) and PLs, either absent.
pub type SampleGenotype = (Option<Vec<Option<usize>>>, Option<Vec<i32>>);

/// One VCF record, reduced to what `getFingerprintFromVc` reads.
#[derive(Debug, Clone)]
pub struct GenotypeRecord {
    pub contig: String,
    pub pos: i32,
    pub filtered: bool,
    /// REF first.
    pub alleles: Vec<String>,
    /// Per sample: the GT allele indices (`None` for a no-call), and the PLs if any.
    pub genotypes: Vec<SampleGenotype>,
}

/// A VCF of genotypes.
#[derive(Debug, Clone)]
pub struct GenotypeFile {
    pub dictionary: Option<Vec<(String, i32)>>,
    pub samples: Vec<String>,
    pub records: Vec<GenotypeRecord>,
}

/// Read a genotypes VCF: the header's contigs and samples, and every record's alleles, filter,
/// GT and PL.
pub fn read_genotype_file(text: &str) -> Result<GenotypeFile, Thrown> {
    let mut dictionary: Vec<(String, i32)> = Vec::new();
    let mut samples: Vec<String> = Vec::new();
    let mut records = Vec::new();
    let mut saw_header = false;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("##contig=<") {
            let body = rest.trim_end_matches('>');
            let mut id = None;
            let mut length = None;
            for kv in body.split(',') {
                if let Some(v) = kv.strip_prefix("ID=") {
                    id = Some(v.to_string());
                } else if let Some(v) = kv.strip_prefix("length=") {
                    length = v.parse::<i32>().ok();
                }
            }
            if let (Some(id), Some(length)) = (id, length) {
                dictionary.push((id, length));
            }
            continue;
        }
        if line.starts_with("##") {
            continue;
        }
        if let Some(rest) = line.strip_prefix('#') {
            let fields: Vec<&str> = rest.split('\t').collect();
            samples = fields.iter().skip(9).map(|s| s.to_string()).collect();
            saw_header = true;
            continue;
        }
        if line.is_empty() || !saw_header {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let mut alleles = vec![f[3].to_string()];
        if f[4] != "." {
            alleles.extend(f[4].split(',').map(str::to_string));
        }
        let filtered = !(f[6] == "PASS" || f[6] == ".");
        let format: Vec<&str> = f.get(8).map_or(Vec::new(), |s| s.split(':').collect());
        let mut genotypes = Vec::new();
        for sample in f.iter().skip(9) {
            let values: Vec<&str> = sample.split(':').collect();
            let mut gt = None;
            let mut pl = None;
            for (k, v) in format.iter().zip(values.iter()) {
                match *k {
                    "GT" => {
                        gt = Some(
                            v.split(['/', '|'])
                                .map(|a| a.parse::<usize>().ok())
                                .collect::<Vec<_>>(),
                        );
                    }
                    "PL" if *v != "." => {
                        pl = Some(v.split(',').filter_map(|x| x.parse().ok()).collect());
                    }
                    _ => {}
                }
            }
            genotypes.push((gt, pl));
        }
        records.push(GenotypeRecord {
            contig: f[0].to_string(),
            pos: f[1].parse().unwrap_or(0),
            filtered,
            alleles,
            genotypes,
        });
    }
    Ok(GenotypeFile {
        dictionary: if dictionary.is_empty() {
            None
        } else {
            Some(dictionary)
        },
        samples,
        records,
    })
}

/// `FingerprintChecker.loadFingerprints` on a VCF with no index, which reads it through.
///
/// The fingerprints are keyed by sample in a `HashMap`, then every header sample that has none
/// gets an empty one.
pub fn load_fingerprints(
    file: &GenotypeFile,
    uri: &str,
    map: &HaplotypeMap,
    specific_sample: Option<&str>,
    genotyping_error_rate: f64,
) -> Result<crate::java_hash_map::JavaHashMap<Fingerprint>, Thrown> {
    let Some(dictionary) = &file.dictionary else {
        return Err(("java.lang.NullPointerException".into(), String::new()));
    };
    check_dictionary(map, dictionary)?;
    let mut fps: crate::java_hash_map::JavaHashMap<Fingerprint> =
        crate::java_hash_map::JavaHashMap::new();
    let mut samples: Option<Vec<String>> = None;
    for record in &file.records {
        if samples.is_none() {
            let s: Vec<String> = match specific_sample {
                Some(s) => vec![s.to_string()],
                // `ctx.getSampleNames()` is a key set of a `HashMap`, so it iterates in hash
                // order; only the insertion into `fingerprints` below sees it.
                None => file.samples.clone(),
            };
            for name in &s {
                fps.put(
                    name,
                    Fingerprint::new(Some(name.clone()), Some(uri.to_string()), None),
                );
            }
            samples = Some(s);
        }
        let Some((block, snp_index)) = map.at(&record.contig, record.pos) else {
            continue;
        };
        let snp = &map.snps[snp_index];
        // `AlleleSubsettingUtils.subsetVCToMatchSnp`, for a biallelic record.
        if record.filtered || record.alleles[0].len() != 1 {
            continue;
        }
        let reference = to_upper(record.alleles[0].as_bytes()[0]);
        let ref_allele = [snp.allele1, snp.allele2]
            .into_iter()
            .find(|b| to_upper(*b) == reference);
        let Some(ref_allele) = ref_allele else {
            continue;
        };
        let other = if snp.allele1 == ref_allele {
            snp.allele2
        } else {
            snp.allele1
        };
        let alt_index = record.alleles[1..]
            .iter()
            .position(|a| a.len() == 1 && to_upper(a.as_bytes()[0]) == to_upper(other));
        let Some(alt_index) = alt_index else {
            continue;
        };
        if record.alleles.len() != 2 {
            // Subsetting a multi-allelic record is not ported.
            let _ = alt_index;
            continue;
        }
        let allele_bases = (
            record.alleles[0].as_bytes()[0],
            record.alleles[1].as_bytes()[0],
        );
        let names: Vec<String> = fps.iter().map(|(k, _)| k.to_string()).collect();
        for sample in names {
            let Some(column) = file.samples.iter().position(|s| *s == sample) else {
                // Thrown inside the loop and caught by the caller as a warning: the record's
                // remaining samples are skipped with it.
                break;
            };
            let (gt, pl) = &record.genotypes[column];
            let fp = fps.get_mut(&sample).expect("put above");
            if let Some(pl) = pl {
                let mut p = Probabilities::genotype_likelihoods(block);
                let gls = [
                    f64::from(pl[0]) / -10.0,
                    f64::from(pl[1]) / -10.0,
                    f64::from(pl[2]) / -10.0,
                ];
                p.add_to_log_likelihoods(snp, allele_bases, gls);
                fp.add(map, p);
            } else {
                let Some(gt) = gt else { continue };
                if gt.iter().all(Option::is_none) || gt.iter().any(Option::is_none) {
                    // `isNoCall()` is every allele a no-call; a half call is not handled here.
                    if gt.iter().all(Option::is_none) {
                        continue;
                    }
                }
                if fp.map.contains_key(&map.blocks[block].key()) {
                    continue;
                }
                let hom = gt.iter().all(|a| *a == gt[0]);
                let first = gt[0].unwrap_or(0);
                let allele = to_upper(record.alleles[first].as_bytes()[0]);
                let half = genotyping_error_rate / 2.0;
                let accuracy = 1.0 - genotyping_error_rate;
                let probs = [
                    if hom && allele == snp.allele1 {
                        accuracy
                    } else {
                        half
                    },
                    if !hom { accuracy } else { half },
                    if hom && allele == snp.allele2 {
                        accuracy
                    } else {
                        half
                    },
                ];
                fp.add(map, Probabilities::genotype(block, snp_index, probs));
            }
        }
    }
    // `computeIfAbsent`, which puts each new entry at the head of its bucket.
    for sample in &file.samples {
        fps.put_front_if_absent(
            sample,
            Fingerprint::new(Some(sample.clone()), Some(uri.to_string()), None),
        );
    }
    Ok(fps)
}

// ---------------------------------------------------------------------------------------------
// Comparing fingerprints.
// ---------------------------------------------------------------------------------------------

/// `LocusResult`.
#[derive(Debug, Clone)]
pub struct LocusResult {
    pub snp: usize,
    pub expected_genotype: DiploidGenotype,
    pub most_likely_genotype: DiploidGenotype,
    pub allele1_count: i32,
    pub allele2_count: i32,
    pub lod_genotype: f64,
}

/// `MatchResults`.
#[derive(Debug, Clone)]
pub struct MatchResults {
    pub sample: Option<String>,
    pub sample_likelihood: f64,
    pub population_likelihood: f64,
    pub lod: f64,
    pub lod_tn: f64,
    pub lod_nt: f64,
    /// Sorted by SNP, duplicates (same contig and position) dropped.
    pub locus_results: Vec<LocusResult>,
}

/// `shiftedLogEvidenceProbabilityUsingGenotypeFrequencies`.
fn shifted_log(likelihoods: &[f64; 3], frequencies: &[f64; 3]) -> f64 {
    log10(sum3(&multiply(likelihoods, frequencies)))
}

/// `FingerprintChecker.calculateMatchResults`.
pub fn calculate_match_results(
    observed: &Fingerprint,
    expected: &Fingerprint,
    map: &HaplotypeMap,
    p_loh: f64,
    locus_info: bool,
    tumor_aware: bool,
) -> Result<MatchResults, Thrown> {
    let mut locus: Vec<LocusResult> = Vec::new();
    let mut no_swap = 0.0;
    let mut swap = 0.0;
    let mut tn = 0.0;
    let mut nt = 0.0;
    for (key, probs2) in &expected.map {
        let Some(probs1) = observed.map.get(key) else {
            continue;
        };
        let l1 = probs1.likelihoods(map);
        let l2 = probs2.likelihoods(map);
        let priors = *probs2.priors(map);
        let post1 = multiply(&l1, &priors);
        let post2 = multiply(&l2, &priors);
        let tumor1 = tumor_likelihoods(l1, p_loh);
        let tumor2 = tumor_likelihoods(l2, p_loh);
        let tumor_post2 = multiply(&tumor2, &priors);
        let snp = probs2.representative_snp(map);
        if locus_info {
            let s = &map.snps[snp];
            if !map.blocks[probs2.block].snps.iter().any(|&i| {
                map.snps[i].name == s.name
                    && map.snps[i].chrom == s.chrom
                    && map.snps[i].pos == s.pos
            }) {
                return Err((
                    "java.lang.IllegalArgumentException".into(),
                    format!(
                        "Snp {} does not belong to haplotype {}.",
                        s.display(),
                        map.blocks[probs2.block].display()
                    ),
                ));
            }
            let expected_genotype = s
                .genotype(probs2.most_likely_index(map))
                .map_err(|m| ("java.lang.IllegalArgumentException".to_string(), m))?;
            let most_likely = s
                .genotype(probs1.most_likely_index(map))
                .map_err(|m| ("java.lang.IllegalArgumentException".to_string(), m))?;
            locus.push(LocusResult {
                snp,
                expected_genotype,
                most_likely_genotype: most_likely,
                allele1_count: probs1.obs_allele1(),
                allele2_count: probs1.obs_allele2(),
                lod_genotype: probs1.lod_most_probable_genotype(map),
            });
        }
        if probs1.has_evidence(map) && probs2.has_evidence(map) {
            no_swap += shifted_log(&l1, &post2);
            swap += shifted_log(&l1, &priors) + shifted_log(&l2, &priors);
            if tumor_aware {
                tn += shifted_log(&tumor1, &post2)
                    - shifted_log(&tumor1, &priors)
                    - shifted_log(&l2, &priors);
                nt += shifted_log(&tumor2, &post1)
                    - shifted_log(&tumor2, &priors)
                    - shifted_log(&l1, &priors);
                let _ = tumor_post2;
            }
        }
    }
    locus.sort_by(|a, b| map.snps[a.snp].compare(&map.snps[b.snp]));
    locus.dedup_by(|a, b| map.snps[a.snp].compare(&map.snps[b.snp]).is_eq());
    Ok(MatchResults {
        sample: expected.sample.clone(),
        sample_likelihood: no_swap,
        population_likelihood: swap,
        lod: no_swap - swap,
        lod_tn: tn,
        lod_nt: nt,
        locus_results: locus,
    })
}

/// `MatchResults.compareTo`: the larger LOD first, then the sample name.
pub fn compare_match_results(a: &MatchResults, b: &MatchResults) -> std::cmp::Ordering {
    if a.lod != b.lod {
        if a.lod > b.lod {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        }
    } else {
        a.sample
            .as_deref()
            .unwrap_or("")
            .cmp(b.sample.as_deref().unwrap_or(""))
    }
}

// ---------------------------------------------------------------------------------------------
// Paths as the reference saw them.
// ---------------------------------------------------------------------------------------------

/// The path a file had where the reference ran.
///
/// A fingerprint's identity hashes its file's URI, and the crosscheck tools write their rows in
/// that hash's order, so a port that hashes the path IT was given orders its rows by a different
/// hash. The coverage harness runs the reference in a container and the port on the host, and
/// passes the mapping between the two in `PICARD_RS_REFERENCE_PATHS` (`host=reference` pairs
/// separated by `;`), which is the inverse of the rewrite it applied to the arguments. Without
/// the variable a path is its own reference path.
pub fn reference_view(path: &str) -> String {
    if let Ok(spec) = std::env::var("PICARD_RS_REFERENCE_PATHS") {
        for pair in spec.split(';') {
            if let Some((host, reference)) = pair.split_once('=') {
                if !host.is_empty() && path.starts_with(host) {
                    return format!("{reference}{}", &path[host.len()..]);
                }
            }
        }
    }
    path.to_string()
}

/// `Path.toAbsolutePath().normalize()` of a string the way `IOUtil.getPath` resolves one.
pub fn absolute_path(path: &str) -> String {
    let joined = if path.starts_with('/') {
        path.to_string()
    } else {
        let cwd = std::env::current_dir()
            .map(|p| p.display().to_string())
            .unwrap_or_default();
        format!("{cwd}/{path}")
    };
    // `normalize`: drop `.` and empty segments, fold `..`.
    let mut parts: Vec<&str> = Vec::new();
    for seg in joined.split('/') {
        match seg {
            "" | "." => {}
            ".." => {
                parts.pop();
            }
            s => parts.push(s),
        }
    }
    format!("/{}", parts.join("/"))
}

/// `path.toUri().toString()`, as the reference would have built it.
pub fn uri_of(path: &str) -> String {
    let abs = reference_view(&absolute_path(path));
    let mut out = String::from("file://");
    for b in abs.bytes() {
        match b {
            b'a'..=b'z'
            | b'A'..=b'Z'
            | b'0'..=b'9'
            | b'/'
            | b'-'
            | b'_'
            | b'.'
            | b'~'
            | b'!'
            | b'$'
            | b'&'
            | b'\''
            | b'('
            | b')'
            | b'*'
            | b'+'
            | b','
            | b';'
            | b'='
            | b':'
            | b'@' => out.push(b as char),
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}

// ---------------------------------------------------------------------------------------------
// Many files.
// ---------------------------------------------------------------------------------------------

/// `CheckFingerprint.fileContainsReads`: by the URI's path, so by extension.
pub fn file_contains_reads(path: &str) -> bool {
    path.ends_with(".bam") || path.ends_with(".sam") || path.ends_with(".cram")
}

/// Whether `getSamReader` would open a queryable reader: `SamFiles.findIndex` looks for `x.bai`
/// (the data extension replaced) or `x.bam.bai`, and the CSI equivalents, beside a BAM; a SAM
/// text file is never queryable, index or not.
fn has_sam_index(path: &str) -> bool {
    if !path.ends_with(".bam") {
        return false;
    }
    let stem = &path[..path.len() - 4];
    [".bai", ".csi"].iter().any(|ext| {
        std::path::Path::new(&format!("{stem}{ext}")).is_file()
            || std::path::Path::new(&format!("{path}{ext}")).is_file()
    })
}

/// `FingerprintChecker.fingerprintVcf`: one fingerprint per sample of the file, keyed by an
/// identity that carries only the sample and the file.
pub fn fingerprint_vcf(
    path: &str,
    map: &HaplotypeMap,
    require_index: bool,
) -> Result<FingerprintMap, Thrown> {
    let uri = uri_of(path);
    if require_index && !std::path::Path::new(&format!("{path}.idx")).exists() {
        return Err((
            "htsjdk.tribble.TribbleException".into(),
            format!(
                "An index is required, but none found with file ending .idx, for input source: {uri}"
            ),
        ));
    }
    let text = std::fs::read_to_string(path)
        .map_err(|e| ("java.io.IOException".to_string(), e.to_string()))?;
    let file = read_genotype_file(&text)?;
    let by_sample = load_fingerprints(&file, &uri, map, None, 0.01)?;
    let mut out: FingerprintMap = JavaMap::new();
    for (sample, fp) in by_sample.iter() {
        let id = IdDetails {
            sample: Some(sample.to_string()),
            file: Some(uri.clone()),
            ..IdDetails::default()
        };
        out.put(id.hash(), id, fp.clone());
    }
    Ok(out)
}

/// `FingerprintChecker.fingerprintFiles` with one thread: every file in order, each file's map
/// `putAll` into a `ConcurrentHashMap` sized for the file count, and the entries handed back in
/// that map's iteration order. A failure is the executor's `Failed to fingerprint`.
pub fn fingerprint_files(
    paths: &[String],
    map: &HaplotypeMap,
    options: &SamOptions,
    require_index: bool,
) -> Result<Vec<(IdDetails, Fingerprint)>, Thrown> {
    let failed = || {
        (
            "picard.PicardException".to_string(),
            "Failed to fingerprint".to_string(),
        )
    };
    let mut all: JavaConcurrentMap<IdDetails, Fingerprint> =
        JavaConcurrentMap::with_capacity(paths.len());
    let mut random = JavaRandom::new(42);
    for path in paths {
        let one = if file_contains_reads(path) {
            if require_index && !has_sam_index(path) {
                return Err(failed());
            }
            let (header, records) = read_reads(path).map_err(|_| failed())?;
            fingerprint_sam(
                &header,
                &records,
                path,
                &uri_of(path),
                map,
                options,
                &mut random,
                &Probabilities::sequence,
            )
            .map_err(|_| failed())?
        } else {
            fingerprint_vcf(path, map, require_index).map_err(|_| failed())?
        };
        if one.is_empty() {
            log(
                "WARN",
                "FingerprintChecker",
                &format!("No fingerprint data was found in file:{path}"),
            );
        }
        all.put_all(one.into_entries());
    }
    Ok(all
        .into_entries()
        .into_iter()
        .map(|(_, k, v)| (k, v))
        .collect())
}

/// A BAM or SAM file, decoded.
pub fn read_reads(path: &str) -> Result<(SamHeader, Vec<BamRecord>), String> {
    let raw = std::fs::read(path).map_err(|e| e.to_string())?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = htsjdk_bam::reader::BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let mut records = Vec::new();
        for r in reader {
            records.push(r.map_err(|e| format!("{e:?}"))?);
        }
        Ok((header, records))
    } else {
        let text = String::from_utf8(raw).map_err(|e| e.to_string())?;
        htsjdk_bam::sam_file::read_sam(&text).map_err(|e| format!("{e:?}"))
    }
}
