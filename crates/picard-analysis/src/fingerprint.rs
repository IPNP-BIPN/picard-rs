//! What every fingerprint tool shares: the haplotype map as `HaplotypeMap` loads it, the per-block
//! evidence (`HaplotypeProbabilities` and its subclasses), a `Fingerprint`, the pileup
//! `FingerprintChecker.fingerprintSamFile` builds one per read group from, the genotypes
//! `getFingerprintFromVc` reads one from, and the comparison `calculateMatchResults` makes.
//!
//! Ported from `picard.fingerprint` and `picard.util.MathUtil` at tag 3.4.0, with the parts of
//! htsjdk 4.2.0's `SamLocusIterator` the pileup goes through.
//!
//! # The arithmetic is the reference's, step for step
//!
//! A block's log-likelihoods are renormalised after every base (`setLogLikelihoods` subtracts the
//! maximum and then the log of the sum), and the posteriors are recomputed from them through
//! `pNormalizeLogProbability`, which unlogs at a bump of 300. Summing the bases first and
//! normalising once gives the same answer to a few ulps and a different file, so the order of
//! operations below is the reference's and not a simplification of it.
//!
//! # Two orders come from hash tables
//!
//! `HaplotypeMap` keeps its blocks in a `HashMap` keyed by anchor name and a block keeps its SNPs in
//! one keyed by SNP name, and `fingerprintSamFile` returns its fingerprints in a `HashMap` keyed by
//! `FingerprintIdDetails`. Where a tool writes in one of those orders, [`java_hash_order`] is it.

#![allow(clippy::needless_range_loop)]

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::java_hash_map::string_hash_code;

/// `MathUtil.MAX_PROB_BELOW_ONE`.
pub const MAX_PROB_BELOW_ONE: f64 = 0.9999999999999999;

/// The order a default `java.util.HashMap` iterates entries put in this order under these hashes:
/// by bucket of the final table (16 doubling past three quarters full), insertion order within a
/// bucket, which `resize` preserves. Keys are assumed distinct; no treeification.
pub fn java_hash_order<T>(items: Vec<(i32, T)>) -> Vec<T> {
    let mut capacity = 16usize;
    while items.len() > capacity * 3 / 4 {
        capacity *= 2;
    }
    let mut keyed: Vec<(usize, T)> = items
        .into_iter()
        .map(|(hash, item)| {
            let h = hash as u32;
            (((h ^ (h >> 16)) as usize) & (capacity - 1), item)
        })
        .collect();
    keyed.sort_by_key(|(bucket, _)| *bucket);
    keyed.into_iter().map(|(_, item)| item).collect()
}

/// `Objects.hashCode` of an optional string.
pub fn hash_or_zero(s: Option<&str>) -> i32 {
    s.map_or(0, string_hash_code)
}

/// `String.compareTo`.
pub fn java_compare(a: &str, b: &str) -> std::cmp::Ordering {
    // UTF-16 code unit order is byte order for ASCII, and every lookup in a block map compares
    // contig names that are almost always equal: the UTF-16 iterators were most of a crosscheck.
    if a == b {
        return std::cmp::Ordering::Equal;
    }
    if a.is_ascii() && b.is_ascii() {
        return a.as_bytes().cmp(b.as_bytes());
    }
    a.encode_utf16().cmp(b.encode_utf16())
}

/// `MathUtil.pNormalizeVector`.
pub fn p_normalize_vector(values: &[f64]) -> Vec<f64> {
    let mut total = 0.0;
    for v in values {
        total += v;
    }
    let min_p = (1.0 - MAX_PROB_BELOW_ONE) / (values.len() - 1) as f64;
    values
        .iter()
        .map(|v| {
            let p = v / total;
            if p > MAX_PROB_BELOW_ONE {
                MAX_PROB_BELOW_ONE
            } else if p < min_p {
                min_p
            } else {
                p
            }
        })
        .collect()
}

/// `MathUtil.max`: the first largest.
fn max(values: &[f64]) -> f64 {
    let mut best = values[0];
    for &v in &values[1..] {
        if v > best {
            best = v;
        }
    }
    best
}

/// `MathUtil.pNormalizeLogProbability`.
pub fn p_normalize_log_probability(values: [f64; 3]) -> [f64; 3] {
    let bump = 300.0 - max(&values);
    let mut tmp = [0.0; 3];
    let mut total = 0.0;
    for i in 0..3 {
        tmp[i] = 10f64.powf(values[i] + bump);
        total += tmp[i];
    }
    let min_p = (1.0 - MAX_PROB_BELOW_ONE) / 2.0;
    for t in &mut tmp {
        *t /= total;
        if *t > MAX_PROB_BELOW_ONE {
            *t = MAX_PROB_BELOW_ONE;
        } else if *t < min_p {
            *t = min_p;
        }
    }
    tmp
}

/// `HaplotypeProbabilitiesUsingLogLikelihoods.setLogLikelihoods`: shift the maximum to zero and
/// divide by the sum.
pub fn normalized_log_likelihoods(ll: [f64; 3]) -> [f64; 3] {
    let m = max(&ll);
    let removed = [ll[0] - m, ll[1] - m, ll[2] - m];
    let mut sum = 0.0;
    for r in removed {
        sum += 10f64.powf(r);
    }
    let shift = -sum.log10();
    [removed[0] + shift, removed[1] + shift, removed[2] + shift]
}

/// `QualityUtil.getErrorProbabilityFromPhredScore`: htsjdk's table, `1 / pow(10, q / 10)`.
pub fn error_probability(quality: u8) -> f64 {
    // htsjdk-rs builds its table again on every call (101 `pow`s): once per quality here.
    static TABLE: std::sync::OnceLock<[f64; 256]> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        std::array::from_fn(|q| {
            htsjdk_bam::quality_util::error_probability_from_phred_score(q as i32).unwrap_or(1.0)
        })
    })[usize::from(quality)]
}

/// `MathUtil.LOG_10_MATH.sum`: the log of a sum of powers, scaled by the largest, with the log
/// taken as `ln(x) / ln(10)` rather than `log10`.
fn log10_math_sum(values: [f64; 3]) -> f64 {
    let scale = max(&values);
    let mut simple = 0.0;
    for v in values {
        simple += 10f64.powf(v - scale);
    }
    simple.ln() / 10f64.ln() + scale
}

/// `DiploidGenotype`, by name: the two bases in enum order, which is alphabetical.
pub fn diploid_genotype(a: u8, b: u8) -> String {
    let (a, b) = (a.to_ascii_uppercase(), b.to_ascii_uppercase());
    let (first, second) = if a <= b { (a, b) } else { (b, a) };
    format!("{}{}", first as char, second as char)
}

/// One row of the map.
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
    /// `Snp.compareTo`: contig name, then position.
    pub fn compare(&self, other: &Snp) -> std::cmp::Ordering {
        java_compare(&self.chrom, &other.chrom).then(self.pos.cmp(&other.pos))
    }

    /// `Snp.hashCode`.
    pub fn hash_code(&self) -> i32 {
        string_hash_code(&self.chrom)
            .wrapping_mul(31)
            .wrapping_add(self.pos)
    }

    /// `getGenotype(DiploidHaplotype.values()[index])`.
    pub fn genotype(&self, index: usize) -> String {
        match index {
            0 => diploid_genotype(self.allele1, self.allele1),
            1 => diploid_genotype(self.allele1, self.allele2),
            _ => diploid_genotype(self.allele2, self.allele2),
        }
    }

    /// `getAlleleString`.
    pub fn allele_string(&self) -> String {
        format!(
            "{}{}",
            self.allele1 as char,
            (self.allele2 as char).to_ascii_lowercase()
        )
    }
}

/// `HaplotypeBlock`.
#[derive(Debug, Clone)]
pub struct Block {
    pub maf: f64,
    pub frequencies: [f64; 3],
    /// `snpsByName`, in its `HashMap` iteration order. Shared with every fingerprint's evidence.
    pub snps: Vec<Arc<Snp>>,
    pub first: Arc<Snp>,
    pub chrom: Arc<str>,
    pub start: i32,
    pub end: i32,
}

/// The key a `TreeMap<HaplotypeBlock, _>` orders by: `HaplotypeBlock.compareTo`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockKey {
    pub chrom: Arc<str>,
    pub start: i32,
    pub end: i32,
}

impl Ord for BlockKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        java_compare(&self.chrom, &other.chrom)
            .then(self.start.cmp(&other.start))
            .then(self.end.cmp(&other.end))
    }
}

impl PartialOrd for BlockKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Block {
    fn new(maf: f64, snp: Snp) -> Block {
        let major = 1.0 - maf;
        Block {
            maf,
            frequencies: [major * major, major * maf * 2.0, maf * maf],
            chrom: Arc::from(snp.chrom.as_str()),
            start: snp.pos,
            end: snp.pos,
            first: Arc::new(snp.clone()),
            snps: vec![Arc::new(snp)],
        }
    }

    /// `addSnp`, after the first.
    fn add(&mut self, snp: Snp) -> Result<(), String> {
        if *self.chrom != *snp.chrom {
            return Err(format!(
                "picard.PicardException: Snp chromosome {} does not agree with chromosome of existing snp(s): {}",
                snp.chrom, self.chrom
            ));
        }
        if snp.pos < self.start {
            self.start = snp.pos;
            self.first = Arc::new(snp.clone());
        }
        if snp.pos > self.end {
            self.end = snp.pos;
        }
        match self.snps.iter().position(|s| s.name == snp.name) {
            Some(at) => self.snps[at] = Arc::new(snp),
            None => self.snps.push(Arc::new(snp)),
        }
        Ok(())
    }

    pub fn key(&self) -> BlockKey {
        BlockKey {
            chrom: self.chrom.clone(),
            start: self.start,
            end: self.end,
        }
    }

    /// `contains`: a SNP of this name at this locus.
    pub fn contains(&self, snp: &Snp) -> bool {
        self.snps
            .iter()
            .any(|s| s.name == snp.name && s.chrom == snp.chrom && s.pos == snp.pos)
    }

    /// `toString`.
    pub fn describe(&self) -> String {
        format!("{}[{}-{}]", self.chrom, self.start, self.end)
    }
}

/// One `@SQ` line of the map's header.
#[derive(Debug, Clone, PartialEq)]
pub struct Sequence {
    pub name: String,
    pub length: i64,
    pub md5: Option<String>,
}

/// `HaplotypeMap`, from the text format.
#[derive(Debug, Clone)]
pub struct HaplotypeMap {
    pub dictionary: Vec<Sequence>,
    /// `getHaplotypes()`: the anchors' `HashMap` order.
    pub blocks: Vec<Block>,
    /// `haplotypesBySnpLocus` and `snpsByPosition`.
    by_locus: std::collections::HashMap<(String, i32), (usize, Arc<Snp>)>,
}

/// The `@SQ` lines of a SAM-format header, in order.
pub fn parse_sq_lines(text: &str) -> Vec<Sequence> {
    text.lines()
        .filter(|l| l.starts_with("@SQ"))
        .map(|l| {
            let mut s = Sequence {
                name: String::new(),
                length: 0,
                md5: None,
            };
            for field in l.split('\t').skip(1) {
                if let Some(v) = field.strip_prefix("SN:") {
                    s.name = v.to_string();
                } else if let Some(v) = field.strip_prefix("LN:") {
                    s.length = v.parse().unwrap_or(0);
                } else if let Some(v) = field.strip_prefix("M5:") {
                    s.md5 = Some(v.to_string());
                }
            }
            s
        })
        .collect()
}

impl HaplotypeMap {
    /// `new HaplotypeMap(file)` on a haplotype database, with the exception it throws as text.
    pub fn load(path: &str) -> Result<HaplotypeMap, String> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| format!("htsjdk.samtools.SAMException: {e}"))?;
        let mut header = String::new();
        let mut rest = Vec::new();
        let mut in_header = true;
        for line in text.lines() {
            if in_header && line.starts_with('@') {
                header.push_str(line);
                header.push('\n');
            } else {
                in_header = false;
                rest.push(line);
            }
        }
        if header.is_empty() {
            let absolute = std::path::absolute(path)
                .map(|p| p.display().to_string())
                .unwrap_or_else(|_| path.to_string());
            return Err(format!(
                "java.lang.IllegalStateException: Haplotype map file must contain header: {absolute}"
            ));
        }
        let dictionary = parse_sq_lines(&header);
        let mut anchors: Vec<(String, Block)> = Vec::new();
        let mut held: Vec<(String, Snp)> = Vec::new();
        // Where each anchor's block is in `anchors`, which keeps the order they were read in.
        let mut anchor_at: foldhash::HashMap<String, usize> = foldhash::HashMap::default();
        for line in rest {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let fields: Vec<&str> = {
                // `String.split` drops trailing empty strings.
                let mut f = fields;
                while f.last() == Some(&"") {
                    f.pop();
                }
                f
            };
            if fields.len() < 6 || fields.len() > 8 {
                return Err(format!(
                    "picard.PicardException: Invalid haplotype map record contains {} fields: {line}",
                    fields.len()
                ));
            }
            let snp = Snp {
                name: fields[2].to_string(),
                chrom: fields[0].to_string(),
                pos: fields[1].trim().parse().map_err(|_| {
                    format!("picard.PicardException: Unparseable int: {}", fields[1])
                })?,
                allele1: fields[3].as_bytes()[0].to_ascii_uppercase(),
                allele2: fields[4].as_bytes()[0].to_ascii_uppercase(),
                maf: fields[5].trim().parse().map_err(|_| {
                    format!("picard.PicardException: Unparseable double: {}", fields[5])
                })?,
            };
            let anchor = fields.get(6).copied();
            match anchor {
                Some(a) if !a.trim().is_empty() && a != snp.name => {
                    held.push((a.to_string(), snp));
                }
                _ => {
                    let block = Block::new(snp.maf, snp.clone());
                    match anchor_at.get(&snp.name) {
                        Some(&at) => anchors[at].1 = block,
                        None => {
                            anchor_at.insert(snp.name.clone(), anchors.len());
                            anchors.push((snp.name.clone(), block));
                        }
                    }
                }
            }
        }
        for (anchor, snp) in held {
            let at = *anchor_at.get(&anchor).ok_or_else(|| {
                format!("picard.PicardException: No haplotype found for anchor snp {anchor}")
            })?;
            anchors[at].1.add(snp)?;
        }
        let blocks: Vec<Block> = java_hash_order(
            anchors
                .into_iter()
                .map(|(name, mut block)| {
                    let snps = std::mem::take(&mut block.snps);
                    block.snps = java_hash_order(
                        snps.into_iter()
                            .map(|s| (string_hash_code(&s.name), s))
                            .collect(),
                    );
                    (string_hash_code(&name), block)
                })
                .collect(),
        );
        let mut by_locus = std::collections::HashMap::new();
        for (index, block) in blocks.iter().enumerate() {
            for snp in &block.snps {
                let key = (snp.chrom.clone(), snp.pos);
                if by_locus.contains_key(&key) {
                    return Err(format!(
                        "java.lang.IllegalStateException: Same snp name cannot be used twice{}:{}",
                        snp.chrom, snp.pos
                    ));
                }
                by_locus.insert(key, (index, snp.clone()));
            }
        }
        Ok(HaplotypeMap {
            dictionary,
            blocks,
            by_locus,
        })
    }

    /// `getHaplotype(chrom, pos)` and `getSnp(chrom, pos)` together.
    pub fn at(&self, chrom: &str, pos: i32) -> Option<(usize, &Arc<Snp>)> {
        self.by_locus
            .get(&(chrom.to_string(), pos))
            .map(|(b, s)| (*b, s))
    }

    /// `getActiveDictionary`: the header's sequences up to the last one a SNP is on.
    pub fn active_dictionary(&self) -> Vec<Sequence> {
        let last = self
            .blocks
            .iter()
            .flat_map(|b| b.snps.iter())
            .filter_map(|s| self.dictionary.iter().position(|d| d.name == s.chrom))
            .max();
        match last {
            Some(i) => self.dictionary[..=i].to_vec(),
            None => self.dictionary.clone(),
        }
    }

    /// Every SNP, by locus, as `uniqued()` leaves the interval list.
    pub fn loci(&self) -> Vec<(String, i32)> {
        let mut loci: Vec<(String, i32)> = self.by_locus.keys().cloned().collect();
        let index = |c: &str| self.dictionary.iter().position(|d| d.name == c);
        loci.sort_by(|a, b| index(&a.0).cmp(&index(&b.0)).then(a.1.cmp(&b.1)));
        loci
    }
}

/// `HaplotypeProbabilities`, as the subclass it is.
#[derive(Debug, Clone)]
pub enum Evidence {
    /// `HaplotypeProbabilitiesFromSequence`: log-likelihoods and the allele counts.
    Sequence {
        ll: [f64; 3],
        obs1: i32,
        obs2: i32,
        other: i32,
    },
    /// `HaplotypeProbabilitiesFromGenotypeLikelihoods`.
    GenotypeLikelihoods { ll: [f64; 3] },
    /// `HaplotypeProbabilitiesFromGenotype`: likelihoods as they are, and the SNP they came from.
    Genotype {
        likelihoods: [f64; 3],
        snp: Arc<Snp>,
    },
    /// `CappedHaplotypeProbabilities`: another's log-likelihoods, floored a distance below their
    /// maximum and renormalised.
    Capped { ll: [f64; 3] },
    /// `HaplotypeProbabilitiesFromContaminatorSequence`: nine models, by the contaminant's
    /// genotype and then the main sample's, kept apart until the likelihoods are read.
    Contaminator {
        map: [[f64; 3]; 3],
        contamination: f64,
        obs1: i32,
        obs2: i32,
        other: i32,
    },
}

/// One block's evidence, with its block's priors and representative SNP to hand.
#[derive(Debug, Clone)]
pub struct Probs {
    pub block: usize,
    pub priors: [f64; 3],
    pub first: Arc<Snp>,
    pub evidence: Evidence,
}

impl Probs {
    pub fn sequence(map: &HaplotypeMap, block: usize) -> Probs {
        Probs::with(
            map,
            block,
            Evidence::Sequence {
                ll: [0.0; 3],
                obs1: 0,
                obs2: 0,
                other: 0,
            },
        )
    }

    pub fn with(map: &HaplotypeMap, block: usize, evidence: Evidence) -> Probs {
        Probs {
            block,
            priors: map.blocks[block].frequencies,
            first: map.blocks[block].first.clone(),
            evidence,
        }
    }

    /// `new CappedHaplotypeProbabilities(probs, cap)`.
    pub fn capped(&self, cap: f64) -> Probs {
        let ll = self.log_likelihoods();
        let m = max(&ll);
        let floored = [
            (ll[0] - m).max(cap),
            (ll[1] - m).max(cap),
            (ll[2] - m).max(cap),
        ];
        Probs {
            block: self.block,
            priors: self.priors,
            first: self.first.clone(),
            evidence: Evidence::Capped {
                ll: normalized_log_likelihoods(floored),
            },
        }
    }

    /// `getLikelihoods`.
    pub fn likelihoods(&self) -> [f64; 3] {
        match &self.evidence {
            Evidence::Sequence { ll, .. }
            | Evidence::GenotypeLikelihoods { ll }
            | Evidence::Capped { ll } => p_normalize_log_probability(*ll),
            Evidence::Contaminator { .. } => p_normalize_log_probability(self.log_likelihoods()),
            Evidence::Genotype { likelihoods, .. } => *likelihoods,
        }
    }

    /// `getLogLikelihoods`.
    pub fn log_likelihoods(&self) -> [f64; 3] {
        match &self.evidence {
            Evidence::Sequence { ll, .. }
            | Evidence::GenotypeLikelihoods { ll }
            | Evidence::Capped { ll } => *ll,
            // `updateLikelihoods`: the main sample's genotype summed out under the priors.
            Evidence::Contaminator { map, .. } => {
                let mut ll = [0.0; 3];
                for c in 0..3 {
                    let mut terms = [0.0; 3];
                    for m in 0..3 {
                        terms[m] = self.priors[m].ln() / 10f64.ln() + map[c][m];
                    }
                    ll[c] = log10_math_sum(terms);
                }
                normalized_log_likelihoods(ll)
            }
            Evidence::Genotype { likelihoods, .. } => likelihoods.map(f64::log10),
        }
    }

    fn uses_logs(&self) -> bool {
        !matches!(self.evidence, Evidence::Genotype { .. })
    }

    fn shifted_log_posterior(&self) -> [f64; 3] {
        let ll = self.log_likelihoods();
        [
            ll[0] + self.priors[0].log10(),
            ll[1] + self.priors[1].log10(),
            ll[2] + self.priors[2].log10(),
        ]
    }

    /// `getPosteriorProbabilities`.
    pub fn posterior_probabilities(&self) -> [f64; 3] {
        if self.uses_logs() {
            p_normalize_log_probability(self.shifted_log_posterior())
        } else {
            let p = p_normalize_vector(&self.posterior_likelihoods());
            [p[0], p[1], p[2]]
        }
    }

    /// `getPosteriorLikelihoods`: likelihoods times priors, unnormalised.
    pub fn posterior_likelihoods(&self) -> [f64; 3] {
        let l = self.likelihoods();
        [
            l[0] * self.priors[0],
            l[1] * self.priors[1],
            l[2] * self.priors[2],
        ]
    }

    /// `getRepresentativeSnp`.
    pub fn representative(&self) -> &Snp {
        match &self.evidence {
            Evidence::Genotype { snp, .. } => snp,
            _ => &self.first,
        }
    }

    pub fn obs1(&self) -> i32 {
        match self.evidence {
            Evidence::Sequence { obs1, .. } | Evidence::Contaminator { obs1, .. } => obs1,
            _ => 0,
        }
    }

    pub fn obs2(&self) -> i32 {
        match self.evidence {
            Evidence::Sequence { obs2, .. } | Evidence::Contaminator { obs2, .. } => obs2,
            _ => 0,
        }
    }

    pub fn total_obs(&self) -> i32 {
        match self.evidence {
            Evidence::Sequence {
                obs1, obs2, other, ..
            }
            | Evidence::Contaminator {
                obs1, obs2, other, ..
            } => obs1 + obs2 + other,
            _ => 0,
        }
    }

    /// `hasEvidence`.
    pub fn has_evidence(&self) -> bool {
        match &self.evidence {
            Evidence::Sequence { ll, obs1, obs2, .. } => {
                ll.iter().any(|d| *d != 0.0) || *obs1 > 0 || *obs2 > 0
            }
            Evidence::GenotypeLikelihoods { ll } | Evidence::Capped { ll } => {
                ll.iter().any(|d| *d != 0.0)
            }
            Evidence::Genotype { .. } => true,
            Evidence::Contaminator { obs1, obs2, .. } => {
                self.log_likelihoods().iter().any(|d| *d != 0.0) || *obs1 > 0 || *obs2 > 0
            }
        }
    }

    /// `getMostLikelyIndex`, which the genotype names come from.
    pub fn most_likely_index(&self) -> usize {
        let p = self.posterior_probabilities();
        if p[0] > p[1] && p[0] > p[2] {
            0
        } else if p[1] > p[2] {
            1
        } else {
            2
        }
    }

    /// `getLodMostProbableGenotype`.
    pub fn lod_most_probable_genotype(&self) -> f64 {
        if self.uses_logs() {
            let mut biggest = -f64::MAX;
            let mut second = biggest;
            for p in self.shifted_log_posterior() {
                if p > biggest {
                    second = biggest;
                    biggest = p;
                } else if p > second {
                    second = p;
                }
            }
            biggest - second
        } else {
            let mut biggest = 0.0;
            let mut second = 0.0;
            for p in self.posterior_probabilities() {
                if p > biggest {
                    second = biggest;
                    biggest = p;
                } else if p > second {
                    second = p;
                }
            }
            f64::log10(biggest) - f64::log10(second)
        }
    }

    /// `HaplotypeProbabilitiesFromSequence.addToProbs`.
    pub fn add_base(&mut self, snp: &Snp, base: u8, quality: u8) {
        let Evidence::Sequence {
            ll,
            obs1,
            obs2,
            other,
        } = &mut self.evidence
        else {
            return;
        };
        let p_error = error_probability(quality);
        if base == snp.allele1 {
            *obs1 += 1;
            for g in 0..3 {
                let p_alt = g as f64 / 2.0;
                ll[g] += ((1.0 - p_alt) * (1.0 - p_error) + p_alt * p_error).log10();
            }
        } else if base == snp.allele2 {
            *obs2 += 1;
            for g in 0..3 {
                let p_alt = 1.0 - g as f64 / 2.0;
                ll[g] += ((1.0 - p_alt) * (1.0 - p_error) + p_alt * p_error).log10();
            }
        } else {
            *other += 1;
        }
        *ll = normalized_log_likelihoods(*ll);
    }

    /// `HaplotypeProbabilitiesFromContaminatorSequence.addToProbs`.
    pub fn add_contaminator_base(&mut self, snp: &Snp, base: u8, quality: u8) {
        let Evidence::Contaminator {
            map,
            contamination,
            obs1,
            obs2,
            other,
        } = &mut self.evidence
        else {
            return;
        };
        let alt = if base == snp.allele1 {
            *obs1 += 1;
            false
        } else if base == snp.allele2 {
            *obs2 += 1;
            true
        } else {
            *other += 1;
            return;
        };
        let p_err = error_probability(quality);
        for c in 0..3 {
            for m in 0..3 {
                let theta = 0.5 * ((1.0 - *contamination) * m as f64 + *contamination * c as f64);
                map[c][m] += ((if alt { theta } else { 1.0 - theta }) * (1.0 - p_err)
                    + (if !alt { theta } else { 1.0 - theta }) * p_err)
                    .log10();
            }
        }
    }

    /// `merge`, for two of the same kind.
    pub fn merge(&mut self, other: &Probs) {
        match (&mut self.evidence, &other.evidence) {
            (
                Evidence::Sequence {
                    ll,
                    obs1,
                    obs2,
                    other: o,
                },
                Evidence::Sequence {
                    ll: ll2,
                    obs1: a,
                    obs2: b,
                    other: c,
                },
            ) => {
                *ll = normalized_log_likelihoods([ll[0] + ll2[0], ll[1] + ll2[1], ll[2] + ll2[2]]);
                *obs1 += a;
                *obs2 += b;
                *o += c;
            }
            (
                Evidence::Contaminator {
                    map,
                    obs1,
                    obs2,
                    other: o,
                    ..
                },
                Evidence::Contaminator {
                    map: m2,
                    obs1: a,
                    obs2: b,
                    other: c,
                    ..
                },
            ) => {
                for g in 0..3 {
                    for h in 0..3 {
                        map[g][h] += m2[g][h];
                    }
                }
                *obs1 += a;
                *obs2 += b;
                *o += c;
            }
            (Evidence::GenotypeLikelihoods { ll }, Evidence::GenotypeLikelihoods { ll: ll2 })
            | (Evidence::Capped { ll }, Evidence::Capped { ll: ll2 }) => {
                *ll = normalized_log_likelihoods([ll[0] + ll2[0], ll[1] + ll2[1], ll[2] + ll2[2]]);
            }
            (
                Evidence::Genotype { likelihoods, .. },
                Evidence::Genotype {
                    likelihoods: l2, ..
                },
            ) => {
                for g in 0..3 {
                    likelihoods[g] *= l2[g];
                }
            }
            _ => {}
        }
    }

    /// `shiftedLogEvidenceProbabilityUsingGenotypeFrequencies`.
    pub fn shifted_log_evidence_using(&self, frequencies: [f64; 3]) -> f64 {
        log_evidence(self.likelihoods(), frequencies)
    }

    /// `shiftedLogEvidenceProbability`.
    pub fn shifted_log_evidence(&self) -> f64 {
        self.shifted_log_evidence_using(self.priors)
    }

    /// `shiftedLogEvidenceProbabilityGivenOtherEvidence`.
    pub fn shifted_log_evidence_given(&self, other: &Probs) -> f64 {
        self.shifted_log_evidence_using(other.posterior_likelihoods())
    }
}

/// `log10(sum(likelihoods * frequencies))`.
fn log_evidence(likelihoods: [f64; 3], frequencies: [f64; 3]) -> f64 {
    let mut sum = 0.0;
    for g in 0..3 {
        sum += likelihoods[g] * frequencies[g];
    }
    sum.log10()
}

/// `HaplotypeProbabilityOfNormalGivenTumor` over some evidence: the likelihoods through the
/// loss-of-heterozygosity transition, and the priors of the block.
pub struct NormalGivenTumor<'a> {
    pub tumor: &'a Probs,
    pub p_loh: f64,
}

impl NormalGivenTumor<'_> {
    pub fn likelihoods(&self) -> [f64; 3] {
        let t = self.tumor.likelihoods();
        let m = [
            [1.0, 0.0, 0.0],
            [self.p_loh / 2.0, 1.0 - self.p_loh, self.p_loh / 2.0],
            [0.0, 0.0, 1.0],
        ];
        let mut n = [0.0; 3];
        for g_n in 0..3 {
            for g_t in 0..3 {
                n[g_n] += t[g_t] * m[g_n][g_t];
            }
        }
        n
    }

    fn posterior_likelihoods(&self) -> [f64; 3] {
        let l = self.likelihoods();
        let p = self.tumor.priors;
        [l[0] * p[0], l[1] * p[1], l[2] * p[2]]
    }

    pub fn shifted_log_evidence(&self) -> f64 {
        log_evidence(self.likelihoods(), self.tumor.priors)
    }

    pub fn shifted_log_evidence_given(&self, other: &Probs) -> f64 {
        log_evidence(self.likelihoods(), other.posterior_likelihoods())
    }
}

/// `Probs::shifted_log_evidence_given` against a tumour-aware other.
pub fn shifted_log_evidence_given_tumor(probs: &Probs, other: &NormalGivenTumor) -> f64 {
    probs.shifted_log_evidence_using(other.posterior_likelihoods())
}

/// `Fingerprint`: a `TreeMap` from block to evidence, with who it is.
#[derive(Debug, Clone)]
pub struct Fingerprint {
    pub sample: Option<String>,
    pub source: Option<String>,
    pub info: Option<String>,
    pub blocks: BTreeMap<BlockKey, Probs>,
}

impl Fingerprint {
    pub fn new(sample: Option<String>, source: Option<String>, info: Option<String>) -> Self {
        Fingerprint {
            sample,
            source,
            info,
            blocks: BTreeMap::new(),
        }
    }

    pub fn add(&mut self, map: &HaplotypeMap, probs: Probs) {
        self.blocks.insert(map.blocks[probs.block].key(), probs);
    }

    /// `merge`.
    pub fn merge(&mut self, other: &Fingerprint) {
        for (key, theirs) in &other.blocks {
            match self.blocks.get_mut(key) {
                Some(mine) => mine.merge(theirs),
                None => {
                    self.blocks.insert(key.clone(), theirs.clone());
                }
            }
        }
    }
}

/// `FingerprintIdDetails`.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct IdDetails {
    pub platform_unit: Option<String>,
    pub run_barcode: Option<String>,
    pub run_lane: Option<i32>,
    pub molecular_barcode: Option<String>,
    pub library: Option<String>,
    pub file: Option<String>,
    pub sample: Option<String>,
    pub group: Option<String>,
}

impl IdDetails {
    /// `new FingerprintIdDetails(platformUnit, file)`.
    pub fn new(platform_unit: Option<&str>, file: &str) -> IdDetails {
        let mut d = IdDetails {
            platform_unit: platform_unit.map(str::to_string),
            run_barcode: Some("?".to_string()),
            run_lane: Some(-1),
            molecular_barcode: Some("?".to_string()),
            file: Some(file.to_string()),
            ..IdDetails::default()
        };
        if let Some(pu) = platform_unit {
            let parts: Vec<&str> = java_split_dot(pu);
            if parts.len() == 3 || parts.len() == 2 {
                d.run_barcode = Some(parts[0].to_string());
                d.molecular_barcode = Some(if parts.len() == 3 {
                    parts[2].to_string()
                } else {
                    String::new()
                });
                if let Ok(lane) = parts[1].parse::<i32>() {
                    d.run_lane = Some(lane);
                }
            }
        }
        d
    }

    /// `hashCode`.
    pub fn hash_code(&self) -> i32 {
        let mut r = hash_or_zero(self.platform_unit.as_deref());
        r = r
            .wrapping_mul(31)
            .wrapping_add(hash_or_zero(self.run_barcode.as_deref()));
        r = r.wrapping_mul(31).wrapping_add(self.run_lane.unwrap_or(0));
        r = r
            .wrapping_mul(31)
            .wrapping_add(hash_or_zero(self.molecular_barcode.as_deref()));
        r = r
            .wrapping_mul(31)
            .wrapping_add(hash_or_zero(self.library.as_deref()));
        r = r
            .wrapping_mul(31)
            .wrapping_add(hash_or_zero(self.file.as_deref()));
        r.wrapping_mul(31)
            .wrapping_add(hash_or_zero(self.sample.as_deref()))
    }

    /// `equals`, which ignores `group`.
    pub fn same(&self, other: &IdDetails) -> bool {
        self.platform_unit == other.platform_unit
            && self.run_barcode == other.run_barcode
            && self.run_lane == other.run_lane
            && self.molecular_barcode == other.molecular_barcode
            && self.library == other.library
            && self.file == other.file
            && self.sample == other.sample
    }

    /// `merge`.
    pub fn merge(&mut self, other: &IdDetails) {
        fn pick<T: Clone + PartialEq>(l: &Option<T>, r: &Option<T>, or: T) -> Option<T> {
            match (l, r) {
                (_, None) => l.clone(),
                (None, _) => r.clone(),
                (Some(a), Some(b)) => Some(if a == b { a.clone() } else { or }),
            }
        }
        let multiple = "<MULTIPLE_VALUES>".to_string();
        self.platform_unit = pick(&self.platform_unit, &other.platform_unit, multiple.clone());
        self.run_barcode = pick(&self.run_barcode, &other.run_barcode, multiple.clone());
        self.run_lane = pick(&self.run_lane, &other.run_lane, i32::MIN);
        self.library = pick(&self.library, &other.library, multiple.clone());
        self.file = pick(&self.file, &other.file, multiple.clone());
        self.sample = pick(&self.sample, &other.sample, multiple.clone());
        self.molecular_barcode = pick(&self.molecular_barcode, &other.molecular_barcode, multiple);
    }
}

/// `String.split("\\.")`: trailing empty strings dropped.
fn java_split_dot(s: &str) -> Vec<&str> {
    let mut parts: Vec<&str> = s.split('.').collect();
    while parts.len() > 1 && parts.last() == Some(&"") {
        parts.pop();
    }
    if parts == [""] {
        return vec![""];
    }
    parts
}

/// The path as the reference saw it: absolute, and with any directory the harness mounted under
/// another name (`PICARD_RS_PATH_MAP`, `host=reference;...`) given that name back, because the
/// reference hashes the path and the iteration order of its tables follows.
pub fn reference_path(path: &str) -> String {
    let absolute = std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string());
    if let Ok(map) = std::env::var("PICARD_RS_PATH_MAP") {
        for pair in map.split(';') {
            if let Some((host, reference)) = pair.split_once('=') {
                if let Some(rest) = absolute.strip_prefix(host) {
                    return format!("{reference}{rest}");
                }
            }
        }
    }
    absolute
}

/// The reverse of `reference_path`: a path the reference was given inside its mounts (read from an
/// input file rather than the command line, which the harness rewrites) as it is on this host.
pub fn host_path(path: &str) -> String {
    if let Ok(map) = std::env::var("PICARD_RS_PATH_MAP") {
        for pair in map.split(';') {
            if let Some((host, reference)) = pair.split_once('=') {
                if let Some(rest) = path.strip_prefix(reference) {
                    if rest.is_empty() || rest.starts_with('/') {
                        return format!("{host}{rest}");
                    }
                }
            }
        }
    }
    path.to_string()
}

/// The `file:` URI `Path.toUri().toString()` gives a local path.
pub fn file_uri(path: &str) -> String {
    format!("file://{}", reference_path(path))
}

/// `LocusResult`.
#[derive(Debug, Clone)]
pub struct LocusResult {
    pub snp: Snp,
    pub expected_genotype: String,
    pub most_likely_genotype: String,
    pub allele1_count: i32,
    pub allele2_count: i32,
    pub lod_genotype: f64,
    pub l_expected_sample: f64,
    pub l_random_sample: f64,
    pub lod_tumor_normal: f64,
    pub lod_normal_tumor: f64,
}

/// `MatchResults`.
#[derive(Debug, Clone)]
pub struct MatchResults {
    pub source: Option<String>,
    pub sample: Option<String>,
    pub sample_likelihood: f64,
    pub population_likelihood: f64,
    pub lod: f64,
    pub lod_tn: f64,
    pub lod_nt: f64,
    /// A `TreeSet` by SNP: sorted, and a second result at one locus dropped.
    pub locus_results: Vec<LocusResult>,
}

/// `FingerprintChecker.calculateMatchResults(observed, expected, pLoH, true, tumorAware)`.
pub fn calculate_match_results(
    map: &HaplotypeMap,
    observed: &Fingerprint,
    expected: &Fingerprint,
    p_loh: f64,
    locus_info: bool,
    tumor_aware: bool,
) -> MatchResults {
    let mut loci: Vec<LocusResult> = Vec::new();
    let (mut no_swap, mut swap, mut tn, mut nt) = (0.0, 0.0, 0.0, 0.0);
    for (key, probs2) in &expected.blocks {
        let Some(probs1) = observed.blocks.get(key) else {
            continue;
        };
        let tumor1 = NormalGivenTumor {
            tumor: probs1,
            p_loh,
        };
        let tumor2 = NormalGivenTumor {
            tumor: probs2,
            p_loh,
        };
        if locus_info {
            let _ = map;
            let snp = probs2.representative().clone();
            loci.push(LocusResult {
                expected_genotype: snp.genotype(probs2.most_likely_index()),
                most_likely_genotype: snp.genotype(probs1.most_likely_index()),
                allele1_count: probs1.obs1(),
                allele2_count: probs1.obs2(),
                lod_genotype: probs1.lod_most_probable_genotype(),
                l_expected_sample: probs1.shifted_log_evidence_given(probs2),
                l_random_sample: probs1.shifted_log_evidence(),
                lod_tumor_normal: if tumor_aware {
                    tumor1.shifted_log_evidence_given(probs2) - tumor1.shifted_log_evidence()
                } else {
                    0.0
                },
                lod_normal_tumor: if tumor_aware {
                    shifted_log_evidence_given_tumor(probs1, &tumor2)
                        - probs1.shifted_log_evidence()
                } else {
                    0.0
                },
                snp: snp.clone(),
            });
        }
        if probs1.has_evidence() && probs2.has_evidence() {
            no_swap += probs1.shifted_log_evidence_given(probs2);
            swap += probs1.shifted_log_evidence() + probs2.shifted_log_evidence();
            if tumor_aware {
                tn += tumor1.shifted_log_evidence_given(probs2)
                    - tumor1.shifted_log_evidence()
                    - probs2.shifted_log_evidence();
                nt += tumor2.shifted_log_evidence_given(probs1)
                    - tumor2.shifted_log_evidence()
                    - probs1.shifted_log_evidence();
            }
        }
    }
    loci.sort_by(|a, b| a.snp.compare(&b.snp));
    loci.dedup_by(|a, b| a.snp.compare(&b.snp).is_eq());
    MatchResults {
        source: expected.source.clone(),
        sample: expected.sample.clone(),
        sample_likelihood: no_swap,
        population_likelihood: swap,
        lod: no_swap - swap,
        lod_tn: tn,
        lod_nt: nt,
        locus_results: loci,
    }
}

/// What [`calculate_match_results`] reads of one block's evidence when it has some, worked out
/// once for a fingerprint that is compared many times.
pub struct BlockTerms<'f> {
    key: &'f BlockKey,
    likelihoods: [f64; 3],
    posterior_likelihoods: [f64; 3],
    shifted_log_evidence: f64,
    /// The same evidence as a tumour's: [`NormalGivenTumor`].
    tumor_likelihoods: [f64; 3],
    tumor_shifted_log_evidence: f64,
}

/// [`BlockTerms`] of every block with evidence, in the fingerprint's block order.
pub fn block_terms(fp: &Fingerprint, p_loh: f64) -> Vec<BlockTerms<'_>> {
    fp.blocks
        .iter()
        .filter(|(_, probs)| probs.has_evidence())
        .map(|(key, probs)| {
            let tumor = NormalGivenTumor {
                tumor: probs,
                p_loh,
            };
            BlockTerms {
                key,
                likelihoods: probs.likelihoods(),
                posterior_likelihoods: probs.posterior_likelihoods(),
                shifted_log_evidence: probs.shifted_log_evidence(),
                tumor_likelihoods: tumor.likelihoods(),
                tumor_shifted_log_evidence: tumor.shifted_log_evidence(),
            }
        })
        .collect()
}

/// [`calculate_match_results`] without locus results, from both sides' [`block_terms`]: the same
/// terms summed in the same order, so the same doubles.
pub fn calculate_match_lods(
    observed: &[BlockTerms],
    expected_fp: &Fingerprint,
    expected: &[BlockTerms],
    tumor_aware: bool,
) -> MatchResults {
    let (mut no_swap, mut swap, mut tn, mut nt) = (0.0, 0.0, 0.0, 0.0);
    let mut at = 0;
    for two in expected {
        while at < observed.len() && observed[at].key < two.key {
            at += 1;
        }
        let Some(one) = observed.get(at).filter(|one| one.key == two.key) else {
            continue;
        };
        no_swap += log_evidence(one.likelihoods, two.posterior_likelihoods);
        swap += one.shifted_log_evidence + two.shifted_log_evidence;
        if tumor_aware {
            tn += log_evidence(one.tumor_likelihoods, two.posterior_likelihoods)
                - one.tumor_shifted_log_evidence
                - two.shifted_log_evidence;
            nt += log_evidence(two.tumor_likelihoods, one.posterior_likelihoods)
                - two.tumor_shifted_log_evidence
                - one.shifted_log_evidence;
        }
    }
    MatchResults {
        source: expected_fp.source.clone(),
        sample: expected_fp.sample.clone(),
        sample_likelihood: no_swap,
        population_likelihood: swap,
        lod: no_swap - swap,
        lod_tn: tn,
        lod_nt: nt,
        locus_results: Vec::new(),
    }
}

/// One base a read shows at a SNP: `SamLocusIterator.RecordAndOffset`, the read held by its place
/// among the reads that showed any.
struct Observed {
    read: u32,
    base: Option<u8>,
    quality: u8,
}

/// `FingerprintChecker.loadFingerprints` on a VCF: one fingerprint per sample, keyed by sample,
/// in the `HashMap` order of the names, with every header sample present even when it has no
/// genotypes (`computeIfAbsent`).
///
/// With an index beside the file the records are read SNP by SNP in `Snp` order, the first record
/// at each position that `subsetVCToMatchSnp` accepts (Picard 3.5.0; 3.4.0 took the first record
/// whatever it was); without one, in file order. A record whose alleles do not match its SNP, or
/// that lacks a sample, is skipped from the sample it failed on, as the caught exception does.
pub fn load_vcf_fingerprints(
    path: &str,
    map: &HaplotypeMap,
    specific_sample: Option<&str>,
    genotyping_error_rate: f64,
) -> Result<Vec<(String, Fingerprint)>, String> {
    let mut vcf = crate::vcf_io::VcfReader::open(path)?;
    if let Err(e) = check_dictionary(&vcf_dictionary(&vcf.file.header), map) {
        // A record the whole-file read refused outranks the dictionary.
        vcf.check()?;
        return Err(e);
    }
    let indexed = std::path::Path::new(&format!("{path}.idx")).exists()
        || std::path::Path::new(&format!("{path}.tbi")).exists();
    let header_samples = vcf.file.header.samples.clone();
    let mut fingerprints: Vec<(String, Fingerprint)> = Vec::new();
    let source = Some(path.to_string());
    let mut started = false;
    let mut start = |fingerprints: &mut Vec<(String, Fingerprint)>| {
        if started {
            return;
        }
        started = true;
        let samples: Vec<String> = match specific_sample {
            Some(s) => vec![s.to_string()],
            None => header_samples.clone(),
        };
        let samples: Vec<String> = java_hash_order(
            samples
                .into_iter()
                .map(|s| (string_hash_code(&s), s))
                .collect(),
        );
        for s in &samples {
            fingerprints.push((
                s.clone(),
                Fingerprint::new(Some(s.clone()), source.clone(), None),
            ));
        }
    };
    if indexed {
        let mut snps: Vec<&Snp> = map
            .blocks
            .iter()
            .flat_map(|b| b.snps.iter().map(|s| &**s))
            .collect();
        snps.sort_by(|a, b| a.compare(b));
        // Each SNP takes the first record in the file that covers its position and matches its
        // alleles: measured on the oracle, a deletion that covers the position before the SNP's
        // own record does not hide it. Per contig, the SNPs as (position, SNP), by position.
        let mut on_contig: foldhash::HashMap<&str, Vec<(i64, usize)>> =
            foldhash::HashMap::default();
        for (at, snp) in snps.iter().enumerate() {
            on_contig
                .entry(snp.chrom.as_str())
                .or_default()
                .push((i64::from(snp.pos), at));
        }
        for on in on_contig.values_mut() {
            on.sort_unstable();
        }
        let mut taken = vec![false; snps.len()];
        // Picard adds the taken records in SNP order, where likelihoods replace a block's evidence
        // and a genotype only fills an empty block: so a block ends with the last likelihoods, or
        // failing those the first genotype. Here the records are added in file order, each kept
        // over what a block holds when it would have come later in SNP order. Per sample and
        // block: whether from likelihoods, and the SNP rank.
        let mut ranks: Vec<Vec<Option<(bool, usize)>>> = Vec::new();
        while let Some(record) = vcf.try_next() {
            let v = record?.variant;
            let Some(on) = on_contig.get(v.contig.as_str()) else {
                continue;
            };
            let from = on.partition_point(|&(pos, _)| pos < v.start);
            // The SNP ranks this record is taken at: the first and the last.
            let mut at: Option<(usize, usize)> = None;
            for &(_, k) in on[from..].iter().take_while(|&&(pos, _)| pos <= v.stop) {
                if !taken[k] && matches_snp(&v, snps[k]) {
                    taken[k] = true;
                    at = Some(at.map_or((k, k), |(lo, hi)| (lo.min(k), hi.max(k))));
                }
            }
            let Some((lo, hi)) = at else {
                continue;
            };
            start(&mut fingerprints);
            ranks.resize_with(fingerprints.len(), || vec![None; map.blocks.len()]);
            for (sample, from_likelihoods, probs) in
                vc_probs(map, &fingerprints, &v, genotyping_error_rate)
            {
                let held = &mut ranks[sample][probs.block];
                let wins = match (*held, from_likelihoods) {
                    (None, _) | (Some((false, _)), true) => true,
                    (Some((true, rank)), true) => hi > rank,
                    (Some((false, rank)), false) => lo < rank,
                    (Some((true, _)), false) => false,
                };
                if wins {
                    *held = Some((from_likelihoods, if from_likelihoods { hi } else { lo }));
                    fingerprints[sample].1.add(map, probs);
                }
            }
        }
    } else {
        while let Some(record) = vcf.try_next() {
            start(&mut fingerprints);
            add_from_vc(
                map,
                &mut fingerprints,
                &record?.variant,
                genotyping_error_rate,
            );
        }
    }
    for s in &header_samples {
        if !fingerprints.iter().any(|(n, _)| n == s) {
            fingerprints.push((
                s.clone(),
                Fingerprint::new(Some(s.clone()), source.clone(), None),
            ));
        }
    }
    Ok(java_hash_order(
        fingerprints
            .into_iter()
            .map(|(n, f)| (string_hash_code(&n), (n, f)))
            .collect(),
    ))
}

/// The `##contig` lines of a VCF header, as a dictionary.
pub fn vcf_dictionary(header: &htsjdk_vcf::header::VcfHeader) -> Vec<Sequence> {
    crate::vcf_io::header_dictionary(header)
        .unwrap_or_default()
        .into_iter()
        .map(|s| Sequence {
            name: s.name,
            length: s.length,
            md5: None,
        })
        .collect()
}

/// `checkDictionaryGoodForFingerprinting`.
pub fn check_dictionary(dictionary: &[Sequence], map: &HaplotypeMap) -> Result<(), String> {
    let active = map.active_dictionary();
    if dictionary.len() < active.len() {
        return Err("htsjdk.samtools.util.SequenceUtil$SequenceListsDifferException: Dictionary on fingerprinted file smaller than that on Haplotype Database!".to_string());
    }
    for (a, b) in active.iter().zip(dictionary) {
        if a.name != b.name || a.length != b.length {
            return Err("picard.PicardException: Dictionary on fingerprinted file does not match dictionary in Haplotype Database.".to_string());
        }
    }
    Ok(())
}

/// `getFingerprintFromVc`.
/// Whether `AlleleSubsettingUtils.subsetVCToMatchSnp` returns a context rather than null: an
/// unfiltered record with a one-base reference that is one of the SNP's alleles, and the SNP's
/// other allele, or `<NON_REF>`, among its alternates.
fn matches_snp(ctx: &htsjdk_vcf::variant::VariantContext, snp: &Snp) -> bool {
    if ctx.is_filtered() || ctx.alleles[0].len() != 1 {
        return false;
    }
    let reference = ctx.alleles[0].base_string().as_bytes()[0].to_ascii_uppercase();
    let Some(ref_allele) = [snp.allele1, snp.allele2]
        .into_iter()
        .find(|b| b.to_ascii_uppercase() == reference)
    else {
        return false;
    };
    let other = if snp.allele1 == ref_allele {
        snp.allele2
    } else {
        snp.allele1
    };
    ctx.alleles[1..].iter().any(|a| {
        // The display string, because a symbolic allele's base string is empty.
        let shown = a.display_string();
        (shown.len() == 1 && shown.as_bytes()[0].eq_ignore_ascii_case(&other))
            || shown == "<NON_REF>"
    })
}

fn add_from_vc(
    map: &HaplotypeMap,
    fingerprints: &mut [(String, Fingerprint)],
    ctx: &htsjdk_vcf::variant::VariantContext,
    error_rate: f64,
) {
    for (at, from_likelihoods, probs) in vc_probs(map, fingerprints, ctx, error_rate) {
        let fp = &mut fingerprints[at].1;
        if from_likelihoods || !fp.blocks.contains_key(&map.blocks[probs.block].key()) {
            fp.add(map, probs);
        }
    }
}

/// What `add_from_vc` adds, sample by sample (by place in `fingerprints`): the evidence, and
/// whether it came from likelihoods, which replace a block's evidence where a genotype only fills
/// a block that has none.
fn vc_probs(
    map: &HaplotypeMap,
    fingerprints: &[(String, Fingerprint)],
    ctx: &htsjdk_vcf::variant::VariantContext,
    error_rate: f64,
) -> Vec<(usize, bool, Probs)> {
    let mut added = Vec::new();
    let Some((block, snp)) = map.at(&ctx.contig, ctx.start as i32) else {
        return added;
    };
    let snp = snp.clone();
    // `subsetVCToMatchSnp`, for a biallelic site.
    if ctx.is_filtered() || ctx.alleles[0].len() != 1 {
        return added;
    }
    let reference = ctx.alleles[0].base_string().as_bytes()[0].to_ascii_uppercase();
    let ref_allele = [snp.allele1, snp.allele2]
        .into_iter()
        .find(|b| b.to_ascii_uppercase() == reference);
    let Some(ref_allele) = ref_allele else {
        return added;
    };
    let other = if snp.allele1 == ref_allele {
        snp.allele2
    } else {
        snp.allele1
    };
    let alts = &ctx.alleles[1..];
    let matching = alts
        .iter()
        .any(|a| a.len() == 1 && a.base_string().as_bytes()[0].eq_ignore_ascii_case(&other));
    if !matching || alts.len() != 1 {
        return added;
    }
    for allele in &ctx.alleles {
        let bases = allele.base_string();
        let b = bases.as_bytes();
        if b.len() > 1 || (b[0] != snp.allele1 && b[0] != snp.allele2) {
            return added;
        }
    }
    let first = ctx.alleles[0].base_string().as_bytes()[0];
    let second = ctx.alleles[1].base_string().as_bytes()[0];
    // The first genotype of each sample, found by name rather than by a scan per sample.
    let mut by_sample: foldhash::HashMap<&str, &_> = foldhash::HashMap::default();
    for g in ctx.genotypes.iter().rev() {
        by_sample.insert(g.sample_name.as_str(), g);
    }
    for (at, (sample, _)) in fingerprints.iter().enumerate() {
        let Some(&genotype) = by_sample.get(sample.as_str()) else {
            return added;
        };
        if let Some(pl) = &genotype.pl {
            if pl.len() != 3 {
                return added;
            }
            let gl = [
                f64::from(pl[0]) / -10.0,
                f64::from(pl[1]) / -10.0,
                f64::from(pl[2]) / -10.0,
            ];
            let ll = if snp.allele1 == first && snp.allele2 == second {
                gl
            } else if snp.allele2 == first && snp.allele1 == second {
                [gl[2], gl[1], gl[0]]
            } else {
                return added;
            };
            let probs = Probs::with(
                map,
                block,
                Evidence::GenotypeLikelihoods {
                    ll: normalized_log_likelihoods(ll),
                },
            );
            added.push((at, true, probs));
        } else {
            if genotype.is_no_call() {
                continue;
            }
            let hom = genotype.is_hom();
            let allele = genotype.alleles[0].base_string().as_bytes()[0].to_ascii_uppercase();
            let half = error_rate / 2.0;
            let accuracy = 1.0 - error_rate;
            let likelihoods = [
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
            added.push((
                at,
                false,
                Probs::with(
                    map,
                    block,
                    Evidence::Genotype {
                        likelihoods,
                        snp: snp.clone(),
                    },
                ),
            ));
        }
    }
    added
}

/// What `fingerprintSamFile` needs to know beyond the file.
pub struct SamOptions {
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub allow_duplicates: bool,
    /// `locusMaxReads`; 0 keeps every read.
    pub locus_max_reads: usize,
    pub strict: bool,
    pub default_sample: String,
}

impl Default for SamOptions {
    fn default() -> Self {
        SamOptions {
            min_mapping_quality: 10,
            min_base_quality: 20,
            allow_duplicates: false,
            locus_max_reads: 0,
            strict: true,
            default_sample: "<UNKNOWN>".to_string(),
        }
    }
}

/// `FingerprintChecker`'s `static Random(42)`, shared by every file the run fingerprints.
pub struct SharedRandom(pub crate::theoretical_sensitivity::JavaRandom);

/// `MathUtil.randomSublist`.
pub fn random_sublist<T: Copy>(list: &[T], n: usize, random: &mut SharedRandom) -> Vec<T> {
    if list.len() <= n {
        return list.to_vec();
    }
    let mut still_needed = n;
    let mut available = list.len();
    let mut short = Vec::with_capacity(n);
    for &item in list {
        if random.0.next_double() < still_needed as f64 / available as f64 {
            short.push(item);
            still_needed -= 1;
        }
        if still_needed == 0 {
            break;
        }
        available -= 1;
    }
    short
}

/// `fingerprintSamFile`: one fingerprint per read group, in the `HashMap` order of their
/// `FingerprintIdDetails`.
pub fn fingerprint_sam_file(
    path: &str,
    map: &HaplotypeMap,
    options: &SamOptions,
    random: &mut SharedRandom,
    make: &dyn Fn(&HaplotypeMap, usize) -> Probs,
    add: &dyn Fn(&mut Probs, &Snp, u8, u8),
) -> Result<Vec<(IdDetails, Fingerprint)>, String> {
    let input = crate::metrics_cli::open_input(path);
    fingerprint_sam_input(input, path, map, options, random, make, add)
}

/// [`fingerprint_sam_file`] on an input already opened from `path`.
pub fn fingerprint_sam_input(
    input: crate::bam_stream::Input,
    path: &str,
    map: &HaplotypeMap,
    options: &SamOptions,
    random: &mut SharedRandom,
    make: &dyn Fn(&HaplotypeMap, usize) -> Probs,
    add: &dyn Fn(&mut Probs, &Snp, u8, u8),
) -> Result<Vec<(IdDetails, Fingerprint)>, String> {
    use htsjdk_bam::tag::{Tag, TagValue};
    let (header, records) = input.into_parts();
    let dictionary: Vec<Sequence> = header
        .sequences
        .iter()
        .map(|s| Sequence {
            name: s.name.clone(),
            length: i64::from(s.length),
            md5: None,
        })
        .collect();
    check_dictionary(&dictionary, map)?;
    let uri = file_uri(path);
    // The read groups, each with its details and fingerprint; `None` is the unknown group.
    let mut groups: Vec<(Option<String>, IdDetails, Fingerprint)> = Vec::new();
    let new_fp = |details: &IdDetails| {
        let mut fp = Fingerprint::new(
            details.sample.clone(),
            Some(path.to_string()),
            details.platform_unit.clone(),
        );
        for b in 0..map.blocks.len() {
            fp.add(map, make(map, b));
        }
        fp
    };
    for rg in &header.read_groups {
        let mut details = IdDetails::new(rg.attributes.get("PU"), &uri);
        details.library = rg.attributes.get("LB").map(str::to_string);
        details.sample = rg.attributes.get("SM").map(str::to_string);
        let fp = new_fp(&details);
        groups.push((Some(rg.id.clone()), details, fp));
    }
    let mut loci = map.loci();
    let order = |c: &str| dictionary.iter().position(|d| d.name == c);
    loci.sort_by(|a, b| order(&a.0).cmp(&order(&b.0)).then(a.1.cmp(&b.1)));
    // Per reference index, the loci on its contig as (position, locus), by position.
    let mut on_contig: std::collections::HashMap<&str, Vec<(i32, usize)>> =
        std::collections::HashMap::new();
    for (at, (contig, pos)) in loci.iter().enumerate() {
        on_contig
            .entry(contig.as_str())
            .or_default()
            .push((*pos, at));
    }
    for on in on_contig.values_mut() {
        on.sort_unstable();
    }
    let none: Vec<(i32, usize)> = Vec::new();
    let by_reference: Vec<&[(i32, usize)]> = header
        .sequences
        .iter()
        .map(|s| on_contig.get(s.name.as_str()).unwrap_or(&none).as_slice())
        .collect();

    // `SamLocusIterator` over the loci with `FingerprintChecker`'s settings, one read at a time:
    // one list per locus, holding each passing base in the order the reads came. Secondary and
    // supplementary reads are dropped, and duplicates too unless `allow_duplicates` (which swaps
    // the filters for `SecondaryAlignmentFilter` alone, so a supplementary read is then KEPT); so
    // are unmapped reads and reads below the mapping-quality cutoff, and bases below the
    // base-quality cutoff. Only the reads that show a base are kept: their name and read group.
    let mut piles: Vec<Vec<Observed>> = (0..loci.len()).map(|_| Vec::new()).collect();
    // shortcut: u32 read numbers, so at most 4G reads that show a base; widen if one ever has more.
    let mut reads: Vec<(Box<str>, Option<usize>)> = Vec::new();
    for r in records {
        let filtered = if options.allow_duplicates {
            r.flags & 0x100 != 0
        } else {
            r.flags & (0x100 | 0x800 | 0x400) != 0
        };
        if filtered || r.flags & 0x4 != 0 || r.mapping_quality < options.min_mapping_quality {
            continue;
        }
        let Some(&on) = usize::try_from(r.reference_index)
            .ok()
            .and_then(|i| by_reference.get(i))
        else {
            continue;
        };
        if on.is_empty() {
            continue;
        }
        let read = reads.len() as u32;
        let mut shown = false;
        let mut position = r.alignment_start;
        let mut offset = 0usize;
        for e in &r.cigar.elements {
            let length = e.length as usize;
            match char::from(e.op.to_char()) {
                'M' | '=' | 'X' => {
                    let from = on.partition_point(|&(p, _)| p < position);
                    for &(p, at) in on[from..]
                        .iter()
                        .take_while(|&&(p, _)| p < position + length as i32)
                    {
                        let step = (p - position) as usize;
                        let q = r.base_qualities.get(offset + step).copied();
                        if r.base_qualities.is_empty()
                            || q.is_some_and(|q| q >= options.min_base_quality)
                        {
                            piles[at].push(Observed {
                                read,
                                base: r.read_bases.get(offset + step).copied(),
                                quality: q.unwrap_or(0),
                            });
                            shown = true;
                        }
                    }
                    offset += length;
                    position += length as i32;
                }
                'I' | 'S' => offset += length,
                'D' | 'N' => position += length as i32,
                _ => {}
            }
        }
        if shown {
            // The header's read group, by its place, which is its place among `groups` too.
            let read_group = match r.tags.get(Tag::new(b"RG")) {
                Some(TagValue::Str(id)) => header.read_groups.iter().position(|g| g.id == *id),
                _ => None,
            };
            reads.push((r.read_name.into_boxed_str(), read_group));
        }
    }
    let mut used: foldhash::HashSet<&str> = foldhash::HashSet::default();
    let mut found = false;
    for ((contig, pos), pile) in loci.iter().zip(&piles) {
        let (block, snp) = map.at(contig, *pos).expect("a map locus");
        let indices: Vec<usize> = (0..pile.len()).collect();
        let chosen = if options.locus_max_reads == 0 {
            indices
        } else {
            random_sublist(&indices, options.locus_max_reads, random)
        };
        for i in chosen {
            let rec = &pile[i];
            let (name, rg) = &reads[rec.read as usize];
            if rg.is_none() && !groups.iter().any(|(g, _, _)| g.is_none()) {
                if options.strict {
                    return Err(format!(
                        "picard.PicardException: Found read with no readgroup: {} in file: {path}",
                        name
                    ));
                }
                let mut details = IdDetails::new(Some("<UNKNOWN>.0.ZZZ"), &uri);
                details.library = Some("<UNKNOWN>".to_string());
                details.sample = Some(options.default_sample.clone());
                let fp = new_fp(&details);
                groups.push((None, details, fp));
            }
            found = true;
            if used.insert(name) {
                let at = match rg {
                    Some(at) => *at,
                    None => groups
                        .iter()
                        .position(|(g, _, _)| g.is_none())
                        .expect("known"),
                };
                let key = map.blocks[block].key();
                let probs = groups[at].2.blocks.get_mut(&key).expect("every block");
                let base = rec.base.expect("a read base at every aligned offset");
                add(probs, snp, base.to_ascii_uppercase(), rec.quality);
            }
        }
    }
    if !found && crate::metrics_cli::sort_order(&header) != "coordinate" {
        return Err(format!("picard.PicardException: Couldn't even find one locus with reads to fingerprint in file {path}, which in addition isn't coordinate-sorted. Please sort the file and try again."));
    }
    // Equal details collapse into one entry, the later fingerprint replacing the earlier.
    let mut unique: Vec<(IdDetails, Fingerprint)> = Vec::new();
    for (_, details, fp) in groups {
        match unique.iter().position(|(d, _)| d.same(&details)) {
            Some(at) => unique[at].1 = fp,
            None => unique.push((details, fp)),
        }
    }
    Ok(java_hash_order(
        unique
            .into_iter()
            .map(|(d, f)| (d.hash_code(), (d, f)))
            .collect(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_hash_order_is_by_bucket_then_insertion() {
        // "b" (98) and "a" (97) land in buckets 2 and 1 of sixteen.
        let order = java_hash_order(vec![
            (string_hash_code("b"), "b"),
            (string_hash_code("a"), "a"),
        ]);
        assert_eq!(order, vec!["a", "b"]);
    }

    #[test]
    fn genotype_names_are_alphabetical() {
        assert_eq!(diploid_genotype(b'G', b'a'), "AG");
        assert_eq!(diploid_genotype(b'T', b'T'), "TT");
    }

    #[test]
    fn normalised_log_likelihoods_sum_to_one() {
        let ll = normalized_log_likelihoods([-1.0, -2.0, -3.0]);
        let total: f64 = ll.iter().map(|l| 10f64.powf(*l)).sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn match_lods_from_block_terms_are_the_match_results() {
        let path = std::env::temp_dir().join(format!("fp-terms-{}.txt", std::process::id()));
        std::fs::write(
            &path,
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
             chr1\t10\ts10\tA\tC\t0.3\nchr1\t50\ts50\tA\tG\t0.2\n\
             chr1\t90\ts90\tC\tT\t0.4\nchr1\t95\ts95\tG\tT\t0.1\n",
        )
        .unwrap();
        let map = HaplotypeMap::load(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).unwrap();
        let fingerprint = |lls: &[Option<[f64; 3]>]| {
            let mut fp = Fingerprint::new(Some("s".into()), Some("f".into()), None);
            for (block, ll) in lls.iter().enumerate() {
                if let Some(ll) = ll {
                    let ll = normalized_log_likelihoods(*ll);
                    fp.add(
                        &map,
                        Probs::with(&map, block, Evidence::GenotypeLikelihoods { ll }),
                    );
                }
            }
            fp
        };
        // A block only one side has, and one with no evidence (all zero).
        let one = fingerprint(&[
            Some([-0.1, -2.0, -5.0]),
            Some([0.0; 3]),
            Some([-3.0, -0.2, -4.0]),
            None,
        ]);
        let two = fingerprint(&[
            Some([-4.0, -0.5, -0.1]),
            Some([-1.0, -1.0, -0.3]),
            Some([-2.5, -0.1, -6.0]),
            Some([-0.2, -1.0, -2.0]),
        ]);
        for tumor_aware in [false, true] {
            let slow = calculate_match_results(&map, &one, &two, 0.5, false, tumor_aware);
            let fast = calculate_match_lods(
                &block_terms(&one, 0.5),
                &two,
                &block_terms(&two, 0.5),
                tumor_aware,
            );
            for (a, b) in [
                (slow.sample_likelihood, fast.sample_likelihood),
                (slow.population_likelihood, fast.population_likelihood),
                (slow.lod, fast.lod),
                (slow.lod_tn, fast.lod_tn),
                (slow.lod_nt, fast.lod_nt),
            ] {
                assert_eq!(a.to_bits(), b.to_bits());
            }
            assert_ne!(fast.lod, 0.0);
        }
    }
}
