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
    htsjdk_bam::quality_util::error_probability_from_phred_score(i32::from(quality)).unwrap_or(1.0)
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
    /// `snpsByName`, in its `HashMap` iteration order.
    pub snps: Vec<Snp>,
    pub first: Snp,
    pub chrom: String,
    pub start: i32,
    pub end: i32,
}

/// The key a `TreeMap<HaplotypeBlock, _>` orders by: `HaplotypeBlock.compareTo`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockKey {
    pub chrom: String,
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
            chrom: snp.chrom.clone(),
            start: snp.pos,
            end: snp.pos,
            first: snp.clone(),
            snps: vec![snp],
        }
    }

    /// `addSnp`, after the first.
    fn add(&mut self, snp: Snp) -> Result<(), String> {
        if self.chrom != snp.chrom {
            return Err(format!(
                "picard.PicardException: Snp chromosome {} does not agree with chromosome of existing snp(s): {}",
                snp.chrom, self.chrom
            ));
        }
        if snp.pos < self.start {
            self.start = snp.pos;
            self.first = snp.clone();
        }
        if snp.pos > self.end {
            self.end = snp.pos;
        }
        match self.snps.iter().position(|s| s.name == snp.name) {
            Some(at) => self.snps[at] = snp,
            None => self.snps.push(snp),
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
    by_locus: std::collections::HashMap<(String, i32), (usize, Snp)>,
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
                    match anchors.iter().position(|(n, _)| *n == snp.name) {
                        Some(at) => anchors[at].1 = block,
                        None => anchors.push((snp.name.clone(), block)),
                    }
                }
            }
        }
        for (anchor, snp) in held {
            let block = anchors
                .iter_mut()
                .find(|(n, _)| *n == anchor)
                .ok_or_else(|| {
                    format!("picard.PicardException: No haplotype found for anchor snp {anchor}")
                })?;
            block.1.add(snp)?;
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
        let mut seen: Vec<(String, i32)> = Vec::new();
        for (index, block) in blocks.iter().enumerate() {
            for snp in &block.snps {
                let key = (snp.chrom.clone(), snp.pos);
                if seen.contains(&key) {
                    return Err(format!(
                        "java.lang.IllegalStateException: Same snp name cannot be used twice{}:{}",
                        snp.chrom, snp.pos
                    ));
                }
                seen.push(key.clone());
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
    pub fn at(&self, chrom: &str, pos: i32) -> Option<(usize, &Snp)> {
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
    Genotype { likelihoods: [f64; 3], snp: Snp },
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
    pub first: Snp,
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
        let snp = probs2.representative().clone();
        if locus_info {
            let _ = map;
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

/// A read as the pileup needs it.
pub struct PileupRead {
    pub name: String,
    pub flags: u16,
    pub reference: Option<String>,
    pub start: i32,
    pub mapping_quality: u8,
    pub cigar: Vec<(usize, char)>,
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
    pub read_group: Option<String>,
}

/// One base a read shows at a SNP: `SamLocusIterator.RecordAndOffset`.
pub struct Observed<'r> {
    pub read: &'r PileupRead,
    pub offset: usize,
}

/// `SamLocusIterator` over the map's loci with `FingerprintChecker`'s settings: one list per
/// locus, in locus order, holding each passing base in the order the reads came.
///
/// Secondary and supplementary reads are dropped, and duplicates too unless `allow_duplicates`
/// (which swaps the filters for `SecondaryAlignmentFilter` alone, so a supplementary read is
/// then KEPT); so are unmapped reads and reads below the mapping-quality cutoff, and bases below
/// the base-quality cutoff.
pub fn pileup<'r>(
    reads: &'r [PileupRead],
    loci: &[(String, i32)],
    min_mapping_quality: u8,
    min_base_quality: u8,
    allow_duplicates: bool,
) -> Vec<Vec<Observed<'r>>> {
    let index: std::collections::HashMap<(&str, i32), usize> = loci
        .iter()
        .enumerate()
        .map(|(i, (c, p))| ((c.as_str(), *p), i))
        .collect();
    let mut piles: Vec<Vec<Observed>> = (0..loci.len()).map(|_| Vec::new()).collect();
    for read in reads {
        let filtered = if allow_duplicates {
            read.flags & 0x100 != 0
        } else {
            read.flags & (0x100 | 0x800 | 0x400) != 0
        };
        if filtered || read.flags & 0x4 != 0 || read.mapping_quality < min_mapping_quality {
            continue;
        }
        let Some(contig) = read.reference.as_deref() else {
            continue;
        };
        let mut position = read.start;
        let mut offset = 0usize;
        for &(length, op) in &read.cigar {
            match op {
                'M' | '=' | 'X' => {
                    for step in 0..length {
                        let q = read.qualities.get(offset + step).copied();
                        if read.qualities.is_empty() || q.is_some_and(|q| q >= min_base_quality) {
                            if let Some(&at) = index.get(&(contig, position + step as i32)) {
                                piles[at].push(Observed {
                                    read,
                                    offset: offset + step,
                                });
                            }
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
    }
    piles
}

/// `FingerprintChecker.loadFingerprints` on a VCF: one fingerprint per sample, keyed by sample,
/// in the `HashMap` order of the names, with every header sample present even when it has no
/// genotypes (`computeIfAbsent`).
///
/// With an index beside the file the records are read SNP by SNP in `Snp` order, the first record
/// at each position; without one, in file order. A record whose alleles do not match its SNP, or
/// that lacks a sample, is skipped from the sample it failed on, as the caught exception does.
pub fn load_vcf_fingerprints(
    path: &str,
    map: &HaplotypeMap,
    specific_sample: Option<&str>,
    genotyping_error_rate: f64,
) -> Result<Vec<(String, Fingerprint)>, String> {
    use crate::vcf_io::read_path;
    let vcf = read_path(path)?;
    check_dictionary(&vcf_dictionary(&vcf.file.header), map)?;
    let indexed = std::path::Path::new(&format!("{path}.idx")).exists()
        || std::path::Path::new(&format!("{path}.tbi")).exists();
    let records: Vec<&htsjdk_vcf::variant::VariantContext> = if indexed {
        let mut snps: Vec<&Snp> = map.blocks.iter().flat_map(|b| b.snps.iter()).collect();
        snps.sort_by(|a, b| a.compare(b));
        snps.iter()
            .filter_map(|snp| {
                vcf.records.iter().map(|r| &r.variant).find(|v| {
                    v.contig == snp.chrom
                        && v.start <= i64::from(snp.pos)
                        && v.stop >= i64::from(snp.pos)
                })
            })
            .collect()
    } else {
        vcf.records.iter().map(|r| &r.variant).collect()
    };
    let mut fingerprints: Vec<(String, Fingerprint)> = Vec::new();
    let source = Some(path.to_string());
    let mut started = false;
    for ctx in records {
        if !started {
            started = true;
            let samples: Vec<String> = match specific_sample {
                Some(s) => vec![s.to_string()],
                None => vcf.file.header.samples.clone(),
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
        }
        add_from_vc(map, &mut fingerprints, ctx, genotyping_error_rate);
    }
    for s in &vcf.file.header.samples {
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
fn add_from_vc(
    map: &HaplotypeMap,
    fingerprints: &mut [(String, Fingerprint)],
    ctx: &htsjdk_vcf::variant::VariantContext,
    error_rate: f64,
) {
    let Some((block, snp)) = map.at(&ctx.contig, ctx.start as i32) else {
        return;
    };
    let snp = snp.clone();
    // `subsetVCToMatchSnp`, for a biallelic site.
    if ctx.is_filtered() || ctx.alleles[0].len() != 1 {
        return;
    }
    let reference = ctx.alleles[0].base_string().as_bytes()[0].to_ascii_uppercase();
    let ref_allele = [snp.allele1, snp.allele2]
        .into_iter()
        .find(|b| b.to_ascii_uppercase() == reference);
    let Some(ref_allele) = ref_allele else { return };
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
        return;
    }
    for allele in &ctx.alleles {
        let bases = allele.base_string();
        let b = bases.as_bytes();
        if b.len() > 1 || (b[0] != snp.allele1 && b[0] != snp.allele2) {
            return;
        }
    }
    let first = ctx.alleles[0].base_string().as_bytes()[0];
    let second = ctx.alleles[1].base_string().as_bytes()[0];
    for (sample, fp) in fingerprints.iter_mut() {
        let Some(genotype) = ctx.genotypes.iter().find(|g| &g.sample_name == sample) else {
            return;
        };
        if let Some(pl) = &genotype.pl {
            if pl.len() != 3 {
                return;
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
                return;
            };
            let probs = Probs::with(
                map,
                block,
                Evidence::GenotypeLikelihoods {
                    ll: normalized_log_likelihoods(ll),
                },
            );
            fp.add(map, probs);
        } else {
            if genotype.is_no_call() {
                continue;
            }
            if fp.blocks.contains_key(&map.blocks[block].key()) {
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
            fp.add(
                map,
                Probs::with(
                    map,
                    block,
                    Evidence::Genotype {
                        likelihoods,
                        snp: snp.clone(),
                    },
                ),
            );
        }
    }
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
    use htsjdk_bam::tag::{Tag, TagValue};
    let (header, records) = crate::metrics_cli::read_input(path);
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
    let reads: Vec<PileupRead> = records
        .iter()
        .map(|r| PileupRead {
            name: r.read_name.clone(),
            flags: r.flags,
            reference: (r.reference_index >= 0)
                .then(|| {
                    header
                        .sequences
                        .get(r.reference_index as usize)
                        .map(|s| s.name.clone())
                })
                .flatten(),
            start: r.alignment_start,
            mapping_quality: r.mapping_quality,
            cigar: r
                .cigar
                .elements
                .iter()
                .map(|e| (e.length as usize, char::from(e.op.to_char())))
                .collect(),
            bases: r.read_bases.clone(),
            qualities: r.base_qualities.clone(),
            read_group: match r.tags.get(Tag::new(b"RG")) {
                Some(TagValue::Str(id)) if header.read_groups.iter().any(|g| g.id == *id) => {
                    Some(id.clone())
                }
                _ => None,
            },
        })
        .collect();
    let mut loci = map.loci();
    let order = |c: &str| dictionary.iter().position(|d| d.name == c);
    loci.sort_by(|a, b| order(&a.0).cmp(&order(&b.0)).then(a.1.cmp(&b.1)));
    let piles = pileup(
        &reads,
        &loci,
        options.min_mapping_quality,
        options.min_base_quality,
        options.allow_duplicates,
    );
    let mut used: std::collections::HashSet<&str> = std::collections::HashSet::new();
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
            let rg = rec.read.read_group.clone();
            if rg.is_none() && !groups.iter().any(|(g, _, _)| g.is_none()) {
                if options.strict {
                    return Err(format!(
                        "picard.PicardException: Found read with no readgroup: {} in file: {path}",
                        rec.read.name
                    ));
                }
                let mut details = IdDetails::new(Some("<UNKNOWN>.0.ZZZ"), &uri);
                details.library = Some("<UNKNOWN>".to_string());
                details.sample = Some(options.default_sample.clone());
                let fp = new_fp(&details);
                groups.push((None, details, fp));
            }
            found = true;
            let name = rec.read.name.as_str();
            if used.insert(name) {
                let at = groups.iter().position(|(g, _, _)| *g == rg).expect("known");
                let key = map.blocks[block].key();
                let probs = groups[at].2.blocks.get_mut(&key).expect("every block");
                let base = rec.read.bases[rec.offset].to_ascii_uppercase();
                let qual = rec.read.qualities.get(rec.offset).copied().unwrap_or(0);
                add(probs, snp, base, qual);
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
}
