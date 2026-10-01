//! `CollectIndependentReplicateMetrics` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.analysis.replicates.CollectIndependentReplicateMetrics.doWork` at tag 3.4.0,
//! with the set classification and the barcode rules of
//! `picard_analysis::collect_independent_replicate_metrics` and the duplicate sets of
//! `picard_analysis::duplicate_set` (htsjdk's `DuplicateSetIterator`, not pre-sorted):
//!
//! * the checks in `doWork`'s order: the index, the two dictionaries, then the sample;
//! * the heterozygous sites: SNPs that pass their filters, whose genotype for the sample is over
//!   `MINIMUM_GQ` and heterozygous, keyed by `QueryInterval` in a `TreeMap`, each with the
//!   genotype's alleles in GT order;
//! * the records a BAM query over those sites returns, filtered (aligned, primary, mapping
//!   quality, paired with a mapped mate) BEFORE the sets are cut;
//! * the walk over the sets, locus by locus, with every quirk of the reference kept: a set whose
//!   read misses the site or is under `MINIMUM_BQ` there is abandoned after its size was counted,
//!   a third allele discards the whole locus, the allele-balance histogram is fed only when a
//!   locus is closed by the next one (not by the end of the run), and the end of the run counts a
//!   three-allele site whenever no locus is pending;
//! * `calculateDerivedFields`, including its `nDuplicateSets - nExactlyDouble - nExactlyDouble`.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use htsjdk_bam::alignment_block::alignment_blocks;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram, MetricBean, MetricsFile, Value};
use picard_analysis::collect_independent_replicate_metrics::{
    classify_set, edit_distance, SetClassification,
};
use picard_analysis::duplicate_set::duplicate_sets;
use picard_analysis::mark_duplicates::{Options, Record as SetRecord};
use picard_analysis::metrics_cli::{absolute, fail, read_input, thrown, Args};

const TOOL: &str = "CollectIndependentReplicateMetrics";

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST: u16 = 0x40;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// Every column, counts then rates: `getFields()` order.
const COLUMNS: [&str; 35] = [
    "nSites",
    "nThreeAllelesSites",
    "nTotalReads",
    "nDuplicateSets",
    "nExactlyTriple",
    "nExactlyDouble",
    "nReadsInBigSets",
    "nDifferentAllelesBiDups",
    "nReferenceAllelesBiDups",
    "nAlternateAllelesBiDups",
    "nDifferentAllelesTriDups",
    "nMismatchingAllelesBiDups",
    "nReferenceAllelesTriDups",
    "nAlternateAllelesTriDups",
    "nMismatchingAllelesTriDups",
    "nReferenceReads",
    "nAlternateReads",
    "nMismatchingUMIsInDiffBiDups",
    "nMatchingUMIsInDiffBiDups",
    "nMismatchingUMIsInSameBiDups",
    "nMatchingUMIsInSameBiDups",
    "nMismatchingUMIsInCoOrientedBiDups",
    "nMismatchingUMIsInContraOrientedBiDups",
    "nBadBarcodes",
    "nGoodBarcodes",
    "biSiteHeterogeneityRate",
    "triSiteHeterogeneityRate",
    "biSiteHomogeneityRate",
    "triSiteHomogeneityRate",
    "independentReplicationRateFromBiDups",
    "independentReplicationRateFromTriDups",
    "pSameUmiInIndependentBiDup",
    "pSameAlleleWhenMismatchingUmi",
    "independentReplicationRateFromUmi",
    "replicationRateFromReplicateSets",
];

/// The counters by name, spelled as the Java fields are so each line reads like its original.
#[derive(Clone, Copy)]
#[allow(non_camel_case_types, clippy::enum_variant_names)]
enum C {
    nSites,
    nThreeAllelesSites,
    nTotalReads,
    nDuplicateSets,
    nExactlyTriple,
    nExactlyDouble,
    nReadsInBigSets,
    nDifferentAllelesBiDups,
    nReferenceAllelesBiDups,
    nAlternateAllelesBiDups,
    nDifferentAllelesTriDups,
    nMismatchingAllelesBiDups,
    nReferenceAllelesTriDups,
    nAlternateAllelesTriDups,
    nMismatchingAllelesTriDups,
    nReferenceReads,
    nAlternateReads,
    nMismatchingUMIsInDiffBiDups,
    nMatchingUMIsInDiffBiDups,
    nMismatchingUMIsInSameBiDups,
    nMatchingUMIsInSameBiDups,
    nMismatchingUMIsInCoOrientedBiDups,
    nMismatchingUMIsInContraOrientedBiDups,
    nBadBarcodes,
    nGoodBarcodes,
}

#[derive(Clone, Default)]
struct Metric {
    counts: [i64; 25],
    rates: [f64; 10],
}

impl Metric {
    fn get(&self, c: C) -> i64 {
        self.counts[c as usize]
    }
    fn add(&mut self, c: C, by: i64) {
        self.counts[c as usize] += by;
    }
    /// `MergeableMetricBase.merge`: every `@MergeByAdding` field summed.
    fn merge(&mut self, other: &Metric) {
        for (mine, theirs) in self.counts.iter_mut().zip(other.counts.iter()) {
            *mine += theirs;
        }
    }
    fn calculate_derived_fields(&mut self) {
        let d = |v: i64| v as f64;
        let bi_het = d(self.get(C::nDifferentAllelesBiDups))
            / d(self.get(C::nDifferentAllelesBiDups)
                + self.get(C::nAlternateAllelesBiDups)
                + self.get(C::nReferenceAllelesBiDups));
        let bi_hom = 1.0 - bi_het;
        let tri_het = d(self.get(C::nDifferentAllelesTriDups))
            / d(self.get(C::nDifferentAllelesTriDups)
                + self.get(C::nAlternateAllelesTriDups)
                + self.get(C::nReferenceAllelesTriDups));
        let tri_hom = 1.0 - tri_het;
        let p_same_umi = d(self.get(C::nMatchingUMIsInDiffBiDups))
            / d(self.get(C::nMismatchingUMIsInDiffBiDups) + self.get(C::nMatchingUMIsInDiffBiDups));
        let p_same_allele = d(self.get(C::nMismatchingUMIsInSameBiDups))
            / d(self.get(C::nMismatchingUMIsInSameBiDups)
                + self.get(C::nMismatchingUMIsInDiffBiDups));
        let from_umi =
            d(self.get(C::nMismatchingUMIsInDiffBiDups)
                + self.get(C::nMismatchingUMIsInSameBiDups))
                / d(self.get(C::nExactlyDouble));
        // The reference subtracts the doubletons twice; kept.
        let big_sets =
            self.get(C::nDuplicateSets) - self.get(C::nExactlyDouble) - self.get(C::nExactlyDouble);
        let from_sets = d(self.get(C::nExactlyDouble)
            + self.get(C::nExactlyTriple) * 2
            + self.get(C::nReadsInBigSets)
            - big_sets)
            / d(self.get(C::nTotalReads));
        self.rates = [
            bi_het,
            tri_het,
            bi_hom,
            tri_hom,
            2.0 * bi_het,
            2.0 * (1.0 - tri_hom.sqrt()),
            p_same_umi,
            p_same_allele,
            from_umi,
            from_sets,
        ];
    }
}

impl MetricBean for Metric {
    fn class_name(&self) -> &str {
        "picard.analysis.replicates.IndependentReplicateMetric"
    }
    fn columns(&self) -> &[&'static str] {
        &COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        self.counts
            .iter()
            .map(|v| Value::Long(*v))
            .chain(self.rates.iter().map(|v| Value::Double(*v)))
            .collect()
    }
}

/// `QueryInterval`: a contig index and a closed range.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct QueryInterval {
    reference: i32,
    start: i32,
    end: i32,
}

impl Ord for QueryInterval {
    /// `QueryInterval.compareTo`.
    fn cmp(&self, other: &Self) -> Ordering {
        let comp = self.reference - other.reference;
        if comp != 0 {
            return comp.cmp(&0);
        }
        let comp = self.start - other.start;
        if comp != 0 {
            return comp.cmp(&0);
        }
        if self.end == other.end {
            Ordering::Equal
        } else if self.end == 0 {
            Ordering::Greater
        } else if other.end == 0 {
            Ordering::Less
        } else {
            (self.end - other.end).cmp(&0)
        }
    }
}

impl PartialOrd for QueryInterval {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl QueryInterval {
    fn overlaps(&self, other: &QueryInterval) -> bool {
        if self.reference != other.reference {
            return false;
        }
        let this_end = if self.end == 0 { i32::MAX } else { self.end };
        let other_end = if other.end == 0 { i32::MAX } else { other.end };
        self.start <= other_end && other.start <= this_end
    }

    /// `isCleanlyBefore`.
    fn cleanly_before(&self, other: &QueryInterval) -> bool {
        !self.overlaps(other) && self.cmp(other) == Ordering::Less
    }
}

/// One VCF data line, as far as the site filters read it.
struct Site {
    contig: String,
    start: i32,
    end: i32,
    filtered: bool,
    snp: bool,
    alleles: Vec<String>,
    /// Per sample: the GT allele indices (`None` for a no-call) and the GQ (`-1` when missing).
    genotypes: Vec<(Vec<Option<usize>>, i32)>,
}

struct Vcf {
    contigs: Vec<(String, i64)>,
    samples: Vec<String>,
    sites: Vec<Site>,
}

fn read_vcf(path: &str) -> Vcf {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{e}")));
    let mut contigs = Vec::new();
    let mut samples = Vec::new();
    let mut sites = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("##contig=<") {
            let rest = rest.trim_end_matches('>');
            let mut id = String::new();
            let mut length = 0i64;
            for field in rest.split(',') {
                if let Some(v) = field.strip_prefix("ID=") {
                    id = v.to_string();
                } else if let Some(v) = field.strip_prefix("length=") {
                    length = v.parse().unwrap_or(0);
                }
            }
            contigs.push((id, length));
            continue;
        }
        if line.starts_with("##") {
            continue;
        }
        if let Some(rest) = line.strip_prefix('#') {
            samples = rest.split('\t').skip(9).map(str::to_string).collect();
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let start: i32 = f[1].parse().unwrap_or(0);
        let reference = f[3].to_string();
        let mut alleles = vec![reference.clone()];
        if f[4] != "." {
            alleles.extend(f[4].split(',').map(str::to_string));
        }
        let snp = alleles.len() > 1
            && alleles
                .iter()
                .all(|a| a.len() == 1 && a.bytes().all(|b| b.is_ascii_alphabetic()) && a != "*");
        let filtered = !(f[6] == "PASS" || f[6] == ".");
        let keys: Vec<&str> = f.get(8).map(|k| k.split(':').collect()).unwrap_or_default();
        let gt_at = keys.iter().position(|k| *k == "GT");
        let gq_at = keys.iter().position(|k| *k == "GQ");
        let genotypes = f
            .iter()
            .skip(9)
            .map(|column| {
                let values: Vec<&str> = column.split(':').collect();
                let gt = gt_at
                    .and_then(|i| values.get(i))
                    .map(|gt| {
                        gt.split(['/', '|'])
                            .map(|a| a.parse::<usize>().ok())
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default();
                let gq = gq_at
                    .and_then(|i| values.get(i))
                    .and_then(|v| v.parse::<i32>().ok())
                    .unwrap_or(-1);
                (gt, gq)
            })
            .collect();
        sites.push(Site {
            contig: f[0].to_string(),
            start,
            end: start + reference.len() as i32 - 1,
            filtered,
            snp,
            alleles,
            genotypes,
        });
    }
    Vcf {
        contigs,
        samples,
        sites,
    }
}

fn string_tag(record: &BamRecord, name: &str) -> Option<String> {
    let bytes = name.as_bytes();
    if bytes.len() != 2 {
        return None;
    }
    match record.tags.get(Tag::new(&[bytes[0], bytes[1]])) {
        Some(TagValue::Str(value)) => Some(value.clone()),
        _ => None,
    }
}

/// `SAMRecord.getReadPositionAtReferencePosition(pos)`: one-based, or 0 off the aligned bases.
fn read_position_at(record: &BamRecord, position: i32) -> i32 {
    for block in alignment_blocks(&record.cigar, record.alignment_start) {
        if position >= block.reference_start && position < block.reference_start + block.length {
            return block.read_start + position - block.reference_start;
        }
    }
    0
}

/// `SamReader.hasIndex()`: a `.bai` beside the BAM, under either of the names htsjdk looks for.
fn has_index(path: &str, raw_is_bam: bool) -> bool {
    if !raw_is_bam {
        return false;
    }
    let p = std::path::Path::new(path);
    let replaced = p.with_extension("bai");
    let appended = format!("{path}.bai");
    replaced.exists() || std::path::Path::new(&appended).exists()
}

fn byte_histogram(bin: &str, value: &str, counts: &BTreeMap<i64, f64>) -> Histogram {
    Histogram {
        bin_label: bin.to_string(),
        value_label: value.to_string(),
        key_class: "java.lang.Byte".to_string(),
        bins: counts.iter().map(|(k, v)| (k.to_string(), *v)).collect(),
    }
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("MO", "MATRIX_OUTPUT"),
        ("V", "VCF"),
        ("GQ", "MINIMUM_GQ"),
        ("MQ", "MINIMUM_MQ"),
        ("BQ", "MINIMUM_BQ"),
        ("ALIAS", "SAMPLE"),
        ("MBQ", "MINIMUM_BARCODE_BQ"),
        ("FUR", "FILTER_UNPAIRED_READS"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let matrix_output = args.get("MATRIX_OUTPUT").map(str::to_string);
    let vcf_path = args.required("VCF");
    let minimum_gq = args.int("MINIMUM_GQ", 90) as i32;
    let minimum_mq = args.int("MINIMUM_MQ", 40) as i32;
    let minimum_bq = args.int("MINIMUM_BQ", 17) as i32;
    let sample_arg = args.get("SAMPLE").map(str::to_string);
    let stop_after = args.int("STOP_AFTER", 0);
    let barcode_tag = args.get("BARCODE_TAG").unwrap_or("RX").to_string();
    let barcode_bq = args.get("BARCODE_BQ").unwrap_or("QX").to_string();
    let minimum_barcode_bq = args.int("MINIMUM_BARCODE_BQ", 30) as i32;
    let filter_unpaired = args.bool("FILTER_UNPAIRED_READS", true);

    // doWork's checks, in its order.
    let raw_is_bam = std::fs::File::open(&input)
        .ok()
        .and_then(|mut f| {
            let mut magic = [0u8; 2];
            std::io::Read::read_exact(&mut f, &mut magic).ok()?;
            Some(magic == [0x1f, 0x8b])
        })
        .unwrap_or(false);
    let (header, records) = read_input(&input);
    if !has_index(&input, raw_is_bam) {
        thrown("picard.PicardException: INPUT file must have an index.");
    }
    let vcf = read_vcf(&vcf_path);
    if !vcf.contigs.is_empty() {
        let ours: Vec<(String, i64)> = header
            .sequences
            .iter()
            .map(|s| (s.name.clone(), i64::from(s.length)))
            .collect();
        if ours != vcf.contigs {
            thrown(&format!(
                "picard.PicardException: Sequence dictionary for ({}) does not match sequence \
                 dictionary for ({})",
                absolute(&input),
                absolute(&vcf_path)
            ));
        }
    }
    let sample = match sample_arg {
        None => {
            if vcf.samples.len() != 1 {
                thrown(&format!(
                    "java.lang.IllegalArgumentException: When sample is null, VCF must have \
                     exactly 1 sample. found {}",
                    vcf.samples.len()
                ));
            }
            vcf.samples[0].clone()
        }
        Some(s) => {
            if !vcf.samples.contains(&s) {
                thrown(&format!(
                    "java.lang.IllegalArgumentException: When sample is not null, VCF must \
                     contain supplied sample. Cannot find sample {s} in vcf."
                ));
            }
            s
        }
    };
    let sample_index = vcf
        .samples
        .iter()
        .position(|s| *s == sample)
        .expect("sample");

    // getQueryIntervalsMap: SNP, passing, GQ, het, in that order; a TreeMap on the interval.
    let mut sites: BTreeMap<QueryInterval, Vec<String>> = BTreeMap::new();
    for site in &vcf.sites {
        if !site.snp || site.filtered {
            continue;
        }
        let (gt, gq) = &site.genotypes[sample_index];
        if *gq < minimum_gq {
            continue;
        }
        let het = gt.len() == 2 && gt.iter().all(Option::is_some) && gt[0] != gt[1];
        if !het {
            continue;
        }
        let reference = vcf
            .contigs
            .iter()
            .position(|(name, _)| *name == site.contig)
            .map(|i| i as i32)
            .unwrap_or_else(|| thrown("java.lang.NullPointerException"));
        let alleles = gt
            .iter()
            .map(|a| site.alleles[a.expect("called")].clone())
            .collect();
        sites.insert(
            QueryInterval {
                reference,
                start: site.start,
                end: site.end,
            },
            alleles,
        );
    }

    // The query, then the filters, then the sets.
    let kept: Vec<&BamRecord> = records
        .iter()
        .filter(|r| {
            if r.reference_index < 0 {
                return false;
            }
            let end = if r.flags & UNMAPPED != 0 {
                r.alignment_start
            } else {
                r.alignment_end()
            };
            sites.keys().any(|q| {
                q.reference == r.reference_index && r.alignment_start <= q.end && end >= q.start
            })
        })
        .filter(|r| r.flags & UNMAPPED == 0)
        .filter(|r| r.flags & (SECONDARY | SUPPLEMENTARY) == 0)
        .filter(|r| i32::from(r.mapping_quality) >= minimum_mq)
        // `CountingPairedFilter`: unpaired, or paired with an unmapped mate.
        .filter(|r| !filter_unpaired || (r.flags & PAIRED != 0 && r.flags & MATE_UNMAPPED == 0))
        .collect();
    let set_records: Vec<SetRecord> = kept
        .iter()
        .map(|record| {
            let group = string_tag(record, "RG").and_then(|id| {
                header
                    .read_groups
                    .iter()
                    .position(|g| g.id == id)
                    .map(|p| (p, &header.read_groups[p]))
            });
            let (library, read_group) = match group {
                Some((p, g)) => (
                    g.attributes
                        .get("LB")
                        .unwrap_or("Unknown Library")
                        .to_string(),
                    p as i32,
                ),
                None => ("Unknown Library".to_string(), -1),
            };
            SetRecord {
                name: record.read_name.clone(),
                flags: record.flags,
                reference_index: record.reference_index,
                alignment_start: record.alignment_start,
                cigar: record.cigar.clone(),
                qualities: record.base_qualities.clone(),
                mate_reference_index: record.mate_reference_index,
                library,
                read_group,
                barcode: None,
                existing_dt: None,
                mate_cigar: match record.tags.get(Tag::new(b"MC")) {
                    Some(TagValue::Str(text)) => htsjdk_bam::text_parse::parse_cigar(text).ok(),
                    _ => None,
                },
                mate_alignment_start: record.mate_alignment_start,
            }
        })
        .collect();
    let sets = duplicate_sets(&set_records, &Options::default());

    let mut confusion: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut confusion_distance: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut metric = Metric::default();
    let mut diff_distances: BTreeMap<i64, f64> = BTreeMap::new();
    let mut same_distances: BTreeMap<i64, f64> = BTreeMap::new();
    let mut allele_balance: BTreeMap<i64, f64> = BTreeMap::new();
    let mut intervals = sites.keys().copied().collect::<Vec<_>>().into_iter();
    let mut query: Option<QueryInterval> = None;
    let mut locus = Metric::default();
    let mut use_locus = true;
    let mut new_locus = false;
    let mut examined = 0i64;

    'set: for set in &sets {
        let reads: Vec<&BamRecord> = set.iter().map(|i| kept[*i]).collect();
        let rep = reads[0];
        let rep_interval = QueryInterval {
            reference: rep.reference_index,
            start: rep.alignment_start,
            end: rep.alignment_end(),
        };
        examined += 1;
        if !use_locus || query.is_some_and(|q| q.cleanly_before(&rep_interval)) {
            if !use_locus {
                metric.add(C::nThreeAllelesSites, 1);
            }
            query = None;
        }
        loop {
            let pending = intervals.len() > 0;
            if !(pending
                && (query.is_none() || query.is_some_and(|q| q.cleanly_before(&rep_interval))))
            {
                break;
            }
            if locus.get(C::nReferenceReads) == 0 || locus.get(C::nAlternateReads) == 0 {
                use_locus = false;
            }
            if use_locus && new_locus {
                metric.merge(&locus);
                let alt = locus.get(C::nAlternateReads) as f64;
                let reference = locus.get(C::nReferenceReads) as f64;
                let balance = (100.0 * (alt + 0.5) / (alt + reference + 1.0) + 0.5).floor() as i64;
                *allele_balance.entry(balance as i8 as i64).or_insert(0.0) += 1.0;
                new_locus = false;
            }
            query = intervals.next();
            locus = Metric::default();
            locus.add(C::nSites, 1);
            use_locus = true;
        }
        new_locus = true;
        let Some(q) = query else {
            break;
        };
        let size = reads.len() as i64;
        locus.add(C::nTotalReads, size);
        if size > 1 {
            locus.add(C::nDuplicateSets, 1);
        }
        if size == 2 {
            locus.add(C::nExactlyDouble, 1);
        } else if size == 3 {
            locus.add(C::nExactlyTriple, 1);
        } else if size > 3 {
            locus.add(C::nReadsInBigSets, size);
        }
        let alleles = &sites[&q];
        let (mut n_ref, mut n_alt, mut n_other) = (0usize, 0usize, 0usize);
        for read in &reads {
            let offset = read_position_at(read, q.start) - 1;
            if offset == -1 {
                continue 'set;
            }
            if i32::from(read.base_qualities[offset as usize]) <= minimum_bq {
                continue 'set;
            }
            let base = read.read_bases[offset as usize].to_ascii_uppercase();
            let base = (base as char).to_string();
            if alleles[0] == base {
                n_ref += 1;
            } else if alleles[1] == base {
                n_alt += 1;
            } else {
                n_other += 1;
                use_locus = false;
            }
        }
        locus.add(C::nAlternateReads, n_alt as i64);
        locus.add(C::nReferenceReads, n_ref as i64);
        if size == 1 || size > 3 {
            continue;
        }
        let classification = classify_set(n_ref, n_alt, n_other);
        if size == 2 {
            let use_barcodes = !reads.iter().any(|r| {
                string_tag(r, &barcode_bq)
                    .unwrap_or_default()
                    .bytes()
                    .any(|q| i32::from(q.wrapping_sub(33) as i8) < minimum_barcode_bq)
            });
            if use_barcodes {
                locus.add(C::nGoodBarcodes, 1);
            } else {
                locus.add(C::nBadBarcodes, 1);
            }
            let barcodes: Vec<String> = reads
                .iter()
                .map(|r| string_tag(r, &barcode_tag).unwrap_or_default())
                .collect();
            let orientations: Vec<bool> = reads
                .iter()
                .map(|r| r.flags & PAIRED == 0 || r.flags & FIRST != 0)
                .collect();
            let multiple_orientations = orientations.iter().any(|o| *o != orientations[0]);
            let distance = edit_distance(&barcodes[0], &barcodes[1]).unwrap_or_else(|| {
                thrown(&format!(
                    "java.lang.IllegalArgumentException: lengths of strings must equal, found \
                     '{}' and '{}'.",
                    barcodes[0], barcodes[1]
                ))
            });
            let distance = i64::from(distance as i8);
            if use_barcodes && distance != 0 {
                if multiple_orientations {
                    locus.add(C::nMismatchingUMIsInContraOrientedBiDups, 1);
                } else {
                    locus.add(C::nMismatchingUMIsInCoOrientedBiDups, 1);
                }
            }
            match classification {
                SetClassification::DifferentAlleles => {
                    locus.add(C::nDifferentAllelesBiDups, 1);
                    if use_barcodes {
                        *diff_distances.entry(distance).or_insert(0.0) += 1.0;
                        if distance == 0 {
                            locus.add(C::nMatchingUMIsInDiffBiDups, 1);
                        } else {
                            locus.add(C::nMismatchingUMIsInDiffBiDups, 1);
                        }
                    }
                }
                SetClassification::MismatchingAllele => {
                    locus.add(C::nMismatchingAllelesBiDups, 1);
                }
                other => {
                    if other == SetClassification::AlternateAllele {
                        locus.add(C::nAlternateAllelesBiDups, 1);
                    } else {
                        locus.add(C::nReferenceAllelesBiDups, 1);
                    }
                    if use_barcodes {
                        *same_distances.entry(distance).or_insert(0.0) += 1.0;
                        let key = (barcodes[0].clone(), barcodes[1].clone());
                        *confusion.entry(key.clone()).or_insert(0.0) += 1.0;
                        confusion_distance.entry(key).or_insert(distance as f64);
                        if distance == 0 {
                            locus.add(C::nMatchingUMIsInSameBiDups, 1);
                        } else {
                            locus.add(C::nMismatchingUMIsInSameBiDups, 1);
                        }
                    }
                }
            }
        }
        if size == 3 {
            locus.add(
                match classification {
                    SetClassification::MismatchingAllele => C::nMismatchingAllelesTriDups,
                    SetClassification::DifferentAlleles => C::nDifferentAllelesTriDups,
                    SetClassification::AlternateAllele => C::nAlternateAllelesTriDups,
                    SetClassification::ReferenceAllele => C::nReferenceAllelesTriDups,
                },
                1,
            );
        }
        if stop_after > 0 && examined > stop_after {
            break;
        }
    }
    if use_locus && new_locus {
        metric.merge(&locus);
    } else {
        metric.add(C::nThreeAllelesSites, 1);
    }

    metric.calculate_derived_fields();
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    file.add_metric(&metric);
    file.histograms.push(byte_histogram(
        "alleleBalance",
        "alleleBalanceCount",
        &allele_balance,
    ));
    file.histograms.push(byte_histogram(
        "editDistance",
        "diffAllelesCount",
        &diff_distances,
    ));
    file.histograms.push(byte_histogram(
        "editDistance",
        "sameAllelesCount",
        &same_distances,
    ));
    write(&output, &file);

    if let Some(matrix) = matrix_output {
        let tuple = |m: &BTreeMap<(String, String), f64>, value: &str| Histogram {
            bin_label: "ConfusionUMI".to_string(),
            value_label: value.to_string(),
            key_class: "htsjdk.samtools.util.ComparableTuple".to_string(),
            bins: m
                .iter()
                .map(|((a, b), v)| (format!("[{a}, {b}]"), *v))
                .collect(),
        };
        let mut file = MetricsFile::new();
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
        file.histograms.push(tuple(&confusion, "Count"));
        file.histograms
            .push(tuple(&confusion_distance, "EditDistance"));
        write(&matrix, &file);
    }
}

fn write(path: &str, file: &MetricsFile) {
    if let Err(e) = std::fs::write(path, file.write()) {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Could not write to file {path}: {e}"
        ));
    }
}
