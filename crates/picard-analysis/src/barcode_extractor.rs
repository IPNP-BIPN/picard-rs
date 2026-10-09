//! `BarcodeExtractor` and the barcode half of `ExtractBarcodesProgram`: which declared barcode a
//! cluster's barcode reads are, under `DistanceMetric`'s three distances, and the `BarcodeMetric`
//! rows the match counts into.
//!
//! Ported from `picard.illumina.BarcodeExtractor`, `ExtractBarcodesProgram`, `DistanceMetric`,
//! `BarcodeMetric` and `picard.util.SingleBarcodeDistanceMetric` at tag 3.4.0, for
//! `ExtractIlluminaBarcodes` and the inline matching of `IlluminaBasecallsToFastq` and
//! `IlluminaBasecallsToSam`.
//!
//! The declared barcodes are walked in the order of a `HashSet` sized to their count, so which of
//! two equidistant barcodes is "best" (and written lower-cased when neither matches) is the
//! bucket order of their hashes. The extractor caches a match by its bases: every declared barcode
//! and the all-`N` one are computed up front as inline matches against perfect qualities, so a
//! cluster whose bases ARE a declared barcode gets that answer, whatever a fresh match would say.

use std::collections::HashMap;

use htsjdk_metrics::file::{MetricBean, Value};

use crate::fingerprint::{java_hash_order, reference_path};
use crate::java_hash_map::string_hash_code;
use crate::metrics_cli::thrown;

const COLUMNS: &[&str] = &[
    "BARCODE",
    "BARCODE_WITHOUT_DELIMITER",
    "BARCODE_NAME",
    "LIBRARY_NAME",
    "READS",
    "PF_READS",
    "PERFECT_MATCHES",
    "PF_PERFECT_MATCHES",
    "ONE_MISMATCH_MATCHES",
    "PF_ONE_MISMATCH_MATCHES",
    "PCT_MATCHES",
    "RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT",
    "PF_PCT_MATCHES",
    "PF_RATIO_THIS_BARCODE_TO_BEST_BARCODE_PCT",
    "PF_NORMALIZED_MATCHES",
];

/// `BarcodeMetric`.
#[derive(Clone, Default)]
pub struct Metric {
    pub barcode: String,
    pub name: Option<String>,
    pub library: Option<String>,
    pub bytes: Vec<Vec<u8>>,
    pub reads: i64,
    pub pf_reads: i64,
    pub perfect: i64,
    pub pf_perfect: i64,
    pub one: i64,
    pub pf_one: i64,
    pub pct: f64,
    pub ratio: f64,
    pub pf_pct: f64,
    pub pf_ratio: f64,
    pub pf_normalized: f64,
}

impl Metric {
    pub fn new(name: Option<String>, library: Option<String>, seqs: &[String]) -> Metric {
        Metric {
            barcode: seqs.join("-"),
            name,
            library,
            bytes: seqs.iter().map(|s| s.as_bytes().to_vec()).collect(),
            ..Metric::default()
        }
    }
}

impl MetricBean for Metric {
    fn class_name(&self) -> &str {
        "picard.illumina.BarcodeMetric"
    }
    fn columns(&self) -> &[&'static str] {
        COLUMNS
    }
    fn values(&self) -> Vec<Value> {
        let text = |s: &Option<String>| s.clone().map_or(Value::Null, Value::Str);
        vec![
            Value::Str(self.barcode.clone()),
            Value::Str(self.barcode.replace('-', "")),
            text(&self.name),
            text(&self.library),
            Value::Long(self.reads),
            Value::Long(self.pf_reads),
            Value::Long(self.perfect),
            Value::Long(self.pf_perfect),
            Value::Long(self.one),
            Value::Long(self.pf_one),
            Value::Double(self.pct),
            Value::Double(self.ratio),
            Value::Double(self.pf_pct),
            Value::Double(self.pf_ratio),
            Value::Double(self.pf_normalized),
        ]
    }
}

#[derive(Clone, Copy, PartialEq)]
pub enum Distance {
    Hamming,
    Lenient,
    Free,
}

/// `BarcodeExtractor.BarcodeMatch`.
#[derive(Clone, Default)]
pub struct Match {
    pub matched: bool,
    pub barcode: String,
    pub mismatches: i32,
    pub to_second: i32,
}

pub fn is_no_call(b: u8) -> bool {
    matches!(b, b'N' | b'n' | b'.')
}

/// `SingleBarcodeDistanceMetric`, for one barcode read.
pub struct Single<'a> {
    pub barcode: &'a [u8],
    pub read: &'a [u8],
    pub quals: Option<&'a [u8]>,
    pub masked: Vec<u8>,
    pub minimum: i32,
    pub max: i32,
}

impl<'a> Single<'a> {
    pub fn new(
        barcode: &'a [u8],
        read: &'a [u8],
        quals: Option<&'a [u8]>,
        minimum: i32,
        max: i32,
    ) -> Single<'a> {
        let mut masked = read.to_vec();
        if let Some(q) = quals {
            if q.iter().any(|v| i32::from(*v) < minimum) {
                for (i, v) in q.iter().enumerate() {
                    if i32::from(*v) < minimum {
                        masked[i] = b'.';
                    }
                }
            }
        }
        Single {
            barcode,
            read,
            quals,
            masked,
            minimum,
            max,
        }
    }

    pub fn hamming(&self) -> i32 {
        let mut n = 0;
        let mut i = 0;
        while i < self.barcode.len() && i < self.read.len() && n <= self.max {
            let r = self.read[i];
            // A different base is a mismatch, and so is an equal one below the quality floor.
            let poor = self.quals.is_some_and(|q| i32::from(q[i]) < self.minimum);
            if !is_no_call(r) && (!self.barcode[i].eq_ignore_ascii_case(&r) || poor) {
                n += 1;
            }
            i += 1;
        }
        n
    }

    pub fn lenient(&self) -> i32 {
        let mut n = 0;
        let mut i = 0;
        while i < self.barcode.len() && i < self.masked.len() && n <= self.max {
            let r = self.masked[i];
            if !is_no_call(r) && !self.barcode[i].eq_ignore_ascii_case(&r) {
                n += 1;
            }
            i += 1;
        }
        n
    }

    pub fn free(&self) -> i32 {
        let n = self.barcode.len();
        if n != self.read.len() {
            thrown(&format!(
                "java.lang.IllegalArgumentException: This version of freeDistance is specifically made for comparing strings of equal length. found {n} and {}.",
                self.read.len()
            ));
        }
        if n == 0 {
            return 0;
        }
        let max = self.max;
        let barcode: Vec<u8> = self.barcode.iter().rev().copied().collect();
        let read: Vec<u8> = self.masked.iter().rev().copied().collect();
        let mut previous = vec![0i32; n + 1];
        let mut cost = vec![i32::MAX; n + 1];
        let boundary = (n as i32).min(max) as usize + 1;
        for p in previous.iter_mut().skip(boundary) {
            *p = i32::MAX;
        }
        for j in 1..=n {
            let r = read[j - 1];
            cost[0] = 0;
            let ji = j as i32;
            let low = 1.max(ji - max) as usize;
            let high = if ji > i32::MAX - max {
                n
            } else {
                (n as i32).min(ji + max) as usize
            };
            if low > 1 {
                cost[low - 1] = i32::MAX - max;
            }
            let mut min_cost = i32::MAX;
            for i in low..=high {
                if barcode[i - 1] == r || is_no_call(r) {
                    cost[i] = previous[i - 1];
                } else {
                    let best = previous[i - 1].min(cost[i - 1]).min(previous[i]);
                    cost[i] = best.wrapping_add(1);
                }
                let to_center = (i as i32 - ji).abs();
                min_cost = min_cost.min(cost[i].wrapping_add(to_center));
            }
            if min_cost > max {
                return max + 1;
            }
            std::mem::swap(&mut previous, &mut cost);
        }
        if previous[n] > max {
            return max + 1;
        }
        previous[n]
    }
}

pub struct Extractor {
    /// The declared barcodes in `HashSet` order.
    pub barcodes: Vec<Vec<Vec<u8>>>,
    pub max_no_calls: i32,
    pub max_mismatches: i32,
    pub min_delta: i32,
    pub minimum: i32,
    pub mode: Distance,
    pub cache: HashMap<Vec<Vec<u8>>, Match>,
}

impl Extractor {
    /// `DistanceMetric.distance`.
    pub fn distance(
        &self,
        barcode: &[Vec<u8>],
        read: &[Vec<u8>],
        quals: Option<&[Vec<u8>]>,
        max: i32,
    ) -> i32 {
        let mut n = 0;
        for j in 0..barcode.len() {
            let single = Single::new(
                &barcode[j],
                &read[j],
                quals.map(|q| q[j].as_slice()),
                self.minimum,
                max - n,
            );
            n += match self.mode {
                Distance::Hamming => single.hamming(),
                Distance::Lenient => single.lenient(),
                Distance::Free => single.free(),
            };
            if n > max {
                return n;
            }
        }
        n
    }

    /// `calculateBarcodeMatch`.
    pub fn calculate(&self, read: &[Vec<u8>], quals: Option<&[Vec<u8>]>, inline: bool) -> Match {
        let mut m = Match::default();
        let mut total = 0;
        let mut no_calls = 0;
        for bc in read {
            total += bc.len() as i32;
            for &b in bc {
                if is_no_call(b) {
                    no_calls += 1;
                }
                if inline && no_calls > self.max_no_calls {
                    m.mismatches = total;
                    m.barcode = String::new();
                    m.matched = false;
                    return m;
                }
            }
        }
        let mut best: Option<String> = None;
        let mut best_n = total + 1;
        let mut second_n = total + 1;
        for barcode in &self.barcodes {
            let n = self.distance(
                barcode,
                read,
                quals,
                self.max_mismatches.min(best_n) + self.min_delta,
            );
            if n < best_n {
                if best.is_some() {
                    second_n = best_n;
                }
                best_n = n;
                best = Some(barcode.iter().map(|b| String::from_utf8_lossy(b)).collect());
            } else if n < second_n {
                second_n = n;
            }
        }
        m.matched = best.is_some()
            && no_calls <= self.max_no_calls
            && best_n <= self.max_mismatches
            && second_n - best_n >= self.min_delta;
        m.mismatches = best_n;
        m.to_second = second_n;
        if m.matched {
            m.barcode = best.unwrap_or_default();
        } else {
            match best {
                Some(b) if !inline && no_calls + best_n < total => m.barcode = b.to_lowercase(),
                _ => {
                    m.mismatches = total;
                    m.barcode = String::new();
                }
            }
        }
        m
    }

    /// `findBestBarcode`.
    pub fn find(&mut self, read: &[Vec<u8>], quals: Option<&[Vec<u8>]>, inline: bool) -> Match {
        let lookup =
            quals.is_none_or(|q| q.iter().flatten().all(|v| i32::from(*v) >= self.minimum));
        if !lookup {
            return self.calculate(read, quals, inline);
        }
        let m = match self.cache.get(read) {
            Some(m) => m.clone(),
            None => self.calculate(read, quals, inline),
        };
        if m.matched {
            self.cache.insert(read.to_vec(), m.clone());
        }
        m
    }
}

/// `TabbedTextFileWithHeaderParser` over the barcode file, then `parseInputFile`.
pub fn parse_input_file(
    path: &str,
    lengths: &[usize],
    messages: &mut Vec<String>,
) -> Vec<(String, Metric)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|_| {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
            reference_path(path)
        ))
    });
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap_or("")
        .split('\t')
        .map(str::to_string)
        .collect();
    // `columnLabels()` is the key set of a `HashMap`.
    let labels: Vec<String> = java_hash_order(
        header
            .iter()
            .map(|h| (string_hash_code(h), h.clone()))
            .collect(),
    );
    let mut valid: Vec<String> = labels
        .into_iter()
        .filter(|name| {
            let upper = name.to_uppercase();
            ["BARCODE_SEQUENCE", "BARCODE"]
                .iter()
                .any(|p| upper.starts_with(p) && !name.eq_ignore_ascii_case("barcode_name"))
        })
        .collect();
    if lengths.len() != valid.len() {
        messages.push(format!(
            "Expected {} valid barcode columns, but found {}",
            lengths.len(),
            valid.join(",")
        ));
    }
    valid.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    let first = valid.first().cloned().unwrap_or_else(|| {
        thrown("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0")
    });
    let numbered = first.rfind('_').filter(|i| {
        let tail = &first[i + 1..];
        tail.len() == 1 && tail.bytes().all(|b| b.is_ascii_digit())
    });
    let sequence_column = match numbered {
        Some(i) => first[..i].to_string(),
        None => first.clone(),
    };
    let index_of = |name: &str| header.iter().position(|h| h == name);
    let mut seen: Vec<String> = Vec::new();
    let mut out: Vec<(String, Metric)> = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let field = |name: &str| -> Option<String> {
            let i = index_of(name).unwrap_or_else(|| {
                thrown(&format!(
                    "java.util.NoSuchElementException: column {name} in {}",
                    reference_path(path)
                ))
            });
            fields.get(i).map(|s| s.to_string())
        };
        let mut seqs = Vec::new();
        for n in 0..lengths.len() {
            let column = if numbered.is_some() {
                format!("{sequence_column}_{}", n + 1)
            } else {
                sequence_column.clone()
            };
            match field(&column) {
                Some(v) => seqs.push(v),
                None => {
                    messages.push(format!("Null barcode in column {column} of row: {line}"));
                    seqs.push(String::new());
                }
            }
        }
        let joined = seqs.concat();
        if !joined.is_empty() && joined.bytes().all(|b| b == b'N' || b == b'n') {
            continue;
        }
        let display = seqs.join("-");
        if seen.contains(&display) {
            messages.push(format!(
                "Barcode {display} specified more than once in {path}"
            ));
        }
        seen.push(display);
        let name = index_of("barcode_name").map_or(Some(String::new()), |_| field("barcode_name"));
        let library =
            index_of("library_name").map_or(Some(String::new()), |_| field("library_name"));
        let metric = Metric::new(name, library, &seqs);
        match out.iter_mut().find(|(k, _)| *k == joined) {
            Some(slot) => slot.1 = metric,
            None => out.push((joined, metric)),
        }
    }
    out
}

/// `HashSet<>(n)` iteration order for keys with these hashes, in insertion order.
pub fn hash_set_order<T>(items: Vec<(i32, T)>, initial: usize) -> Vec<T> {
    let mut capacity = initial.max(1).next_power_of_two();
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

/// `ExtractBarcodesProgram.finalizeMetrics`.
pub fn finalize(metrics: &mut [Metric], no_match: &mut Metric) {
    let mut total = no_match.reads;
    let mut total_pf = no_match.pf_reads;
    let mut assigned = 0;
    for m in metrics.iter() {
        total += m.reads;
        total_pf += m.pf_reads;
        assigned += m.pf_reads;
    }
    if total > 0 {
        no_match.pct = no_match.reads as f64 / total as f64;
        let mut best = 0.0;
        for m in metrics.iter_mut() {
            m.pct = m.reads as f64 / total as f64;
            if m.pct > best {
                best = m.pct;
            }
        }
        if best > 0.0 {
            no_match.ratio = no_match.pct / best;
            for m in metrics.iter_mut() {
                m.ratio = m.pct / best;
            }
        }
    }
    if total_pf > 0 {
        no_match.pf_pct = no_match.pf_reads as f64 / total_pf as f64;
        let mut best = 0.0;
        for m in metrics.iter_mut() {
            m.pf_pct = m.pf_reads as f64 / total_pf as f64;
            if m.pf_pct > best {
                best = m.pf_pct;
            }
        }
        if best > 0.0 {
            no_match.pf_ratio = no_match.pf_pct / best;
            for m in metrics.iter_mut() {
                m.pf_ratio = m.pf_pct / best;
            }
        }
    }
    if assigned > 0 {
        let mean = assigned as f64 / metrics.len() as f64;
        for m in metrics.iter_mut() {
            m.pf_normalized = m.pf_reads as f64 / mean;
        }
    }
}

impl Extractor {
    /// `new BarcodeExtractor(...)`: the declared barcodes in `HashSet` order, and the cache seeded
    /// with each of them and the all-`N` barcode as inline matches against qualities of sixty.
    pub fn new(
        declared: &[Metric],
        lengths: &[usize],
        options: (i32, i32, i32, i32),
        mode: Distance,
    ) -> Extractor {
        let (max_no_calls, max_mismatches, min_delta, minimum) = options;
        let barcodes = hash_set_order(
            declared
                .iter()
                .map(|m| {
                    let mut h: i32 = 0;
                    for b in m.bytes.iter().flatten() {
                        h = h.wrapping_mul(31).wrapping_add(i32::from(*b as i8));
                    }
                    (h, m.bytes.clone())
                })
                .collect(),
            declared.len(),
        );
        let mut extractor = Extractor {
            barcodes,
            max_no_calls,
            max_mismatches,
            min_delta,
            minimum,
            mode,
            cache: HashMap::new(),
        };
        let perfect: Vec<Vec<u8>> = lengths.iter().map(|l| vec![60u8; *l]).collect();
        let no_match: Vec<Vec<u8>> = lengths.iter().map(|l| vec![b'N'; *l]).collect();
        for bytes in declared
            .iter()
            .map(|m| &m.bytes)
            .chain(std::iter::once(&no_match))
        {
            let found = extractor.calculate(bytes, Some(&perfect), true);
            extractor.cache.insert(bytes.clone(), found);
        }
        extractor
    }
}

impl Distance {
    /// The `DISTANCE_MODE` argument; HAMMING when absent.
    pub fn parse(value: Option<&str>) -> Distance {
        match value.unwrap_or("HAMMING") {
            "LENIENT_HAMMING" => Distance::Lenient,
            "FREE" => Distance::Free,
            _ => Distance::Hamming,
        }
    }
}
