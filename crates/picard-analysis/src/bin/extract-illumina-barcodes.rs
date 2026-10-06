//! `ExtractIlluminaBarcodes` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.ExtractIlluminaBarcodes`, `ExtractBarcodesProgram`, `BarcodeExtractor`,
//! `DistanceMetric` and `SingleBarcodeDistanceMetric` at tag 3.4.0 over
//! `picard_analysis::illumina_reader`.
//!
//! The declared barcodes are walked in the order of a `HashSet` sized to their count, so which of
//! two equidistant barcodes is "best" (and written lower-cased when neither matches) is the
//! bucket order of their hashes. The extractor caches a match by its bases: every declared barcode
//! and the all-`N` one are computed up front as inline matches against perfect qualities, so a
//! cluster whose bases ARE a declared barcode gets that answer, whatever a fresh match would say.
//!
//! `BARCODE_FILE` replaces `INPUT_PARAMS_FILE` before it is read, and `MINIMUM_QUALITY` is read in
//! a field initialiser before the parser assigns it, so neither changes anything. `COMPRESS_OUTPUTS`
//! gzips through the JDK's deflater and is not ported.

use std::collections::HashMap;

use htsjdk_metrics::file::{MetricBean, MetricsFile, Value};
use picard_analysis::fingerprint::{java_hash_order, reference_path};
use picard_analysis::illumina_files::{parse_read_structure, SegmentKind};
use picard_analysis::illumina_reader::Run;
use picard_analysis::java_hash_map::string_hash_code;
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};

const TOOL: &str = "ExtractIlluminaBarcodes";

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
struct Metric {
    barcode: String,
    name: Option<String>,
    library: Option<String>,
    bytes: Vec<Vec<u8>>,
    reads: i64,
    pf_reads: i64,
    perfect: i64,
    pf_perfect: i64,
    one: i64,
    pf_one: i64,
    pct: f64,
    ratio: f64,
    pf_pct: f64,
    pf_ratio: f64,
    pf_normalized: f64,
}

impl Metric {
    fn new(name: Option<String>, library: Option<String>, seqs: &[String]) -> Metric {
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
enum Distance {
    Hamming,
    Lenient,
    Free,
}

/// `BarcodeExtractor.BarcodeMatch`.
#[derive(Clone, Default)]
struct Match {
    matched: bool,
    barcode: String,
    mismatches: i32,
    to_second: i32,
}

fn is_no_call(b: u8) -> bool {
    matches!(b, b'N' | b'n' | b'.')
}

/// `SingleBarcodeDistanceMetric`, for one barcode read.
struct Single<'a> {
    barcode: &'a [u8],
    read: &'a [u8],
    quals: Option<&'a [u8]>,
    masked: Vec<u8>,
    minimum: i32,
    max: i32,
}

impl<'a> Single<'a> {
    fn new(
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

    fn hamming(&self) -> i32 {
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

    fn lenient(&self) -> i32 {
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

    fn free(&self) -> i32 {
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

struct Extractor {
    /// The declared barcodes in `HashSet` order.
    barcodes: Vec<Vec<Vec<u8>>>,
    max_no_calls: i32,
    max_mismatches: i32,
    min_delta: i32,
    minimum: i32,
    mode: Distance,
    cache: HashMap<Vec<Vec<u8>>, Match>,
}

impl Extractor {
    /// `DistanceMetric.distance`.
    fn distance(
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
    fn calculate(&self, read: &[Vec<u8>], quals: Option<&[Vec<u8>]>, inline: bool) -> Match {
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
    fn find(&mut self, read: &[Vec<u8>], quals: Option<&[Vec<u8>]>) -> Match {
        let lookup =
            quals.is_none_or(|q| q.iter().flatten().all(|v| i32::from(*v) >= self.minimum));
        if !lookup {
            return self.calculate(read, quals, false);
        }
        let m = match self.cache.get(read) {
            Some(m) => m.clone(),
            None => self.calculate(read, quals, false),
        };
        if m.matched {
            self.cache.insert(read.to_vec(), m.clone());
        }
        m
    }
}

/// `TabbedTextFileWithHeaderParser` over the barcode file, then `parseInputFile`.
fn parse_input_file(
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
fn hash_set_order<T>(items: Vec<(i32, T)>, initial: usize) -> Vec<T> {
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
fn finalize(metrics: &mut [Metric], no_match: &mut Metric) {
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

fn main() {
    let args = Args::from_env(&[
        ("B", "BASECALLS_DIR"),
        ("L", "LANE"),
        ("M", "METRICS_FILE"),
        ("RS", "READ_STRUCTURE"),
    ]);
    let basecalls = std::path::PathBuf::from(args.required("BASECALLS_DIR"));
    let metrics_file = args.required("METRICS_FILE");
    let output_dir = args
        .get("OUTPUT_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| basecalls.clone());
    let lanes: Vec<i32> = args
        .collection("LANE", &[])
        .iter()
        .filter_map(|v| v.parse().ok())
        .collect();
    let structure = parse_read_structure(&args.required("READ_STRUCTURE"))
        .unwrap_or_else(|| thrown("picard.PicardException: Read structure could not be parsed"));
    let barcode_file = args.get("BARCODE_FILE").map(str::to_string);
    let barcode_args = args.collection("BARCODE", &[]);
    if args.bool("COMPRESS_OUTPUTS", false) {
        // The reference gzips through the JDK's deflater, whose bytes are not the tool's logic;
        // the array holds the argument false.
        thrown("java.lang.UnsupportedOperationException: COMPRESS_OUTPUTS is not ported");
    }
    let lengths: Vec<usize> = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::Barcode)
        .map(|s| s.cycles)
        .collect();

    // `customCommandLineValidation`: this class's messages first, then its parent's.
    let mut messages = Vec::new();
    let mut declared: Vec<(String, Metric)> = Vec::new();
    if barcode_file.is_none() {
        let mut seen: Vec<&String> = Vec::new();
        for barcode in &barcode_args {
            if seen.contains(&barcode) {
                messages.push(format!("Barcode {barcode} specified more than once."));
            }
            seen.push(barcode);
            let mut seqs = Vec::new();
            let mut at = 0;
            for &l in &lengths {
                let piece = barcode.get(at..at + l).unwrap_or_else(|| {
                    thrown(&format!(
                        "java.lang.StringIndexOutOfBoundsException: begin {at}, end {}, length {}",
                        at + l,
                        barcode.len()
                    ))
                });
                seqs.push(piece.to_string());
                at += l;
            }
            let metric = Metric::new(None, None, &seqs);
            match declared.iter_mut().find(|(k, _)| k == barcode) {
                Some(slot) => slot.1 = metric,
                None => declared.push((barcode.clone(), metric)),
            }
        }
    }
    let mut parent = Vec::new();
    if let Some(path) = &barcode_file {
        declared = parse_input_file(path, &lengths, &mut parent);
        if declared.is_empty() {
            parent.push("No barcodes have been specified.".to_string());
        }
    }
    if (barcode_file.is_some() || !barcode_args.is_empty()) && declared.is_empty() {
        messages.push("No barcodes have been specified.".to_string());
    }
    messages.extend(parent);
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }

    let minimum = args.int("MINIMUM_BASE_QUALITY", 0) as i32;
    let mode = match args.get("DISTANCE_MODE").unwrap_or("HAMMING") {
        "LENIENT_HAMMING" => Distance::Lenient,
        "FREE" => Distance::Free,
        _ => Distance::Hamming,
    };
    let no_match_seqs: Vec<String> = lengths.iter().map(|l| "N".repeat(*l)).collect();
    let mut no_match = Metric::new(None, None, &no_match_seqs);
    let ordered: Vec<Vec<Vec<u8>>> = hash_set_order(
        declared
            .iter()
            .map(|(_, m)| {
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
        barcodes: ordered,
        max_no_calls: args.int("MAX_NO_CALLS", 2) as i32,
        max_mismatches: args.int("MAX_MISMATCHES", 1) as i32,
        min_delta: args.int("MIN_MISMATCH_DELTA", 1) as i32,
        minimum,
        mode,
        cache: HashMap::new(),
    };
    let perfect: Vec<Vec<u8>> = lengths.iter().map(|l| vec![60u8; *l]).collect();
    for (_, m) in &declared {
        let found = extractor.calculate(&m.bytes, Some(&perfect), true);
        extractor.cache.insert(m.bytes.clone(), found);
    }
    let found = extractor.calculate(&no_match.bytes, Some(&perfect), true);
    extractor.cache.insert(no_match.bytes.clone(), found);

    // The provider reads every cycle the structure does not skip; the barcode reads are the `B`
    // segments among them.
    let mut reads: Vec<Vec<i32>> = Vec::new();
    let mut barcode_reads: Vec<usize> = Vec::new();
    let mut cycle = 1;
    for s in &structure {
        let cycles: Vec<i32> = (cycle..cycle + s.cycles as i32).collect();
        cycle += s.cycles as i32;
        if s.kind == SegmentKind::Skip {
            continue;
        }
        if s.kind == SegmentKind::Barcode {
            barcode_reads.push(reads.len());
        }
        reads.push(cycles);
    }
    let mut metrics: Vec<Metric> = declared.iter().map(|(_, m)| m.clone()).collect();
    for &lane in &lanes {
        let run = Run::new(&basecalls, lane);
        let tiles = run.available_tiles().unwrap_or_else(|e| thrown(&e));
        for tile in tiles {
            let clusters = run.clusters(tile, &reads, 2).unwrap_or_else(|e| thrown(&e));
            let mut text = String::new();
            for c in &clusters {
                let read: Vec<Vec<u8>> = barcode_reads
                    .iter()
                    .map(|i| c.reads[*i].0.clone())
                    .collect();
                let quals: Vec<Vec<u8>> = barcode_reads
                    .iter()
                    .map(|i| c.reads[*i].1.clone())
                    .collect();
                let m = extractor.find(&read, (minimum > 0).then_some(quals.as_slice()));
                // `updateMetrics`.
                let target = if m.matched {
                    let at = declared
                        .iter()
                        .position(|(k, _)| *k == m.barcode)
                        .unwrap_or_else(|| thrown("java.lang.NullPointerException"));
                    &mut metrics[at]
                } else {
                    &mut no_match
                };
                target.reads += 1;
                if c.pf {
                    target.pf_reads += 1;
                }
                if m.matched {
                    if m.mismatches == 0 {
                        target.perfect += 1;
                        if c.pf {
                            target.pf_perfect += 1;
                        }
                    } else if m.mismatches == 1 {
                        target.one += 1;
                        if c.pf {
                            target.pf_one += 1;
                        }
                    }
                }
                for r in &read {
                    text.push_str(&String::from_utf8_lossy(r));
                }
                text.push_str(&format!(
                    "\t{}\t{}\t{}\t{}\n",
                    if m.matched { "Y" } else { "N" },
                    m.barcode,
                    m.mismatches,
                    m.to_second
                ));
            }
            let path = output_dir.join(format!("s_{lane}_{tile:04}_barcode.txt"));
            let written = std::fs::write(&path, text);
            if let Err(e) = written {
                thrown(&format!("htsjdk.samtools.SAMException: {e}"));
            }
        }
    }
    finalize(&mut metrics, &mut no_match);
    let mut file = MetricsFile::new();
    file.add_header(&format!("{TOOL} <command line>"));
    file.add_header("Started on: <timestamp>");
    for m in &metrics {
        file.add_metric(m);
    }
    file.add_metric(&no_match);
    if let Err(e) = std::fs::write(&metrics_file, file.write()) {
        thrown(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
