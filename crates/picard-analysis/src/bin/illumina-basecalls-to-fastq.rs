//! `IlluminaBasecallsToFastq` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.IlluminaBasecallsToFastq` with `BasecallsConverter` (sorted and
//! unsorted), the EAMSS filter, `Casava18ReadNameEncoder` and `IlluminaReadNameEncoder`,
//! `TrimmingUtil.findQualityTrimPoint` and `AdapterMarker` at tag 3.4.0, for per-tile runs.
//!
//! Every writer is opened, and so every file created, from `MULTIPLEX_PARAMS` before a cluster is
//! read. A cluster's matched barcode comes from the per-tile barcode file, or, with
//! `MATCH_BARCODES_INLINE`, from matching its barcode reads against the parameters' barcodes; only
//! then is `METRICS_FILE` written, and its unmatched row is the program's own, which the converter
//! never counts into. A template read cut by `TRIMMING_QUALITY` or an adapter to fewer than
//! `MIN_TRIMMED_LENGTH` bases is cut there instead, `Arrays.copyOfRange` padding it with zero bytes.

use std::io::Write;

use htsjdk_metrics::file::MetricsFile;
use picard_analysis::barcode_extractor::{finalize, parse_input_file, Distance, Extractor, Metric};
use picard_analysis::basecalls_converter::{unexpected, ClusterData, Converter, Options};
use picard_analysis::fingerprint::{host_path, java_hash_order, reference_path};
use picard_analysis::illumina_files::{parse_read_structure, SegmentKind};
use picard_analysis::java_hash_map::string_hash_code;
use picard_analysis::mark_illumina_adapters::{
    adapter_pair, find_index_of_clip_sequence, reverse_complement,
    substring_and_remove_trailing_ns, NO_MATCH,
};
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};

const TOOL: &str = "IlluminaBasecallsToFastq";

/// One truncated adapter pair: its three prime and its five prime, both in read order.
struct AdapterPair {
    three: Vec<u8>,
    five_read_order: Vec<u8>,
}

/// `new AdapterMarker(30, adapters)`: truncate each pair, collapsing those that became the same.
fn adapter_marker(pairs: &[(String, String)]) -> Vec<AdapterPair> {
    let mut out: Vec<AdapterPair> = Vec::new();
    for (three, five_read_order) in pairs {
        let candidate = AdapterPair {
            three: substring_and_remove_trailing_ns(three, 30).into_bytes(),
            five_read_order: substring_and_remove_trailing_ns(five_read_order, 30).into_bytes(),
        };
        if !out
            .iter()
            .any(|a| a.three == candidate.three && a.five_read_order == candidate.five_read_order)
        {
            out.push(candidate);
        }
    }
    out
}

/// `TrimmingUtil.findQualityTrimPoint`.
fn quality_trim_point(quals: &[u8], trim: i32) -> usize {
    let length = quals.len();
    if trim < 1 || length == 0 {
        return 0;
    }
    let (mut score, mut max_score, mut point) = (0i32, 0i32, length);
    for i in (0..length).rev() {
        score += trim - i32::from(quals[i]);
        if score < 0 {
            break;
        }
        if score > max_score {
            max_score = score;
            point = i;
        }
    }
    point
}

/// `Arrays.copyOfRange(a, 0, n)`: cut, or padded with zeros past the end.
fn copy_of_range(a: &[u8], n: usize) -> Vec<u8> {
    let mut out = a[..n.min(a.len())].to_vec();
    out.resize(n, 0);
    out
}

/// A writer's files: one per template, sample barcode and molecular barcode read.
struct FastqWriter {
    template: Vec<Out>,
    sample: Vec<Out>,
    molecular: Vec<Out>,
}

/// A FASTQ file, plain or BGZF (`BlockCompressedOutputStream`, which ends with its empty block).
enum Out {
    Plain(std::io::BufWriter<std::fs::File>),
    Bgzf(htsjdk_bgzf::BgzfWriter<std::fs::File>),
}

impl Out {
    fn write_all(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        match self {
            Out::Plain(w) => w.write_all(bytes),
            Out::Bgzf(w) => w.write_all(bytes),
        }
    }

    fn close(&mut self) -> std::io::Result<()> {
        match self {
            Out::Plain(w) => w.flush(),
            Out::Bgzf(w) => w.finish(),
        }
    }
}

fn open(path: &std::path::Path, compress: bool) -> Out {
    let file = std::fs::File::create(path).unwrap_or_else(|e| {
        thrown(&format!(
            "htsjdk.samtools.util.RuntimeIOException: Error opening file: file://{}: {e}",
            reference_path(&path.display().to_string())
        ))
    });
    if compress {
        Out::Bgzf(htsjdk_bgzf::BgzfWriter::new(file))
    } else {
        Out::Plain(std::io::BufWriter::new(file))
    }
}

struct Settings {
    name_format: String,
    machine: Option<String>,
    run_barcode: Option<String>,
    flowcell: Option<String>,
    trimming: Option<i32>,
    min_trimmed: usize,
    adapters: Vec<AdapterPair>,
    lane: i32,
}

impl Settings {
    fn short_name(&self, c: &ClusterData) -> String {
        let run = self.run_barcode.as_deref().unwrap_or("null");
        if self.name_format == "ILLUMINA" {
            format!("{run}:{}:{}:{}:{}", self.lane, c.tile, c.x, c.y)
        } else {
            format!(
                "{}:{run}:{}:{}:{}:{}:{}",
                self.machine.as_deref().unwrap_or("null"),
                self.flowcell.as_deref().unwrap_or("null"),
                self.lane,
                c.tile,
                c.x,
                c.y
            )
        }
    }

    fn read_name(&self, c: &ClusterData, pair: Option<usize>) -> String {
        let short = self.short_name(c);
        if self.name_format == "ILLUMINA" {
            match pair {
                Some(p) => format!("{short}/{p}"),
                None => short,
            }
        } else {
            format!(
                "{short} {}:{}:0:{}",
                pair.map(|p| p.to_string()).unwrap_or_default(),
                if c.pf { 'N' } else { 'Y' },
                c.barcode.as_deref().unwrap_or("")
            )
        }
    }

    /// `ClusterToFastqWriter.write`.
    fn write(&self, w: &mut FastqWriter, c: &ClusterData) {
        let templates = w.template.len();
        let molecular = w.molecular.len();
        let (mut t, mut b, mut m) = (0usize, 0usize, 0usize);
        for (kind, bases, quals) in &c.reads {
            let (out, name, bases, quals) = match kind {
                SegmentKind::Template => {
                    t += 1;
                    let name = self.read_name(c, (templates > 1).then_some(t));
                    let (bases, quals) = self.trim(bases, quals, t);
                    (&mut w.template[t - 1], name, bases, quals)
                }
                SegmentKind::Barcode => {
                    b += 1;
                    (
                        &mut w.sample[b - 1],
                        self.read_name(c, None),
                        bases.clone(),
                        quals.clone(),
                    )
                }
                _ => {
                    m += 1;
                    (
                        &mut w.molecular[m - 1],
                        self.read_name(c, (molecular > 1).then_some(m)),
                        bases.clone(),
                        quals.clone(),
                    )
                }
            };
            let fastq: Vec<u8> = quals.iter().map(|q| q.wrapping_add(33)).collect();
            let mut record = Vec::with_capacity(name.len() + 2 * bases.len() + 6);
            record.push(b'@');
            record.extend_from_slice(name.as_bytes());
            record.push(b'\n');
            record.extend_from_slice(&bases);
            record.extend_from_slice(b"\n+\n");
            record.extend_from_slice(&fastq);
            record.push(b'\n');
            if let Err(e) = out.write_all(&record) {
                thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}"));
            }
        }
    }

    /// `ClusterToFastqWriter.trimRead`.
    fn trim(&self, bases: &[u8], quals: &[u8], template: usize) -> (Vec<u8>, Vec<u8>) {
        let mut bases = bases.to_vec();
        let mut quals = quals.to_vec();
        if let Some(trim) = self.trimming {
            let index = quality_trim_point(&quals, trim).max(self.min_trimmed);
            quals = copy_of_range(&quals, index);
            bases = copy_of_range(&bases, index);
        }
        if !self.adapters.is_empty() {
            // `ClippingUtility.findAdapterPairAndIndexForSingleRead`.
            for adapter in &self.adapters {
                let sequence = match template {
                    1 => &adapter.three,
                    2 => &adapter.five_read_order,
                    _ => thrown("picard.PicardException: Read template index must be 1 or 2"),
                };
                let index = find_index_of_clip_sequence(&bases, sequence, 12, 0.1);
                if index != NO_MATCH {
                    let index = (index as usize).max(self.min_trimmed);
                    quals = copy_of_range(&quals, index);
                    bases = copy_of_range(&bases, index);
                    break;
                }
            }
        }
        (bases, quals)
    }
}

fn main() {
    let args = Args::from_env(&[
        ("B", "BASECALLS_DIR"),
        ("L", "LANE"),
        ("M", "METRICS_FILE"),
        ("RS", "READ_STRUCTURE"),
        ("O", "OUTPUT_PREFIX"),
    ]);
    let basecalls = std::path::PathBuf::from(args.required("BASECALLS_DIR"));
    let barcodes_dir = args.get("BARCODES_DIR").map(std::path::PathBuf::from);
    let multiplex = args.get("MULTIPLEX_PARAMS").map(str::to_string);
    let output_prefix = args.get("OUTPUT_PREFIX").map(str::to_string);
    let input_params = args
        .get("INPUT_PARAMS_FILE")
        .map(str::to_string)
        .or_else(|| multiplex.clone());
    let name_format = args
        .get("READ_NAME_FORMAT")
        .unwrap_or("CASAVA_1_8")
        .to_string();
    let machine = args.get("MACHINE_NAME").map(str::to_string);
    let flowcell = args.get("FLOWCELL_BARCODE").map(str::to_string);
    let five = args.get("FIVE_PRIME_ADAPTER").map(str::to_string);
    let three = args.get("THREE_PRIME_ADAPTER").map(str::to_string);
    let structure = parse_read_structure(&args.required("READ_STRUCTURE"))
        .unwrap_or_else(|| thrown("picard.PicardException: Read structure could not be parsed"));
    let lengths: Vec<usize> = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::Barcode)
        .map(|s| s.cycles)
        .collect();

    // `customCommandLineValidation`, then the parent's.
    let mut errors = Vec::new();
    if name_format == "CASAVA_1_8" && machine.is_none() {
        errors.push(
            "MACHINE_NAME is required when using Casava1.8-style read name headers.".to_string(),
        );
    }
    if name_format == "CASAVA_1_8" && flowcell.is_none() {
        errors.push(
            "FLOWCELL_BARCODE is required when using Casava1.8-style read name headers."
                .to_string(),
        );
    }
    if five.is_none() != three.is_none() {
        errors.push(
            "THREE_PRIME_ADAPTER and FIVE_PRIME_ADAPTER must either both be null or both be set."
                .to_string(),
        );
    }
    let mut declared: Vec<(String, Metric)> = Vec::new();
    if let Some(path) = &input_params {
        declared = parse_input_file(path, &lengths, &mut errors);
        if declared.is_empty() {
            errors.push("No barcodes have been specified.".to_string());
        }
    }
    if !errors.is_empty() {
        refuse_validation(TOOL, &errors);
    }

    // `initialize`.
    let mut adapter_pairs: Vec<(String, String)> = args
        .collection("ADAPTERS_TO_CHECK", &[])
        .iter()
        .filter_map(|n| adapter_pair(n))
        .map(|a| (a.three_prime_read_order, a.five_prime_read_order))
        .collect();
    if let (Some(f), Some(t)) = (&five, &three) {
        adapter_pairs.push((t.clone(), reverse_complement(f)));
    }
    let lane_list: Vec<i32> = args
        .collection("LANE", &[])
        .iter()
        .filter_map(|v| v.parse().ok())
        .collect();
    let compress = args.bool("COMPRESS_OUTPUTS", false);
    let suffix = if compress { "fastq.gz" } else { "fastq" };
    let count = |k: SegmentKind| structure.iter().filter(|s| s.kind == k).count();
    let build_writer = |prefix: &str| -> FastqWriter {
        let host = std::path::PathBuf::from(host_path(prefix));
        let dir = std::path::absolute(&host)
            .ok()
            .and_then(|p| p.parent().map(std::path::Path::to_path_buf))
            .unwrap_or_default();
        let name = host
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let files = |format: &dyn Fn(usize) -> String, n: usize| -> Vec<Out> {
            (1..=n)
                .map(|i| open(&dir.join(format(i)), compress))
                .collect()
        };
        FastqWriter {
            template: files(
                &|i| format!("{name}.{i}.{suffix}"),
                count(SegmentKind::Template),
            ),
            sample: files(
                &|i| format!("{name}.barcode_{i}.{suffix}"),
                count(SegmentKind::Barcode),
            ),
            molecular: files(
                &|i| format!("{name}.index_{i}.{suffix}"),
                count(SegmentKind::MolecularIndex),
            ),
        }
    };
    let mut writers: Vec<(Option<String>, FastqWriter)> = Vec::new();
    let demultiplex;
    if let Some(prefix) = &output_prefix {
        writers.push((None, build_writer(prefix)));
        demultiplex = false;
    } else {
        let path = multiplex.clone().unwrap_or_default();
        let text = std::fs::read_to_string(&path).unwrap_or_else(|_| {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                reference_path(&path)
            ))
        });
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let labels: Vec<String> = (1..=lengths.len())
            .map(|i| format!("BARCODE_{i}"))
            .collect();
        let mut expected: Vec<String> = vec!["OUTPUT_PREFIX".to_string()];
        expected.extend(labels.iter().cloned());
        let missing: Vec<String> = java_hash_order(
            expected
                .into_iter()
                .filter(|e| !header.contains(&e.as_str()))
                .map(|e| (string_hash_code(&e), e))
                .collect(),
        );
        if !missing.is_empty() {
            thrown(&format!(
                "picard.PicardException: MULTIPLEX_PARAMS file {} is missing the following columns: {}.",
                reference_path(&path),
                missing.join(", ")
            ));
        }
        let column = |name: &str| header.iter().position(|h| *h == name);
        let rows: Vec<Vec<&str>> = lines
            .filter(|l| !l.is_empty())
            .map(|l| l.split('\t').collect())
            .collect();
        let mut seen: Vec<Option<String>> = Vec::new();
        for row in &rows {
            let field = |name: &str| {
                column(name)
                    .and_then(|i| row.get(i).copied())
                    .unwrap_or("")
                    .to_string()
            };
            let values: Vec<String> = labels.iter().map(|l| field(l)).collect();
            let key = if labels.is_empty() || values.iter().any(|v| v == "N") {
                None
            } else {
                Some(values.concat())
            };
            if seen.contains(&key) {
                thrown(&format!(
                    "picard.PicardException: Row for barcode {} appears more than once in MULTIPLEX_PARAMS file {path}",
                    key.as_deref().unwrap_or("null")
                ));
            }
            seen.push(key.clone());
            let writer = build_writer(&field("OUTPUT_PREFIX"));
            match writers.iter_mut().find(|(k, _)| *k == key) {
                Some(slot) => slot.1 = writer,
                None => writers.push((key, writer)),
            }
        }
        if seen.is_empty() {
            thrown(&format!(
                "picard.PicardException: MULTIPLEX_PARAMS file {path} does have any data rows."
            ));
        }
        demultiplex = true;
    }
    let inline = args.bool("MATCH_BARCODES_INLINE", false) && demultiplex;
    let metrics: Vec<Metric> = declared.iter().map(|(_, m)| m.clone()).collect();
    let no_match_seqs: Vec<String> = lengths.iter().map(|l| "N".repeat(*l)).collect();
    let mut no_match = Metric::new(None, None, &no_match_seqs);
    let extractor = inline.then(|| {
        Extractor::new(
            &metrics,
            &lengths,
            (
                args.int("MAX_NO_CALLS", 2) as i32,
                args.int("MAX_MISMATCHES", 1) as i32,
                args.int("MIN_MISMATCH_DELTA", 1) as i32,
                args.int("MINIMUM_BASE_QUALITY", 0) as i32,
            ),
            Distance::parse(args.get("DISTANCE_MODE")),
        )
    });
    let mut converter = Converter::new(Options {
        basecalls: &basecalls,
        lanes: &lane_list,
        structure: &structure,
        demultiplex,
        barcodes_dir: (!inline).then(|| barcodes_dir.clone().unwrap_or_else(|| basecalls.clone())),
        first_tile: args.get("FIRST_TILE").and_then(|v| v.parse().ok()),
        tile_limit: args.get("TILE_LIMIT").and_then(|v| v.parse().ok()),
        include_non_pf: args.bool("INCLUDE_NON_PF_READS", true),
        apply_eamss: args.bool("APPLY_EAMSS_FILTER", true),
        extractor,
        metrics,
    });

    let mut settings = Settings {
        name_format,
        machine,
        run_barcode: args.get("RUN_BARCODE").map(str::to_string),
        flowcell,
        trimming: args.get("TRIMMING_QUALITY").and_then(|v| v.parse().ok()),
        min_trimmed: args.int("MIN_TRIMMED_LENGTH", 20).max(0) as usize,
        adapters: adapter_marker(&adapter_pairs),
        lane: 0,
    };
    let ignore_unexpected = args.bool("IGNORE_UNEXPECTED_BARCODES", false);
    let sort = args.bool("SORT", true);
    for tile in converter.tiles.clone() {
        // Clusters by barcode, each in the order it was read.
        let mut by_barcode: Vec<(Option<String>, Vec<ClusterData>)> = Vec::new();
        for data in converter.tile(tile) {
            let key = data.barcode.clone();
            if !writers.iter().any(|(k, _)| *k == key) && !ignore_unexpected {
                unexpected(&key, sort);
            }
            match by_barcode.iter_mut().find(|(k, _)| *k == key) {
                Some((_, v)) => v.push(data),
                None => by_barcode.push((key, vec![data])),
            }
        }
        for (key, clusters) in by_barcode.iter_mut() {
            if sort {
                clusters.sort_by_cached_key(|c| {
                    settings.lane = c.lane;
                    settings.short_name(c)
                });
            }
            if let Some((_, w)) = writers.iter_mut().find(|(k, _)| k == key) {
                for c in clusters.iter() {
                    settings.lane = c.lane;
                    settings.write(w, c);
                }
            }
        }
    }
    for (_, w) in writers.iter_mut() {
        for out in w
            .template
            .iter_mut()
            .chain(w.sample.iter_mut())
            .chain(w.molecular.iter_mut())
        {
            if let Err(e) = out.close() {
                thrown(&format!("htsjdk.samtools.util.RuntimeIOException: {e}"));
            }
        }
    }
    drop(writers);
    if let (Some(path), true) = (args.get("METRICS_FILE"), inline) {
        let mut metrics = converter.metrics;
        finalize(&mut metrics, &mut no_match);
        let mut file = MetricsFile::new();
        file.add_header(&format!("{TOOL} <command line>"));
        file.add_header("Started on: <timestamp>");
        for m in &metrics {
            file.add_metric(m);
        }
        file.add_metric(&no_match);
        if let Err(e) = std::fs::write(path, file.write()) {
            thrown(&format!("htsjdk.samtools.SAMException: {e}"));
        }
    }
}
