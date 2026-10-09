//! `IlluminaBasecallsToSam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.IlluminaBasecallsToSam` and `ClusterDataToSamConverter` at tag 3.4.0
//! over `picard_analysis::basecalls_converter`, for per-tile runs.
//!
//! One unmapped record per template read, named `<run>:<lane>:<tile>:<x>:<y>`; two templates are a
//! pair. The filter's verdict is the vendor-check flag, a read of nothing but `A`s and no-calls
//! gets `XN:i:1`, and the barcode is written as `BC` (with `QT` when asked) by the population
//! strategy. Each `LIBRARY_PARAMS` row is an output with its own read group per lane, whose
//! attributes are the sample, the library, then the platform, platform unit, centre and run date in
//! that order, then the row's own two-letter columns. The adapters mark `XT`.

use htsjdk_bam::header::{ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use htsjdk_metrics::file::MetricsFile;
use picard_analysis::barcode_extractor::{finalize, parse_input_file, Distance, Extractor, Metric};
use picard_analysis::basecalls_converter::{unexpected, ClusterData, Converter, Options};
use picard_analysis::fingerprint::{host_path, java_hash_order, reference_path};
use picard_analysis::illumina_files::{parse_read_structure, SegmentKind};
use picard_analysis::java_hash_map::string_hash_code;
use picard_analysis::mark_illumina_adapters::{
    adapter_pair, reverse_complement, substring_and_remove_trailing_ns, trim_paired_reads,
    trim_single_read, Adapter, MAX_ERROR_RATE, MAX_PE_ERROR_RATE, MIN_MATCH_BASES,
    MIN_MATCH_PE_BASES,
};
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};

const TOOL: &str = "IlluminaBasecallsToSam";

/// One output: its file, its header, and the records it has been given.
struct Output {
    path: std::path::PathBuf,
    header: SamHeader,
    records: Vec<BamRecord>,
}

/// `new AdapterMarker(30, adapters)`: truncate each pair, collapsing those that became the same.
fn adapter_marker(pairs: &[Adapter]) -> Vec<Adapter> {
    let mut out: Vec<Adapter> = Vec::new();
    for pair in pairs {
        let candidate = Adapter {
            name: format!("truncated {}", pair.name),
            three_prime_read_order: substring_and_remove_trailing_ns(
                &pair.three_prime_read_order,
                30,
            ),
            five_prime_read_order: substring_and_remove_trailing_ns(
                &pair.five_prime_read_order,
                30,
            ),
        };
        match out.iter_mut().find(|a| {
            a.three_prime_read_order == candidate.three_prime_read_order
                && a.five_prime_read_order == candidate.five_prime_read_order
        }) {
            Some(a) => a.name = format!("{}|{}", a.name, pair.name),
            None => out.push(candidate),
        }
    }
    out
}

/// A two-letter tag given on the command line.
fn tag(name: &str) -> Tag {
    match name.as_bytes() {
        [a, b] => Tag::new(&[*a, *b]),
        _ => thrown(&format!(
            "java.lang.IllegalArgumentException: SAM tags must be two characters: {name}"
        )),
    }
}

/// `Iso8601Date.toString()` of a date given as `yyyy/MM/dd`, at midnight UTC.
fn iso8601(date: &str) -> String {
    let parts: Vec<&str> = date.split(['/', '-']).collect();
    match parts.as_slice() {
        [y, m, d] => format!("{y:0>4}-{m:0>2}-{d:0>2}T00:00:00+0000"),
        _ => date.to_string(),
    }
}

fn main() {
    let args = Args::from_env(&[
        ("B", "BASECALLS_DIR"),
        ("L", "LANE"),
        ("M", "METRICS_FILE"),
        ("RS", "READ_STRUCTURE"),
        ("O", "OUTPUT"),
        ("ALIAS", "SAMPLE_ALIAS"),
        ("RG", "READ_GROUP_ID"),
        ("LIB", "LIBRARY_NAME"),
    ]);
    let basecalls = std::path::PathBuf::from(args.required("BASECALLS_DIR"));
    let barcodes_dir = args.get("BARCODES_DIR").map(std::path::PathBuf::from);
    let library_params = args
        .get("BARCODE_PARAMS")
        .or(args.get("LIBRARY_PARAMS"))
        .map(str::to_string);
    let output = args.get("OUTPUT").map(str::to_string);
    let run_barcode = args.required("RUN_BARCODE");
    let structure = parse_read_structure(&args.required("READ_STRUCTURE"))
        .unwrap_or_else(|| thrown("picard.PicardException: Read structure could not be parsed"));
    let lengths: Vec<usize> = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::Barcode)
        .map(|s| s.cycles)
        .collect();
    let molecular = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::MolecularIndex)
        .count();
    let five = args.get("FIVE_PRIME_ADAPTER").map(str::to_string);
    let three = args.get("THREE_PRIME_ADAPTER").map(str::to_string);
    let tag_per_index = args.collection("TAG_PER_MOLECULAR_INDEX", &[]);
    let lanes: Vec<i32> = args
        .collection("LANE", &[])
        .iter()
        .filter_map(|v| v.parse().ok())
        .collect();

    // `customCommandLineValidation`, then the parent's.
    let mut messages = Vec::new();
    if !lengths.is_empty() && library_params.is_none() {
        messages.push("BARCODE_PARAMS or LIBRARY_PARAMS is missing.  If READ_STRUCTURE contains a B (barcode) then either LIBRARY_PARAMS or BARCODE_PARAMS(deprecated) must be provided!".to_string());
    }
    let read_group_id = args
        .get("READ_GROUP_ID")
        .map(str::to_string)
        .unwrap_or_else(|| run_barcode.chars().take(5).collect());
    if !tag_per_index.is_empty() && tag_per_index.len() != molecular {
        messages.push("The number of tags given in TAG_PER_MOLECULAR_INDEX does not match the number of molecular indexes in READ_STRUCTURE".to_string());
    }
    if five.is_none() != three.is_none() {
        messages.push(
            "THREE_PRIME_ADAPTER and FIVE_PRIME_ADAPTER must either both be null or both be set."
                .to_string(),
        );
    }
    let mut first_tile = args.get("FIRST_TILE").and_then(|v| v.parse().ok());
    let mut tile_limit = args.get("TILE_LIMIT").and_then(|v| v.parse().ok());
    let single_tile: Option<i32> = args.get("PROCESS_SINGLE_TILE").and_then(|v| v.parse().ok());
    if let Some(t) = single_tile {
        tile_limit = Some(1);
        first_tile = Some(t);
    }
    let mut declared: Vec<(String, Metric)> = Vec::new();
    if let Some(path) = &library_params {
        declared = parse_input_file(path, &lengths, &mut messages);
        if declared.is_empty() {
            messages.push("No barcodes have been specified.".to_string());
        }
    }
    if !messages.is_empty() {
        refuse_validation(TOOL, &messages);
    }

    // `initialize`.
    let sort = args.bool("SORT", true);
    let platform = if args.given("PLATFORM") {
        args.get("PLATFORM").map(str::to_string)
    } else {
        Some("ILLUMINA".to_string())
    };
    let center = args.get("SEQUENCING_CENTER").map(str::to_string);
    let run_date = args.get("RUN_START_DATE").map(iso8601);
    let include_bc = args.bool("INCLUDE_BC_IN_RG_TAG", false);
    let lane_string = lanes
        .iter()
        .map(|l| l.to_string())
        .collect::<Vec<_>>()
        .join(",");
    // `buildSamHeaderParameters`, a `LinkedHashMap`.
    let header_parameters = |barcodes: Option<&[String]>| -> Vec<(String, Option<String>)> {
        let mut params = Vec::new();
        let mut unit = format!("{run_barcode}.{lane_string}");
        if let Some(b) = barcodes {
            let joined = b.join("-");
            unit = format!("{unit}.{joined}");
            if include_bc {
                params.push(("BC".to_string(), Some(joined)));
            }
        }
        if let Some(p) = &platform {
            params.push(("PL".to_string(), Some(p.clone())));
        }
        params.push(("PU".to_string(), Some(unit)));
        if let Some(c) = &center {
            params.push(("CN".to_string(), Some(c.clone())));
        }
        params.push(("DT".to_string(), run_date.clone()));
        params
    };
    let build = |path: &str,
                 sample: Option<&str>,
                 library: Option<&str>,
                 params: &[(String, Option<String>)]|
     -> Output {
        let mut header = SamHeader::new();
        header.set_sort_order(if sort { "queryname" } else { "unsorted" });
        for lane in &lanes {
            let mut rg = ReadGroup::new(&format!("{read_group_id}.{lane}"));
            if let Some(s) = sample {
                rg.attributes.set("SM", s);
            }
            if let Some(l) = library {
                rg.attributes.set("LB", l);
            }
            for (k, v) in params {
                if let Some(v) = v {
                    rg.attributes.set(k, v);
                }
            }
            header.read_groups.push(rg);
        }
        let host = std::path::PathBuf::from(host_path(path));
        if let Err(e) = std::fs::File::create(&host) {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Cannot write file {}: {e}",
                reference_path(&host.display().to_string())
            ));
        }
        Output {
            path: host,
            header,
            records: Vec::new(),
        }
    };
    let mut writers: Vec<(Option<String>, Output)> = Vec::new();
    if let Some(out) = &output {
        let params = header_parameters(None);
        writers.push((
            None,
            build(
                out,
                args.get("SAMPLE_ALIAS"),
                args.get("LIBRARY_NAME"),
                &params,
            ),
        ));
    } else {
        let path = library_params.clone().unwrap_or_default();
        let text = std::fs::read_to_string(&path).unwrap_or_else(|_| {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                reference_path(&path)
            ))
        });
        let mut lines = text.lines();
        let header: Vec<String> = lines
            .next()
            .unwrap_or("")
            .split('\t')
            .map(str::to_string)
            .collect();
        let has = |c: &str| header.iter().any(|h| h == c);
        let mut labels: Vec<String> = Vec::new();
        if lengths.len() == 1 {
            if has("BARCODE") {
                labels.push("BARCODE".to_string());
            } else if has("BARCODE_1") {
                labels.push("BARCODE_1".to_string());
            } else {
                thrown(&format!("picard.PicardException: LIBRARY_PARAMS(BARCODE_PARAMS) file {path} does not have column BARCODE or BARCODE_1."));
            }
        } else {
            labels = (1..=lengths.len())
                .map(|i| format!("BARCODE_{i}"))
                .collect();
        }
        let mut expected: Vec<String> = ["OUTPUT", "SAMPLE_ALIAS", "LIBRARY_NAME"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        expected.extend(labels.iter().cloned());
        let missing: Vec<String> = java_hash_order(
            expected
                .iter()
                .filter(|e| !has(e))
                .map(|e| (string_hash_code(e), e.clone()))
                .collect(),
        );
        if !missing.is_empty() {
            thrown(&format!(
                "picard.PicardException: LIBRARY_PARAMS file {} is missing the following columns: {}.",
                reference_path(&path),
                missing.join(", ")
            ));
        }
        // The remaining columns, a `HashSet`, are read group tags.
        let tags: Vec<String> = java_hash_order(
            header
                .iter()
                .filter(|h| !expected.contains(h))
                .map(|h| (string_hash_code(h), h.clone()))
                .collect(),
        );
        let forbidden: Vec<String> = java_hash_order(
            header_parameters(None)
                .into_iter()
                .map(|(k, _)| k)
                .filter(|k| tags.contains(k))
                .map(|k| (string_hash_code(&k), k))
                .collect(),
        );
        if !forbidden.is_empty() {
            thrown(&format!(
                "picard.PicardException: Illegal ReadGroup tags in library params(barcode params) file({}) Offending headers = {}",
                reference_path(&path),
                forbidden.join(", ")
            ));
        }
        for column in &tags {
            if column.chars().count() > 2 {
                thrown(&format!("picard.PicardException: Column label ({column}) unrecognized.  Library params(barcode params) can only contain the columns (OUTPUT, LIBRARY_NAME, SAMPLE_ALIAS, BARCODE, BARCODE_<X> where X is a positive integer) OR two letter RG tags!"));
            }
        }
        let column = |name: &str| header.iter().position(|h| h == name);
        for line in lines.filter(|l| !l.is_empty()) {
            let row: Vec<&str> = line.split('\t').collect();
            let field = |name: &str| -> Option<String> {
                column(name).and_then(|i| row.get(i)).map(|s| s.to_string())
            };
            let values: Option<Vec<String>> = (!labels.is_empty()).then(|| {
                labels
                    .iter()
                    .map(|l| field(l).unwrap_or_default())
                    .collect()
            });
            let key = match &values {
                Some(v) if !v.iter().any(|x| x == "N") => Some(v.concat()),
                _ => None,
            };
            if writers.iter().any(|(k, _)| *k == key) {
                thrown(&format!("picard.PicardException: Row for barcode {} appears more than once in LIBRARY_PARAMS or BARCODE_PARAMS file {path}", key.as_deref().unwrap_or("null")));
            }
            let mut params = header_parameters(values.as_deref());
            for t in &tags {
                params.push((t.clone(), field(t)));
            }
            let mut out = field("OUTPUT").unwrap_or_default();
            if let Some(t) = single_tile {
                let p = std::path::Path::new(&out);
                let name = p
                    .file_name()
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_default();
                out = p
                    .parent()
                    .map(|d| d.join(format!("{t}.{name}")).display().to_string())
                    .unwrap_or(out.clone());
            }
            let sample = field("SAMPLE_ALIAS");
            let library = field("LIBRARY_NAME");
            writers.push((
                key,
                build(&out, sample.as_deref(), library.as_deref(), &params),
            ));
        }
        if writers.is_empty() {
            thrown(&format!("picard.PicardException: LIBRARY_PARAMS(BARCODE_PARAMS) file {path} does have any data rows."));
        }
    }

    let mut adapter_list: Vec<Adapter> = args
        .collection(
            "ADAPTERS_TO_CHECK",
            &["INDEXED", "DUAL_INDEXED", "NEXTERA_V2", "FLUIDIGM"],
        )
        .iter()
        .filter_map(|n| adapter_pair(n))
        .collect();
    if let (Some(f), Some(t)) = (&five, &three) {
        adapter_list.push(Adapter {
            name: "Custom adapter pair".to_string(),
            three_prime_read_order: t.clone(),
            five_prime_read_order: reverse_complement(f),
        });
    }
    let adapters = adapter_marker(&adapter_list);
    let demultiplex = !lengths.is_empty();
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
        lanes: &lanes,
        structure: &structure,
        demultiplex,
        barcodes_dir: (!inline).then(|| barcodes_dir.clone().unwrap_or_else(|| basecalls.clone())),
        first_tile,
        tile_limit,
        include_non_pf: args.bool("INCLUDE_NON_PF_READS", true),
        apply_eamss: args.bool("APPLY_EAMSS_FILTER", true),
        extractor,
        metrics,
    });

    let strategy = args
        .get("BARCODE_POPULATION_STRATEGY")
        .unwrap_or("ORPHANS_ONLY")
        .to_string();
    let with_quality = args.bool("INCLUDE_BARCODE_QUALITY", false);
    let index_tag = args.get("MOLECULAR_INDEX_TAG").unwrap_or("RX").to_string();
    let index_quality_tag = args
        .get("MOLECULAR_INDEX_BASE_QUALITY_TAG")
        .unwrap_or("QX")
        .to_string();
    let ignore_unexpected = args.bool("IGNORE_UNEXPECTED_BARCODES", false);
    let templates = structure
        .iter()
        .filter(|s| s.kind == SegmentKind::Template)
        .count();
    let paired = templates == 2;
    // `convertClusterToOutputRecord`.
    let convert = |c: &ClusterData| -> Vec<BamRecord> {
        let name = format!("{run_barcode}:{}:{}:{}:{}", c.lane, c.tile, c.x, c.y);
        let sample: Vec<&(SegmentKind, Vec<u8>, Vec<u8>)> = c.of(SegmentKind::Barcode);
        let seqs: Vec<String> = sample
            .iter()
            .map(|r| String::from_utf8_lossy(&r.1).into_owned())
            .collect();
        let unmatched = (!sample.is_empty()
            && match strategy.as_str() {
                "ALWAYS" => true,
                "INEXACT_MATCH" => Some(seqs.concat()) != c.barcode,
                _ => c.barcode.is_none(),
            })
        .then(|| seqs.join("-").replace('.', "N"));
        let barcode_quality = (unmatched.is_some() && with_quality).then(|| {
            sample
                .iter()
                .map(|r| {
                    String::from_utf8_lossy(&r.2.iter().map(|q| q + 33).collect::<Vec<_>>())
                        .into_owned()
                })
                .collect::<Vec<_>>()
                .join("~")
        });
        let indexes: Vec<(String, String)> = c
            .of(SegmentKind::MolecularIndex)
            .iter()
            .map(|r| {
                (
                    String::from_utf8_lossy(&r.1).replace('.', "N"),
                    String::from_utf8_lossy(&r.2.iter().map(|q| q + 33).collect::<Vec<_>>())
                        .into_owned(),
                )
            })
            .collect();
        let mut records: Vec<BamRecord> = c
            .of(SegmentKind::Template)
            .iter()
            .take(if paired { 2 } else { 1 })
            .enumerate()
            .map(|(i, r)| {
                let mut flags = 0x4u16;
                if paired {
                    flags |= 0x1 | 0x8 | if i == 0 { 0x40 } else { 0x80 };
                }
                if !c.pf {
                    flags |= 0x200;
                }
                let mut rec = BamRecord {
                    read_name: name.clone(),
                    flags,
                    read_bases: r.1.clone(),
                    base_qualities: r.2.clone(),
                    ..Default::default()
                };
                if r.1
                    .iter()
                    .all(|b| matches!(b, b'A' | b'a' | b'N' | b'n' | b'.'))
                {
                    rec.tags.insert(Tag::new(b"XN"), TagValue::Int(1));
                }
                rec.tags.insert(
                    Tag::new(b"RG"),
                    TagValue::Str(format!("{read_group_id}.{}", c.lane)),
                );
                if let Some(bc) = &unmatched {
                    rec.tags.insert(Tag::new(b"BC"), TagValue::Str(bc.clone()));
                    if let Some(q) = &barcode_quality {
                        rec.tags.insert(Tag::new(b"QT"), TagValue::Str(q.clone()));
                    }
                }
                if !indexes.is_empty() {
                    let join = |f: &dyn Fn(&(String, String)) -> String| {
                        indexes.iter().map(f).collect::<Vec<_>>().join("-")
                    };
                    if !index_tag.is_empty() {
                        rec.tags
                            .insert(tag(&index_tag), TagValue::Str(join(&|x| x.0.clone())));
                    }
                    if !index_quality_tag.is_empty() {
                        rec.tags.insert(
                            tag(&index_quality_tag),
                            TagValue::Str(join(&|x| x.1.clone())),
                        );
                    }
                    for (t, (bases, _)) in tag_per_index.iter().zip(&indexes) {
                        rec.tags.insert(tag(t), TagValue::Str(bases.clone()));
                    }
                }
                rec
            })
            .collect();
        if !adapters.is_empty() {
            let xt = Tag::new(b"XT");
            if paired {
                let trim = trim_paired_reads(
                    &records[0].read_bases,
                    &records[1].read_bases,
                    &adapters,
                    MIN_MATCH_PE_BASES,
                    MAX_PE_ERROR_RATE,
                );
                if let Some(t) = trim.first_tag {
                    records[0].tags.insert(xt, TagValue::Int(i64::from(t)));
                }
                if let Some(t) = trim.second_tag {
                    records[1].tags.insert(xt, TagValue::Int(i64::from(t)));
                }
            } else if let Some((_, index)) = trim_single_read(
                &records[0].read_bases,
                &adapters,
                MIN_MATCH_BASES,
                MAX_ERROR_RATE,
            ) {
                records[0]
                    .tags
                    .insert(xt, TagValue::Int(i64::from(index + 1)));
            }
        }
        records
    };

    for tile in converter.tiles.clone() {
        let mut by_barcode: Vec<(Option<String>, Vec<Vec<BamRecord>>)> = Vec::new();
        for data in converter.tile(tile) {
            let key = data.barcode.clone();
            if !writers.iter().any(|(k, _)| *k == key) && !ignore_unexpected {
                unexpected(&key, sort);
            }
            let records = convert(&data);
            match by_barcode.iter_mut().find(|(k, _)| *k == key) {
                Some((_, v)) => v.push(records),
                None => by_barcode.push((key, vec![records])),
            }
        }
        for (key, clusters) in by_barcode.iter_mut() {
            if sort {
                clusters.sort_by(|a, b| a[0].read_name.cmp(&b[0].read_name));
            }
            if let Some((_, w)) = writers.iter_mut().find(|(k, _)| k == key) {
                for c in clusters.drain(..) {
                    w.records.extend(c);
                }
            }
        }
    }
    for (_, w) in &writers {
        let is_sam = w.path.extension().is_some_and(|e| e == "sam");
        let bytes = if is_sam {
            write_sam(&w.header, &w.records)
                .unwrap_or_else(|| {
                    thrown("htsjdk.samtools.SAMException: record cannot be written as SAM")
                })
                .into_bytes()
        } else {
            let mut writer = BamWriter::new(Vec::new(), &w.header)
                .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.SAMException: {e}")));
            for r in &w.records {
                if let Err(e) = writer.write(r) {
                    thrown(&format!("htsjdk.samtools.SAMException: {e:?}"));
                }
            }
            writer
                .finish()
                .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.SAMException: {e}")))
        };
        if let Err(e) = std::fs::write(&w.path, bytes) {
            thrown(&format!("htsjdk.samtools.SAMException: {e}"));
        }
    }
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
