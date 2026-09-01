//! `MergeSamFiles` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.MergeSamFiles.doWork` at tag 3.4.0. The merge itself, header records
//! included, is `picard_analysis::merge_sam_files`; what this adds is the three decisions `doWork`
//! makes around it, each of which the array varies.
//!
//! **Whether anything is sorted.** `doWork` computes `matchedSortOrders` over the inputs, and then:
//!
//! ```java
//! if (matchedSortOrders || SORT_ORDER == unsorted || ASSUME_SORTED || INTERVALS != null) {
//!     headerMergerSortOrder = SORT_ORDER; presorted = true;
//! } else {
//!     headerMergerSortOrder = unsorted;  presorted = false;
//! }
//! ```
//!
//! A presorted writer keeps the merge order. So `ASSUME_SORTED=true` does not mean "the inputs are
//! sorted", it means *nothing sorts them*: a coordinate-sorted input written under
//! `SORT_ORDER=queryname` comes out in coordinate order with a queryname header, and a port that
//! sorted anyway would be wrong on exactly those rows.
//!
//! **What the header says.** The merged header is built by `SamFileHeaderMerger`, which sets the
//! group order to `none`, and `doWork` then overwrites `SO` with `SORT_ORDER` after the merge. The
//! `SO` is therefore what was asked for, whatever the records are actually in.
//!
//! **Whether the run happens at all.** `customCommandLineValidation` refuses
//! `CREATE_INDEX` with any `SORT_ORDER` but `coordinate`, before `doWork` runs, and Barclay prints
//! that as an argument error rather than as an exception.
//!
//! **`INTERVALS`, which is not optional in practice.** The covering array holds it at the corpus's
//! interval list, so every row takes the interval path, and that path has two rules of its own.
//! `doWork` refuses an input it cannot query (`Merging with interval but file is not indexed`),
//! which is every fixture but the indexed BAM, and `SamRecordIntervalIteratorFactory` then reads
//! through `queryOverlapping`, so the output is the mapped records overlapping a uniqued interval,
//! in coordinate order, with the unmapped ones gone. `presorted` is true whenever `INTERVALS` is
//! set, so nothing re-sorts them.
//!
//! `USE_THREADING` only turns on the writer's async I/O, which changes nothing about the bytes;
//! the port accepts it so a row naming it measures the tool.

use std::io::{Read, Write};

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::writer::BamWriter;
use picard_analysis::merge_sam_files::merge_parsed;
use picard_analysis::sort_sam::{
    mate_cigar_refusal_in_order, MateCigarCheck, RequestedOrder, SortOrder,
};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST_OF_PAIR: u16 = 0x40;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// Every value given for a repeatable argument, in command-line order.
fn args_all(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

fn flag(args: &[String], key: &str, default: bool) -> Result<bool, String> {
    match arg(args, key).as_deref() {
        None => Ok(default),
        Some("true") => Ok(true),
        Some("false") => Ok(false),
        Some(other) => Err(format!("unknown value for {key}{other}")),
    }
}

fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
    }
}

fn read_one(
    path: &str,
    stringency: ValidationStringency,
) -> Result<(SamHeader, Vec<BamRecord>), String> {
    let mut raw = Vec::new();
    std::fs::File::open(path)
        .map_err(|e| e.to_string())?
        .read_to_end(&mut raw)
        .map_err(|e| e.to_string())?;
    if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        Ok((header, records))
    } else {
        let text = String::from_utf8(raw).map_err(|e| e.to_string())?;
        read_sam_with(&text, stringency).map_err(|e| format!("{e:?}"))
    }
}

/// One interval of a Picard interval list: 1-based, inclusive.
struct Interval {
    contig: String,
    start: i32,
    end: i32,
}

/// `IntervalList.fromFile(INTERVALS).uniqued()`: the intervals, coordinate-ordered per contig with
/// overlapping and abutting ones merged. Header lines are skipped; a data line is
/// `contig<TAB>start<TAB>end<TAB>strand<TAB>name`.
fn read_intervals(path: &str) -> Result<Vec<Interval>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut intervals: Vec<Interval> = Vec::new();
    for line in text.lines() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        intervals.push(Interval {
            contig: fields[0].to_string(),
            start: fields[1]
                .parse()
                .map_err(|_| format!("interval start: {line}"))?,
            end: fields[2]
                .parse()
                .map_err(|_| format!("interval end: {line}"))?,
        });
    }
    intervals.sort_by(|a, b| {
        (a.contig.as_str(), a.start, a.end).cmp(&(b.contig.as_str(), b.start, b.end))
    });
    let mut uniqued: Vec<Interval> = Vec::new();
    for interval in intervals {
        match uniqued.last_mut() {
            Some(last) if last.contig == interval.contig && interval.start <= last.end + 1 => {
                last.end = last.end.max(interval.end);
            }
            _ => uniqued.push(interval),
        }
    }
    Ok(uniqued)
}

/// Whether htsjdk would say this input `hasIndex()`: a BAM with a `.bai` beside it.
fn has_index(path: &str) -> bool {
    let stem = path.strip_suffix(".bam");
    match stem {
        Some(stem) => {
            std::path::Path::new(&format!("{stem}.bai")).exists()
                || std::path::Path::new(&format!("{path}.bai")).exists()
        }
        None => false,
    }
}

fn reference_span(cigar: &htsjdk_bam::cigar::Cigar) -> i32 {
    cigar
        .elements
        .iter()
        .filter(|element| element.op.consumes_reference_bases())
        .map(|element| element.length as i32)
        .sum()
}

/// What the duplicate-order comparator asks a record for, and the refusal message names.
fn check_of(record: &BamRecord, header: &SamHeader) -> MateCigarCheck {
    let reference = if record.reference_index >= 0 {
        header
            .sequences
            .get(record.reference_index as usize)
            .map(|sequence| sequence.name.clone())
    } else {
        None
    };
    MateCigarCheck {
        read_name: record.read_name.clone(),
        paired: record.flags & READ_PAIRED != 0,
        first_of_pair: record.flags & FIRST_OF_PAIR != 0,
        unmapped: record.flags & READ_UNMAPPED != 0,
        mate_unmapped: record.flags & MATE_UNMAPPED != 0,
        read_length: record.read_bases.len(),
        reference,
        start: record.alignment_start,
        end: record.alignment_start + reference_span(&record.cigar) - 1,
        mate_cigar: match record.tags.get(htsjdk_bam::tag::Tag::new(b"MC")) {
            Some(htsjdk_bam::tag::TagValue::Str(text)) => Some(text.to_string()),
            _ => None,
        },
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let mut inputs = args_all(&args, "INPUT=");
    inputs.extend(args_all(&args, "I="));
    if inputs.is_empty() {
        return Err("INPUT= is required".into());
    }
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let sort_order = arg(&args, "SORT_ORDER=")
        .or_else(|| arg(&args, "SO="))
        .unwrap_or_else(|| "coordinate".to_string());
    let assume_sorted = flag(&args, "ASSUME_SORTED=", false)?;
    let merge_dictionaries = flag(&args, "MERGE_SEQUENCE_DICTIONARIES=", false)?;
    // Accepted and ignored: it selects the writer's async I/O path, not its bytes.
    let _use_threading = flag(&args, "USE_THREADING=", false)?;
    let create_index = flag(&args, "CREATE_INDEX=", false)?;
    let stringency = match arg(&args, "VALIDATION_STRINGENCY=").as_deref() {
        None | Some("STRICT") => ValidationStringency::Strict,
        Some("LENIENT") => ValidationStringency::Lenient,
        Some("SILENT") => ValidationStringency::Silent,
        Some(other) => return Err(format!("unknown VALIDATION_STRINGENCY: {other}").into()),
    };

    // `customCommandLineValidation`, which runs before `doWork` and is reported by Barclay as an
    // argument error rather than as an uncaught exception.
    if create_index && sort_order != "coordinate" {
        eprintln!("Can't CREATE_INDEX unless SORT_ORDER is coordinate");
        std::process::exit(1);
    }

    let intervals = match arg(&args, "INTERVALS=").or_else(|| arg(&args, "RGN=")) {
        Some(path) => Some(read_intervals(&path)?),
        None => None,
    };

    // `doWork` opens each input and, with INTERVALS set, refuses the ones it cannot query. The
    // check is per input and in input order, so the message names the first such file.
    if intervals.is_some() {
        for path in &inputs {
            if !has_index(path) {
                eprintln!(
                    "Exception in thread \"main\" picard.PicardException: \
                     Merging with interval but file is not indexed: {path}"
                );
                std::process::exit(1);
            }
        }
    }

    let mut parsed: Vec<(SamHeader, Vec<BamRecord>)> = inputs
        .iter()
        .map(|path| read_one(path, stringency))
        .collect::<Result<_, _>>()?;

    // `queryOverlapping`: the mapped records that overlap a uniqued interval, in the file's own
    // coordinate order. An unmapped record has no interval to overlap and is dropped.
    if let Some(intervals) = &intervals {
        for (header, records) in parsed.iter_mut() {
            records.retain(|record| {
                // "Unmapped reads are discarded" means *unplaced* ones. A read flagged unmapped
                // but given its mate's position is still at a coordinate, so the index returns it,
                // and with no cigar it spans the single base it starts on. Five of this corpus's
                // records are exactly that, and dropping them was a five-record difference.
                if record.reference_index < 0 || record.alignment_start < 1 {
                    return false;
                }
                let Some(contig) = header.sequences.get(record.reference_index as usize) else {
                    return false;
                };
                let start = record.alignment_start;
                let end = start + reference_span(&record.cigar).max(1) - 1;
                intervals.iter().any(|interval| {
                    interval.contig == contig.name && start <= interval.end && end >= interval.start
                })
            });
        }
    }

    // `matchedSortOrders`: every input already claims the order that was asked for.
    let matched = parsed
        .iter()
        .all(|(header, _)| header.attributes.get("SO").unwrap_or("unsorted") == sort_order);
    let presorted = matched || sort_order == "unsorted" || assume_sorted || intervals.is_some();

    // A presorted writer keeps the merge order; only the other branch sorts, and then it sorts by
    // the order the writer was given.
    // The duplicate order refuses before anything is written, and it refuses whether or not the
    // writer was presorted: `SamFileHeaderMerger` is given SORT_ORDER on the presorted branch, so
    // it is `MergingSamRecordIterator`'s own comparator that asks for the mate cigar, and the
    // writer's `SortingCollection` that asks for it on the other branch. Either way the message is
    // the same and no output exists.
    if sort_order == "duplicate" {
        let checks: Vec<MateCigarCheck> = parsed
            .iter()
            .flat_map(|(header, records)| {
                records.iter().map(move |record| check_of(record, header))
            })
            .collect();
        // File order, not the SortingCollection's: the merging iterator asks each stream's
        // records in turn, so the first record that needs a mate cigar is the one named.
        let order: Vec<usize> = (0..checks.len()).collect();
        if let Some(message) = mate_cigar_refusal_in_order(&checks, &order) {
            eprintln!("Exception in thread \"main\" htsjdk.samtools.SAMException: {message}");
            std::process::exit(1);
        }
        eprintln!("sorting by duplicate order is a separate surface");
        std::process::exit(1);
    }

    let sort_with = if presorted {
        None
    } else {
        match RequestedOrder::parse(&sort_order) {
            Some(RequestedOrder::Coordinate) => Some(SortOrder::Coordinate),
            Some(RequestedOrder::Queryname) => Some(SortOrder::Queryname),
            // The duplicate order never reaches here: it is refused above, before the writer.
            Some(RequestedOrder::Duplicate) => unreachable!("duplicate is refused before the sort"),
            None => {
                // `unknown`: `getComparatorInstance()` is null and `SortingCollection` falls back
                // to natural ordering, which casts a record to `Comparable`.
                eprintln!(
                    "Exception in thread \"main\" java.lang.ClassCastException: class \
                     htsjdk.samtools.BAMRecord cannot be cast to class java.lang.Comparable \
                     (htsjdk.samtools.BAMRecord is in unnamed module of loader 'app'; \
                     java.lang.Comparable is in module java.base of loader 'bootstrap')"
                );
                std::process::exit(1);
            }
        }
    };

    let (header, records) = merge_parsed(parsed, &sort_order, merge_dictionaries, sort_with)
        .map_err(|e| format!("{e:?}"))?;

    if output.ends_with(".sam") {
        let sam = write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?;
        let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
        out.write_all(sam.as_bytes())?;
        out.flush()?;
        return Ok(());
    }

    let mut writer = BamWriter::new(Vec::new(), &header)?;
    for record in &records {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    let bam = writer.finish()?;
    std::fs::write(&output, &bam)?;
    if create_index && header.attributes.get("SO") == Some("coordinate") {
        let index = build_bam_index(&bam).map_err(|e| format!("{e:?}"))?;
        std::fs::write(index_path(&output), index)?;
    }
    Ok(())
}
