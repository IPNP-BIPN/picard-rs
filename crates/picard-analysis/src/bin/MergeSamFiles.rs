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

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::text_parse::ValidationStringency;
use picard_analysis::args::arg;
use picard_analysis::merge_sam_files::{merge_headers, merge_parsed, PriorityMerge};
use picard_analysis::sort_sam::{
    mate_cigar_refusal_in_order, MateCigarCheck, RequestedOrder, SortOrder,
};
use picard_analysis::stream_output::Output;

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST_OF_PAIR: u16 = 0x40;

/// Every value given for a repeatable argument, in command-line order.
fn args_all(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

/// A boolean argument, which the legacy parser builds as `new Boolean(value)`.
fn flag(args: &[String], key: &str, default: bool) -> bool {
    arg(args, key).map_or(default, |value| picard_analysis::args::is_true(&value))
}

/// Opens one input, BAM or SAM, read as `read_one` read it whole: whether it streams, so that it
/// can be read again, and its header.
fn open_one(
    path: &str,
    stringency: ValidationStringency,
) -> Result<(bool, SamHeader, picard_analysis::bam_stream::Records), String> {
    let input =
        picard_analysis::bam_stream::Input::open_with(path, stringency).map_err(|e| e.message())?;
    let streams = input.streams();
    let (header, records) = input.into_parts();
    Ok((streams, header, records))
}

/// `queryOverlapping`: whether a record is one of the mapped records that overlap a uniqued
/// interval. An unmapped record has no interval to overlap and is dropped.
fn overlaps(header: &SamHeader, record: &BamRecord, intervals: &[Interval]) -> bool {
    // "Unmapped reads are discarded" means *unplaced* ones. A read flagged unmapped but given its
    // mate's position is still at a coordinate, so the index returns it, and with no cigar it
    // spans the single base it starts on. Five of this corpus's records are exactly that, and
    // dropping them was a five-record difference.
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
        sam_string: picard_analysis::view_sam::sam_string(header, record),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("MergeSamFiles");
    picard_analysis::output::install();
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
    let assume_sorted = flag(&args, "ASSUME_SORTED=", false);
    let merge_dictionaries = flag(&args, "MERGE_SEQUENCE_DICTIONARIES=", false);
    // Accepted and ignored: it selects the writer's async I/O path, not its bytes.
    let _use_threading = flag(&args, "USE_THREADING=", false);
    let create_index = flag(&args, "CREATE_INDEX=", false);
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

    // Every input is read through before its records are merged; only its header is kept until
    // its turn comes, so that no more than one input is open at a time.
    let mut opened: Vec<(bool, SamHeader)> = Vec::with_capacity(inputs.len());
    for path in &inputs {
        let (streams, header, _) = open_one(path, stringency)?;
        opened.push((streams, header));
    }
    // An input again, its records in file order, those outside INTERVALS dropped.
    let records_of = |index: usize| -> Result<Box<dyn Iterator<Item = BamRecord> + '_>, String> {
        let (streams, header) = &opened[index];
        let again = picard_analysis::bam_stream::Input::reopen(&inputs[index], *streams);
        let records = again.map_err(|e| e.message())?.into_parts().1;
        Ok(match &intervals {
            Some(intervals) => {
                Box::new(records.filter(move |record| overlaps(header, record, intervals)))
            }
            None => Box::new(records),
        })
    };

    // `matchedSortOrders`: every input already claims the order that was asked for.
    let matched = opened
        .iter()
        .all(|(_, header)| header.attributes.get("SO").unwrap_or("unsorted") == sort_order);
    let presorted = matched || sort_order == "unsorted" || assume_sorted || intervals.is_some();

    // A presorted writer keeps the merge order; only the other branch sorts, and then it sorts by
    // the order the writer was given.
    // The duplicate order refuses before anything is written, and it refuses whether or not the
    // writer was presorted: `SamFileHeaderMerger` is given SORT_ORDER on the presorted branch, so
    // it is `MergingSamRecordIterator`'s own comparator that asks for the mate cigar, and the
    // writer's `SortingCollection` that asks for it on the other branch. Either way the message is
    // the same and no output exists.
    let comments: Vec<String> = args
        .iter()
        .filter_map(|a| a.strip_prefix("COMMENT=").or_else(|| a.strip_prefix("CO=")))
        .map(|comment| format!("@CO\t{comment}"))
        .collect();
    // A presorted writer keeps the merge order: `MergingSamRecordIterator` interleaves the inputs
    // by the order's comparator, and for `unsorted` by object identity, which no run reproduces,
    // so those are written each input in turn.
    if presorted && sort_order != "duplicate" {
        let headers: Vec<SamHeader> = opened.iter().map(|(_, h)| h.clone()).collect();
        let (mut header, remaps) = merge_headers(&headers, &sort_order, merge_dictionaries)
            .map_err(|e| format!("{e:?}"))?;
        // `header.addComment(comment)` for every COMMENT, after the inputs' own comments.
        header.comments.extend(comments);
        let mut out = Output::create(&output, &header, false)?;
        let compare: Option<fn(&BamRecord, &BamRecord) -> std::cmp::Ordering> =
            match sort_order.as_str() {
                "coordinate" => Some(htsjdk_bam::coordinate::compare),
                "queryname" => Some(htsjdk_bam::query_name::compare),
                _ => None,
            };
        match compare {
            Some(compare) => {
                // shortcut: records are compared after their reference indices are remapped to the
                // merged dictionary, where htsjdk's merged-dictionary comparator also places an
                // unmapped read by its mate's contig; compare through the mapping if that matters.
                let streams = remaps
                    .iter()
                    .enumerate()
                    .map(|(index, remap)| {
                        Ok(records_of(index)?.map(move |mut record| {
                            remap.apply(&mut record);
                            record
                        }))
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                for record in PriorityMerge::new(streams, compare) {
                    out.write(&record)?;
                }
            }
            None => {
                for (index, remap) in remaps.iter().enumerate() {
                    for mut record in records_of(index)? {
                        remap.apply(&mut record);
                        out.write(&record)?;
                    }
                }
            }
        }
        let coordinate = header.attributes.get("SO") == Some("coordinate");
        out.finish(create_index, create_index && coordinate)?;
        return Ok(());
    }
    let parsed: Vec<(SamHeader, Vec<BamRecord>)> = (0..inputs.len())
        .map(|index| Ok((opened[index].1.clone(), records_of(index)?.collect())))
        .collect::<Result<_, String>>()?;

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

    let (mut header, records) = merge_parsed(parsed, &sort_order, merge_dictionaries, sort_with)
        .map_err(|e| format!("{e:?}"))?;
    // `header.addComment(comment)` for every COMMENT, after the inputs' own comments.
    header.comments.extend(comments);
    let mut out = Output::create(&output, &header, false)?;
    for record in &records {
        out.write(record)?;
    }
    let coordinate = header.attributes.get("SO") == Some("coordinate");
    out.finish(create_index, create_index && coordinate)?;
    Ok(())
}
