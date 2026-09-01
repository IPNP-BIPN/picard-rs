//! `AddOrReplaceReadGroups` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.AddOrReplaceReadGroups.doWork` at tag 3.4.0. The reheader and the `RG` stamp
//! are `picard_analysis::add_or_replace_read_groups`; what this adds is the one line of `doWork`
//! that decides whether the output is sorted:
//!
//! ```java
//! final SAMFileWriter outWriter = new SAMFileWriterFactory().makeWriter(outHeader,
//!         outHeader.getSortOrder() == inHeader.getSortOrder(), OUTPUT, REFERENCE_SEQUENCE);
//! ```
//!
//! `presorted` is not a property of the records: it is whether `SORT_ORDER` asked for the order the
//! input already claims. So the same `SORT_ORDER=coordinate` writes the input's own order for a
//! coordinate-sorted file and re-sorts a queryname-sorted one, and `SORT_ORDER` unset is always
//! presorted because the output header keeps the input's `SO`.
//!
//! `SAMFileWriterImpl.setHeader` then decides what "not presorted" means:
//!
//! * `unsorted` never sorts, presorted or not: `addAlignment` writes straight through.
//! * any other order, not presorted, builds a `SortingCollection` on that order's comparator.
//! * `duplicate` is that comparator asking each record for its mate's unclipped end, which it reads
//!   from the `MC` tag. No fixture here carries one, so those rows are refusals rather than
//!   outputs, and the message names the record the sort asked first.
//!
//! `unknown` is the row worth having. `SortOrder.unknown.getComparatorInstance()` is null, and a
//! `SortingCollection` built without a comparator falls back to natural ordering: it casts the
//! first record to `Comparable`, which no `SAMRecord` implements. So the reference does not write
//! an output for that combination at all, it dies with a `ClassCastException` naming the record
//! class its reader produced (`BAMRecord` for a BAM, `SAMRecord` for SAM text). Reproducing a
//! crash is reproducing the tool.
//!
//! `CREATE_INDEX` follows `initializeBAMWriter`: an index only for a coordinate-sorted BAM, and
//! silence otherwise. `VALIDATION_STRINGENCY` reaches the reader, and the corpus is valid under all
//! three.

use std::io::{Read, Write};

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::writer::BamWriter;
use picard_analysis::add_or_replace_read_groups::{replace_and_stamp, Options};
use picard_analysis::sort_sam::{mate_cigar_refusal, MateCigarCheck, RequestedOrder};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST_OF_PAIR: u16 = 0x40;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
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
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let opts = Options {
        rgid: arg(&args, "RGID=")
            .or_else(|| arg(&args, "ID="))
            .unwrap_or_else(|| "1".to_string()),
        rglb: arg(&args, "RGLB=")
            .or_else(|| arg(&args, "LB="))
            .ok_or("RGLB= is required")?,
        rgpl: arg(&args, "RGPL=")
            .or_else(|| arg(&args, "PL="))
            .ok_or("RGPL= is required")?,
        rgpu: arg(&args, "RGPU=")
            .or_else(|| arg(&args, "PU="))
            .ok_or("RGPU= is required")?,
        rgsm: arg(&args, "RGSM=")
            .or_else(|| arg(&args, "SM="))
            .ok_or("RGSM= is required")?,
    };
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(false);
    let stringency = match arg(&args, "VALIDATION_STRINGENCY=").as_deref() {
        None | Some("STRICT") => ValidationStringency::Strict,
        Some("LENIENT") => ValidationStringency::Lenient,
        Some("SILENT") => ValidationStringency::Silent,
        Some(other) => return Err(format!("unknown VALIDATION_STRINGENCY: {other}").into()),
    };

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    // Which reader htsjdk would have used, because a message below names its record class.
    let bam_input = raw.starts_with(&[0x1f, 0x8b]);
    let (mut header, mut records): (SamHeader, Vec<BamRecord>) = if bam_input {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let text = reader.header.text.clone();
        let decoded = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<Vec<_>, _>>()?;
        (text, decoded)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam_with(&text, stringency).map_err(|e| format!("{e:?}"))?
    };

    // The input's own `SO`, read before the header is rewritten: it is half of the presorted test.
    let input_order = header
        .attributes
        .get("SO")
        .unwrap_or("unsorted")
        .to_string();
    replace_and_stamp(&mut header, &mut records, &opts);

    let requested = arg(&args, "SORT_ORDER=").or_else(|| arg(&args, "SO="));
    if let Some(order) = requested.as_deref() {
        header.set_sort_order(order);
    }
    let output_order = requested.unwrap_or_else(|| input_order.clone());
    let presorted = output_order == input_order;

    // `unsorted` is written straight through whatever `presorted` says; every other order sorts
    // when the requested order is not the one the input already claims.
    if !presorted && output_order != "unsorted" {
        match RequestedOrder::parse(&output_order) {
            Some(RequestedOrder::Duplicate) => {
                let checks: Vec<MateCigarCheck> = records
                    .iter()
                    .map(|record| check_of(record, &header))
                    .collect();
                if let Some(message) = mate_cigar_refusal(&checks) {
                    eprintln!(
                        "Exception in thread \"main\" htsjdk.samtools.SAMException: {message}"
                    );
                    std::process::exit(1);
                }
                eprintln!("sorting by duplicate order is a separate surface");
                std::process::exit(1);
            }
            Some(order) => {
                let sort = order
                    .sort_order()
                    .expect("coordinate and queryname both sort");
                // Stable, so records equal under the comparator keep input order (decision 0021).
                records.sort_by(sort.comparator());
            }
            None => {
                // `unknown`. `SortOrder.unknown.getComparatorInstance()` is null, and
                // `SortingCollection` with no comparator falls back to natural ordering, so it
                // casts the first record to `Comparable` and fails. The reference therefore does
                // not write an output for this row at all: it dies with a ClassCastException whose
                // text names the record class the reader produced, `BAMRecord` for a BAM and
                // `SAMRecord` for SAM text. Reproducing a crash is reproducing the tool, so the
                // port names the same class its own input implies.
                let record_class = if bam_input { "BAMRecord" } else { "SAMRecord" };
                eprintln!(
                    "Exception in thread \"main\" java.lang.ClassCastException: class \
                     htsjdk.samtools.{record_class} cannot be cast to class java.lang.Comparable \
                     (htsjdk.samtools.{record_class} is in unnamed module of loader 'app'; \
                     java.lang.Comparable is in module java.base of loader 'bootstrap')"
                );
                std::process::exit(1);
            }
        }
    }

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
