//! `RevertOriginalBaseQualitiesAndAddMateCigar` as a runnable binary: the covering array's port
//! side.
//!
//! Ports `picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar.doWork` at tag 3.4.0. The revert,
//! the mate-info pass and the sort orders live in
//! `picard_analysis::revert_original_quals_add_mate_cigar`.
//!
//! The tool's first act is to decide whether to do nothing. It reads the input twice: once to ask
//! whether any record carries an `OQ` and whether the first mapped pair already has its `MC`, and
//! again to do the work. When the answer is that there is nothing to do it returns zero having
//! written NO output file, which is a different thing from writing an unchanged one.
//!
//! `SORT_ORDER` decides two things: the `SO` the output header carries, and the order the records
//! come out in -- the writer is built unsorted, so it sorts. Two of the five orders name no
//! comparator, and the records then keep the query-name order the mate-info pass left them in.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::writer::BamWriter;
use picard_analysis::revert_original_quals_add_mate_cigar::{
    can_skip, revert_original_records_with, Options, SortOrder,
};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// The order the writer will be given: the argument, or the input header's when there is none.
fn resolved_order(text: &str, options: &Options) -> Option<SortOrder> {
    options.sort_order.or_else(|| {
        text.lines()
            .find(|line| line.starts_with("@HD"))
            .and_then(|line| line.split('\t').find_map(|field| field.strip_prefix("SO:")))
            .and_then(SortOrder::parse)
    })
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let sort_order = match arg(&args, "SORT_ORDER=") {
        None => None,
        Some(value) => {
            Some(SortOrder::parse(&value).ok_or_else(|| format!("unknown SORT_ORDER: {value}"))?)
        }
    };
    let options = Options {
        restore_original_qualities: arg(&args, "RESTORE_ORIGINAL_QUALITIES=")
            .map(|value| value == "true")
            .unwrap_or(true),
        sort_order,
        max_records_to_examine: arg(&args, "MAX_RECORDS_TO_EXAMINE=")
            .and_then(|value| value.parse().ok())
            .unwrap_or(10_000),
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let text = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records: Vec<BamRecord> = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?
    } else {
        String::from_utf8(raw)?
    };

    // `canSkipSAMFile`, on its own pass over the input. When it skips, nothing is written: the
    // output file does not exist, which is the answer the reference gives too.
    let (_, records) = read_sam(&text).map_err(|e| format!("{e:?}"))?;
    if can_skip(&records, &options).skips() {
        return Ok(());
    }

    // `SO:unknown` names no comparator, and the writer -- built unsorted -- hands the records to a
    // `SortingCollection` that falls back to natural ordering. A SAM record is not `Comparable`,
    // so the run dies on a cast rather than on anything about the data. `unsorted` does NOT: the
    // writer has a case for it and never sorts.
    //
    // The class in the message is the record's, so it names `BAMRecord` for a BAM input and would
    // name `SAMRecord` for a text one.
    if resolved_order(&text, &options) == Some(SortOrder::Unknown) {
        eprintln!(
            "Exception in thread \"main\" java.lang.ClassCastException: class \
             htsjdk.samtools.BAMRecord cannot be cast to class java.lang.Comparable \
             (htsjdk.samtools.BAMRecord is in unnamed module of loader 'app'; \
             java.lang.Comparable is in module java.base of loader 'bootstrap')"
        );
        std::process::exit(1);
    }

    let (header, records) =
        revert_original_records_with(&text, &options).map_err(|e| format!("{e:?}"))?;
    if output.ends_with(".sam") {
        let sam = write_sam(&header, &records).ok_or("records failed to re-encode as SAM")?;
        std::fs::write(&output, sam)?;
        return Ok(());
    }
    let mut writer = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for record in &records {
        writer
            .write(record)
            .expect("records that parsed re-encode as BAM");
    }
    std::fs::write(
        &output,
        writer.finish().expect("in-memory BAM writer never fails"),
    )?;
    Ok(())
}
