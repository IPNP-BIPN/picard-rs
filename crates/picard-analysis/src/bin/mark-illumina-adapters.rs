//! `MarkIlluminaAdapters` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.illumina.MarkIlluminaAdapters.doWork` at tag 3.4.0. The search itself, the
//! truncation and the paired rule live in `picard_analysis::mark_illumina_adapters`.
//!
//! Three things about this tool are easy to get wrong from its name.
//!
//! It MARKS: the record's bases are written out untouched and the answer is the `XT` tag, a
//! one-based position. It clears any `XT` the input carried before it searches, so a file marked
//! twice carries the second run's answer rather than the union of the two.
//!
//! It marks a PAIR as a pair. The paired path finds one index for both reads and writes it onto
//! both, so a pair where only one read carries an adapter comes back with two tags at the same
//! position -- or with none, when the one that matched does not survive the stricter re-check.
//!
//! And it never searches with the adapters it was given: `AdapterMarker` truncates them to
//! `ADAPTER_TRUNCATION_LENGTH` first, which is what makes `INDEXED` and `DUAL_INDEXED` the same
//! adapter.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;
use htsjdk_metrics::file::{Histogram, MetricsFile};
use picard_analysis::mark_illumina_adapters::{
    read_bases_in_read_order, trim_paired_reads, trim_single_read, Adapter, ADAPTER_PAIRS,
    DEFAULT_ADAPTER_LENGTH, MAX_ERROR_RATE, MAX_PE_ERROR_RATE, MIN_MATCH_BASES, MIN_MATCH_PE_BASES,
};

const READ_PAIRED_FLAG: u16 = 0x1;
const READ_STRAND_FLAG: u16 = 0x10;
const FIRST_OF_PAIR_FLAG: u16 = 0x40;
const SECOND_OF_PAIR_FLAG: u16 = 0x80;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// Every value given for a repeatable argument, in the order the command line gives them.
fn args_of(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

/// The reference throws rather than exiting, so its handler prints the class before the message.
fn throw(message: &str) -> ! {
    eprintln!("Exception in thread \"main\" picard.PicardException: {message}");
    std::process::exit(1);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    // OUTPUT is optional: "If output is not specified, just the metrics are generated".
    let output = arg(&args, "OUTPUT=").or_else(|| arg(&args, "O="));
    let metrics = arg(&args, "METRICS=")
        .or_else(|| arg(&args, "M="))
        .ok_or("METRICS= is required")?;

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let number = |key: &str, default: usize| -> Result<usize, Box<dyn std::error::Error>> {
        match arg(&args, key) {
            None => Ok(default),
            Some(value) => Ok(value.parse()?),
        }
    };
    let rate = |key: &str, default: f64| -> Result<f64, Box<dyn std::error::Error>> {
        match arg(&args, key) {
            None => Ok(default),
            Some(value) => Ok(value.parse()?),
        }
    };
    let min_match_se = number("MIN_MATCH_BASES_SE=", MIN_MATCH_BASES)?;
    let min_match_pe = number("MIN_MATCH_BASES_PE=", MIN_MATCH_PE_BASES)?;
    let max_error_se = rate("MAX_ERROR_RATE_SE=", MAX_ERROR_RATE)?;
    let max_error_pe = rate("MAX_ERROR_RATE_PE=", MAX_PE_ERROR_RATE)?;
    let truncation = number("ADAPTER_TRUNCATION_LENGTH=", DEFAULT_ADAPTER_LENGTH)?;

    // `customCommandLineValidation`: the two custom sequences are a pair, and half a pair is a
    // command line error rather than a run with a default filled in.
    let five_prime = arg(&args, "FIVE_PRIME_ADAPTER=");
    let three_prime = arg(&args, "THREE_PRIME_ADAPTER=");
    if five_prime.is_some() != three_prime.is_some() {
        eprintln!(
            "THREE_PRIME_ADAPTER and FIVE_PRIME_ADAPTER must either both be null or both be set."
        );
        std::process::exit(1);
    }

    // Picard builds its Barclay parser with `APPEND_TO_COLLECTIONS`, so a list argument the caller
    // names is ADDED to its declared default rather than replacing it. Naming one adapter pair
    // therefore searches FOUR: the three defaults, and then the named one. It is the opposite of
    // GATK, whose parser clears the default the moment the argument is named at all, and it is
    // visible in the output: a row that asks for `INDEXED` alone still marks a `PAIRED_END`
    // adapter, because `PAIRED_END` is still in the list.
    //
    // A value of `null` clears what has been accumulated, which is the documented way to search
    // only the pairs the caller named.
    let mut names: Vec<String> = ["INDEXED", "DUAL_INDEXED", "PAIRED_END"]
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    for value in args_of(&args, "ADAPTERS=") {
        if value == "null" {
            names.clear();
        } else {
            names.push(value);
        }
    }
    for name in &names {
        if !ADAPTER_PAIRS.iter().any(|(known, _, _)| known == name) {
            eprintln!("Argument ADAPTERS has bad value: '{name}' is not a valid value");
            std::process::exit(1);
        }
    }
    let borrowed: Vec<&str> = names.iter().map(String::as_str).collect();
    let mut adapters =
        picard_analysis::mark_illumina_adapters::truncated_adapters(&borrowed, truncation);
    if let (Some(five), Some(three)) = (five_prime.as_deref(), three_prime.as_deref()) {
        // The custom pair is appended after the named ones and truncated with them, so it is
        // searched last and can be collapsed into a named pair that truncates to the same bytes.
        let custom = Adapter {
            name: "truncated CustomAdapterPair".to_string(),
            three_prime_read_order:
                picard_analysis::mark_illumina_adapters::substring_and_remove_trailing_ns(
                    three, truncation,
                ),
            five_prime_read_order:
                picard_analysis::mark_illumina_adapters::substring_and_remove_trailing_ns(
                    &picard_analysis::mark_illumina_adapters::reverse_complement(five),
                    truncation,
                ),
        };
        match adapters.iter_mut().find(|existing| {
            existing.three_prime_read_order == custom.three_prime_read_order
                && existing.five_prime_read_order == custom.five_prime_read_order
        }) {
            Some(existing) => existing.name = format!("{}|CustomAdapterPair", existing.name),
            None => adapters.push(custom),
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, mut records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };
    let queryname = header.attributes.get("SO") == Some("queryname");

    let xt = Tag::new(b"XT");
    let mut index = 0usize;
    while index < records.len() {
        let paired = records[index].flags & READ_PAIRED_FLAG != 0;
        // `rec.getReadPairedFlag() && iterator.hasNext() ? iterator.next() : null`: the second
        // record is consumed BEFORE anything is checked, so a paired record at the end of the file
        // is a missing mate rather than a single-end read.
        let mate = if paired && index + 1 < records.len() {
            Some(index + 1)
        } else {
            None
        };
        records[index].tags.remove(xt);

        if paired {
            // The order is asserted only once a paired read is seen, so a single-end file in any
            // order runs.
            if !queryname {
                throw("Input file must be sorted by queryname");
            }
            let Some(second) = mate else {
                throw(&format!(
                    "Missing mate pair for paired read: {}",
                    records[index].read_name
                ));
            };
            records[second].tags.remove(xt);
            if records[index].read_name != records[second].read_name {
                throw(&format!(
                    "Adjacent reads expected to be mate-pairs have different names: {}, {}",
                    records[index].read_name, records[second].read_name
                ));
            }
            let (first_of_pair, second_of_pair) = if records[index].flags & FIRST_OF_PAIR_FLAG != 0
                && records[second].flags & SECOND_OF_PAIR_FLAG != 0
            {
                (index, second)
            } else if records[index].flags & SECOND_OF_PAIR_FLAG != 0
                && records[second].flags & FIRST_OF_PAIR_FLAG != 0
            {
                (second, index)
            } else {
                throw(&format!(
                    "Two reads with same name but not correctly marked as 1st/2nd of pair: {}",
                    records[index].read_name
                ));
            };
            let first_bases = read_bases_in_read_order(
                &records[first_of_pair].read_bases,
                records[first_of_pair].flags & READ_STRAND_FLAG != 0,
            );
            let second_bases = read_bases_in_read_order(
                &records[second_of_pair].read_bases,
                records[second_of_pair].flags & READ_STRAND_FLAG != 0,
            );
            let trim = trim_paired_reads(
                &first_bases,
                &second_bases,
                &adapters,
                min_match_pe,
                max_error_pe,
            );
            if let Some(tag) = trim.first_tag {
                records[first_of_pair]
                    .tags
                    .insert(xt, TagValue::Int(tag as i64));
            }
            if let Some(tag) = trim.second_tag {
                records[second_of_pair]
                    .tags
                    .insert(xt, TagValue::Int(tag as i64));
            }
        } else {
            let bases = read_bases_in_read_order(
                &records[index].read_bases,
                records[index].flags & READ_STRAND_FLAG != 0,
            );
            if let Some((_, found)) =
                trim_single_read(&bases, &adapters, min_match_se, max_error_se)
            {
                records[index]
                    .tags
                    .insert(xt, TagValue::Int((found + 1) as i64));
            }
        }
        index += if mate.is_some() { 2 } else { 1 };
    }

    // `histo.increment(r.getReadLength() - clip + 1)`: the bases a read would lose if it were
    // clipped at the tag, counted for every marked record and for no other.
    let mut bins: std::collections::BTreeMap<i64, f64> = std::collections::BTreeMap::new();
    for record in &records {
        if let Some(TagValue::Int(tag)) = record.tags.get(xt) {
            let clipped = record.read_bases.len() as i64 - tag + 1;
            *bins.entry(clipped).or_insert(0.0) += 1.0;
        }
    }

    // `CREATE_INDEX` is accepted and writes nothing: an index needs a coordinate-sorted file, and
    // this tool refuses anything but queryname order as soon as it sees a paired read. htsjdk logs
    // that it is not creating one and writes the file anyway, which is a line on stderr rather
    // than a difference in the answer.
    if let Some(output) = output {
        let mut writer =
            BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
        for record in &records {
            writer
                .write(record)
                .expect("records that parsed re-encode as BAM");
        }
        std::fs::write(
            &output,
            writer
                .finish()
                .expect("in-memory BAM writer never fails to finish"),
        )?;
    }

    let mut file = MetricsFile::new();
    file.add_header("MarkIlluminaAdapters <command line>");
    file.add_header("Started on: <timestamp>");
    file.histograms.push(Histogram {
        bin_label: "clipped_bases".to_string(),
        value_label: "read_count".to_string(),
        key_class: "java.lang.Integer".to_string(),
        bins: bins
            .into_iter()
            .map(|(key, value)| (key.to_string(), value))
            .collect(),
    });
    std::fs::write(&metrics, file.write())?;
    Ok(())
}
