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

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_metrics::file::{Histogram, MetricsFile};
use picard_analysis::args::arg;
use picard_analysis::mark_illumina_adapters::{
    read_bases_in_read_order, trim_paired_reads, trim_single_read, Adapter, ADAPTER_PAIRS,
    DEFAULT_ADAPTER_LENGTH, MAX_ERROR_RATE, MAX_PE_ERROR_RATE, MIN_MATCH_BASES, MIN_MATCH_PE_BASES,
};

const READ_PAIRED_FLAG: u16 = 0x1;
const READ_STRAND_FLAG: u16 = 0x10;
const FIRST_OF_PAIR_FLAG: u16 = 0x40;
const SECOND_OF_PAIR_FLAG: u16 = 0x80;

/// Every value given for a repeatable argument, in the order the command line gives them.
fn args_of(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

/// The reference throws rather than exiting, so its handler prints the class before the message.
/// `AdapterMarker`'s pruning: once `threshold` adapter hits have been tallied, the list searched
/// from then on is the `keep` adapters seen most often, plus any tied with the last of them.
///
/// The ranking goes through a `TreeMap` keyed by COUNT, filled from a `HashMap` of adapter to count,
/// so two adapters seen equally often collide and the one the hash map yields last survives. The
/// reads already written keep their tags whatever the pruning decides; the pair that crossed the
/// threshold has not been written yet, and loses its tag if its adapter was not kept.
struct Pruning {
    threshold: i64,
    keep: usize,
    reached: bool,
    seen: i64,
    /// `(adapter, count)` in the order each adapter was first tallied.
    counts: Vec<(Adapter, i64)>,
}

impl Pruning {
    fn new(threshold: i64, keep: usize) -> Self {
        Pruning {
            threshold,
            keep,
            reached: false,
            seen: 0,
            counts: Vec::new(),
        }
    }

    /// `TruncatedAdapterPair.hashCode`: the five prime sequence (not in read order) and then the
    /// three prime one, through `String.hashCode`.
    fn hash(adapter: &Adapter) -> i32 {
        let string = |s: &str| {
            s.encode_utf16()
                .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(c as i32))
        };
        let five = picard_analysis::mark_illumina_adapters::reverse_complement(
            &adapter.five_prime_read_order,
        );
        string(&five)
            .wrapping_mul(31)
            .wrapping_add(string(&adapter.three_prime_read_order))
    }

    /// Tallies a hit for `adapters[position]` and prunes when the threshold is crossed. Returns
    /// whether the reads that hit it keep their tags.
    fn tally(&mut self, adapters: &mut Vec<Adapter>, position: usize) -> bool {
        if self.reached || self.threshold < 1 {
            return true;
        }
        let found = adapters[position].clone();
        let same = |a: &Adapter, b: &Adapter| {
            a.three_prime_read_order == b.three_prime_read_order
                && a.five_prime_read_order == b.five_prime_read_order
        };
        match self.counts.iter_mut().find(|(a, _)| same(a, &found)) {
            Some(entry) => entry.1 += 1,
            None => self.counts.push((found.clone(), 1)),
        }
        self.seen += 1;
        if self.seen < self.threshold {
            return true;
        }
        // A `HashMap` this small never resizes: its order is the bucket, then insertion.
        let mut ordered: Vec<(usize, &(Adapter, i64))> = self
            .counts
            .iter()
            .map(|entry| {
                let h = Self::hash(&entry.0) as u32;
                (((h ^ (h >> 16)) & 15) as usize, entry)
            })
            .collect();
        ordered.sort_by_key(|(bucket, _)| *bucket);
        let mut by_count: std::collections::BTreeMap<std::cmp::Reverse<i64>, Adapter> =
            std::collections::BTreeMap::new();
        for (_, (adapter, count)) in ordered {
            by_count.insert(std::cmp::Reverse(*count), adapter.clone());
        }
        let mut best: Vec<Adapter> = Vec::new();
        let mut last = i64::MAX;
        for (std::cmp::Reverse(count), adapter) in by_count {
            if best.len() >= self.keep {
                if count == last {
                    best.push(adapter);
                } else {
                    break;
                }
            } else {
                last = count;
                best.push(adapter);
            }
        }
        self.reached = true;
        let kept = best.iter().any(|a| same(a, &found));
        *adapters = best;
        kept
    }
}

fn throw(message: &str) -> ! {
    picard_analysis::stream_output::on_thrown();
    eprintln!("Exception in thread \"main\" picard.PicardException: {message}");
    std::process::exit(1);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("MarkIlluminaAdapters");
    picard_analysis::output::install();
    let args: Vec<String> = picard_analysis::args::process().into_iter().collect();
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
    let num_adapters_to_keep: i64 = arg(&args, "NUM_ADAPTERS_TO_KEEP=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(1);
    let prune_after: i64 = arg(&args, "PRUNE_ADAPTER_LIST_AFTER_THIS_MANY_ADAPTERS_SEEN=")
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(100);

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

    let (header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();
    let queryname = header.attributes.get("SO") == Some("queryname");

    // `AdapterMarker.setNumAdaptersToKeep`, called while the marker is built, before any read.
    if num_adapters_to_keep <= 0 {
        eprintln!(
            "Exception in thread \"main\" java.lang.IllegalArgumentException: numAdaptersToKeep should be positive: {num_adapters_to_keep}"
        );
        std::process::exit(1);
    }
    let mut pruning = Pruning::new(prune_after, num_adapters_to_keep as usize);

    // `CREATE_INDEX` is accepted and writes nothing: an index needs a coordinate-sorted file, and
    // this tool refuses anything but queryname order as soon as it sees a paired read. htsjdk logs
    // that it is not creating one and writes the file anyway, which is a line on stderr rather
    // than a difference in the answer. The output is a BAM whatever its name, opened before the
    // first record as `makeWriter` opens it, so a run that throws leaves it empty.
    let mut out = output
        .as_deref()
        .map(|output| picard_analysis::stream_output::Output::bam(output, &header, true))
        .transpose()?;
    // `histo.increment(r.getReadLength() - clip + 1)`: the bases a read would lose if it were
    // clipped at the tag, counted for every marked record and for no other.
    let mut bins: std::collections::BTreeMap<i64, f64> = std::collections::BTreeMap::new();
    let xt = Tag::new(b"XT");
    let mut emit = |record: &BamRecord| -> Result<(), Box<dyn std::error::Error>> {
        if let Some(TagValue::Int(tag)) = record.tags.get(xt) {
            let clipped = record.read_bases.len() as i64 - tag + 1;
            *bins.entry(clipped).or_insert(0.0) += 1.0;
        }
        match &mut out {
            Some(out) => out.write(record),
            None => Ok(()),
        }
    };
    let mut records = records;
    while let Some(mut record) = records.next() {
        let paired = record.flags & READ_PAIRED_FLAG != 0;
        // `rec.getReadPairedFlag() && iterator.hasNext() ? iterator.next() : null`: the second
        // record is consumed BEFORE anything is checked, so a paired record at the end of the file
        // is a missing mate rather than a single-end read.
        let mate = if paired { records.next() } else { None };
        record.tags.remove(xt);

        if paired {
            // The order is asserted only once a paired read is seen, so a single-end file in any
            // order runs.
            if !queryname {
                throw("Input file must be sorted by queryname");
            }
            let Some(mut second) = mate else {
                throw(&format!(
                    "Missing mate pair for paired read: {}",
                    record.read_name
                ));
            };
            second.tags.remove(xt);
            if record.read_name != second.read_name {
                throw(&format!(
                    "Adjacent reads expected to be mate-pairs have different names: {}, {}",
                    record.read_name, second.read_name
                ));
            }
            let (first_of_pair, second_of_pair) = if record.flags & FIRST_OF_PAIR_FLAG != 0
                && second.flags & SECOND_OF_PAIR_FLAG != 0
            {
                (&mut record, &mut second)
            } else if record.flags & SECOND_OF_PAIR_FLAG != 0
                && second.flags & FIRST_OF_PAIR_FLAG != 0
            {
                (&mut second, &mut record)
            } else {
                throw(&format!(
                    "Two reads with same name but not correctly marked as 1st/2nd of pair: {}",
                    record.read_name
                ));
            };
            let first_bases = read_bases_in_read_order(
                &first_of_pair.read_bases,
                first_of_pair.flags & READ_STRAND_FLAG != 0,
            );
            let second_bases = read_bases_in_read_order(
                &second_of_pair.read_bases,
                second_of_pair.flags & READ_STRAND_FLAG != 0,
            );
            let trim = trim_paired_reads(
                &first_bases,
                &second_bases,
                &adapters,
                min_match_pe,
                max_error_pe,
            );
            if let Some(tag) = trim.first_tag {
                first_of_pair.tags.insert(xt, TagValue::Int(tag as i64));
            }
            if let Some(tag) = trim.second_tag {
                second_of_pair.tags.insert(xt, TagValue::Int(tag as i64));
            }
            if let Some(position) = trim.matched {
                if !pruning.tally(&mut adapters, position) {
                    first_of_pair.tags.remove(xt);
                    second_of_pair.tags.remove(xt);
                }
            }
            emit(&record)?;
            emit(&second)?;
        } else {
            let bases =
                read_bases_in_read_order(&record.read_bases, record.flags & READ_STRAND_FLAG != 0);
            if let Some((position, found)) =
                trim_single_read(&bases, &adapters, min_match_se, max_error_se)
            {
                record.tags.insert(xt, TagValue::Int((found + 1) as i64));
                if !pruning.tally(&mut adapters, position) {
                    record.tags.remove(xt);
                }
            }
            emit(&record)?;
        }
    }
    if let Some(out) = out {
        // `cram_output::write`: a CRAM takes `CREATE_INDEX` from the command line.
        out.finish(picard_analysis::cram_output::create_index_argument(), false)?;
    }

    let mut file = MetricsFile::new();
    for header in picard_analysis::command_line::default_headers("MarkIlluminaAdapters") {
        file.add_header(&header);
    }
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
