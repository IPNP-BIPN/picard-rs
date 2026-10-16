//! `BaitDesigner` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.util.BaitDesigner.doWork` at tag 3.4.0. The design itself, one target at a time
//! under one of three strategies, is `picard_analysis::bait_designer`; this is everything around
//! it: padding the targets, merging the ones that are close enough, walking the reference one
//! contig at a time (and refusing to go back), the statistics, and the six kinds of file the tool
//! writes into `OUTPUT_DIRECTORY`: the targets and baits as interval lists, the parameters file
//! (which the reference writes by reflection over its own fields, USAGE strings included), the
//! design FASTA, and the pools, as FASTA and as Agilent tables.

use std::io::Write;

use htsjdk_bam::interval::{Interval, IntervalList};
use picard_analysis::bait_designer::{
    design_target, masked_base_count, prepare_targets, reverse_complement, Bait, Options, Strategy,
    Target,
};
use picard_analysis::metrics_cli::{
    absolute, fail, refuse_validation, thrown, Args, ReferenceWalker,
};

const USAGE_SUMMARY: &str = "Designs oligonucleotide baits for hybrid selection reactions.";
const USAGE_DETAILS: &str = "<p>This tool is used to design custom bait sets for hybrid selection experiments. The following files are input into BaitDesigner: a (TARGET) interval list indicating the sequences of interest, e.g. exons with their respective coordinates, a reference sequence, and a unique identifier string (DESIGN_NAME). </p><p>The tool will output interval_list files of both bait and target sequences as well as the actual bait sequences in FastA format. At least two baits are output for each target sequence, with greater numbers for larger intervals. Although the default values for both bait size  (120 bases) nd offsets (80 bases) are suitable for most applications, these values can be customized. Offsets represent the distance between sequential baits on a contiguous stretch of target DNA sequence. </p><p>The tool will also output a pooled set of 55,000 (default) oligonucleotides representing all of the baits redundantly. This redundancy achieves a uniform concentration of oligonucleotides for synthesis by a vendor as well as equal numbersof each bait to prevent bias during the hybrid selection reaction. </p><h4>Usage example:</h4><pre>java -jar picard.jar BaitDesigner \\<br />      TARGET=targets.interval_list \\<br />      DESIGN_NAME=new_baits \\<br />      R=reference_sequence.fasta </pre> <hr />";

/// `Double.toString`: plain between 10^-3 and 10^7, computerized scientific notation outside.
fn java_double_to_string(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let magnitude = value.abs();
    if magnitude == 0.0 {
        return format!("{sign}0.0");
    }
    let scientific = format!("{magnitude:e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if (1e-3..1e7).contains(&magnitude) {
        if exponent >= 0 {
            let whole_len = exponent as usize + 1;
            let mut padded = digits.clone();
            while padded.len() < whole_len {
                padded.push('0');
            }
            let (whole, fraction) = padded.split_at(whole_len);
            let fraction = if fraction.is_empty() { "0" } else { fraction };
            format!("{sign}{whole}.{fraction}")
        } else {
            let zeros = "0".repeat((-exponent - 1) as usize);
            format!("{sign}0.{zeros}{digits}")
        }
    } else {
        let (first, rest) = digits.split_at(1);
        let rest = if rest.is_empty() { "0" } else { rest };
        format!("{sign}{first}.{rest}E{exponent}")
    }
}

/// The header of an interval list as `SAMTextHeaderCodec` writes it back: `@HD` with the CURRENT
/// version first and the other attributes after it, every `@SQ` as `SN`, `LN` and then the rest.
type Tags = Vec<(String, String)>;

struct Header {
    hd: Tags,
    sequences: Vec<(String, i32, Tags)>,
    other: Vec<String>,
}

impl Header {
    fn parse(text: &str) -> Header {
        let mut header = Header {
            hd: Vec::new(),
            sequences: Vec::new(),
            other: Vec::new(),
        };
        for line in text.lines().filter(|l| l.starts_with('@')) {
            let fields: Vec<&str> = line.split('\t').collect();
            let tags: Vec<(String, String)> = fields[1..]
                .iter()
                .filter_map(|f| f.split_once(':'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect();
            match fields[0] {
                "@HD" => header.hd = tags,
                "@SQ" => {
                    let name = tags.iter().find(|(k, _)| k == "SN").map(|(_, v)| v.clone());
                    let length = tags
                        .iter()
                        .find(|(k, _)| k == "LN")
                        .and_then(|(_, v)| v.parse().ok());
                    if let (Some(name), Some(length)) = (name, length) {
                        let rest = tags
                            .into_iter()
                            .filter(|(k, _)| k != "SN" && k != "LN")
                            .collect();
                        header.sequences.push((name, length, rest));
                    }
                }
                _ => header.other.push(line.to_string()),
            }
        }
        header
    }

    fn encode(&self) -> String {
        let mut out = String::from("@HD\tVN:1.6");
        for (k, v) in self.hd.iter().filter(|(k, _)| k != "VN") {
            out.push_str(&format!("\t{k}:{v}"));
        }
        out.push('\n');
        for (name, length, rest) in &self.sequences {
            out.push_str(&format!("@SQ\tSN:{name}\tLN:{length}"));
            for (k, v) in rest {
                out.push_str(&format!("\t{k}:{v}"));
            }
            out.push('\n');
        }
        for line in &self.other {
            out.push_str(line);
            out.push('\n');
        }
        out
    }

    fn length_of(&self, contig: &str) -> i32 {
        self.sequences
            .iter()
            .find(|(n, _, _)| n == contig)
            .map(|(_, l, _)| *l)
            .unwrap_or(i32::MAX)
    }

    fn names(&self) -> Vec<String> {
        self.sequences.iter().map(|(n, _, _)| n.clone()).collect()
    }
}

/// The contigs of a FASTA, in file order, with their bases as stored (case included).
fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read(path).unwrap_or_else(|e| fail(&format!("{e}")));
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.split(|b| *b == b'\n') {
        if let Some(name) = line.strip_prefix(b">") {
            let name = String::from_utf8_lossy(name);
            records.push((
                name.split_whitespace().next().unwrap_or("").to_string(),
                Vec::new(),
            ));
        } else if let Some(last) = records.last_mut() {
            last.1
                .extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
        }
    }
    records
}

fn make_target(i: &Interval) -> Target {
    Target {
        contig: i.contig.clone(),
        start: i.start,
        end: i.end,
        negative_strand: i.negative_strand,
        name: i.name.clone().unwrap_or_else(|| "null".to_string()),
    }
}

fn as_interval(bait: &Bait) -> Interval {
    Interval::with_strand_and_name(
        &bait.contig,
        bait.start,
        bait.end,
        bait.negative_strand,
        Some(&bait.name),
    )
}

/// `getBaitSequence`: the primers around the bases as stored, complemented for a reverse copy.
fn bait_sequence(bait: &Bait, options: &Options, reverse: bool) -> String {
    let mut sequence = options.left_primer.clone().into_bytes();
    sequence.extend_from_slice(&bait.bases);
    sequence.extend_from_slice(options.right_primer.as_bytes());
    if reverse {
        sequence = reverse_complement(&sequence);
    }
    String::from_utf8_lossy(&sequence).into_owned()
}

fn write_file(path: &std::path::Path, text: &str) {
    if let Err(e) = std::fs::write(path, text) {
        thrown(&format!(
            "picard.PicardException: Error writing {}: {e}",
            path.display()
        ));
    }
}

fn main() {
    let args = Args::from_env(&[
        ("T", "TARGETS"),
        ("O", "OUTPUT_DIRECTORY"),
        ("R", "REFERENCE_SEQUENCE"),
    ]);
    let targets_path = args.required("TARGETS");
    let design_name = args.required("DESIGN_NAME");
    let reference_path = args.required("REFERENCE_SEQUENCE");

    let strategy = match args.get("DESIGN_STRATEGY").unwrap_or("FixedOffset") {
        "CenteredConstrained" => Strategy::CenteredConstrained,
        "FixedOffset" => Strategy::FixedOffset,
        "Simple" => Strategy::Simple,
        other => fail(&format!(
            "Argument 'DESIGN_STRATEGY' cannot be set to '{other}': invalid value"
        )),
    };
    let options = Options {
        strategy,
        bait_size: args.int("BAIT_SIZE", 120) as i32,
        bait_offset: args.int("BAIT_OFFSET", 80) as i32,
        minimum_baits_per_target: args.int("MINIMUM_BAITS_PER_TARGET", 2) as i32,
        padding: 0,
        merge_nearby_targets: args.bool("MERGE_NEARBY_TARGETS", true),
        design_on_target_strand: args.bool("DESIGN_ON_TARGET_STRAND", false),
        left_primer: args
            .get("LEFT_PRIMER")
            .unwrap_or("ATCGCACCAGCGTGT")
            .to_string(),
        right_primer: args
            .get("RIGHT_PRIMER")
            .unwrap_or("CACTGCGGCTCCTCA")
            .to_string(),
        pool_size: args.int("POOL_SIZE", 55_000).max(0) as usize,
        fill_pools: args.bool("FILL_POOLS", true),
        repeat_tolerance: args.int("REPEAT_TOLERANCE", 50) as i32,
        design_name: design_name.clone(),
    };
    let padding = args.int("PADDING", 0) as i32;
    let output_agilent_files = args.bool("OUTPUT_AGILENT_FILES", true);
    let output_directory = args
        .get("OUTPUT_DIRECTORY")
        .map(str::to_string)
        .unwrap_or_else(|| design_name.clone());

    // `customCommandLineValidation`: a primer is bases and nothing else.
    let mut errors = Vec::new();
    let valid = |p: &str| {
        p.bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
    };
    if !valid(&options.left_primer) {
        errors.push(format!(
            "Left primer {} is not a valid primer sequence.",
            options.left_primer
        ));
    }
    if !valid(&options.right_primer) {
        errors.push(format!(
            "Right primer {} is not a valid primer sequence.",
            options.right_primer
        ));
    }
    if !errors.is_empty() {
        refuse_validation("BaitDesigner", &errors);
    }

    // `IOUtil.assertFileIsReadable` on each, then the directory.
    for path in [&targets_path, &reference_path] {
        if !std::path::Path::new(path).is_file() {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                absolute(path)
            ));
        }
    }
    let directory = std::path::Path::new(&output_directory);
    let _ = std::fs::create_dir_all(directory);

    let targets_text =
        std::fs::read_to_string(&targets_path).unwrap_or_else(|e| fail(&format!("{e}")));
    let header = Header::parse(&targets_text);
    let original = IntervalList::parse_body(header.names(), &targets_text)
        .unwrap_or_else(|e| thrown(&format!("htsjdk.samtools.SAMException: {e:?}")));

    // Padding, then merging, over the targets in the order the file lists them.
    let padded: Vec<Target> = original
        .intervals
        .iter()
        .map(|i| {
            let mut target = make_target(i);
            target.start = (i.start - padding).max(1);
            target.end = (i.end + padding).min(header.length_of(&i.contig));
            target
        })
        .collect();
    let mut merge_options = options.clone();
    merge_options.padding = 0;
    let targets = prepare_targets(&padded, i32::MAX, &merge_options);

    let reference = read_fasta(&reference_path);
    let mut walker = ReferenceWalker::default();
    let mut baits: Vec<Bait> = Vec::new();
    for target in &targets {
        let index = header
            .sequences
            .iter()
            .position(|(n, _, _)| *n == target.contig)
            .map(|i| i as i32)
            .unwrap_or(-1);
        walker.get(index);
        let bases = &reference
            .get(index.max(0) as usize)
            .unwrap_or_else(|| thrown("java.lang.NullPointerException"))
            .1;
        for bait in design_target(target, bases, &options) {
            if bait.length() != options.bait_size {
                thrown(&format!(
                    "picard.PicardException: Bait designed at wrong length: Bait{{name={}, bases={}}}",
                    bait.name,
                    String::from_utf8_lossy(&bait.bases)
                ));
            }
            if masked_base_count(&bait.bases) <= options.repeat_tolerance {
                baits.push(bait);
            }
        }
    }

    // `calculateStatistics`.
    let target_intervals: Vec<Interval> = targets
        .iter()
        .map(|t| {
            Interval::with_strand_and_name(
                &t.contig,
                t.start,
                t.end,
                t.negative_strand,
                Some(&t.name),
            )
        })
        .collect();
    let bait_intervals: Vec<Interval> = baits.iter().map(as_interval).collect();
    let unique_base_count = |intervals: &[Interval]| -> i64 {
        let list = IntervalList {
            dictionary: header.names(),
            intervals: intervals.to_vec(),
        };
        list.uniqued(true)
            .intervals
            .iter()
            .map(|i| i64::from(i.end - i.start + 1))
            .sum()
    };
    let target_territory = unique_base_count(&target_intervals) as i32;
    let target_count = target_intervals.len();
    let bait_territory = unique_base_count(&bait_intervals) as i32;
    let bait_count = bait_intervals.len();
    let design_efficiency = f64::from(target_territory) / f64::from(bait_territory);
    let mut zero_bait_targets = 0;
    let mut intersection: i64 = 0;
    for target in &target_intervals {
        let overlaps: Vec<&Interval> = bait_intervals
            .iter()
            .filter(|b| b.contig == target.contig && b.start <= target.end && b.end >= target.start)
            .collect();
        if overlaps.is_empty() {
            zero_bait_targets += 1;
        } else {
            for bait in overlaps {
                intersection +=
                    i64::from(target.end.min(bait.end) - target.start.max(bait.start) + 1);
            }
        }
    }

    // The files.
    let mut targets_file = header.encode();
    for interval in &original.intervals {
        targets_file.push_str(&interval.to_file_line());
        targets_file.push('\n');
    }
    write_file(
        &directory.join(format!("{design_name}.targets.interval_list")),
        &targets_file,
    );
    let mut baits_file = header.encode();
    for interval in &bait_intervals {
        baits_file.push_str(&interval.to_file_line());
        baits_file.push('\n');
    }
    write_file(
        &directory.join(format!("{design_name}.baits.interval_list")),
        &baits_file,
    );

    let strategy_name = match options.strategy {
        Strategy::CenteredConstrained => "CenteredConstrained",
        Strategy::FixedOffset => "FixedOffset",
        Strategy::Simple => "Simple",
    };
    let mut parameters = String::new();
    for (name, value) in [
        ("USAGE_SUMMARY", USAGE_SUMMARY.to_string()),
        ("USAGE_DETAILS", USAGE_DETAILS.to_string()),
        ("TARGETS", targets_path.clone()),
        ("DESIGN_NAME", design_name.clone()),
        ("LEFT_PRIMER", options.left_primer.clone()),
        ("RIGHT_PRIMER", options.right_primer.clone()),
        ("DESIGN_STRATEGY", strategy_name.to_string()),
        ("BAIT_SIZE", options.bait_size.to_string()),
        (
            "MINIMUM_BAITS_PER_TARGET",
            options.minimum_baits_per_target.to_string(),
        ),
        ("BAIT_OFFSET", options.bait_offset.to_string()),
        ("PADDING", padding.to_string()),
        ("REPEAT_TOLERANCE", options.repeat_tolerance.to_string()),
        ("POOL_SIZE", options.pool_size.to_string()),
        ("FILL_POOLS", options.fill_pools.to_string()),
        (
            "DESIGN_ON_TARGET_STRAND",
            options.design_on_target_strand.to_string(),
        ),
        (
            "MERGE_NEARBY_TARGETS",
            options.merge_nearby_targets.to_string(),
        ),
        ("OUTPUT_AGILENT_FILES", output_agilent_files.to_string()),
        ("OUTPUT_DIRECTORY", output_directory.clone()),
        ("TARGET_TERRITORY", target_territory.to_string()),
        ("TARGET_COUNT", target_count.to_string()),
        ("BAIT_TERRITORY", bait_territory.to_string()),
        ("BAIT_COUNT", bait_count.to_string()),
        (
            "BAIT_TARGET_TERRITORY_INTERSECTION",
            (intersection as i32).to_string(),
        ),
        ("ZERO_BAIT_TARGETS", zero_bait_targets.to_string()),
        (
            "DESIGN_EFFICIENCY",
            java_double_to_string(design_efficiency),
        ),
    ] {
        parameters.push_str(&format!("{name}={value}\n"));
    }
    write_file(
        &directory.join(format!("{design_name}.design_parameters.txt")),
        &parameters,
    );

    let fasta_entry = |bait: &Bait, reverse: bool| {
        format!(
            ">{}\n{}\n",
            bait.name,
            bait_sequence(bait, &options, reverse)
        )
    };
    let mut design_fasta = String::new();
    for bait in &baits {
        design_fasta.push_str(&fasta_entry(bait, false));
    }
    write_file(
        &directory.join(format!("{design_name}.design.fasta")),
        &design_fasta,
    );

    if options.pool_size > 0 && !baits.is_empty() {
        let copies = if options.fill_pools && baits.len() < options.pool_size {
            options.pool_size / baits.len()
        } else {
            1
        };
        let prefix = format!("{}_", &design_name[..design_name.len().min(8)]);
        let mut written = 0usize;
        let mut next_pool = 0usize;
        let mut fasta = String::new();
        let mut agilent = String::new();
        let mut current: Option<String> = None;
        let flush = |current: &Option<String>, fasta: &str, agilent: &str| {
            if let Some(stem) = current {
                write_file(&directory.join(format!("{stem}fasta")), fasta);
                if output_agilent_files {
                    write_file(&directory.join(format!("{stem}txt")), agilent);
                }
            }
        };
        for copy in 0..copies {
            let reverse = copy % 2 == 1;
            let mut bait_id = 1;
            for bait in &baits {
                if written.is_multiple_of(options.pool_size) {
                    flush(&current, &fasta, &agilent);
                    fasta.clear();
                    agilent.clear();
                    current = Some(format!("{design_name}.pool{next_pool}.design."));
                    next_pool += 1;
                }
                written += 1;
                fasta.push_str(&fasta_entry(bait, reverse));
                if output_agilent_files {
                    agilent.push_str(&format!(
                        "{prefix}{bait_id:06}\t{}\n",
                        bait_sequence(bait, &options, reverse).to_uppercase()
                    ));
                    bait_id += 1;
                }
            }
        }
        flush(&current, &fasta, &agilent);
    }
    let _ = std::io::stdout().flush();
}
