//! The command line and the input plumbing shared by the metrics binaries that take a BAM, a
//! reference and a handful of masks: `CollectSequencingArtifactMetrics`, `CollectOxoGMetrics`,
//! `CollectGcBiasMetrics` and `CollectRrbsMetrics`.
//!
//! What is here is what every one of them does before its own collector sees a read:
//!
//! * the arguments, as the `NAME=value` pairs the harness hands a port, under either the long or
//!   the short name, with Barclay's `null` meaning "not given" and a collection appending;
//! * the input, BAM or SAM, told apart by its first two bytes the way htsjdk does;
//! * `SinglePassSamProgram.makeItSo`'s sort-order refusal and `ReferenceSequenceFileWalker`'s
//!   refusal to rewind;
//! * `IntervalListReferenceSequenceMask` and `DbSnpBitSetUtil`, which both tools that take
//!   `INTERVALS` and `DB_SNP` consult per base.
//!
//! Ported from Picard 3.4.0 (`SinglePassSamProgram`, `CommandLineProgram.parseArgs`,
//! `DbSnpBitSetUtil`, `ByIntervalListVariantContextIterator`) and htsjdk 4.2.0
//! (`ReferenceSequenceFileWalker`, `IntervalListReferenceSequenceMask`).

use std::collections::HashMap;
use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::interval::IntervalList;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::tag::{Tag, TagValue};

/// The parsed command line: every `NAME=value` pair, names resolved to their long form.
pub struct Args {
    pairs: Vec<(String, String)>,
}

impl Args {
    /// `std::env::args()` with the short names in `aliases` (short, long) mapped to long ones.
    pub fn from_env(aliases: &[(&str, &str)]) -> Self {
        let mut pairs = Vec::new();
        for raw in std::env::args().skip(1) {
            let raw = raw.trim_start_matches('-').to_string();
            if let Some((name, value)) = raw.split_once('=') {
                let long = aliases
                    .iter()
                    .find(|(short, _)| *short == name)
                    .map(|(_, long)| long.to_string())
                    .unwrap_or_else(|| name.to_string());
                pairs.push((long, value.to_string()));
            }
        }
        Args { pairs }
    }

    /// The last value given for a scalar argument; `null` reads as not given.
    pub fn get(&self, name: &str) -> Option<&str> {
        self.pairs
            .iter()
            .rev()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_str())
            .filter(|v| *v != "null")
    }

    /// Whether the argument was given at all, `null` included: for a scalar whose default is not
    /// null, where `null` and absence differ.
    pub fn given(&self, name: &str) -> bool {
        self.pairs.iter().any(|(n, _)| n == name)
    }

    /// Every value of a collection argument. Barclay appends each one; `null` empties it.
    pub fn all(&self, name: &str) -> Vec<String> {
        let mut out = Vec::new();
        for (n, v) in &self.pairs {
            if n == name {
                if v == "null" {
                    out.clear();
                } else {
                    out.push(v.clone());
                }
            }
        }
        out
    }

    /// A collection argument with a default: each value is APPENDED to the default, the way
    /// Barclay treats a collection field that was initialised, and `null` empties it first.
    pub fn collection(&self, name: &str, default: &[&str]) -> Vec<String> {
        let mut out: Vec<String> = default.iter().map(|s| s.to_string()).collect();
        for (n, v) in &self.pairs {
            if n == name {
                if v == "null" {
                    out.clear();
                } else {
                    out.push(v.clone());
                }
            }
        }
        out
    }

    pub fn required(&self, name: &str) -> String {
        self.get(name)
            .map(str::to_string)
            .unwrap_or_else(|| fail(&format!("Argument '{name}' is required")))
    }

    pub fn int(&self, name: &str, default: i64) -> i64 {
        match self.get(name) {
            None => default,
            Some(v) => v.parse().unwrap_or_else(|_| {
                fail(&format!(
                    "Argument '{name}' cannot be set to '{v}': it is not a number"
                ))
            }),
        }
    }

    pub fn double(&self, name: &str, default: f64) -> f64 {
        match self.get(name) {
            None => default,
            Some(v) => v.parse().unwrap_or_else(|_| {
                fail(&format!(
                    "Argument '{name}' cannot be set to '{v}': it is not a number"
                ))
            }),
        }
    }

    pub fn bool(&self, name: &str, default: bool) -> bool {
        match self.get(name) {
            None => default,
            Some(v) if v.eq_ignore_ascii_case("true") => true,
            Some(v) if v.eq_ignore_ascii_case("false") => false,
            Some(v) => fail(&format!("Argument '{name}' cannot be set to '{v}'")),
        }
    }
}

/// A refusal printed the way the harness reads one, and the run ended at one.
pub fn fail(message: &str) -> ! {
    eprintln!("{message}");
    std::process::exit(1);
}

/// An uncaught Java throwable: `Exception in thread "main" <class>: <message>`.
pub fn thrown(message: &str) -> ! {
    eprintln!("Exception in thread \"main\" {message}");
    std::process::exit(1);
}

/// `customCommandLineValidation`'s messages: the usage, then one message per line.
pub fn refuse_validation(tool: &str, messages: &[String]) -> ! {
    eprintln!("USAGE: {tool} [arguments]\n");
    for message in messages {
        eprintln!("{message}");
    }
    std::process::exit(1);
}

/// The whole input, decoded.
pub fn read_input(path: &str) -> (SamHeader, Vec<BamRecord>) {
    let mut raw = Vec::new();
    if let Err(e) = std::fs::File::open(path).and_then(|mut f| f.read_to_end(&mut raw)) {
        fail(&format!("{e}"));
    }
    if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).unwrap_or_else(|e| fail(&format!("{e:?}")));
        let reader = BamReader::new(&plain).unwrap_or_else(|e| fail(&format!("{e:?}")));
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.unwrap_or_else(|e| fail(&format!("{e:?}"))))
            .collect();
        (header, records)
    } else {
        let text = String::from_utf8(raw).unwrap_or_else(|e| fail(&format!("{e}")));
        read_sam(&text).unwrap_or_else(|e| fail(&format!("{e:?}")))
    }
}

/// `SAMFileHeader.getSortOrder().name()`: an absent or unrecognised `SO` reads as unsorted.
pub fn sort_order(header: &SamHeader) -> &str {
    match header.attributes.get("SO") {
        Some(so @ ("coordinate" | "queryname" | "duplicate" | "unsorted" | "unknown")) => so,
        _ => "unsorted",
    }
}

/// `SinglePassSamProgram.makeItSo`'s sort check: a warning when `ASSUME_SORTED`, else a refusal.
pub fn check_coordinate_sorted(input: &str, header: &SamHeader, assume_sorted: bool) {
    let sort = sort_order(header);
    if sort != "coordinate" && !assume_sorted {
        thrown(&format!(
            "picard.PicardException: File {} should be coordinate sorted but the header says the \
             sort order is {sort}. If you believe the file to be coordinate sorted you may pass \
             ASSUME_SORTED=true",
            absolute(input)
        ));
    }
}

/// `File.getAbsolutePath()`.
pub fn absolute(path: &str) -> String {
    if path.starts_with('/') {
        return path.to_string();
    }
    let cwd = std::env::current_dir()
        .map(|p| p.display().to_string())
        .unwrap_or_default();
    format!("{cwd}/{path}")
}

/// `ReferenceSequenceFileWalker.get`: the walker keeps the contig it last served and refuses a
/// request for an earlier one.
#[derive(Default)]
pub struct ReferenceWalker {
    current: Option<i32>,
}

impl ReferenceWalker {
    pub fn get(&mut self, index: i32) {
        if let Some(current) = self.current {
            if index < current {
                thrown(&format!(
                    "htsjdk.samtools.SAMException: Requesting earlier reference sequence: {index} < {current}"
                ));
            }
        }
        self.current = Some(index);
    }
}

/// The record's read group, from its `RG` tag, if the header declares it.
pub fn read_group<'h>(
    header: &'h SamHeader,
    record: &BamRecord,
) -> Option<&'h htsjdk_bam::header::ReadGroup> {
    let id = match record.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(id)) => id,
        _ => return None,
    };
    header.read_groups.iter().find(|g| g.id == *id)
}

/// `SAMRecord.getOriginalBaseQualities()`: the `OQ` string, FASTQ-decoded.
pub fn original_qualities(record: &BamRecord) -> Option<Vec<u8>> {
    match record.tags.get(Tag::new(b"OQ")) {
        Some(TagValue::Str(text)) => Some(text.bytes().map(|b| b.wrapping_sub(33)).collect()),
        _ => None,
    }
}

/// The `@SQ` names of an interval list's own header, in order.
fn interval_list_dictionary(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| l.starts_with("@SQ"))
        .filter_map(|l| {
            l.split('\t')
                .find_map(|f| f.strip_prefix("SN:"))
                .map(str::to_string)
        })
        .collect()
}

/// `IntervalList.fromFile(file).uniqued()`.
pub fn read_interval_list(path: &str) -> IntervalList {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{e}")));
    let dictionary = interval_list_dictionary(&text);
    IntervalList::parse_body(dictionary, &text)
        .unwrap_or_else(|e| fail(&format!("{e:?}")))
        .uniqued(true)
}

/// `IntervalListReferenceSequenceMask`: which positions of each contig the intervals cover, keyed
/// by the contig's index in the interval list's own dictionary.
pub struct IntervalMask {
    by_index: HashMap<i32, Vec<(i32, i32)>>,
}

impl IntervalMask {
    pub fn new(list: &IntervalList) -> Self {
        let mut by_index: HashMap<i32, Vec<(i32, i32)>> = HashMap::new();
        for interval in &list.intervals {
            let index = list
                .dictionary
                .iter()
                .position(|n| *n == interval.contig)
                .map(|i| i as i32)
                .unwrap_or(-1);
            by_index
                .entry(index)
                .or_default()
                .push((interval.start, interval.end));
        }
        IntervalMask { by_index }
    }

    pub fn get(&self, sequence_index: i32, position: i32) -> bool {
        self.by_index
            .get(&sequence_index)
            .is_some_and(|spans| spans.iter().any(|&(s, e)| s <= position && position <= e))
    }
}

/// `DbSnpBitSetUtil` with no variant types named, which marks every variant's whole span.
pub struct DbSnpMask {
    sites: HashMap<String, Vec<(i32, i32)>>,
}

impl DbSnpMask {
    /// Every record of the VCF, or with `intervals` only those overlapping one of them, which is
    /// what `ByIntervalListVariantContextIterator`'s per-interval queries return between them.
    pub fn load(path: &str, intervals: Option<&IntervalList>) -> Self {
        let text = std::fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{e}")));
        let mut sites: HashMap<String, Vec<(i32, i32)>> = HashMap::new();
        for line in text.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                continue;
            }
            let start: i32 = f[1].parse().unwrap_or(0);
            let end = start + f[3].len() as i32 - 1;
            if let Some(list) = intervals {
                let overlaps = list
                    .intervals
                    .iter()
                    .any(|iv| iv.contig == f[0] && iv.start <= end && start <= iv.end);
                if !overlaps {
                    continue;
                }
            }
            sites
                .entry(f[0].to_string())
                .or_default()
                .push((start, end));
        }
        DbSnpMask { sites }
    }

    pub fn is_db_snp_site(&self, contig: &str, position: i32) -> bool {
        self.sites
            .get(contig)
            .is_some_and(|spans| spans.iter().any(|&(s, e)| s <= position && position <= e))
    }
}
