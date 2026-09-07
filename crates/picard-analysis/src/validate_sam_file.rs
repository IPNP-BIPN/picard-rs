//! `ValidateSamFile`.
//!
//! Ports `picard.sam.ValidateSamFile` + `htsjdk.samtools.SamFileValidator` at tags 3.4.0 / 4.2.0,
//! for **both output modes** (VERBOSE, the default, and SUMMARY), SAM text input, and an optional
//! reference. In VERBOSE mode the tool prints one raw [`SAMValidationError.toString`] line per error
//! as it is found (no timestamp, no banner); in SUMMARY mode it prints a `MetricsFile` histogram of
//! error-type counts (rendered by [`htsjdk_metrics`], with no command-line banner because
//! `SamFileValidator` adds none). Either mode prints `No errors found` when the file is clean, so
//! the whole output is compared raw.
//!
//! Covered so far are the header checks and the per-record checks that need no reference:
//!
//! - Header: version present/acceptable, empty read groups, missing/invalid read-group platform,
//!   and the empty-sequence-dictionary tracking (`REF_SEQ_TOO_LONG_FOR_BAI` is wired but never
//!   exercised by the corpus).
//! - Per record: a record missing its read group, a mapped read missing its `NM` tag (presence
//!   only, since there is no reference), the `QUAL == *` warning, and the first mapped read under
//!   an empty dictionary.
//! - Mate pairs (`SamFileValidator.validateMateFields` / `PairEndInfo.validateMates`): a paired
//!   read whose mate never arrives (`MATE_NOT_FOUND`), and the mate-field mismatches
//!   `MISMATCH_MATE_ALIGNMENT_START`, `MISMATCH_FLAG_MATE_NEG_STRAND`, `MISMATCH_MATE_REF_INDEX`,
//!   `MISMATCH_FLAG_MATE_UNMAPPED`, `MISMATCH_MATE_CIGAR_STRING` (via the `MC` tag), and
//!   `MATES_ARE_SAME_END`. Pairs are matched on one reference, as htsjdk's coordinate-sorted map
//!   does; cross-reference pairing and multi-leftover ordering are deferred.
//! - `SAMRecord.isValid()`, the reference-free / dictionary-independent subset, emitted first per
//!   record in htsjdk's own order: the unpaired-read flag checks (`INVALID_FLAG_PROPER_PAIR`,
//!   `INVALID_FLAG_MATE_UNMAPPED`, `INVALID_FLAG_MATE_NEG_STRAND`, `INVALID_FLAG_FIRST_OF_PAIR`,
//!   `INVALID_FLAG_SECOND_OF_PAIR`), the unmapped-read checks (`INVALID_FLAG_NOT_PRIM_ALIGNMENT`,
//!   `INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT`, `INVALID_MAPPING_QUALITY`), the mapped-read
//!   `INVALID_CIGAR`, and `READ_GROUP_NOT_FOUND`.
//!
//! - The reference-based `NM` **value** check (`INVALID_TAG_NM`), when a `REFERENCE_SEQUENCE` is
//!   given: each mapped read's `NM` tag is compared against the value recomputed from the reference
//!   by [`calculate_md_and_nm`], reusing the same primitive as `SetNmMdAndUqTags`.
//! - The sort-order check (`SAMSortOrderChecker` / `RECORD_OUT_OF_ORDER`): each record is compared
//!   against the previous one under the header's comparator for `SO:coordinate`
//!   (`SAMRecordCoordinateComparator.fileOrderCompare`) and `SO:queryname` (`String.compareTo`).
//!   `unsorted` / `unknown` / a missing `SO` skip the check; `duplicate` is deferred.
//!
//! Deferred to follow-up slices (each entangled with a wider htsjdk surface): the rest of
//! `SAMRecord.isValid()` (the unpaired mate-reference checks, the paired branch's
//! reference/position checks, `INVALID_INSERT_SIZE`, the mapped read's empty-dictionary /
//! missing-reference-name checks, and `isValidReferenceIndexAndPosition`), the `duplicate`
//! sort order, the secondary base calls (`E2`/`U2`) and duplicate/`CG` tag checks, the header parser's own
//! validation errors (`HEADER_RECORD_MISSING_REQUIRED_TAG`, `MISSING_VERSION_NUMBER` when `VN` is
//! absent), BAM/CRAM input, and the quality-format detector.
//!
//! The quality-format detector is safe to defer here: `ValidateSamFile` calls
//! `generateBestGuess(SAM, Standard)`, the expected-quality branch with `checkExpected = false`,
//! which excludes Standard (Phred) only if some observed `QUAL` byte falls outside ASCII `[33,126]`.
//! Every spec-conformant SAM `QUAL` character is in `[33,126]`, so Standard always survives and the
//! tool never emits `INVALID_QUALITY_FORMAT` for conformant input.

use std::collections::{BTreeMap, HashMap};

use htsjdk_bam::fasta::{read_fasta, FastaError};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::md_nm::calculate_md_and_nm;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::tag::Tag;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

use crate::java_hash_map::JavaHashMap;

/// Why `ValidateSamFile` could not run: the input failed to parse, or the reference FASTA did.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValidateError {
    Parse(ParseError),
    Fasta(String),
}

impl From<ParseError> for ValidateError {
    fn from(e: ParseError) -> Self {
        ValidateError::Parse(e)
    }
}

impl From<FastaError> for ValidateError {
    fn from(e: FastaError) -> Self {
        ValidateError::Fasta(format!("{e:?}"))
    }
}

const READ_PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// `SAMFileHeader.ACCEPTABLE_VERSIONS`.
const ACCEPTABLE_VERSIONS: [&str; 5] = ["1.0", "1.3", "1.4", "1.5", "1.6"];

/// `SAMReadGroupRecord.PlatformValue` (compared case-insensitively, as htsjdk uppercases first).
const PLATFORM_VALUES: [&str; 14] = [
    "BGI",
    "CAPILLARY",
    "DNBSEQ",
    "ELEMENT",
    "HELICOS",
    "ILLUMINA",
    "IONTORRENT",
    "LS454",
    "ONT",
    "OTHER",
    "PACBIO",
    "SINGULAR",
    "SOLID",
    "ULTIMA",
];

/// `GenomicIndexUtil.BIN_GENOMIC_SPAN` (512 MiB): a reference longer than this cannot be BAI-indexed.
const BIN_GENOMIC_SPAN: i64 = 512 * 1024 * 1024;

const WARNING: &str = "WARNING";
const ERROR: &str = "ERROR";

/// Accumulates a validation run. In verbose mode `out` gets one `SAMValidationError.toString` line
/// per error; in summary mode nothing is printed per error, but `counts` (keyed by the histogram
/// string `SEVERITY:TYPE`, ordered as htsjdk's `TreeMap` orders it) drives the summary histogram.
/// `count` is `errorsByType.getCount()`, which decides whether `No errors found` is printed.
struct Report {
    out: String,
    count: usize,
    verbose: bool,
    counts: BTreeMap<String, f64>,
    /// `errorsToIgnore`: types dropped before they are counted at all.
    ignore: Vec<String>,
    /// `ignoreWarnings`: a warning is dropped before it is counted, exactly as an ignored type is.
    ignore_warnings: bool,
    /// `maxVerboseOutput`, 100 by default. In verbose mode, reaching it throws
    /// `MaxOutputExceededException`, which ends the whole validation.
    max_output: usize,
    exceeded: bool,
    errors: usize,
    warnings: usize,
}

impl Report {
    fn new(verbose: bool) -> Self {
        Report {
            out: String::new(),
            count: 0,
            verbose,
            counts: BTreeMap::new(),
            ignore: Vec::new(),
            ignore_warnings: false,
            max_output: 100,
            exceeded: false,
            errors: 0,
            warnings: 0,
        }
    }

    /// `errorsByType.increment` plus, in verbose mode, `out.println(error)`: bumps the count and the
    /// per-type bin, and renders one `SAMValidationError.toString` (there is never a file source
    /// here) when verbose.
    fn add(
        &mut self,
        severity: &str,
        ty: &str,
        record_number: Option<i64>,
        read_name: Option<&str>,
        message: &str,
    ) {
        // `addError`: an ignored type never reaches the counters, and neither does a warning when
        // warnings are ignored. Once the verbose cap has been reached htsjdk is unwinding through
        // `MaxOutputExceededException`, so nothing further is recorded at all.
        if self.exceeded || self.ignore.iter().any(|t| t == ty) {
            return;
        }
        if severity == WARNING {
            if self.ignore_warnings {
                return;
            }
            self.warnings += 1;
        } else {
            self.errors += 1;
        }
        self.count += 1;
        *self.counts.entry(format!("{severity}:{ty}")).or_insert(0.0) += 1.0;
        if !self.verbose {
            return;
        }
        self.out.push_str(severity);
        self.out.push_str("::");
        self.out.push_str(ty);
        self.out.push(':');
        if let Some(n) = record_number {
            if n > 0 {
                self.out.push_str(&format!("Record {n}, "));
            }
        }
        if let Some(r) = read_name {
            self.out.push_str(&format!("Read name {r}, "));
        }
        self.out.push_str(message);
        self.out.push('\n');
        if self.count >= self.max_output {
            self.exceeded = true;
        }
    }
}

/// `SamFileValidator.validateHeader`, restricted to the reference-free checks in scope.
fn validate_header(header: &SamHeader, rep: &mut Report) {
    // Version: getVersion() is the @HD VN attribute. (A missing VN is a parser-level error here and
    // is deferred; if present it must be one of the acceptable versions.)
    match header.attributes.get("VN") {
        None => rep.add(
            ERROR,
            "MISSING_VERSION_NUMBER",
            None,
            None,
            "Header has no version number",
        ),
        Some(v) if !ACCEPTABLE_VERSIONS.contains(&v) => rep.add(
            ERROR,
            "INVALID_VERSION_NUMBER",
            None,
            None,
            &format!(
                "Header version: {v} does not match any of the acceptable versions: {}",
                ACCEPTABLE_VERSIONS.join(", ")
            ),
        ),
        Some(_) => {}
    }

    // Sequence dictionary: an empty one only arms a warning that fires on the first mapped read.
    if !header.sequences.is_empty() {
        let long: Vec<&str> = header
            .sequences
            .iter()
            .filter(|s| s.length as i64 > BIN_GENOMIC_SPAN)
            .map(|s| s.name.as_str())
            .collect();
        if !long.is_empty() {
            rep.add(
                WARNING,
                "REF_SEQ_TOO_LONG_FOR_BAI",
                None,
                None,
                &format!(
                    "Reference sequences are too long for BAI indexing: {}",
                    long.join(", ")
                ),
            );
        }
    }

    if header.read_groups.is_empty() {
        rep.add(
            ERROR,
            "MISSING_READ_GROUP",
            None,
            None,
            "Read groups is empty",
        );
    }

    // Read groups: duplicate id, then missing / invalid platform.
    let mut seen: Vec<&str> = Vec::new();
    for rg in &header.read_groups {
        let id = rg.id.as_str();
        if seen.contains(&id) {
            rep.add(
                ERROR,
                "DUPLICATE_READ_GROUP_ID",
                None,
                None,
                &format!("Duplicate read group id: {id}"),
            );
        } else {
            seen.push(id);
        }

        match rg.attributes.get("PL") {
            None | Some("") => rep.add(
                ERROR,
                "MISSING_PLATFORM_VALUE",
                None,
                Some(id),
                "A platform (PL) attribute was not found for read group ",
            ),
            Some(pl) if !PLATFORM_VALUES.contains(&pl.to_ascii_uppercase().as_str()) => rep.add(
                ERROR,
                "INVALID_PLATFORM_VALUE",
                None,
                Some(id),
                &format!(
                    "The platform (PL) attribute ({pl}) + was not one of the valid values for read group "
                ),
            ),
            Some(_) => {}
        }
    }
}

/// `SamFileValidator.validateReadGroup`: the record's read group is unknown if it has no `RG` tag or
/// the tag's id is not in the header.
fn read_group_present(header: &SamHeader, rec: &BamRecord) -> bool {
    match rec.tags.get(Tag::new(b"RG")) {
        Some(htsjdk_bam::tag::TagValue::Str(id)) => {
            header.read_groups.iter().any(|rg| rg.id == *id)
        }
        _ => false,
    }
}

/// `SamFileValidator.PairEndInfo`: the per-read view kept while waiting to meet a read's mate,
/// carrying both the read's own fields and what the read asserts about its mate.
struct PairEndInfo {
    read_alignment_start: i32,
    read_reference_index: i32,
    read_neg_strand: bool,
    read_unmapped: bool,
    read_cigar: String,
    mate_alignment_start: i32,
    mate_reference_index: i32,
    mate_neg_strand: bool,
    mate_unmapped: bool,
    mate_cigar: Option<String>,
    first_of_pair: bool,
    record_number: i64,
}

impl PairEndInfo {
    fn new(rec: &BamRecord, record_number: i64) -> Self {
        let mate_cigar = match rec.tags.get(Tag::new(b"MC")) {
            Some(htsjdk_bam::tag::TagValue::Str(s)) => Some(s.clone()),
            _ => None,
        };
        PairEndInfo {
            read_alignment_start: rec.alignment_start,
            read_reference_index: rec.reference_index,
            read_neg_strand: rec.flags & READ_REVERSE != 0,
            read_unmapped: rec.flags & READ_UNMAPPED != 0,
            read_cigar: rec.cigar.to_text(),
            mate_alignment_start: rec.mate_alignment_start,
            mate_reference_index: rec.mate_reference_index,
            mate_neg_strand: rec.flags & MATE_REVERSE != 0,
            mate_unmapped: rec.flags & MATE_UNMAPPED != 0,
            mate_cigar,
            first_of_pair: rec.flags & FIRST_OF_PAIR != 0,
            record_number,
        }
    }
}

/// `PairEndInfo.validateMateFields(end1, end2)`: the mate fields `end1` asserts must agree with
/// `end2`'s own fields. All errors carry `end1`'s record number.
fn validate_mate_fields(end1: &PairEndInfo, end2: &PairEndInfo, read_name: &str, rep: &mut Report) {
    let rn = Some(end1.record_number);
    if end1.mate_alignment_start != end2.read_alignment_start {
        rep.add(
            ERROR,
            "MISMATCH_MATE_ALIGNMENT_START",
            rn,
            Some(read_name),
            "Mate alignment does not match alignment start of mate",
        );
    }
    if end1.mate_neg_strand != end2.read_neg_strand {
        rep.add(
            ERROR,
            "MISMATCH_FLAG_MATE_NEG_STRAND",
            rn,
            Some(read_name),
            "Mate negative strand flag does not match read negative strand flag of mate",
        );
    }
    if end1.mate_reference_index != end2.read_reference_index {
        rep.add(
            ERROR,
            "MISMATCH_MATE_REF_INDEX",
            rn,
            Some(read_name),
            "Mate reference index (MRNM) does not match reference index of mate",
        );
    }
    if end1.mate_unmapped != end2.read_unmapped {
        rep.add(
            ERROR,
            "MISMATCH_FLAG_MATE_UNMAPPED",
            rn,
            Some(read_name),
            "Mate unmapped flag does not match read unmapped flag of mate",
        );
    }
    if let Some(mc) = &end1.mate_cigar {
        if mc != &end2.read_cigar {
            rep.add(
                ERROR,
                "MISMATCH_MATE_CIGAR_STRING",
                rn,
                Some(read_name),
                "Mate CIGAR string does not match CIGAR string of mate",
            );
        }
    }
}

/// `PairEndInfo.validateMates`: both directions, then the both-marked-same-end check (reported once,
/// against the first-seen read's record number).
fn validate_mates(first: &PairEndInfo, second: &PairEndInfo, read_name: &str, rep: &mut Report) {
    validate_mate_fields(first, second, read_name, rep);
    validate_mate_fields(second, first, read_name, rep);
    if first.first_of_pair == second.first_of_pair {
        let which = if first.first_of_pair {
            "first"
        } else {
            "second"
        };
        rep.add(
            ERROR,
            "MATES_ARE_SAME_END",
            Some(first.record_number),
            Some(read_name),
            &format!("Both mates are marked as {which} of pair"),
        );
    }
}

/// `SAMRecord.isValid`, restricted to the reference-free, dictionary-independent flag / mapping /
/// CIGAR / read-group checks, emitted in htsjdk's own order. Every error carries the record number
/// (`SamFileValidator` calls `setRecordNumber` on each). The mate-reference checks for unpaired
/// reads, the paired branch's reference/position checks, `INVALID_INSERT_SIZE`, the mapped read's
/// empty-dictionary / missing-reference-name checks, and `isValidReferenceIndexAndPosition` are
/// deferred (each needs the mate reference, the insert-size bound, or the sequence dictionary).
fn is_valid_record(header: &SamHeader, rec: &BamRecord, record_number: i64, rep: &mut Report) {
    let rn = Some(record_number);
    let name = Some(rec.read_name.as_str());
    let paired = rec.flags & READ_PAIRED != 0;
    let unmapped = rec.flags & READ_UNMAPPED != 0;

    if !paired {
        if rec.flags & PROPER_PAIR != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_PROPER_PAIR",
                rn,
                name,
                "Proper pair flag should not be set for unpaired read.",
            );
        }
        if rec.flags & MATE_UNMAPPED != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_MATE_UNMAPPED",
                rn,
                name,
                "Mate unmapped flag should not be set for unpaired read.",
            );
        }
        if rec.flags & MATE_REVERSE != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_MATE_NEG_STRAND",
                rn,
                name,
                "Mate negative strand flag should not be set for unpaired read.",
            );
        }
        if rec.flags & FIRST_OF_PAIR != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_FIRST_OF_PAIR",
                rn,
                name,
                "First of pair flag should not be set for unpaired read.",
            );
        }
        if rec.flags & SECOND_OF_PAIR != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_SECOND_OF_PAIR",
                rn,
                name,
                "Second of pair flag should not be set for unpaired read.",
            );
        }
    }

    if unmapped {
        if rec.flags & SECONDARY != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_NOT_PRIM_ALIGNMENT",
                rn,
                name,
                "Secondary alignment flag should not be set for unmapped read.",
            );
        }
        if rec.flags & SUPPLEMENTARY != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT",
                rn,
                name,
                "Supplementary alignment flag should not be set for unmapped read.",
            );
        }
        if rec.mapping_quality != 0 {
            rep.add(
                ERROR,
                "INVALID_MAPPING_QUALITY",
                rn,
                name,
                "MAPQ should be 0 for unmapped read.",
            );
        }
    } else if rec.cigar.elements.is_empty() {
        // (MAPQ >= 256 is unreachable: the field is a single byte.)
        rep.add(
            ERROR,
            "INVALID_CIGAR",
            rn,
            name,
            "CIGAR should have > zero elements for mapped read.",
        );
    }

    // The RG ID, when present, must resolve in the header.
    if let Some(htsjdk_bam::tag::TagValue::Str(id)) = rec.tags.get(Tag::new(b"RG")) {
        if !header.read_groups.iter().any(|rg| rg.id == *id) {
            rep.add(
                ERROR,
                "READ_GROUP_NOT_FOUND",
                rn,
                name,
                &format!("RG ID on SAMRecord not found in header: {id}"),
            );
        }
    }
}

/// The header sort orders that carry a comparator (`SortOrder.getComparatorInstance`). `unsorted`,
/// `unknown`, a missing `SO`, and `duplicate` (deferred) get no order check.
enum SortOrder {
    Coordinate,
    Queryname,
    Unchecked,
}

impl SortOrder {
    fn of(header: &SamHeader) -> Self {
        match header.attributes.get("SO") {
            Some("coordinate") => SortOrder::Coordinate,
            Some("queryname") => SortOrder::Queryname,
            _ => SortOrder::Unchecked,
        }
    }

    fn name(&self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
            SortOrder::Unchecked => "unsorted",
        }
    }

    /// `SAMRecordCoordinateComparator` / `SAMRecordQueryNameComparator` `fileOrderCompare`. `prev` is
    /// in order iff this is `<= 0`.
    fn file_order_compare(&self, prev: &BamRecord, rec: &BamRecord) -> i32 {
        match self {
            SortOrder::Coordinate => {
                let (r1, r2) = (prev.reference_index, rec.reference_index);
                if r1 == -1 {
                    return if r2 == -1 { 0 } else { 1 };
                }
                if r2 == -1 {
                    return -1;
                }
                if r1 != r2 {
                    return r1 - r2;
                }
                prev.alignment_start - rec.alignment_start
            }
            // compareReadNames is String.compareTo, i.e. UTF-16 code-unit order, which equals Rust's
            // byte order for the ASCII read names in practice.
            SortOrder::Queryname => match prev.read_name.cmp(&rec.read_name) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            },
            SortOrder::Unchecked => 0,
        }
    }
}

/// The tool's arguments, with Picard's defaults.
#[derive(Debug, Clone)]
pub struct Options {
    /// `MODE`: VERBOSE prints one line per error, SUMMARY prints the histogram of counts.
    pub verbose: bool,
    /// `IGNORE`: error types dropped before they are counted.
    pub ignore: Vec<String>,
    /// `IGNORE_WARNINGS`: drop every warning, which also removes it from the exit code.
    pub ignore_warnings: bool,
    /// `SKIP_MATE_VALIDATION`: skip the mate checks and the unmatched-pair pass.
    pub skip_mate_validation: bool,
    /// `MAX_OUTPUT`, 100: in verbose mode, the error that reaches this count is the last one.
    pub max_output: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            verbose: true,
            ignore: Vec::new(),
            ignore_warnings: false,
            skip_mate_validation: false,
            max_output: 100,
        }
    }
}

/// What a validation found: its report, and the counts the exit code is derived from.
pub struct Validation {
    pub report: String,
    pub errors: usize,
    pub warnings: usize,
}

impl Validation {
    /// `ValidateSamFile.ReturnTypes`: 0 clean, 1 warnings only, 2 errors and warnings, 3 errors
    /// only. It is the tool's answer rather than a failure, which is why the row records it.
    pub fn exit_code(&self) -> i32 {
        match (self.errors, self.warnings) {
            (0, 0) => 0,
            (0, _) => 1,
            (_, 0) => 3,
            _ => 2,
        }
    }
}

/// `ValidateSamFile` over already-parsed records, with the tool's arguments.
pub fn validate_records(
    header: &SamHeader,
    records: &[BamRecord],
    fasta: Option<&[u8]>,
    options: &Options,
) -> Result<Validation, ValidateError> {
    validate_inner(header, records, fasta, options)
}

/// `ValidateSamFile MODE=VERBOSE`, SAM text input, no reference: the raw verbose report.
pub fn validate_sam_file_verbose(input_sam: &str) -> Result<String, ValidateError> {
    validate(input_sam, None, true)
}

/// `ValidateSamFile MODE=VERBOSE REFERENCE_SEQUENCE=<fasta>`: as above, plus the reference-based
/// `NM` value check (`INVALID_TAG_NM`). `fasta` is the `REFERENCE_SEQUENCE` bytes.
pub fn validate_sam_file_verbose_with_reference(
    input_sam: &str,
    fasta: &[u8],
) -> Result<String, ValidateError> {
    validate(input_sam, Some(fasta), true)
}

/// `ValidateSamFile MODE=SUMMARY`, SAM text input, no reference: the summary histogram (or
/// `No errors found`).
pub fn validate_sam_file_summary(input_sam: &str) -> Result<String, ValidateError> {
    validate(input_sam, None, false)
}

/// `ValidateSamFile MODE=SUMMARY REFERENCE_SEQUENCE=<fasta>`: the summary histogram, with the
/// reference-based `NM` value check included in the counts.
pub fn validate_sam_file_summary_with_reference(
    input_sam: &str,
    fasta: &[u8],
) -> Result<String, ValidateError> {
    validate(input_sam, Some(fasta), false)
}

/// `SamFileValidator.validateSamFileVerbose` / `validateSamFileSummary`. When `fasta` is present,
/// each mapped read's `NM` tag is checked against the value recomputed from the reference
/// (`INVALID_TAG_NM`); otherwise only its presence is checked (`MISSING_TAG_NM`), exactly as htsjdk
/// skips the value check without a reference. `verbose` selects the per-error report vs the summary
/// histogram.
fn validate(input_sam: &str, fasta: Option<&[u8]>, verbose: bool) -> Result<String, ValidateError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Silent)?;
    let options = Options {
        verbose,
        ..Options::default()
    };
    Ok(validate_inner(&header, &records, fasta, &options)?.report)
}

/// The validation itself, over parsed records.
fn validate_inner(
    header: &SamHeader,
    records: &[BamRecord],
    fasta: Option<&[u8]>,
    options: &Options,
) -> Result<Validation, ValidateError> {
    let verbose = options.verbose;

    // The reference bases by contig name, resolved once, for the `NM` value check.
    let contigs = match fasta {
        Some(bytes) => read_fasta(bytes)?,
        None => Vec::new(),
    };
    let ref_by_name: HashMap<&str, &[u8]> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.bases.as_slice()))
        .collect();
    let have_reference = fasta.is_some();

    let mut rep = Report::new(verbose);
    rep.ignore = options.ignore.clone();
    rep.ignore_warnings = options.ignore_warnings;
    rep.max_output = options.max_output;
    validate_header(header, &mut rep);

    // Armed by an empty dictionary, disarmed by the first mapped read it reports on.
    let mut dict_empty_pending = header.sequences.is_empty();

    // Sort-order checking: the comparator from the header's SO, and the previous record seen.
    let sort_order = SortOrder::of(header);
    let mut prev_record: Option<&BamRecord> = None;

    // Reads awaiting their mate. htsjdk keeps them in an `InMemoryPairEndInfoMap`, which is a
    // `java.util.HashMap<String, PairEndInfo>` keyed on the read name alone -- its `remove` takes a
    // reference index and ignores it -- and reports the leftovers by ITERATING that map. So the
    // order of the `MATE_NOT_FOUND` lines is Java's bucket order, which is why this is a
    // [`JavaHashMap`] and not a vector: the names are the same either way and the file is not.
    let mut pending: JavaHashMap<PairEndInfo> = JavaHashMap::new();

    for (i, rec) in records.iter().enumerate() {
        let record_number = (i + 1) as i64;
        let unmapped = rec.flags & READ_UNMAPPED != 0;

        // isValid(): the per-record flag / mapping / CIGAR / read-group checks, emitted first.
        is_valid_record(header, rec, record_number, &mut rep);

        // validateMateFields: only for paired, primary reads, and only when the caller has not
        // asked for the mate checks to be skipped. (The MC-as-valid-cigar check is deferred; the
        // corpus MC tags are valid cigars.)
        if !options.skip_mate_validation
            && rec.flags & READ_PAIRED != 0
            && rec.flags & (SECONDARY | SUPPLEMENTARY) == 0
        {
            match pending.remove(&rec.read_name) {
                Some(first) => {
                    let second = PairEndInfo::new(rec, record_number);
                    validate_mates(&first, &second, &rec.read_name, &mut rep);
                }
                None => pending.put(&rec.read_name, PairEndInfo::new(rec, record_number)),
            }
        }

        // validateSortOrder: compare against the previous record under the header's comparator.
        if let Some(prev) = prev_record {
            if sort_order.file_order_compare(prev, rec) > 0 {
                rep.add(
                    ERROR,
                    "RECORD_OUT_OF_ORDER",
                    Some(record_number),
                    Some(&rec.read_name),
                    &format!(
                        "The record is out of [{}] order, prior read name [{}], prior coodinates [{}:{}]",
                        sort_order.name(),
                        prev.read_name,
                        prev.reference_index,
                        prev.alignment_start,
                    ),
                );
            }
        }
        prev_record = Some(rec);

        // validateReadGroup
        if !read_group_present(header, rec) {
            rep.add(
                WARNING,
                "RECORD_MISSING_READ_GROUP",
                None,
                Some(&rec.read_name),
                "A record is missing a read group",
            );
        }

        // validateNmTag: the tag must be present (MISSING_TAG_NM) and, when a reference is given,
        // must equal the value recomputed from the reference (INVALID_TAG_NM).
        if !unmapped {
            match rec.tags.get(Tag::new(b"NM")) {
                None => rep.add(
                    WARNING,
                    "MISSING_TAG_NM",
                    Some(record_number),
                    Some(&rec.read_name),
                    "NM tag (nucleotide differences) is missing",
                ),
                Some(htsjdk_bam::tag::TagValue::Int(tag_nm)) if have_reference => {
                    let name = &header.sequences[rec.reference_index as usize].name;
                    if let Some(ref_bases) = ref_by_name.get(name.as_str()) {
                        let (_, actual) = calculate_md_and_nm(
                            rec.alignment_start,
                            &rec.cigar,
                            &rec.read_bases,
                            ref_bases,
                        );
                        if *tag_nm != actual as i64 {
                            rep.add(
                                ERROR,
                                "INVALID_TAG_NM",
                                Some(record_number),
                                Some(&rec.read_name),
                                &format!(
                                    "NM tag (nucleotide differences) in file [{tag_nm}] does not match reality [{actual}]"
                                ),
                            );
                        }
                    }
                }
                Some(_) => {}
            }
        }

        // Empty dictionary reported once, on the first mapped read.
        if dict_empty_pending && !unmapped {
            rep.add(
                ERROR,
                "MISSING_SEQUENCE_DICTIONARY",
                None,
                None,
                "Sequence dictionary is empty",
            );
            dict_empty_pending = false;
        }

        // QUAL == '*' (no stored qualities).
        if rec.base_qualities.is_empty() {
            rep.add(
                WARNING,
                "QUALITY_NOT_STORED",
                Some(record_number),
                Some(&rec.read_name),
                "QUAL field is set to * (unspecified quality scores), this is allowed by the SAM \
                 specification but many tools expect reads to include qualities ",
            );
        }
    }

    // validateUnmatchedPairs: reads marked paired whose mate never arrived. Skipped with the mate
    // checks, because `pending` is only filled when they run.
    for (name, _) in pending.iter() {
        rep.add(
            ERROR,
            "MATE_NOT_FOUND",
            None,
            Some(name),
            "Mate not found for paired read",
        );
    }

    // A clean file prints `No errors found` in either mode. In verbose mode the per-error lines are
    // already in `out`; in summary mode the errors, if any, become the histogram.
    let (errors, warnings) = (rep.errors, rep.warnings);
    if rep.count == 0 {
        return Ok(Validation {
            report: "No errors found\n".to_string(),
            errors,
            warnings,
        });
    }
    if verbose {
        let mut report = rep.out;
        if rep.exceeded {
            // The catch in `validateSamFileVerbose`, which is where the run ends.
            report.push_str(&format!(
                "Maximum output of [{}] errors reached.\n",
                rep.max_output
            ));
        }
        return Ok(Validation {
            report,
            errors,
            warnings,
        });
    }
    let mut metrics = htsjdk_metrics::file::MetricsFile::new();
    metrics.histograms.push(htsjdk_metrics::file::Histogram {
        bin_label: "Error Type".to_string(),
        value_label: "Count".to_string(),
        key_class: "java.lang.String".to_string(),
        bins: rep.counts.into_iter().collect(),
    });
    Ok(Validation {
        report: metrics.write(),
        errors,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_clean_file_reports_no_errors() {
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n";
        assert_eq!(validate_sam_file_verbose(sam).unwrap(), "No errors found\n");
    }

    #[test]
    fn an_unpaired_read_with_a_pairing_flag_is_invalid() {
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t2\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n";
        assert_eq!(
            validate_sam_file_verbose(sam).unwrap(),
            "ERROR::INVALID_FLAG_PROPER_PAIR:Record 1, Read name a, Proper pair flag should not be set for unpaired read.\n"
        );
    }

    #[test]
    fn a_record_with_an_unknown_read_group_is_reported() {
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:other\tNM:i:0\n";
        assert_eq!(
            validate_sam_file_verbose(sam).unwrap(),
            "ERROR::READ_GROUP_NOT_FOUND:Record 1, Read name a, RG ID on SAMRecord not found in header: other\n\
             WARNING::RECORD_MISSING_READ_GROUP:Read name a, A record is missing a read group\n"
        );
    }

    #[test]
    fn summary_mode_prints_a_histogram_of_error_types() {
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t2\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n";
        // One ERROR: the proper-pair flag on an unpaired read.
        let out = validate_sam_file_summary(sam).unwrap();
        assert_eq!(
            out,
            "\n\n## HISTOGRAM\tjava.lang.String\nError Type\tCount\n\
             ERROR:INVALID_FLAG_PROPER_PAIR\t1\n\n"
        );
        // A clean file prints the same line in either mode.
        let clean = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n";
        assert_eq!(
            validate_sam_file_summary(clean).unwrap(),
            "No errors found\n"
        );
    }

    #[test]
    fn a_coordinate_out_of_order_record_is_reported() {
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t0\tchr1\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n\
            b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n";
        assert_eq!(
            validate_sam_file_verbose(sam).unwrap(),
            "ERROR::RECORD_OUT_OF_ORDER:Record 2, Read name b, The record is out of [coordinate] order, prior read name [a], prior coodinates [0:50]\n"
        );
    }

    #[test]
    fn a_wrong_nm_tag_is_reported_against_the_reference() {
        let fasta: &[u8] = b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:40\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACCTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\n";
        assert_eq!(
            validate_sam_file_verbose_with_reference(sam, fasta).unwrap(),
            "ERROR::INVALID_TAG_NM:Record 1, Read name a, NM tag (nucleotide differences) in file [0] does not match reality [1]\n"
        );
        // Without a reference the value check is skipped and the correct-presence read is clean.
        assert_eq!(validate_sam_file_verbose(sam).unwrap(), "No errors found\n");
    }

    #[test]
    fn an_absent_mate_is_reported() {
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            p\t99\tchr1\t10\t60\t4M\t=\t20\t14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n";
        assert_eq!(
            validate_sam_file_verbose(sam).unwrap(),
            "ERROR::MATE_NOT_FOUND:Read name p, Mate not found for paired read\n"
        );
    }

    #[test]
    fn a_mapped_read_without_nm_warns() {
        let sam = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
            @RG\tID:rg1\tSM:s\tPL:illumina\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        assert_eq!(
            validate_sam_file_verbose(sam).unwrap(),
            "WARNING::MISSING_TAG_NM:Record 1, Read name a, NM tag (nucleotide differences) is missing\n"
        );
    }
}
