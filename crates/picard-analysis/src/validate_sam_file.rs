//! `ValidateSamFile`.
//!
//! Ports `picard.sam.ValidateSamFile` + `htsjdk.samtools.SamFileValidator` at tags 3.4.0 / 4.2.0,
//! for the **default VERBOSE mode**, SAM text input, and no reference. In VERBOSE mode the tool
//! prints one raw [`SAMValidationError.toString`] line per error as it is found (no timestamp, no
//! banner), and prints `No errors found` when the file is clean, so the whole output is compared
//! raw.
//!
//! This first slice covers the header checks and the per-record checks that need neither the
//! reference nor cross-record state:
//!
//! - Header: version present/acceptable, empty read groups, missing/invalid read-group platform,
//!   and the empty-sequence-dictionary tracking (`REF_SEQ_TOO_LONG_FOR_BAI` is wired but never
//!   exercised by the corpus).
//! - Per record: a record missing its read group, a mapped read missing its `NM` tag (presence
//!   only, since there is no reference), the `QUAL == *` warning, and the first mapped read under
//!   an empty dictionary.
//!
//! Deferred to follow-up slices (each entangled with a wider htsjdk surface): everything from
//! `SAMRecord.isValid()` (flag/CIGAR/mapping-quality/reference-index checks, and the
//! `READ_GROUP_NOT_FOUND` message it raises), the mate-pair validation
//! (`MATE_NOT_FOUND`/`MISMATCH_MATE_*`/`MATES_ARE_SAME_END`), the sort-order checker
//! (`RECORD_OUT_OF_ORDER`), the reference-based `NM` **value** check (`INVALID_TAG_NM`), the
//! secondary base calls (`E2`/`U2`) and duplicate/`CG` tag checks, the header parser's own
//! validation errors (`HEADER_RECORD_MISSING_REQUIRED_TAG`, `MISSING_VERSION_NUMBER` when `VN` is
//! absent), the SUMMARY mode histogram, BAM/CRAM input, and the quality-format detector.
//!
//! The quality-format detector is safe to defer here: `ValidateSamFile` calls
//! `generateBestGuess(SAM, Standard)`, the expected-quality branch with `checkExpected = false`,
//! which excludes Standard (Phred) only if some observed `QUAL` byte falls outside ASCII `[33,126]`.
//! Every spec-conformant SAM `QUAL` character is in `[33,126]`, so Standard always survives and the
//! tool never emits `INVALID_QUALITY_FORMAT` for conformant input.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::tag::Tag;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

const READ_UNMAPPED: u16 = 0x4;

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

/// Accumulates the verbose report and the running error+warning count (`errorsByType.getCount()`),
/// which decides whether `No errors found` is printed.
struct Report {
    out: String,
    count: usize,
}

impl Report {
    fn new() -> Self {
        Report {
            out: String::new(),
            count: 0,
        }
    }

    /// `out.println(error)` in verbose mode plus `errorsByType.increment`: renders one
    /// `SAMValidationError.toString` (there is never a file source here) and bumps the count.
    fn add(
        &mut self,
        severity: &str,
        ty: &str,
        record_number: Option<i64>,
        read_name: Option<&str>,
        message: &str,
    ) {
        self.count += 1;
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

/// `ValidateSamFile MODE=VERBOSE`, SAM text input, no reference: the raw verbose report.
pub fn validate_sam_file_verbose(input_sam: &str) -> Result<String, ParseError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Silent)?;

    let mut rep = Report::new();
    validate_header(&header, &mut rep);

    // Armed by an empty dictionary, disarmed by the first mapped read it reports on.
    let mut dict_empty_pending = header.sequences.is_empty();

    for (i, rec) in records.iter().enumerate() {
        let record_number = (i + 1) as i64;
        let unmapped = rec.flags & READ_UNMAPPED != 0;

        // validateReadGroup
        if !read_group_present(&header, rec) {
            rep.add(
                WARNING,
                "RECORD_MISSING_READ_GROUP",
                None,
                Some(&rec.read_name),
                "A record is missing a read group",
            );
        }

        // validateNmTag (presence only; the reference-based value check is deferred).
        if !unmapped && rec.tags.get(Tag::new(b"NM")).is_none() {
            rep.add(
                WARNING,
                "MISSING_TAG_NM",
                Some(record_number),
                Some(&rec.read_name),
                "NM tag (nucleotide differences) is missing",
            );
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

    if rep.count == 0 {
        rep.out.push_str("No errors found\n");
    }
    Ok(rep.out)
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
