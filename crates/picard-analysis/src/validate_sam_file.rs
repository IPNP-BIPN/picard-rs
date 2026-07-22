//! `ValidateSamFile`.
//!
//! Ports `picard.sam.ValidateSamFile` + `htsjdk.samtools.SamFileValidator` at tags 3.4.0 / 4.2.0,
//! for the **default VERBOSE mode**, SAM text input, and no reference. In VERBOSE mode the tool
//! prints one raw [`SAMValidationError.toString`] line per error as it is found (no timestamp, no
//! banner), and prints `No errors found` when the file is clean, so the whole output is compared
//! raw.
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
//!
//! Deferred to follow-up slices (each entangled with a wider htsjdk surface): everything from
//! `SAMRecord.isValid()` (flag/CIGAR/mapping-quality/reference-index checks, and the
//! `READ_GROUP_NOT_FOUND` message it raises), the sort-order checker
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

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
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

/// `ValidateSamFile MODE=VERBOSE`, SAM text input, no reference: the raw verbose report.
pub fn validate_sam_file_verbose(input_sam: &str) -> Result<String, ParseError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Silent)?;

    let mut rep = Report::new();
    validate_header(&header, &mut rep);

    // Armed by an empty dictionary, disarmed by the first mapped read it reports on.
    let mut dict_empty_pending = header.sequences.is_empty();

    // Reads awaiting their mate. Keyed, as htsjdk's coordinate-sorted map is, by (reference bucket,
    // read name): a read is stored under the reference index it claims for its mate, and matched
    // when a later read on that reference arrives. A linear vector keeps a deterministic order for
    // the leftover `MATE_NOT_FOUND` pass. (Cross-reference pairing and multi-leftover ordering are
    // deferred; the covered corpus keeps every pair on one reference with at most one leftover.)
    let mut pending: Vec<(i32, String, PairEndInfo)> = Vec::new();

    for (i, rec) in records.iter().enumerate() {
        let record_number = (i + 1) as i64;
        let unmapped = rec.flags & READ_UNMAPPED != 0;

        // validateMateFields: only for paired, primary reads. (The MC-as-valid-cigar check is
        // deferred; the corpus MC tags are valid cigars.)
        if rec.flags & READ_PAIRED != 0 && rec.flags & (SECONDARY | SUPPLEMENTARY) == 0 {
            let found = pending.iter().position(|(bucket, name, _)| {
                *bucket == rec.reference_index && *name == rec.read_name
            });
            if let Some(pos) = found {
                let (_, _, first) = pending.remove(pos);
                let second = PairEndInfo::new(rec, record_number);
                validate_mates(&first, &second, &rec.read_name, &mut rep);
            } else {
                pending.push((
                    rec.mate_reference_index,
                    rec.read_name.clone(),
                    PairEndInfo::new(rec, record_number),
                ));
            }
        }

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

    // validateUnmatchedPairs: reads marked paired whose mate never arrived.
    for (_, name, _) in &pending {
        rep.add(
            ERROR,
            "MATE_NOT_FOUND",
            None,
            Some(name),
            "Mate not found for paired read",
        );
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
