//! `SplitSamByLibrary`.
//!
//! Ports `picard.sam.SplitSamByLibrary.doWork` at tag 3.4.0: split one SAM/BAM into one output per
//! library (the `LB` of each `@RG`), routing every record by its read group's library, plus an
//! `unknown` output for records whose read group has no library (or that carry no known read group).
//! Each output's header is the input header with its `@RG` block filtered to that library's read
//! groups (the `unknown` header keeps the read groups that have no `LB`); everything else in the
//! header is unchanged, no `@PG` is added, and records are written `presorted=true`, so each output
//! is comparable raw. If no `@RG` in the header declares a library the tool errors
//! (`NO_LIBRARIES_SPECIFIED_IN_HEADER`).
//!
//! Output file base names are `IOUtil.makeFileNameSafe(library)` (and `unknown`); the shards are
//! returned in output-name order. `REFERENCE_SEQUENCE`/CRAM are out of scope.

use std::collections::HashMap;

use htsjdk_bam::header::{ReadGroup, SamHeader};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::tag::Tag;
use htsjdk_bam::tag::TagValue;
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

/// Why `SplitSamByLibrary` could not run.
#[derive(Debug)]
pub enum SplitLibError {
    Parse(ParseError),
    /// `NO_LIBRARIES_SPECIFIED_IN_HEADER` (exit code 2): no `@RG` declares a library.
    NoLibrariesInHeader,
}

impl From<ParseError> for SplitLibError {
    fn from(e: ParseError) -> Self {
        SplitLibError::Parse(e)
    }
}

/// `IOUtil.makeFileNameSafe`: trim, then replace each whitespace or reserved punctuation character
/// with `_`.
fn make_file_name_safe(s: &str) -> String {
    const RESERVED: &str = "!\"#$%&'()*/:;<=>?@[]\\^`{|}~";
    s.trim()
        .chars()
        .map(|c| {
            if c.is_whitespace() || RESERVED.contains(c) {
                '_'
            } else {
                c
            }
        })
        .collect()
}

/// The library of a read group, i.e. its `LB` attribute.
fn library_of(rg: &ReadGroup) -> Option<&str> {
    rg.attributes.get("LB")
}

/// A record's read group id, from its `RG` tag.
fn read_group_id(rec: &BamRecord) -> Option<&str> {
    match rec.tags.get(Tag::new(b"RG")) {
        Some(TagValue::Str(id)) => Some(id.as_str()),
        _ => None,
    }
}

/// `SplitSamByLibrary.doWork`: the outputs as `(file base name, SAM text)`, in output-name order.
pub fn split_sam_by_library(input_sam: &str) -> Result<Vec<(String, String)>, SplitLibError> {
    let (header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    // Group read groups by library, preserving header order; those with no LB seed the unknown header.
    let mut lib_order: Vec<String> = Vec::new();
    let mut lib_to_rgs: HashMap<String, Vec<ReadGroup>> = HashMap::new();
    let mut unknown_rgs: Vec<ReadGroup> = Vec::new();
    for rg in &header.read_groups {
        match library_of(rg) {
            Some(lib) => {
                let lib = lib.to_string();
                if !lib_to_rgs.contains_key(&lib) {
                    lib_order.push(lib.clone());
                }
                lib_to_rgs.entry(lib).or_default().push(rg.clone());
            }
            None => unknown_rgs.push(rg.clone()),
        }
    }
    if lib_to_rgs.is_empty() {
        return Err(SplitLibError::NoLibrariesInHeader);
    }

    // id -> library, to route records.
    let mut id_to_lib: HashMap<&str, &str> = HashMap::new();
    for rg in &header.read_groups {
        if let Some(lib) = library_of(rg) {
            id_to_lib.insert(rg.id.as_str(), lib);
        }
    }

    // Route each record into its library's bucket, or unknown.
    let mut lib_records: HashMap<String, Vec<&BamRecord>> = HashMap::new();
    let mut unknown_records: Vec<&BamRecord> = Vec::new();
    for rec in &records {
        match read_group_id(rec).and_then(|id| id_to_lib.get(id)) {
            Some(lib) => lib_records.entry((*lib).to_string()).or_default().push(rec),
            None => unknown_records.push(rec),
        }
    }

    // Build the outputs. A library output is always written (even with no records); the unknown
    // output only when at least one record routed there.
    let mut out: Vec<(String, String)> = Vec::new();
    for lib in &lib_order {
        let mut h = header.clone();
        h.read_groups = lib_to_rgs[lib].clone();
        let recs = lib_records.get(lib).cloned().unwrap_or_default();
        out.push((make_file_name_safe(lib), render(&h, &recs)));
    }
    if !unknown_records.is_empty() {
        let mut h = header.clone();
        h.read_groups = unknown_rgs.clone();
        out.push(("unknown".to_string(), render(&h, &unknown_records)));
    }

    out.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(out)
}

fn render(header: &SamHeader, recs: &[&BamRecord]) -> String {
    let owned: Vec<BamRecord> = recs.iter().map(|r| (*r).clone()).collect();
    write_sam(header, &owned).expect("records re-encode as SAM")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn routes_records_by_library_and_filters_read_groups() {
        let input = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n\
            @RG\tID:rg1\tLB:libA\tSM:s\n@RG\tID:rg2\tLB:libB\tSM:s\n\
            a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
            b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n";
        let out = split_sam_by_library(input).unwrap();
        assert_eq!(out.len(), 2);
        let (names, _): (Vec<_>, Vec<_>) = out.iter().cloned().unzip();
        assert_eq!(names, vec!["libA".to_string(), "libB".to_string()]);
        // libA output carries only rg1 and record a.
        let liba = &out[0].1;
        assert!(liba.contains("@RG\tID:rg1\tLB:libA"));
        assert!(!liba.contains("ID:rg2"));
        assert!(liba.contains("a\t0\tchr1\t1\t") && !liba.contains("b\t0\tchr1"));
    }

    #[test]
    fn records_without_a_library_go_to_unknown() {
        let input = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n\
            @RG\tID:rg1\tLB:libA\tSM:s\n@RG\tID:rg2\tSM:s\n\
            a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
            b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n\
            c\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let out = split_sam_by_library(input).unwrap();
        let names: Vec<_> = out.iter().map(|(n, _)| n.clone()).collect();
        assert_eq!(names, vec!["libA".to_string(), "unknown".to_string()]);
        // unknown keeps rg2 (no LB) in its header and records b and c.
        let unknown = &out[1].1;
        assert!(unknown.contains("@RG\tID:rg2") && !unknown.contains("ID:rg1"));
        assert!(unknown.contains("b\t0\tchr1\t2\t") && unknown.contains("c\t0\tchr1\t3\t"));
    }

    #[test]
    fn no_library_in_header_is_an_error() {
        let input = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s\n\
            a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        assert!(matches!(
            split_sam_by_library(input),
            Err(SplitLibError::NoLibrariesInHeader)
        ));
    }
}
