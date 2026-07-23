//! `MergeSamFiles`.
//!
//! Ports `picard.sam.MergeSamFiles.doWork` at tag 3.4.0 for **inputs that share a sequence
//! dictionary**: merge several SAM files into one, sorted by `SORT_ORDER`. The tool builds the output
//! header with `SamFileHeaderMerger` and merges the records with a `MergingSamRecordIterator` (a k-way
//! merge on the `SORT_ORDER` comparator), and it adds **no** `@PG` (only optional `@CO`, default
//! none), so the output compares raw.
//!
//! The merged header always gains the one change `SamFileHeaderMerger` makes
//! (`SamFileHeaderMerger.java:185`): the group order is set to `none`. htsjdk builds the merged `@HD`
//! fresh, so its attributes come out in insertion order `VN`, `GO`, `SO` (not the input's `VN`, `SO`).
//!
//! `@RG` and `@PG` records are merged by `SamFileHeaderMerger.mergeHeaderRecords`
//! (`SamFileHeaderMerger.java:415`): records are binned by ID in file order, and within an ID by
//! identical content. An ID seen only with one content is kept as-is; an ID reused with **different**
//! content collides, and every content after the first is **renamed** `ID.<n>` where `n` is a
//! base36 render (`positiveFourDigitBase36Str`, digits `0-9a-z`) of a per-merge counter. The `@RG` and
//! `@PG` results are each finally **sorted by ID** (`RECORD_ID_COMPARATOR`, `String.compareTo`). When a
//! record's read group (`RG:Z`) or program (`PG:Z`) was renamed, the reading iterator rewrites that tag
//! to the new ID; this port remaps each input's records the same way before merging them. `@CO` is
//! every input's comments concatenated (the merger does not dedupe comments). Records from all inputs
//! are concatenated and **stably** sorted by the `SORT_ORDER` comparator; the coordinate/queryname
//! comparators fully order any two distinct records, so this reproduces the k-way merge except for
//! wholly-identical records, where the stable sort keeps input order (first file first), as the merge
//! does.
//!
//! The sequence dictionaries either match (`SequenceUtil.assertSequenceListsEqual`) or, with
//! `MERGE_SEQUENCE_DICTIONARIES`, are unioned into the name superset (`mergeSequences`); when unioned,
//! each input's records have their `reference_index`/`mate_reference_index` remapped by name, as
//! `MergingSamRecordIterator.setHeader` re-resolves them. A mismatch without the flag, or a shared
//! sequence ordered incompatibly across inputs, is an error.
//!
//! Scope of this slice: `@PG` records carry no `PP` chain. Deferred to a further slice: the `@PG`
//! previous-program (`PP`) tree merge, which rewrites `PP` pointers through the same rename table
//! across processing passes; it is reported as an error here rather than merged.

use std::collections::{HashMap, HashSet};

use htsjdk_bam::header::{ProgramRecord, ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::ParseError;

use crate::sort_sam::SortOrder;

/// Why a merge could not run.
#[derive(Debug)]
pub enum MergeError {
    Parse(ParseError),
    /// One input declares two `@RG` records with the same ID. `mergeReadGroups` throws
    /// `contains more than one RG with the same id`.
    DuplicateReadGroupId(String),
    /// The same for a `@PG` ID.
    DuplicateProgramId(String),
    /// A `@PG` record carries a `PP` (previous-program) attribute. The `PP`-chain tree merge is a
    /// separate slice.
    ProgramChainNotSupported,
    /// The inputs' sequence dictionaries differ and `MERGE_SEQUENCE_DICTIONARIES` was not set.
    SequenceDictionaryMismatch,
    /// `MERGE_SEQUENCE_DICTIONARIES` was set but two dictionaries order a shared sequence
    /// incompatibly, so `mergeSequences` cannot find a consistent superset.
    SequenceDictionaryOrderConflict,
}

impl From<ParseError> for MergeError {
    fn from(e: ParseError) -> Self {
        MergeError::Parse(e)
    }
}

/// `positiveFourDigitBase36Str`: render `n` in base 36 with digits `0-9a-z`, no leading zeros
/// (`0` -> `"0"`, `1` -> `"1"`, `36` -> `"10"`). Despite the name it is not zero-padded.
fn base36(mut n: u32) -> String {
    const DIGITS: &[u8; 36] = b"0123456789abcdefghijklmnopqrstuvwxyz";
    if n == 0 {
        return "0".to_string();
    }
    let mut buf = Vec::new();
    while n > 0 {
        buf.push(DIGITS[(n % 36) as usize]);
        n /= 36;
    }
    buf.reverse();
    String::from_utf8(buf).expect("base36 digits are ASCII")
}

/// One merged header record (`@RG` or `@PG`) reduced to the two things the merge needs: its ID and a
/// content key that is equal exactly when `SamFileHeaderMerger` treats two records as identical.
trait IdRecord: Clone {
    fn id(&self) -> &str;
    fn with_id(&self, id: &str) -> Self;
    /// True when the record carries data the deferred `PP`-chain slice would have to handle.
    fn has_program_chain(&self) -> bool {
        false
    }
}

impl IdRecord for ReadGroup {
    fn id(&self) -> &str {
        &self.id
    }
    fn with_id(&self, id: &str) -> Self {
        let mut c = self.clone();
        c.id = id.to_string();
        c
    }
}

impl IdRecord for ProgramRecord {
    fn id(&self) -> &str {
        &self.id
    }
    fn with_id(&self, id: &str) -> Self {
        let mut c = self.clone();
        c.id = id.to_string();
        c
    }
    fn has_program_chain(&self) -> bool {
        self.attributes.get("PP").is_some()
    }
}

/// The result of merging one record kind: the merged records, and per input index the old-ID ->
/// new-ID map that `MergingSamRecordIterator` applies to that input's records.
type MergedRecords<T> = (Vec<T>, Vec<HashMap<String, String>>);

/// `SamFileHeaderMerger.mergeHeaderRecords` for one record kind (`@RG` or `@PG`), for the no-`PP` case:
/// bin the per-input records by ID in file order, then by identical content; keep the first content of
/// each ID under its own ID and rename every later content `ID.<base36(counter)>`; sort the result by
/// ID (`RECORD_ID_COMPARATOR`). Returns the merged records and, per input index, the old-ID -> new-ID
/// map that `MergingSamRecordIterator` applies to that input's records (identity entries included, as
/// htsjdk records them).
fn merge_header_records<T: IdRecord + PartialEq>(
    per_input: &[Vec<T>],
    duplicate_err: impl Fn(String) -> MergeError,
) -> Result<MergedRecords<T>, MergeError> {
    // idToRecord: insertion-ordered map ID -> [(content, [file indices])]. Both levels are
    // LinkedHashMap in htsjdk, so file order is deterministic and reproduced by a Vec.
    let mut id_order: Vec<String> = Vec::new();
    let mut groups: HashMap<String, Vec<(T, Vec<usize>)>> = HashMap::new();

    for (fi, records) in per_input.iter().enumerate() {
        let mut seen_in_file: HashSet<&str> = HashSet::new();
        for record in records {
            let id = record.id();
            if !seen_in_file.insert(id) {
                return Err(duplicate_err(id.to_string()));
            }
            if !groups.contains_key(id) {
                id_order.push(id.to_string());
                groups.insert(id.to_string(), Vec::new());
            }
            let bin = groups.get_mut(id).unwrap();
            match bin.iter_mut().find(|(content, _)| content == record) {
                Some((_, files)) => files.push(fi),
                None => bin.push((record.clone(), vec![fi])),
            }
        }
    }

    // Resolve collisions by remapping the second and later contents of each ID.
    let mut ids_taken: HashSet<String> = HashSet::new();
    let mut counter: u32 = 0;
    let mut translations: Vec<HashMap<String, String>> = vec![HashMap::new(); per_input.len()];
    let mut result: Vec<T> = Vec::new();

    for id in &id_order {
        for (content, files) in &groups[id] {
            let new_id = if !ids_taken.contains(id) {
                // Don't remap the first record with this ID.
                ids_taken.insert(id.clone());
                counter += 1;
                id.clone()
            } else {
                let mut candidate;
                loop {
                    candidate = format!("{}.{}", id, base36(counter));
                    counter += 1;
                    if !ids_taken.contains(&candidate) {
                        break;
                    }
                }
                ids_taken.insert(candidate.clone());
                candidate
            };
            for &fi in files {
                translations[fi].insert(id.clone(), new_id.clone());
            }
            result.push(content.with_id(&new_id));
        }
    }

    // RECORD_ID_COMPARATOR: String.compareTo, i.e. lexicographic by code unit; for the ASCII IDs here
    // this is Rust's str ordering.
    result.sort_by(|a, b| a.id().cmp(b.id()));
    Ok((result, translations))
}

/// Apply a read-group/program ID translation to one record's `RG:Z` or `PG:Z` tag, as
/// `MergingSamRecordIterator` does while reading.
fn remap_tag(record: &mut BamRecord, tag: &[u8; 2], translation: &HashMap<String, String>) {
    let new_value = match record.tags.get(Tag::new(tag)) {
        Some(TagValue::Str(value)) => translation.get(value).filter(|n| *n != value).cloned(),
        _ => None,
    };
    if let Some(value) = new_value {
        record.tags.insert(Tag::new(tag), TagValue::Str(value));
    }
}

/// `SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH`.
const UNKNOWN_SEQUENCE_LENGTH: i32 = 0;

/// `SAMSequenceRecord.isSameSequence`, for two records compared **at the same dictionary index** (as
/// `assertSequenceListsEqual` does): equal length (unless either is unknown), equal `M5` when both
/// carry one (compared as a hex number, so case and leading zeros do not matter), else equal name.
///
/// The alternative-name (`AN`) cross-match that htsjdk falls back to is not reproduced; a dictionary
/// that relies on `AN` to be considered equal is out of this slice's scope.
fn is_same_sequence(a: &SequenceRecord, b: &SequenceRecord) -> bool {
    if a.length != UNKNOWN_SEQUENCE_LENGTH
        && b.length != UNKNOWN_SEQUENCE_LENGTH
        && a.length != b.length
    {
        return false;
    }
    match (a.attributes.get("M5"), b.attributes.get("M5")) {
        (Some(m5a), Some(m5b)) => m5a
            .trim_start_matches('0')
            .eq_ignore_ascii_case(m5b.trim_start_matches('0')),
        _ => a.name == b.name,
    }
}

/// `SequenceUtil.assertSequenceListsEqual` as a predicate: same length, and every record the same
/// sequence at its index.
fn dictionaries_equal(a: &[SequenceRecord], b: &[SequenceRecord]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| is_same_sequence(x, y))
}

/// `SamFileHeaderMerger.mergeSequences`: fold `from` into `into`, finding the sequence-name superset.
/// New names are held and inserted just before the next shared name (or appended); a shared name that
/// appears out of order relative to an earlier shared one is unmergeable.
fn merge_sequences(
    into: &[SequenceRecord],
    from: &[SequenceRecord],
) -> Result<Vec<SequenceRecord>, MergeError> {
    let mut result: Vec<SequenceRecord> = into.to_vec();
    let mut holder: Vec<SequenceRecord> = Vec::new();
    let mut prevloc: i64 = -1;
    for seq in from {
        match result.iter().position(|s| s.name == seq.name) {
            None => holder.push(seq.clone()),
            Some(loc) => {
                if prevloc > loc as i64 {
                    return Err(MergeError::SequenceDictionaryOrderConflict);
                }
                let held = std::mem::take(&mut holder);
                let inserted = held.len();
                for (k, h) in held.into_iter().enumerate() {
                    result.insert(loc + k, h);
                }
                prevloc = (loc + inserted) as i64;
            }
        }
    }
    result.extend(holder);
    Ok(result)
}

/// `mergeSequenceDictionaries`: fold every input's dictionary together in file order, starting empty.
fn merge_sequence_dictionaries(headers: &[SamHeader]) -> Result<Vec<SequenceRecord>, MergeError> {
    let mut merged: Vec<SequenceRecord> = Vec::new();
    for header in headers {
        merged = merge_sequences(&merged, &header.sequences)?;
    }
    Ok(merged)
}

/// `createSequenceMapping` for one input: old reference index -> merged reference index, by name.
fn sequence_index_mapping(input: &[SequenceRecord], merged: &[SequenceRecord]) -> Vec<i32> {
    input
        .iter()
        .map(|s| {
            merged
                .iter()
                .position(|m| m.name == s.name)
                .map(|i| i as i32)
                .expect("every input sequence is present in the merged superset")
        })
        .collect()
}

/// Apply an index mapping to one reference index, leaving the unaligned sentinel (`-1`) untouched.
fn remap_reference_index(index: &mut i32, mapping: &[i32]) {
    if *index >= 0 {
        *index = mapping[*index as usize];
    }
}

/// The merge itself: the merged header and the merged, sorted records. Shared by the SAM and BAM
/// renderers so the header construction and record ordering cannot drift. `merge_dictionaries` is
/// `MERGE_SEQUENCE_DICTIONARIES` (default false): when the inputs' dictionaries are not equal it
/// unions them, otherwise a mismatch is an error.
fn merge_records(
    inputs: &[&str],
    order: SortOrder,
    merge_dictionaries: bool,
) -> Result<(SamHeader, Vec<BamRecord>), MergeError> {
    let mut headers: Vec<SamHeader> = Vec::with_capacity(inputs.len());
    let mut per_input_records: Vec<Vec<BamRecord>> = Vec::with_capacity(inputs.len());
    let mut comments: Vec<String> = Vec::new();
    for input in inputs {
        let (header, records) = read_sam(input)?;
        comments.extend(header.comments.iter().cloned());
        headers.push(header);
        per_input_records.push(records);
    }

    // Resolve the output dictionary. `getSequenceDictionary` requires every input to be equal (by
    // `SequenceUtil.assertSequenceListsEqual`); if they differ, `MERGE_SEQUENCE_DICTIONARIES` either
    // unions them (`mergeSequenceDictionaries`) or the tool errors.
    let first: &[SequenceRecord] = headers
        .first()
        .map(|h| h.sequences.as_slice())
        .unwrap_or(&[]);
    let dictionaries_equal = headers
        .iter()
        .all(|h| dictionaries_equal(first, &h.sequences));
    let sequences: Vec<SequenceRecord> = if dictionaries_equal {
        first.to_vec()
    } else if merge_dictionaries {
        merge_sequence_dictionaries(&headers)?
    } else {
        return Err(MergeError::SequenceDictionaryMismatch);
    };

    // When the dictionary was unioned, each input's records point into its own (smaller/differently
    // ordered) dictionary, so `MergingSamRecordIterator`'s `record.setHeader(merged)` re-resolves the
    // reference indices by name. Reproduce that: remap `reference_index` and `mate_reference_index`
    // from each input's dictionary to the merged one. Equal dictionaries need no remap.
    if !dictionaries_equal {
        for (fi, records) in per_input_records.iter_mut().enumerate() {
            let mapping = sequence_index_mapping(&headers[fi].sequences, &sequences);
            for record in records.iter_mut() {
                remap_reference_index(&mut record.reference_index, &mapping);
                remap_reference_index(&mut record.mate_reference_index, &mapping);
            }
        }
    }

    // Merge @RG and @PG, collecting the per-input ID translations.
    let rg_per_input: Vec<Vec<ReadGroup>> = headers.iter().map(|h| h.read_groups.clone()).collect();
    let (read_groups, rg_translation) =
        merge_header_records(&rg_per_input, MergeError::DuplicateReadGroupId)?;

    if headers
        .iter()
        .any(|h| h.programs.iter().any(ProgramRecord::has_program_chain))
    {
        return Err(MergeError::ProgramChainNotSupported);
    }
    let pg_per_input: Vec<Vec<ProgramRecord>> =
        headers.iter().map(|h| h.programs.clone()).collect();
    let (programs, pg_translation) =
        merge_header_records(&pg_per_input, MergeError::DuplicateProgramId)?;

    // Rewrite each input's records' RG/PG tags through that input's translation, then concatenate.
    for (fi, records) in per_input_records.iter_mut().enumerate() {
        for record in records.iter_mut() {
            remap_tag(record, b"RG", &rg_translation[fi]);
            remap_tag(record, b"PG", &pg_translation[fi]);
        }
    }
    let mut records: Vec<BamRecord> = per_input_records.into_iter().flatten().collect();

    let mut header = SamHeader::new(); // @HD VN:<current>
    header.set_group_order("none"); // VN, GO
    header.set_sort_order(order.name()); // VN, GO, SO
    header.sequences = sequences;
    header.read_groups = read_groups;
    header.programs = programs;
    header.comments = comments;

    // A stable sort so wholly-identical records keep input order (decision 0021).
    records.sort_by(order.comparator());
    Ok((header, records))
}

/// `MergeSamFiles.doWork`: the merged, sorted SAM. `merge_dictionaries` is
/// `MERGE_SEQUENCE_DICTIONARIES` (default false).
pub fn merge_sam_files(
    inputs: &[&str],
    order: SortOrder,
    merge_dictionaries: bool,
) -> Result<String, MergeError> {
    let (header, records) = merge_records(inputs, order, merge_dictionaries)?;
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// The same merge with **BAM** output, written through htsjdk-rs's byte-identical `BamWriter`.
/// Byte-identical to Picard with `USE_JDK_DEFLATER=true` (the merged records are those the SAM path
/// already reproduces).
pub fn merge_sam_files_to_bam(
    inputs: &[&str],
    order: SortOrder,
    merge_dictionaries: bool,
) -> Result<Vec<u8>, MergeError> {
    use htsjdk_bam::writer::BamWriter;

    let (header, records) = merge_records(inputs, order, merge_dictionaries)?;
    let mut w = BamWriter::new(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        w.write(rec).expect("record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: &str =
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s\tLB:lib1\n";

    fn rec(name: &str, start: i32) -> String {
        format!("{name}\t0\tchr1\t{start}\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n")
    }

    #[test]
    fn base36_matches_positive_four_digit_base36_str() {
        assert_eq!(base36(0), "0");
        assert_eq!(base36(1), "1");
        assert_eq!(base36(2), "2");
        assert_eq!(base36(35), "z");
        assert_eq!(base36(36), "10");
    }

    #[test]
    fn merges_and_sorts_under_a_group_ordered_header() {
        let a = format!("{H}{}{}", rec("a", 10), rec("c", 30));
        let b = format!("{H}{}{}", rec("b", 20), rec("d", 40));
        let out = merge_sam_files(&[&a, &b], SortOrder::Coordinate, false).unwrap();
        assert_eq!(
            out,
            "@HD\tVN:1.6\tGO:none\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:1000\n\
             @RG\tID:rg1\tSM:s\tLB:lib1\n\
             a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
             b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
             c\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n\
             d\t0\tchr1\t40\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
        );
    }

    #[test]
    fn identical_read_groups_are_not_duplicated() {
        let a = format!("{H}{}", rec("a", 10));
        let b = format!("{H}{}", rec("b", 20));
        let out = merge_sam_files(&[&a, &b], SortOrder::Coordinate, false).unwrap();
        assert_eq!(out.matches("@RG\t").count(), 1);
    }

    #[test]
    fn bam_output_decodes_to_the_same_as_the_sam_output() {
        use htsjdk_bam::reader::BamReader;
        let a = format!("{H}{}{}", rec("a", 10), rec("c", 30));
        let b = format!("{H}{}{}", rec("b", 20), rec("d", 40));
        let sam = merge_sam_files(&[&a, &b], SortOrder::Coordinate, false).unwrap();
        let bam = merge_sam_files_to_bam(&[&a, &b], SortOrder::Coordinate, false).unwrap();
        let decoded = htsjdk_bgzf::decompress_all(&bam).unwrap();
        let reader = BamReader::new(&decoded).unwrap();
        let header = reader.header.text.clone();
        let recs: Vec<_> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(write_sam(&header, &recs).unwrap(), sam);
    }

    #[test]
    fn distinct_read_groups_are_unioned() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s1\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg2\tSM:s2\n\
            b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n";
        let out = merge_sam_files(&[a, b], SortOrder::Coordinate, false).unwrap();
        assert!(out.contains("@RG\tID:rg1\tSM:s1"));
        assert!(out.contains("@RG\tID:rg2\tSM:s2"));
    }

    #[test]
    fn a_colliding_read_group_id_is_renamed_and_records_are_remapped() {
        // Both inputs declare ID:rg1 with different content -> rg1 (file 0) and rg1.1 (file 1). The
        // second file's records' RG:Z is rewritten to rg1.1.
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s1\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s2\n\
            b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        let out = merge_sam_files(&[a, b], SortOrder::Coordinate, false).unwrap();
        assert!(out.contains("@RG\tID:rg1\tSM:s1"), "{out}");
        assert!(out.contains("@RG\tID:rg1.1\tSM:s2"), "{out}");
        // rg1 sorts before rg1.1, so the @RG order is rg1 then rg1.1.
        assert!(
            out.find("ID:rg1\t").unwrap() < out.find("ID:rg1.1").unwrap(),
            "{out}"
        );
        // Record a keeps RG:Z:rg1, record b becomes RG:Z:rg1.1.
        let a_line = out.lines().find(|l| l.starts_with("a\t")).unwrap();
        let b_line = out.lines().find(|l| l.starts_with("b\t")).unwrap();
        assert!(a_line.ends_with("RG:Z:rg1"), "{a_line}");
        assert!(b_line.ends_with("RG:Z:rg1.1"), "{b_line}");
    }

    #[test]
    fn three_colliding_read_groups_get_counter_suffixes() {
        let mk = |sm: &str, name: &str| {
            format!(
                "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:{sm}\n\
                 {name}\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            )
        };
        let a = mk("s1", "a");
        let b = mk("s2", "b");
        let c = mk("s3", "c");
        let out = merge_sam_files(&[&a, &b, &c], SortOrder::Coordinate, false).unwrap();
        assert!(out.contains("@RG\tID:rg1\tSM:s1"), "{out}");
        assert!(out.contains("@RG\tID:rg1.1\tSM:s2"), "{out}");
        assert!(out.contains("@RG\tID:rg1.2\tSM:s3"), "{out}");
        let c_line = out.lines().find(|l| l.starts_with("c\t")).unwrap();
        assert!(c_line.ends_with("RG:Z:rg1.2"), "{c_line}");
    }

    #[test]
    fn a_duplicate_read_group_within_one_input_is_an_error() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            @RG\tID:rg1\tSM:s1\n@RG\tID:rg1\tSM:s2\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
        assert!(matches!(
            merge_sam_files(&[a], SortOrder::Coordinate, false),
            Err(MergeError::DuplicateReadGroupId(_))
        ));
    }

    #[test]
    fn a_mismatched_sequence_dictionary_is_an_error() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:2000\n\
            b\t0\tchr2\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        assert!(matches!(
            merge_sam_files(&[a, b], SortOrder::Coordinate, false),
            Err(MergeError::SequenceDictionaryMismatch)
        ));
    }

    #[test]
    fn a_program_chain_is_not_yet_supported() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            @PG\tID:p1\tPN:a\n@PG\tID:p2\tPN:b\tPP:p1\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        assert!(matches!(
            merge_sam_files(&[a], SortOrder::Coordinate, false),
            Err(MergeError::ProgramChainNotSupported)
        ));
    }

    #[test]
    fn dictionaries_are_unioned_into_the_superset_and_records_remapped() {
        // File 0 knows chr1, chr2; file 1 knows chr2, chr3. Superset is chr1, chr2, chr3. File 1's
        // records point into its own dictionary (chr2=0, chr3=1) and must be remapped (chr2=1,
        // chr3=2) so their RNAME renders correctly and they sort into the merged order.
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
            b\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let c = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chr3\tLN:3000\n\
            c\t0\tchr2\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
            d\t0\tchr3\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let out = merge_sam_files(&[a, c], SortOrder::Coordinate, true).unwrap();
        assert_eq!(
            out,
            "@HD\tVN:1.6\tGO:none\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:1000\n\
             @SQ\tSN:chr2\tLN:2000\n\
             @SQ\tSN:chr3\tLN:3000\n\
             a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
             b\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
             c\t0\tchr2\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
             d\t0\tchr3\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        );
    }

    #[test]
    fn an_incompatible_dictionary_order_is_an_error() {
        // chr1 then chr2 vs chr2 then chr1: the shared sequences are ordered oppositely, so no
        // consistent superset exists.
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chr1\tLN:1000\n\
            b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        assert!(matches!(
            merge_sam_files(&[a, b], SortOrder::Coordinate, true),
            Err(MergeError::SequenceDictionaryOrderConflict)
        ));
    }

    #[test]
    fn equal_dictionaries_need_no_flag_and_no_remap() {
        // Same dictionary in both: the union path is not taken even with the flag off.
        let a = format!("{H}{}", rec("a", 10));
        let b = format!("{H}{}", rec("b", 20));
        let off = merge_sam_files(&[&a, &b], SortOrder::Coordinate, false).unwrap();
        let on = merge_sam_files(&[&a, &b], SortOrder::Coordinate, true).unwrap();
        assert_eq!(off, on);
    }
}
