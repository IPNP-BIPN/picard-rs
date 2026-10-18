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
//! `@PG` records additionally honour the `PP` (previous-program) chain (`mergeProgramGroups`): tree
//! roots are merged first, then each level of children whose `PP` points at a just-merged record from
//! the same input, with `PP` rewritten through the accumulating rename table (so a child chaining to a
//! renamed parent collides and is itself renamed) and the base36 counter shared across every level. A
//! `PP` pointing at a program group that never appears is an error.
//!
//! The sequence dictionaries either match (`SequenceUtil.assertSequenceListsEqual`) or, with
//! `MERGE_SEQUENCE_DICTIONARIES`, are unioned into the name superset (`mergeSequences`); when unioned,
//! each input's records have their `reference_index`/`mate_reference_index` remapped by name, as
//! `MergingSamRecordIterator.setHeader` re-resolves them. A mismatch without the flag, or a shared
//! sequence ordered incompatibly across inputs, is an error.
//!
//! This completes the `SamFileHeaderMerger` surface `MergeSamFiles` exercises.

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
    /// A `@PG` record's `PP` (previous-program) attribute points at a program group that does not
    /// exist, so `mergeProgramGroups` cannot place it in the tree (`program groups weren't
    /// processed`).
    ProgramChainBroken,
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
pub fn base36(mut n: u32) -> String {
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

/// `SAMFileHeader.PgIdGenerator.getNonCollidingId` over a fresh generator: the id itself when no
/// `@PG` of the header has it, else the id with `.` and the first base-36 counter, from the number
/// of ids already taken, that is free.
pub fn non_colliding_pg_id(header: &htsjdk_bam::header::SamHeader, id: &str) -> String {
    let taken: std::collections::BTreeSet<&str> =
        header.programs.iter().map(|pg| pg.id.as_str()).collect();
    if !taken.contains(id) {
        return id.to_string();
    }
    let mut counter = taken.len() as u32;
    loop {
        let candidate = format!("{id}.{}", base36(counter));
        counter += 1;
        if !taken.contains(candidate.as_str()) {
            return candidate;
        }
    }
}

/// One merged header record (`@RG` or `@PG`) reduced to the two things the merge needs: its ID and a
/// content key that is equal exactly when `SamFileHeaderMerger` treats two records as identical.
trait IdRecord: Clone {
    fn id(&self) -> &str;
    fn with_id(&self, id: &str) -> Self;
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
}

/// The result of merging one record kind: the merged records, and per input index the old-ID ->
/// new-ID map that `MergingSamRecordIterator` applies to that input's records.
type MergedRecords<T> = (Vec<T>, Vec<HashMap<String, String>>);

/// State threaded across `mergeHeaderRecords` calls: the merged records so far, the ID set (`@PG`
/// resolves its tree one level per call, all sharing this state), the base36 counter, and the
/// per-input old-ID -> new-ID translation.
struct MergeState<T> {
    ids_taken: HashSet<String>,
    counter: u32,
    translations: Vec<HashMap<String, String>>,
    result: Vec<T>,
}

impl<T: IdRecord + PartialEq> MergeState<T> {
    fn new(n_inputs: usize) -> Self {
        MergeState {
            ids_taken: HashSet::new(),
            counter: 0,
            translations: vec![HashMap::new(); n_inputs],
            result: Vec::new(),
        }
    }

    /// `SamFileHeaderMerger.mergeHeaderRecords` for one batch of `(file index, record)`: bin by ID in
    /// batch order, then by identical content; keep the first content of each ID under its own ID and
    /// rename every later content `ID.<base36(counter)>`; append to `result`, recording each mapping.
    fn resolve_batch(&mut self, batch: &[(usize, T)]) {
        // idToRecord: insertion-ordered map ID -> [(content, [file indices])]. Both levels are
        // LinkedHashMap in htsjdk, so batch order is deterministic and reproduced by a Vec.
        let mut id_order: Vec<String> = Vec::new();
        let mut groups: HashMap<String, Vec<(T, Vec<usize>)>> = HashMap::new();
        for (fi, record) in batch {
            let id = record.id();
            if !groups.contains_key(id) {
                id_order.push(id.to_string());
                groups.insert(id.to_string(), Vec::new());
            }
            let bin = groups.get_mut(id).unwrap();
            match bin.iter_mut().find(|(content, _)| content == record) {
                Some((_, files)) => files.push(*fi),
                None => bin.push((record.clone(), vec![*fi])),
            }
        }

        for id in &id_order {
            for (content, files) in &groups[id] {
                let new_id = if !self.ids_taken.contains(id) {
                    // Don't remap the first record with this ID.
                    self.ids_taken.insert(id.clone());
                    self.counter += 1;
                    id.clone()
                } else {
                    let mut candidate;
                    loop {
                        candidate = format!("{}.{}", id, base36(self.counter));
                        self.counter += 1;
                        if !self.ids_taken.contains(&candidate) {
                            break;
                        }
                    }
                    self.ids_taken.insert(candidate.clone());
                    candidate
                };
                for &fi in files {
                    self.translations[fi].insert(id.clone(), new_id.clone());
                }
                self.result.push(content.with_id(&new_id));
            }
        }
    }
}

/// No two records in one input may share an ID (`mergeReadGroups`/`mergeProgramGroups` throw
/// `contains more than one RG/PG with the same id`).
fn assert_unique_ids<T: IdRecord>(
    per_input: &[Vec<T>],
    duplicate_err: impl Fn(String) -> MergeError,
) -> Result<(), MergeError> {
    for records in per_input {
        let mut seen: HashSet<&str> = HashSet::new();
        for r in records {
            if !seen.insert(r.id()) {
                return Err(duplicate_err(r.id().to_string()));
            }
        }
    }
    Ok(())
}

/// The whole batch flattened to `(file index, record)` in file order.
fn flatten<T: Clone>(per_input: &[Vec<T>]) -> Vec<(usize, T)> {
    per_input
        .iter()
        .enumerate()
        .flat_map(|(fi, rs)| rs.iter().cloned().map(move |r| (fi, r)))
        .collect()
}

/// `mergeReadGroups`: one `mergeHeaderRecords` pass over every `@RG`, then sort by ID.
fn merge_read_groups(per_input: &[Vec<ReadGroup>]) -> Result<MergedRecords<ReadGroup>, MergeError> {
    assert_unique_ids(per_input, MergeError::DuplicateReadGroupId)?;
    let mut state = MergeState::new(per_input.len());
    state.resolve_batch(&flatten(per_input));
    state.result.sort_by(|a, b| a.id.cmp(&b.id));
    Ok((state.result, state.translations))
}

/// `translateIds` for one `@PG` record: rename its ID (always) and, when `translate_pp`, its `PP`
/// pointer, through this record's input translation table. A record that needs no change is returned
/// unchanged.
fn translate_program(
    record: &ProgramRecord,
    translation: &HashMap<String, String>,
    translate_pp: bool,
) -> ProgramRecord {
    let mut out = record.clone();
    if let Some(new_id) = translation.get(&record.id) {
        if new_id != &record.id {
            out.id = new_id.clone();
        }
    }
    if translate_pp {
        if let Some(pp) = record.attributes.get("PP") {
            if let Some(new_pp) = translation.get(pp) {
                if new_pp != pp {
                    out.attributes.set("PP", new_pp);
                }
            }
        }
    }
    out
}

/// One level of `@PG` records to process, each tagged with its input's index.
type Batch = Vec<(usize, ProgramRecord)>;

/// `mergeProgramGroups`: merge `@PG` records honouring the `PP` (previous-program) chain. Process the
/// tree roots (no `PP`) first, then each level of children whose `PP` points at a just-processed
/// record from the same input, rewriting `PP` through the accumulating rename table so a child that
/// now chains to a renamed parent collides (and is itself renamed) exactly as htsjdk does. The base36
/// counter is shared across every level. Finally sort by ID.
fn merge_programs(
    per_input: &[Vec<ProgramRecord>],
) -> Result<MergedRecords<ProgramRecord>, MergeError> {
    assert_unique_ids(per_input, MergeError::DuplicateProgramId)?;
    let mut state = MergeState::new(per_input.len());

    // Split off the tree roots (no PP) as the first level; the rest wait for their parent.
    let (mut current, mut left): (Batch, Batch) = flatten(per_input)
        .into_iter()
        .partition(|(_, r)| r.attributes.get("PP").is_none());

    while !current.is_empty() {
        state.resolve_batch(&current);

        // Reflect this level's renames: current's own IDs, and every leftover's PP pointer.
        current = current
            .iter()
            .map(|(fi, r)| (*fi, translate_program(r, &state.translations[*fi], false)))
            .collect();
        left = left
            .iter()
            .map(|(fi, r)| (*fi, translate_program(r, &state.translations[*fi], true)))
            .collect();

        // The next level is every leftover whose PP now names a just-processed record in its file.
        let mut next = Vec::new();
        let mut still_left = Vec::new();
        for (fi, r) in left {
            let matches = r
                .attributes
                .get("PP")
                .is_some_and(|pp| current.iter().any(|(cfi, cr)| *cfi == fi && cr.id == pp));
            if matches {
                next.push((fi, r));
            } else {
                still_left.push((fi, r));
            }
        }
        left = still_left;
        current = next;
    }

    if !left.is_empty() {
        return Err(MergeError::ProgramChainBroken);
    }
    state.result.sort_by(|a, b| a.id.cmp(&b.id));
    Ok((state.result, state.translations))
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
    let mut parsed: Vec<(SamHeader, Vec<BamRecord>)> = Vec::with_capacity(inputs.len());
    for input in inputs {
        parsed.push(read_sam(input)?);
    }
    merge_parsed(parsed, order.name(), merge_dictionaries, Some(order))
}

/// The merge over inputs that are already parsed, which is what the runnable binary has: its
/// inputs are BAMs as often as they are SAM text.
///
/// Two things `doWork` decides and this cannot are passed in rather than derived. `sort_order_name`
/// is the `SO` the merged header carries, which is `SORT_ORDER` verbatim and may name an order
/// this port does not sort by (`unsorted`, `unknown`). `sort_with` is whether the *writer* would
/// have sorted: `makeWriter` is told the stream is presorted when every input's `SO` already equals
/// `SORT_ORDER`, or when `ASSUME_SORTED` says to believe it does, and a presorted writer keeps the
/// merge order rather than imposing one.
pub fn merge_parsed(
    parsed: Vec<(SamHeader, Vec<BamRecord>)>,
    sort_order_name: &str,
    merge_dictionaries: bool,
    sort_with: Option<SortOrder>,
) -> Result<(SamHeader, Vec<BamRecord>), MergeError> {
    let mut headers: Vec<SamHeader> = Vec::with_capacity(parsed.len());
    let mut per_input_records: Vec<Vec<BamRecord>> = Vec::with_capacity(parsed.len());
    let mut comments: Vec<String> = Vec::new();
    for (header, records) in parsed {
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
    let (read_groups, rg_translation) = merge_read_groups(&rg_per_input)?;

    let pg_per_input: Vec<Vec<ProgramRecord>> =
        headers.iter().map(|h| h.programs.clone()).collect();
    let (programs, pg_translation) = merge_programs(&pg_per_input)?;

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
    header.set_sort_order(sort_order_name); // VN, GO, SO
    header.sequences = sequences;
    header.read_groups = read_groups;
    header.programs = programs;
    header.comments = comments;

    // A stable sort so wholly-identical records keep input order (decision 0021). Skipped
    // entirely when the writer was presorted, because then nothing reordered the merge.
    if let Some(order) = sort_with {
        records.sort_by(order.comparator());
    }
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
    let (header, records) = merge_records(inputs, order, merge_dictionaries)?;
    let mut w =
        crate::output::bam_writer(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        w.write(rec).expect("record re-encodes as BAM");
    }
    Ok(w.finish().expect("finish never fails on a Vec"))
}

/// `MergingSamRecordIterator` over inputs already read: a `java.util.PriorityQueue` of the inputs
/// keyed by each one's next record, polled one record at a time and offered back while it has
/// more. The queue's sift order is reproduced rather than a stable sort, because two heads the
/// comparator calls equal leave in the heap's order, not in file order.
pub fn priority_merge(
    per_input: Vec<Vec<BamRecord>>,
    compare: fn(&BamRecord, &BamRecord) -> std::cmp::Ordering,
) -> Vec<BamRecord> {
    use std::cmp::Ordering;
    let total = per_input.iter().map(Vec::len).sum();
    let mut rest: Vec<std::vec::IntoIter<BamRecord>> =
        per_input.into_iter().map(Vec::into_iter).collect();
    let mut heads: Vec<Option<BamRecord>> = rest.iter_mut().map(Iterator::next).collect();
    let order = |heads: &[Option<BamRecord>], a: usize, b: usize| -> Ordering {
        let head = |i: usize| heads[i].as_ref().expect("a queued input has a record");
        compare(head(a), head(b))
    };
    // `offer`: siftUp.
    let offer = |queue: &mut Vec<usize>, heads: &[Option<BamRecord>], x: usize| {
        let mut k = queue.len();
        queue.push(x);
        while k > 0 {
            let parent = (k - 1) >> 1;
            let e = queue[parent];
            if order(heads, x, e) != Ordering::Less {
                break;
            }
            queue[k] = e;
            k = parent;
        }
        queue[k] = x;
    };
    let mut queue: Vec<usize> = Vec::with_capacity(heads.len());
    for index in 0..heads.len() {
        if heads[index].is_some() {
            offer(&mut queue, &heads, index);
        }
    }
    let mut out = Vec::with_capacity(total);
    while !queue.is_empty() {
        // `poll`: the head leaves, and the last element sifts down from the root.
        let head = queue[0];
        let last = queue.pop().expect("the queue is not empty");
        let n = queue.len();
        if n > 0 {
            let mut k = 0;
            while k < n >> 1 {
                let mut child = 2 * k + 1;
                let right = child + 1;
                if right < n && order(&heads, queue[child], queue[right]) == Ordering::Greater {
                    child = right;
                }
                if order(&heads, last, queue[child]) != Ordering::Greater {
                    break;
                }
                queue[k] = queue[child];
                k = child;
            }
            queue[k] = last;
        }
        let next = rest[head].next();
        out.push(std::mem::replace(&mut heads[head], next).expect("a queued input has a record"));
        if heads[head].is_some() {
            offer(&mut queue, &heads, head);
        }
    }
    out
}

/// What `MergingSamRecordIterator` hands out for each input after `SamFileHeaderMerger` has merged
/// their headers without merging dictionaries: the merged header, and every input's records with
/// their `RG` and `PG` tags renamed where the merge renamed the group.
pub fn merge_headers_for_iteration(
    parsed: Vec<(SamHeader, Vec<BamRecord>)>,
    sort_order_name: &str,
) -> Result<(SamHeader, Vec<Vec<BamRecord>>), MergeError> {
    let lengths: Vec<usize> = parsed.iter().map(|(_, records)| records.len()).collect();
    let (header, records) = merge_parsed(parsed, sort_order_name, false, None)?;
    let mut records = records.into_iter();
    let per_input = lengths
        .iter()
        .map(|&n| records.by_ref().take(n).collect())
        .collect();
    Ok((header, per_input))
}

/// The merged header, and the records or what reading the first of them throws.
pub type OpenedInputs = (SamHeader, Result<Vec<BamRecord>, String>);

/// `AbstractMarkDuplicatesCommandLineProgram.openInputs`: the first header takes the assumed sort
/// order when there is one, a single input is read as it is, and several go through a
/// `SamFileHeaderMerger` in the first header's order (dictionaries not merged) and a
/// `MergingSamRecordIterator` that checks every input's order unless one was assumed. Each `Err`
/// is what the JVM prints after `Exception in thread "main" `: the outer one is thrown while the
/// inputs are opened, the inner one only once the first record is read.
pub fn open_inputs(
    mut inputs: Vec<(SamHeader, Vec<BamRecord>)>,
    assumed: Option<&str>,
) -> Result<OpenedInputs, String> {
    if let Some(order) = assumed {
        inputs[0].0.set_sort_order(order);
    }
    if inputs.len() == 1 {
        let (header, records) = inputs.pop().expect("one input");
        return Ok((header, Ok(records)));
    }
    let order_of = |h: &SamHeader| h.attributes.get("SO").unwrap_or("unsorted").to_string();
    let sort_order = order_of(&inputs[0].0);
    let compatible = assumed.is_some()
        || sort_order == "unsorted"
        || inputs.iter().all(|(h, _)| order_of(h) == sort_order);
    let headers: Vec<SamHeader> = inputs.iter().map(|(h, _)| h.clone()).collect();
    let duplicate = |kind: &str, id: &str| -> String {
        let header = headers
            .iter()
            .find(|h| {
                let ids: Vec<&str> = if kind == "RG" {
                    h.read_groups.iter().map(|g| g.id.as_str()).collect()
                } else {
                    h.programs.iter().map(|p| p.id.as_str()).collect()
                };
                ids.iter().filter(|&&i| i == id).count() > 1
            })
            .map(|h| {
                let attributes: Vec<String> = h
                    .attributes
                    .iter()
                    .map(|(k, v)| format!("{k}={v}"))
                    .collect();
                format!("SAMFileHeader{{{}}}", attributes.join(", "))
            })
            .unwrap_or_default();
        format!(
            "htsjdk.samtools.SAMException: Input file: {header} contains more than one {kind} \
             with the same id ({id})"
        )
    };
    let (header, per_input) = match merge_headers_for_iteration(inputs, &sort_order) {
        Ok(merged) => merged,
        Err(MergeError::DuplicateReadGroupId(id)) => return Err(duplicate("RG", &id)),
        Err(MergeError::DuplicateProgramId(id)) => return Err(duplicate("PG", &id)),
        Err(MergeError::ProgramChainBroken) => {
            return Err(
                "htsjdk.samtools.SAMException: program groups weren't processed. Do \
                        their PP ids point to existing PGs? "
                    .to_string(),
            )
        }
        Err(other) => return Err(format!("htsjdk.samtools.SAMException: {other:?}")),
    };
    if !compatible {
        return Err(
            "htsjdk.samtools.SAMException: Files are not compatible with sort order".into(),
        );
    }
    // The queue compares the first two non-empty inputs' heads as soon as it is read, the
    // second's first.
    let heads: Vec<&BamRecord> = per_input.iter().filter_map(|r| r.first()).take(2).collect();
    let records = match sort_order.as_str() {
        "coordinate" => Ok(priority_merge(per_input, htsjdk_bam::coordinate::compare)),
        "queryname" => Ok(priority_merge(per_input, htsjdk_bam::query_name::compare)),
        // `SortOrder.unknown` has no comparator, and `ComparableSamRecordIterator.compareTo` asks it
        // for its class.
        "unknown" if heads.len() == 2 => Err("java.lang.NullPointerException: Cannot invoke \
                                              \"Object.getClass()\" because \"this.comparator\" \
                                              is null"
            .to_string()),
        // `SAMRecordDuplicateComparator` reads every mapped pair's mate coordinate from its MC.
        "duplicate" if heads.len() == 2 => {
            let lacks_mate_cigar = |r: &BamRecord| {
                r.flags & 0x1 != 0
                    && r.flags & 0x4 == 0
                    && r.flags & 0x8 == 0
                    && r.tags.get(Tag::new(b"MC")).is_none()
            };
            if let Some(r) = [heads[1], heads[0]]
                .into_iter()
                .find(|r| lacks_mate_cigar(r))
            {
                Err(format!(
                    "htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found: {}",
                    crate::view_sam::sam_string(&header, r).trim_end()
                ))
            } else if per_input.iter().flatten().any(lacks_mate_cigar) {
                // shortcut: a later record without MC throws where the duplicate order meets it,
                // which needs that comparator ported; every caller refuses this order afterwards.
                Err("merging inputs in duplicate order is not ported".to_string())
            } else {
                Ok(per_input.concat())
            }
        }
        // shortcut: any other order merges by identity hash in htsjdk, which no caller reads on
        // (every MarkDuplicates tool refuses that order next); port it if one ever accepts it.
        _ => Ok(per_input.concat()),
    };
    Ok((header, records))
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
    fn open_inputs_merges_by_coordinate_and_refuses_a_mixed_order() {
        let a = read_sam(&format!("{H}{}{}", rec("a", 10), rec("c", 30))).unwrap();
        let b = read_sam(&format!("{H}{}{}", rec("b", 20), rec("d", 40))).unwrap();
        let (header, records) = open_inputs(vec![a.clone(), b.clone()], None).unwrap();
        let records = records.unwrap();
        assert_eq!(header.attributes.get("GO"), Some("none"));
        let names: Vec<&str> = records.iter().map(|r| r.read_name.as_str()).collect();
        assert_eq!(names, ["a", "b", "c", "d"]);

        let queryname = read_sam(&format!(
            "{}{}",
            H.replace("SO:coordinate", "SO:queryname"),
            rec("e", 5)
        ))
        .unwrap();
        let refused = open_inputs(vec![a.clone(), queryname.clone()], None).unwrap_err();
        assert!(
            refused.ends_with("Files are not compatible with sort order"),
            "{refused}"
        );
        // An assumed order is not checked against the files.
        assert!(open_inputs(vec![a.clone(), queryname], Some("coordinate")).is_ok());
        // `unknown` has no comparator: opening succeeds, reading the first record throws.
        let (_, unread) = open_inputs(vec![a, b], Some("unknown")).unwrap();
        assert!(unread
            .unwrap_err()
            .starts_with("java.lang.NullPointerException"));
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
    fn an_identical_program_chain_is_deduped() {
        // Both files carry the same p1 (root) and p2 (PP:p1); the whole chain dedupes to one of each.
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            @PG\tID:p1\tPN:a\n@PG\tID:p2\tPN:b\tPP:p1\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            @PG\tID:p1\tPN:a\n@PG\tID:p2\tPN:b\tPP:p1\n\
            b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2\n";
        let out = merge_sam_files(&[a, b], SortOrder::Coordinate, false).unwrap();
        assert_eq!(out.matches("@PG\t").count(), 2, "{out}");
        assert!(out.contains("@PG\tID:p1\tPN:a\n"), "{out}");
        assert!(out.contains("@PG\tID:p2\tPN:b\tPP:p1\n"), "{out}");
    }

    #[test]
    fn a_colliding_program_chain_renames_the_child_and_rechains_its_pp() {
        // p1 differs (VN) between files -> p1, p1.1. Each file's p2 (PP:p1) rechains to its own
        // parent, so file 1's p2 becomes PP:p1.1; the two p2 now differ and collide. The shared
        // counter is at 2 after the roots, so the child rename uses base36(3) = "3": p2.3.
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            @PG\tID:p1\tPN:a\tVN:1\n@PG\tID:p2\tPN:b\tPP:p1\n\
            r0\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2\n";
        let b = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            @PG\tID:p1\tPN:a\tVN:2\n@PG\tID:p2\tPN:b\tPP:p1\n\
            r1\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2\n";
        let out = merge_sam_files(&[a, b], SortOrder::Coordinate, false).unwrap();
        assert_eq!(
            out,
            "@HD\tVN:1.6\tGO:none\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:1000\n\
             @PG\tID:p1\tPN:a\tVN:1\n\
             @PG\tID:p1.1\tPN:a\tVN:2\n\
             @PG\tID:p2\tPN:b\tPP:p1\n\
             @PG\tID:p2.3\tPN:b\tPP:p1.1\n\
             r0\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2\n\
             r1\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2.3\n"
        );
    }

    #[test]
    fn a_program_pointing_at_a_missing_parent_is_an_error() {
        let a = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            @PG\tID:p2\tPN:b\tPP:nope\n\
            a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        assert!(matches!(
            merge_sam_files(&[a], SortOrder::Coordinate, false),
            Err(MergeError::ProgramChainBroken)
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
