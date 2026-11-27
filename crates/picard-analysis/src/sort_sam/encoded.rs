//! `SortSam`'s `SortingCollection`, holding each record as `BamWriter` would encode it rather than
//! decoded.
//!
//! The comparators read the fixed fields and the name out of those bytes, and the bytes go to the
//! runs and to the output as they are, so a record is decoded once and encoded once, where the
//! decoded collection in the parent module encodes it again for every run it goes to and for the
//! output. The other tools that sort keep that collection, with its comparators on decoded records.

use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::binary_heap::PeekMut;
use std::collections::{BinaryHeap, VecDeque};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use htsjdk_bam::bin;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::index::{BamIndexer, Chunk};
use htsjdk_bam::record::{BamRecord, EncodeError};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::WriteError;
use htsjdk_bgzf::vfp::make_file_pointer;
use htsjdk_bgzf::{BgzfWriter, DEFAULT_UNCOMPRESSED_BLOCK_SIZE};
use rayon::prelude::*;

use super::{SortOrder, MAX_OPEN_RUNS};

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

/// A record as `BAMRecordCodec` encodes it, read where the comparators and the index look, with
/// the two things its fixed fields do not carry: the `HI` tag and the alignment end.
#[derive(Clone, Copy)]
pub struct Encoded<'a> {
    /// The record's bytes, its length prefix included.
    pub bytes: &'a [u8],
    hit: Option<i64>,
    end: i32,
}

impl Encoded<'_> {
    fn int(&self, at: usize) -> i32 {
        i32::from_le_bytes(self.bytes[at..at + 4].try_into().unwrap())
    }
    fn reference_index(&self) -> i32 {
        self.int(4)
    }
    fn alignment_start(&self) -> i32 {
        self.int(8) + 1
    }
    fn mapping_quality(&self) -> u8 {
        self.bytes[13]
    }
    fn flags(&self) -> u16 {
        u16::from_le_bytes([self.bytes[18], self.bytes[19]])
    }
    fn has(&self, bit: u16) -> bool {
        self.flags() & bit != 0
    }
    /// The name's bytes, which order as the decoded name does: a BAM name is read a byte to a
    /// character.
    // shortcut: a SAM name with a character above U+00FF encodes lossily and may sort differently
    // from the decoded name; read names are ASCII in practice.
    fn read_name(&self) -> &[u8] {
        &self.bytes[36..36 + (self.bytes[12] as usize).saturating_sub(1)]
    }
}

/// `SAMRecordCoordinateComparator`, as [`coordinate::compare`] ports it, on encoded records.
fn coordinate_compare(a: &Encoded, b: &Encoded) -> Ordering {
    // `fileOrderCompare`: no reference sorts after every reference, and two such are equal.
    let (r1, r2) = (a.reference_index(), b.reference_index());
    let cmp = match (r1 == -1, r2 == -1) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => r1
            .cmp(&r2)
            .then_with(|| a.alignment_start().cmp(&b.alignment_start())),
    };
    if cmp != Ordering::Equal {
        return cmp;
    }
    let field = |r: &Encoded, at: usize| r.int(at);
    a.has(READ_NEGATIVE_STRAND)
        .cmp(&b.has(READ_NEGATIVE_STRAND))
        .then_with(|| a.read_name().cmp(b.read_name()))
        .then_with(|| a.flags().cmp(&b.flags()))
        .then_with(|| a.mapping_quality().cmp(&b.mapping_quality()))
        .then_with(|| field(a, 24).cmp(&field(b, 24)))
        .then_with(|| field(a, 28).cmp(&field(b, 28)))
        .then_with(|| field(a, 32).cmp(&field(b, 32)))
}

/// `SAMRecordQueryNameComparator`, as [`query_name::compare`] ports it, on encoded records.
fn queryname_compare(a: &Encoded, b: &Encoded) -> Ordering {
    let cmp = a.read_name().cmp(b.read_name());
    if cmp != Ordering::Equal {
        return cmp;
    }
    let (a_paired, b_paired) = (a.has(READ_PAIRED), b.has(READ_PAIRED));
    if a_paired || b_paired {
        if !a_paired {
            return Ordering::Greater;
        }
        if !b_paired {
            return Ordering::Less;
        }
        if a.has(FIRST_OF_PAIR) && b.has(SECOND_OF_PAIR) {
            return Ordering::Less;
        }
        if a.has(SECOND_OF_PAIR) && b.has(FIRST_OF_PAIR) {
            return Ordering::Greater;
        }
    }
    for bit in [
        READ_NEGATIVE_STRAND,
        SECONDARY_ALIGNMENT,
        SUPPLEMENTARY_ALIGNMENT,
    ] {
        let cmp = a.has(bit).cmp(&b.has(bit));
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    match (a.hit, b.hit) {
        (Some(h1), Some(h2)) => h1.cmp(&h2),
        (Some(_), None) => Ordering::Greater,
        (None, Some(_)) => Ordering::Less,
        (None, None) => Ordering::Equal,
    }
}

type Compare = fn(&Encoded, &Encoded) -> Ordering;

impl SortOrder {
    /// [`Self::comparator`] on records as a [`SortingCollection`] holds them.
    fn encoded_comparator(self) -> Compare {
        match self {
            SortOrder::Coordinate => coordinate_compare,
            SortOrder::Queryname => queryname_compare,
        }
    }
}

/// Where one record of the batch in RAM starts, in which of its arenas, with what its bytes do not
/// carry.
struct Entry {
    arena: u32,
    offset: usize,
    hit: Option<i64>,
    end: i32,
}

/// Appends `record` to `arena` as `BamWriter::write` encodes it, with a bin of 0 where `bin_zero`
/// says its reference is too long for one; on an error `arena` is as it was.
fn push_encoded(
    record: &BamRecord,
    bin_zero: bool,
    arena: &mut Vec<u8>,
) -> Result<Entry, EncodeError> {
    let offset = arena.len();
    let encoded = if bin_zero {
        record.encode_with_bin_into(0, arena)
    } else {
        record.encode_into(arena)
    };
    if let Err(error) = encoded {
        arena.truncate(offset);
        return Err(error);
    }
    let hit = match record.tags.get(Tag::new(b"HI")) {
        Some(TagValue::Int(v)) => Some(*v),
        _ => None,
    };
    Ok(Entry {
        arena: 0,
        offset,
        hit,
        end: record.alignment_end(),
    })
}

/// Records encoded as [`SortingCollection::add`] encodes them, on any thread, for
/// [`SortingCollection::add_batch`] to add in order.
#[derive(Default)]
pub struct EncodedBatch {
    bytes: Vec<u8>,
    entries: Vec<Entry>,
}

impl EncodedBatch {
    /// A batch with room for `bytes` of records: about what the bytes they were decoded from take.
    pub fn with_capacity(bytes: usize) -> Self {
        EncodedBatch {
            bytes: Vec::with_capacity(bytes),
            entries: Vec::new(),
        }
    }

    pub fn push(&mut self, record: &BamRecord, bin_zero: bool) -> Result<(), EncodeError> {
        let entry = push_encoded(record, bin_zero, &mut self.bytes)?;
        self.entries.push(entry);
        Ok(())
    }

    /// The records from the `at`th on, moved to a batch of their own.
    fn split_off(&mut self, at: usize) -> EncodedBatch {
        let from = self.entries[at].offset;
        let entries = self.entries.split_off(at).into_iter();
        EncodedBatch {
            bytes: self.bytes.split_off(from),
            entries: entries
                .map(|entry| Entry {
                    offset: entry.offset - from,
                    ..entry
                })
                .collect(),
        }
    }
}

/// htsjdk's `SortingCollection` of records, holding each as `BAMRecordCodec` encodes it rather
/// than decoded: the bytes are written to the runs and to the output as they are, so a record is
/// encoded once.
pub struct SortingCollection {
    compare: Compare,
    /// Whether the batch in RAM may be sorted on every core.
    parallel: bool,
    max_records_in_ram: usize,
    tmp_dir: PathBuf,
    /// The bytes of the batch in RAM: one arena that [`Self::add`] appends to, or one per
    /// [`EncodedBatch`] added.
    arenas: Vec<Vec<u8>>,
    entries: Vec<Entry>,
    runs: Vec<File>,
    /// The last batch spilled, still being written on its own thread while the next one fills.
    writing: Option<std::thread::JoinHandle<io::Result<File>>>,
}

impl SortingCollection {
    pub fn new(order: SortOrder, max_records_in_ram: usize, tmp_dir: PathBuf) -> Self {
        SortingCollection {
            compare: order.encoded_comparator(),
            // The coordinate comparator is a total preorder, so every stable sort orders as the
            // serial one does; the queryname one is not transitive (see the tests), and keeps the
            // serial sort's answer.
            parallel: order == SortOrder::Coordinate,
            max_records_in_ram: max_records_in_ram.max(1),
            tmp_dir,
            arenas: Vec::new(),
            entries: Vec::new(),
            runs: Vec::new(),
            writing: None,
        }
    }

    /// `add`: spill first when the buffer is already full, as htsjdk does. The record is encoded
    /// as `BamWriter::write` encodes it, with a bin of 0 where `bin_zero` says its reference is
    /// too long for one.
    pub fn add(&mut self, record: &BamRecord, bin_zero: bool) -> Result<(), WriteError> {
        if self.entries.len() == self.max_records_in_ram {
            self.spill().map_err(WriteError::Io)?;
        }
        if self.arenas.is_empty() {
            self.arenas.push(Vec::new());
        }
        let arena = self.arenas.len() - 1;
        let entry =
            push_encoded(record, bin_zero, &mut self.arenas[arena]).map_err(WriteError::Encode)?;
        self.entries.push(Entry {
            arena: arena as u32,
            ..entry
        });
        Ok(())
    }

    /// [`Self::add`] for each record of `batch`, in its order, its bytes kept where they are
    /// unless a spill falls inside it.
    pub fn add_batch(&mut self, mut batch: EncodedBatch) -> io::Result<()> {
        while !batch.entries.is_empty() {
            if self.entries.len() == self.max_records_in_ram {
                self.spill()?;
            }
            let room = self.max_records_in_ram - self.entries.len();
            let rest = (batch.entries.len() > room).then(|| batch.split_off(room));
            let arena = self.arenas.len() as u32;
            self.arenas.push(batch.bytes);
            self.entries.extend(
                batch
                    .entries
                    .into_iter()
                    .map(|entry| Entry { arena, ..entry }),
            );
            batch = rest.unwrap_or_default();
        }
        Ok(())
    }

    /// How many runs have gone to disk.
    pub fn spills(&self) -> usize {
        self.runs.len() + self.writing.is_some() as usize
    }

    /// The batch in RAM, sorted.
    fn sort(&mut self) {
        sort_entries(&mut self.entries, &self.arenas, self.compare, self.parallel);
    }

    /// Sorts the batch in RAM and writes it to a run, both on a thread of its own while the next
    /// batch fills.
    fn spill(&mut self) -> io::Result<()> {
        self.wait_for_run()?;
        let mut run = RunWriter::create(&self.tmp_dir)?;
        // shortcut: an error writing this run is answered at the next spill or at the end, after
        // more records are read; the serial write answered it before reading on.
        let arenas = std::mem::take(&mut self.arenas);
        let records = self.entries.len();
        let mut entries = std::mem::replace(&mut self.entries, Vec::with_capacity(records));
        let (compare, parallel) = (self.compare, self.parallel);
        self.writing = Some(std::thread::spawn(move || {
            sort_entries(&mut entries, &arenas, compare, parallel);
            for entry in &entries {
                run.push(view(&arenas, entry))?;
            }
            run.finish()
        }));
        if self.runs.len() + 1 == MAX_OPEN_RUNS {
            self.wait_for_run()?;
            let runs = std::mem::take(&mut self.runs).into_iter().map(Run::file);
            let mut merged = Merge::new(self.compare, self.parallel, runs.collect())?;
            let mut run = RunWriter::create(&self.tmp_dir)?;
            merged.for_each(|record| run.push(record))?;
            self.runs.push(run.finish()?);
        }
        Ok(())
    }

    /// Adds the run being written to the runs, once it is written.
    fn wait_for_run(&mut self) -> io::Result<()> {
        if let Some(writing) = self.writing.take() {
            let run = writing.join().expect("the run writer does not panic")?;
            self.runs.push(run);
        }
        Ok(())
    }

    /// `doneAdding` and `iterator`: hands `emit` the records in order, from RAM when nothing
    /// spilled.
    pub fn for_each_sorted<E: From<io::Error>>(
        mut self,
        mut emit: impl FnMut(Encoded) -> Result<(), E>,
    ) -> Result<(), E> {
        if self.spills() == 0 {
            self.sort();
            for entry in &self.entries {
                emit(view(&self.arenas, entry))?;
            }
            return Ok(());
        }
        // The last batch is merged from RAM, as the run it would have been spilled to: the last.
        self.sort();
        self.wait_for_run()?;
        let mut runs: Vec<Run> = self.runs.into_iter().map(Run::file).collect();
        if !self.entries.is_empty() {
            runs.push(Run::Ram {
                arenas: self.arenas,
                entries: self.entries.into_iter(),
            });
        }
        Merge::new(self.compare, self.parallel, runs)?.for_each(emit)
    }
}

/// Sorts a batch's entries. Stable, so equal records keep input order (decision 0021).
fn sort_entries(entries: &mut [Entry], arenas: &[Vec<u8>], compare: Compare, parallel: bool) {
    let by = |a: &Entry, b: &Entry| compare(&view(arenas, a), &view(arenas, b));
    if parallel {
        entries.par_sort_by(by);
    } else {
        entries.sort_by(by);
    }
}

fn view<'a>(arenas: &'a [Vec<u8>], entry: &Entry) -> Encoded<'a> {
    let arena = &arenas[entry.arena as usize];
    let size = i32::from_le_bytes(arena[entry.offset..entry.offset + 4].try_into().unwrap());
    Encoded {
        bytes: &arena[entry.offset..entry.offset + 4 + size as usize],
        hit: entry.hit,
        end: entry.end,
    }
}

/// How much of a run is written or read at a time: a run is a whole batch, and the default 8 KiB
/// made a spill wait on its own system calls.
const RUN_BUFFER: usize = 1 << 20;

/// A fresh file in `dir`, unlinked at once so nothing is left behind, that records are written to
/// each as whether it has an `HI`, its end, the `HI` if it has one, how much of its name the record
/// before it has, and its bytes without that much of the name. Records are not otherwise compressed:
/// LZ4 took a third off a run for 7% more CPU over a whole sort, where a name's prefix costs a compare.
struct RunWriter {
    out: BufWriter<File>,
    name: Vec<u8>,
}

/// Where a record's name starts in its bytes, and the byte giving its length, NUL included.
const NAME_AT: usize = 36;
const NAME_LENGTH_AT: usize = 12;

impl RunWriter {
    fn create(dir: &Path) -> io::Result<RunWriter> {
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        // Not the decoded collection's name: each counts its own runs, and both run in one test
        // process.
        let path = dir.join(format!(
            "sortingcollection.encoded.{}.{n}.tmp",
            std::process::id()
        ));
        let file = File::options()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)?;
        // shortcut: unlinking an open file keeps it readable on Unix only; Windows is not a target.
        std::fs::remove_file(&path)?;
        Ok(RunWriter {
            out: BufWriter::with_capacity(RUN_BUFFER, file),
            name: Vec::new(),
        })
    }

    fn push(&mut self, record: Encoded) -> io::Result<()> {
        self.out.write_all(&[record.hit.is_some() as u8])?;
        self.out.write_all(&record.end.to_le_bytes())?;
        if let Some(hit) = record.hit {
            self.out.write_all(&hit.to_le_bytes())?;
        }
        let name_end = NAME_AT + record.bytes[NAME_LENGTH_AT] as usize;
        let name = &record.bytes[NAME_AT..name_end];
        let shared = shared_prefix(name, &self.name);
        self.out.write_all(&[shared as u8])?;
        self.out.write_all(&record.bytes[..NAME_AT])?;
        self.out.write_all(&record.bytes[NAME_AT + shared..])?;
        self.name.clear();
        self.name.extend_from_slice(name);
        Ok(())
    }

    /// The run, rewound for reading.
    fn finish(self) -> io::Result<File> {
        let mut file = self.out.into_inner().map_err(|e| e.into_error())?;
        file.seek(SeekFrom::Start(0))?;
        Ok(file)
    }
}

/// A record read back from a run.
#[derive(Default)]
struct Owned {
    bytes: Vec<u8>,
    hit: Option<i64>,
    end: i32,
}

impl Owned {
    fn view(&self) -> Encoded<'_> {
        Encoded {
            bytes: &self.bytes,
            hit: self.hit,
            end: self.end,
        }
    }
}

/// How many leading bytes `a` and `b` have in common, eight at a time: this is asked once a record.
fn shared_prefix(a: &[u8], b: &[u8]) -> usize {
    let mut shared = 0;
    for (x, y) in a.as_chunks::<8>().0.iter().zip(b.as_chunks::<8>().0) {
        let difference = u64::from_le_bytes(*x) ^ u64::from_le_bytes(*y);
        if difference != 0 {
            return shared + (difference.trailing_zeros() / 8) as usize;
        }
        shared += 8;
    }
    let tails = a[shared..].iter().zip(&b[shared..]);
    shared + tails.take_while(|(x, y)| x == y).count()
}

/// One sorted run, read back record by record.
enum Run {
    File {
        reader: BufReader<File>,
        /// The name of the record read last.
        name: Vec<u8>,
    },
    /// The last batch, sorted in RAM and never spilled.
    Ram {
        arenas: Vec<Vec<u8>>,
        entries: std::vec::IntoIter<Entry>,
    },
}

impl Run {
    fn file(file: File) -> Run {
        Run::File {
            reader: BufReader::with_capacity(RUN_BUFFER, file),
            name: Vec::new(),
        }
    }

    /// Reads the next record over `record`; `false` when the run has ended.
    fn next_into(&mut self, record: &mut Owned) -> io::Result<bool> {
        let (reader, name) = match self {
            Run::File { reader, name } => (reader, name),
            Run::Ram { arenas, entries } => {
                let Some(entry) = entries.next() else {
                    return Ok(false);
                };
                record.bytes.clear();
                record.bytes.extend_from_slice(view(arenas, &entry).bytes);
                record.hit = entry.hit;
                record.end = entry.end;
                return Ok(true);
            }
        };
        if reader.fill_buf()?.is_empty() {
            return Ok(false);
        }
        let mut fixed = [0u8; 5];
        reader.read_exact(&mut fixed)?;
        record.end = i32::from_le_bytes(fixed[1..5].try_into().unwrap());
        record.hit = if fixed[0] != 0 {
            let mut hit = [0u8; 8];
            reader.read_exact(&mut hit)?;
            Some(i64::from_le_bytes(hit))
        } else {
            None
        };
        let mut shared = [0u8; 1];
        reader.read_exact(&mut shared)?;
        let shared = shared[0] as usize;
        let mut length = [0u8; 4];
        reader.read_exact(&mut length)?;
        let size = i32::from_le_bytes(length) as usize;
        record.bytes.clear();
        record.bytes.extend_from_slice(&length);
        record.bytes.resize(4 + size, 0);
        reader.read_exact(&mut record.bytes[4..NAME_AT])?;
        record.bytes[NAME_AT..NAME_AT + shared].copy_from_slice(&name[..shared]);
        reader.read_exact(&mut record.bytes[NAME_AT + shared..])?;
        name.clear();
        let name_end = NAME_AT + record.bytes[NAME_LENGTH_AT] as usize;
        name.extend_from_slice(&record.bytes[NAME_AT..name_end]);
        Ok(true)
    }
}

/// The next record of one run, ordered so the heap's top is the least record, and of equal
/// records the one from the earliest run: `MergingIterator`'s tie-break.
struct Head {
    record: Owned,
    run: usize,
    compare: Compare,
}

impl Ord for Head {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.compare)(&other.record.view(), &self.record.view()).then(other.run.cmp(&self.run))
    }
}

impl PartialOrd for Head {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Head {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Head {}

/// A k-way merge of sorted runs.
struct Merge {
    compare: Compare,
    total_preorder: bool,
    runs: Vec<Run>,
    heap: BinaryHeap<Head>,
}

impl Merge {
    fn new(compare: Compare, total_preorder: bool, runs: Vec<Run>) -> io::Result<Merge> {
        let mut merge = Merge {
            compare,
            total_preorder,
            runs,
            heap: BinaryHeap::new(),
        };
        for run in 0..merge.runs.len() {
            merge.refill(run)?;
        }
        Ok(merge)
    }

    fn refill(&mut self, run: usize) -> io::Result<()> {
        let mut record = Owned::default();
        if self.runs[run].next_into(&mut record)? {
            self.heap.push(Head {
                record,
                run,
                compare: self.compare,
            });
        }
        Ok(())
    }

    /// Hands `emit` every record in order, reading each run's next record over the buffer of the
    /// one it replaces.
    fn for_each<E: From<io::Error>>(
        &mut self,
        mut emit: impl FnMut(Encoded) -> Result<(), E>,
    ) -> Result<(), E> {
        if self.total_preorder {
            // Sifted down once, half the comparisons of a pop and a push. Only a total preorder
            // gives the same order whichever way the heap is walked.
            while let Some(mut head) = self.heap.peek_mut() {
                emit(head.record.view())?;
                let run = head.run;
                if !self.runs[run].next_into(&mut head.record)? {
                    PeekMut::pop(head);
                }
            }
        } else {
            while let Some(mut head) = self.heap.pop() {
                emit(head.record.view())?;
                if self.runs[head.run].next_into(&mut head.record)? {
                    self.heap.push(head);
                }
            }
        }
        Ok(())
    }
}

/// The header of the clean BAM at `input`, sort order set, and its records in a collection,
/// decoded, checked and encoded on every core a batch at a time and added in file order: what
/// adding them one by one as [`crate::bam_stream::Input`] hands them over makes. `None` where
/// any record is one that read answers for: a stream that does not end cleanly, or a record that
/// does not decode, does not encode, or does not pass validation in silence.
pub fn read_in_parallel(
    input: &str,
    sort_order: SortOrder,
    max_records_in_ram: usize,
    tmp_dir: &Path,
) -> io::Result<Option<(SamHeader, SortingCollection)>> {
    let Some((mut header, validator, mut batches)) = crate::bam_stream::batches(input) else {
        return Ok(None);
    };
    let lengths: Vec<i32> = header.sequences.iter().map(|s| s.length).collect();
    let encode = |raw: &Vec<u8>| {
        let mut batch = EncodedBatch::with_capacity(raw.len());
        let mut at = 0;
        while let Ok(Some((record, used))) = BamRecord::decode(&raw[at..]) {
            at += used;
            if !validator.passes(&record) {
                return None;
            }
            batch
                .push(&record, bin_zero(&lengths, record.reference_index))
                .ok()?;
        }
        // A batch holds whole records only, so anything left did not decode.
        (at == raw.len()).then_some(batch)
    };
    let mut sorting = SortingCollection::new(sort_order, max_records_in_ram, tmp_dir.into());
    let threads = rayon::current_num_threads();
    // The stream is cut on a thread of its own, a round or two ahead of the decoding: cutting is
    // the one serial step, so it never waits for the cores.
    let (rounds, received) = std::sync::mpsc::sync_channel(2);
    std::thread::spawn(move || loop {
        // Four batches a thread, so the cores that finish first take the rest.
        let round = batches.next_round(threads * 4);
        let last = !matches!(&round, Some(raw) if !raw.is_empty());
        if rounds.send(round).is_err() || last {
            break;
        }
    });
    while let Ok(Some(raw)) = received.recv() {
        if raw.is_empty() {
            header.set_sort_order(sort_order.name());
            return Ok(Some((header, sorting)));
        }
        let Some(encoded) = raw.par_iter().map(encode).collect::<Option<Vec<_>>>() else {
            return Ok(None);
        };
        for batch in encoded {
            sorting.add_batch(batch)?;
        }
    }
    Ok(None)
}

/// `BAMFileWriter` for records a [`SortingCollection`] encoded: the header is what
/// [`crate::output::bam_writer`] writes, and each record's bytes and index entry are what its
/// `write` makes of the record they were encoded from.
pub struct EncodedBamWriter<W: Write> {
    bgzf: BgzfWriter<Logged<W>>,
    reference_lengths: Vec<i32>,
    indexer: Option<BamIndexer>,
    /// Where each block reached the file, read off the stream under the writer.
    addresses: Rc<RefCell<Addresses>>,
    /// Inflated bytes written so far, header included.
    written: u64,
    /// Records written whose index entry waits for the address of a block still being deflated.
    pending: VecDeque<Pending>,
}

/// What [`BamIndexer::process`] takes for one record, its chunk as inflated positions.
struct Pending {
    reference: i32,
    alignment_start: i32,
    end: i32,
    index_bin: i32,
    unmapped: bool,
    from: u64,
    to: u64,
}

impl<W: Write> EncodedBamWriter<W> {
    pub fn new(inner: W, header: &SamHeader, index: bool) -> io::Result<Self> {
        // The header alone, through the writer every tool writes with, inflated back.
        let block = crate::output::bam_writer(Vec::new(), header)?.finish()?;
        let text =
            crate::bgzf::decompress_owned(block).map_err(|e| io::Error::other(format!("{e:?}")))?;
        let addresses = Rc::new(RefCell::new(Addresses::default()));
        let mut bgzf = crate::output::bgzf_writer(Logged {
            inner,
            addresses: Rc::clone(&addresses),
        });
        bgzf.write_all(&text)?;
        let reference_lengths: Vec<i32> = header.sequences.iter().map(|s| s.length).collect();
        Ok(EncodedBamWriter {
            bgzf,
            indexer: index.then(|| BamIndexer::new(reference_lengths.clone())),
            reference_lengths,
            addresses,
            written: text.len() as u64,
            pending: VecDeque::new(),
        })
    }

    /// `BamWriter::write`, from the record's bytes.
    pub fn write(&mut self, record: Encoded) -> io::Result<()> {
        self.bgzf.write_all(record.bytes)?;
        let from = self.written;
        self.written += record.bytes.len() as u64;
        if self.indexer.is_none() {
            return Ok(());
        }
        // `file_pointer` before and after the record would wait for every block in flight, so
        // the chunk is kept as inflated positions until the blocks they fall in are written.
        let (reference, alignment_start) = (record.reference_index(), record.alignment_start());
        let index_bin = if bin_zero(&self.reference_lengths, reference)
            || alignment_start == bin::NO_ALIGNMENT_START
        {
            0
        } else {
            bin::compute_indexing_bin(alignment_start, record.end).unwrap_or(0)
        };
        self.pending.push_back(Pending {
            reference,
            alignment_start,
            end: record.end,
            index_bin,
            unmapped: record.has(READ_UNMAPPED),
            from,
            to: self.written,
        });
        self.replay();
        Ok(())
    }

    /// Hands the indexer, in order, every waiting record whose blocks have been written.
    fn replay(&mut self) {
        let Some(indexer) = &mut self.indexer else {
            return;
        };
        let mut addresses = self.addresses.borrow_mut();
        while let Some(next) = self.pending.front() {
            let (Some(start), Some(end)) =
                (addresses.pointer(next.from), addresses.pointer(next.to))
            else {
                break;
            };
            indexer.process(
                next.reference,
                next.alignment_start,
                next.end,
                next.index_bin,
                next.unmapped,
                Chunk { start, end },
            );
            self.pending.pop_front();
        }
        addresses.forget_before(self.pending.front().map_or(self.written, |p| p.from));
    }

    /// Closes the stream with its terminator block, and the BAI when one was built.
    pub fn finish(mut self) -> io::Result<(W, Option<Vec<u8>>)> {
        self.bgzf.finish()?;
        self.replay();
        assert!(self.pending.is_empty(), "every block is written");
        let index = self.indexer.map(BamIndexer::finish);
        Ok((self.bgzf.into_inner()?.inner, index))
    }
}

/// The stream under a [`BgzfWriter`], noting where each block starts as its header goes by.
struct Logged<W> {
    inner: W,
    addresses: Rc<RefCell<Addresses>>,
}

impl<W: Write> Write for Logged<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(bytes)?;
        self.addresses.borrow_mut().saw(&bytes[..n]);
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// The compressed address of each block from `first` on, read from the block headers: blocks
/// are cut every [`DEFAULT_UNCOMPRESSED_BLOCK_SIZE`] inflated bytes, nothing flushing between, so
/// an inflated position names its block, and the block's address is what `file_pointer` gives.
struct Addresses {
    /// The address of block `first + i`; the last is that of the block whose header is next.
    known: VecDeque<u64>,
    first: u64,
    /// Compressed bytes seen, where the block they are in starts, and its header so far.
    at: u64,
    block: u64,
    header: [u8; BGZF_HEADER],
}

/// A BGZF block header's length; the block's own length, less one, is in its last two bytes.
const BGZF_HEADER: usize = 18;

impl Default for Addresses {
    fn default() -> Self {
        Addresses {
            known: VecDeque::from([0]),
            first: 0,
            at: 0,
            block: 0,
            header: [0; BGZF_HEADER],
        }
    }
}

impl Addresses {
    fn saw(&mut self, mut bytes: &[u8]) {
        while !bytes.is_empty() {
            let into = (self.at - self.block) as usize;
            let n = if into < BGZF_HEADER {
                let n = (BGZF_HEADER - into).min(bytes.len());
                self.header[into..into + n].copy_from_slice(&bytes[..n]);
                if into + n == BGZF_HEADER {
                    let length = u16::from_le_bytes([self.header[16], self.header[17]]) as u64 + 1;
                    self.known.push_back(self.block + length);
                }
                n
            } else {
                let next = *self.known.back().expect("never empty");
                ((next - self.at) as usize).min(bytes.len())
            };
            self.at += n as u64;
            bytes = &bytes[n..];
            if Some(&self.at) == self.known.back() && self.at - self.block > BGZF_HEADER as u64 {
                self.block = self.at;
            }
        }
    }

    /// `file_pointer` at inflated position `at`, once the address of its block is known.
    fn pointer(&self, at: u64) -> Option<u64> {
        let block = at / DEFAULT_UNCOMPRESSED_BLOCK_SIZE as u64;
        let address = *self
            .known
            .get(usize::try_from(block.checked_sub(self.first)?).ok()?)?;
        let offset = (at % DEFAULT_UNCOMPRESSED_BLOCK_SIZE as u64) as u32;
        Some(make_file_pointer(address, offset).expect("a BAM's addresses fit"))
    }

    /// Drops the addresses of blocks before the one inflated position `at` falls in.
    fn forget_before(&mut self, at: u64) {
        let block = at / DEFAULT_UNCOMPRESSED_BLOCK_SIZE as u64;
        while self.first < block && self.known.len() > 1 {
            self.known.pop_front();
            self.first += 1;
        }
    }
}

/// Whether `BamWriter::write` writes a record on this reference with a bin of 0: one too long for
/// the BAI bin field.
pub fn bin_zero(reference_lengths: &[i32], reference_index: i32) -> bool {
    usize::try_from(reference_index)
        .ok()
        .and_then(|i| reference_lengths.get(i))
        .is_some_and(|&length| length > bin::BIN_GENOMIC_SPAN)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sort_sam::DEFAULT_MAX_RECORDS_IN_RAM;
    use htsjdk_bam::sam_file::read_sam;

    /// Sorts through a collection holding `max` records at most.
    fn sorted_names(sam: &str, order: SortOrder, max: usize) -> (Vec<Vec<u8>>, usize) {
        let (_, records) = read_sam(sam).unwrap();
        let mut sorting = SortingCollection::new(order, max, std::env::temp_dir());
        for record in &records {
            sorting.add(record, false).unwrap();
        }
        let spills = sorting.spills();
        let mut out = Vec::new();
        sorting
            .for_each_sorted(|r| {
                out.push(r.bytes.to_vec());
                Ok::<(), io::Error>(())
            })
            .unwrap();
        (out, spills)
    }

    /// The records sorted decoded, by htsjdk's comparator, as the collection used to hold them.
    fn sorted_decoded(sam: &str, order: SortOrder) -> Vec<Vec<u8>> {
        let (_, mut records) = read_sam(sam).unwrap();
        records.sort_by(order.comparator());
        records.iter().map(|r| r.encode().unwrap()).collect()
    }

    #[test]
    fn a_run_gives_back_what_it_was_given_for_less_than_its_names() {
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n");
        // Long shared prefixes, a name that is a prefix of the next, one sharing nothing, the same
        // name twice, and the unnamed read.
        for name in [
            "A00123:456:HXXXXDSXY:1:1101:1000:1000",
            "A00123:456:HXXXXDSXY:1:1101:1000:1001",
            "A00123:456:HXXXXDSXY:1:1101:1000:1001",
            "A00123",
            "A00123:456",
            "z",
            "*",
            "A00123:456:HXXXXDSXY:1:1102:2000:1000",
        ] {
            sam.push_str(&format!("{name}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"));
        }
        let (_, records) = read_sam(&sam).unwrap();
        let encoded: Vec<Vec<u8>> = records.iter().map(|r| r.encode().unwrap()).collect();
        let mut run = RunWriter::create(&std::env::temp_dir()).unwrap();
        for (i, bytes) in encoded.iter().enumerate() {
            run.push(Encoded {
                bytes,
                hit: (i % 2 == 0).then_some(i as i64),
                end: i as i32,
            })
            .unwrap();
        }
        let file = run.finish().unwrap();
        let everything: usize = encoded.iter().map(Vec::len).sum();
        assert!((file.metadata().unwrap().len() as usize) < everything);
        let mut run = Run::file(file);
        let mut record = Owned::default();
        for (i, bytes) in encoded.iter().enumerate() {
            assert!(run.next_into(&mut record).unwrap());
            assert_eq!(&record.bytes, bytes);
            assert_eq!(record.hit, (i % 2 == 0).then_some(i as i64));
            assert_eq!(record.end, i as i32);
        }
        assert!(!run.next_into(&mut record).unwrap());
    }

    #[test]
    fn the_shared_prefix_counts_what_comparing_bytes_counts() {
        let long = b"A00123:456:HXXXXDSXY:1:1101:1000:1000\0";
        for length in 0..=long.len() {
            for at in 0..length {
                let mut other = long[..length].to_vec();
                other[at] ^= 0x40;
                assert_eq!(shared_prefix(&long[..length], &other), at);
            }
            assert_eq!(shared_prefix(&long[..length], long), length);
            assert_eq!(shared_prefix(long, &long[..length]), length);
        }
    }

    #[test]
    fn a_spilled_sort_writes_what_the_in_memory_sort_writes() {
        // Equal under both comparators except for the tag, so only stability orders them.
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n");
        for i in 0..300 {
            let name = ["b", "a", "c"][i % 3];
            let pos = [500, 100, 300, 100][i % 4];
            sam.push_str(&format!(
                "{name}\t0\tchr1\t{pos}\t60\t4M\t*\t0\t0\tACGT\tIIII\tXI:i:{i}\n"
            ));
        }
        for order in [SortOrder::Coordinate, SortOrder::Queryname] {
            let (whole, none) = sorted_names(&sam, order, DEFAULT_MAX_RECORDS_IN_RAM);
            assert_eq!(none, 0);
            assert_eq!(whole, sorted_decoded(&sam, order), "{order:?}");
            for max in [2, 7, 299, 300] {
                let (spilled, spills) = sorted_names(&sam, order, max);
                assert_eq!(spills, 300usize.div_ceil(max) - 1, "max {max}");
                assert_eq!(spilled, whole, "{order:?} with MAX_RECORDS_IN_RAM={max}");
            }
        }
    }

    #[test]
    fn the_encoded_comparators_order_as_htsjdks_do() {
        // Every field either comparator reads, varied so that each of them decides some pairs.
        let mut sam = String::from(
            "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n",
        );
        let flags = [
            0,
            1,
            65,
            129,
            83,
            163,
            16,
            256,
            2048,
            4,
            0x41 | 0x80,
            113,
            177,
        ];
        for i in 0..600usize {
            let name = ["q", "p", "q2", "p"][i % 4];
            let flag = flags[i % flags.len()];
            let (rname, pos) = match i % 5 {
                0 => ("*", 0),
                1 => ("chr2", 7 + i % 3),
                _ => ("chr1", 7 + i % 2),
            };
            let cigar = if rname == "*" { "*" } else { "4M" };
            let mref = match ["=", "chr2", "*"][i % 3] {
                "=" if rname == "*" => "*",
                mref => mref,
            };
            let mpos = if mref == "*" { 0 } else { 1 + i % 4 };
            // Consistent as the SAM parser wants: no reference is unmapped, no mate reference is a
            // mate unmapped.
            let flag = flag
                | if rname == "*" { 0x4 } else { 0 }
                | if mref == "*" && flag & 1 != 0 { 0x8 } else { 0 };
            let mq = if rname == "*" { 0 } else { [60, 0, 37][i % 3] };
            let tlen = [0, 50, -50][i % 3];
            let hi = match i % 7 {
                0 => String::new(),
                k => format!("\tHI:i:{}", k % 3),
            };
            sam.push_str(&format!(
                "{name}\t{flag}\t{rname}\t{pos}\t{mq}\t{cigar}\t{mref}\t{mpos}\t{tlen}\tACGT\tIIII\tXI:i:{i}{hi}\n"
            ));
        }
        for order in [SortOrder::Coordinate, SortOrder::Queryname] {
            let whole = sorted_names(&sam, order, DEFAULT_MAX_RECORDS_IN_RAM).0;
            assert_eq!(whole, sorted_decoded(&sam, order), "{order:?}");
        }
        // The queryname comparator is not transitive over these flags (a record flagged both first
        // and second of pair is less and greater than its mate), so merged runs order as the merge
        // decides; the coordinate one is a total preorder, and spilling changes nothing.
        let spilled = sorted_names(&sam, SortOrder::Coordinate, 13).0;
        assert_eq!(spilled, sorted_decoded(&sam, SortOrder::Coordinate));
    }

    /// A BAM of `n` valid records, many of them equal under both comparators, written to a file
    /// of its own in the temporary directory.
    fn bam_of(n: usize, name: &str) -> std::path::PathBuf {
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n");
        let bases: String = (0..100)
            .map(|k| ['A', 'C', 'G', 'T'][(k * 7 + k / 3) % 4])
            .collect();
        let quals = "I".repeat(100);
        for i in 0..n {
            let name = format!("q{}", (i * 7919) % 97);
            let pos = 1 + (i * 104_729) % 5000;
            let hi = if i % 5 == 0 {
                format!("\tHI:i:{}", i % 3)
            } else {
                String::new()
            };
            let line = match i % 6 {
                0 => format!("{name}\t99\tchr1\t{pos}\t60\t100M\t=\t{}\t250", pos + 150),
                1 => format!("{name}\t147\tchr1\t{pos}\t60\t100M\t=\t{}\t-250", pos + 150),
                2 => format!("{name}\t16\tchr1\t{pos}\t{}\t100M\t*\t0\t0", i % 61),
                3 => format!("{name}\t4\t*\t0\t0\t*\t*\t0\t0"),
                4 => format!("{name}\t4\tchr1\t{pos}\t0\t*\t*\t0\t0"),
                _ => format!("{name}\t0\tchr1\t{pos}\t37\t60M40S\t*\t0\t0"),
            };
            sam.push_str(&format!("{line}\t{bases}\t{quals}\tXI:i:{i}{hi}\n"));
        }
        let (header, records) = read_sam(&sam).unwrap();
        let mut writer = crate::output::bam_writer(Vec::new(), &header).unwrap();
        for record in &records {
            writer.write(record).unwrap();
        }
        let path = std::env::temp_dir().join(format!("{name}.{}.bam", std::process::id()));
        std::fs::write(&path, writer.finish().unwrap()).unwrap();
        path
    }

    /// Every record's bytes, in the order the collection hands them over.
    fn emitted(sorting: SortingCollection) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        sorting
            .for_each_sorted(|r| {
                out.push(r.bytes.to_vec());
                Ok::<(), io::Error>(())
            })
            .unwrap();
        out
    }

    #[test]
    fn the_parallel_read_collects_what_adding_each_record_collects() {
        // Enough records for many batches, and so for spills that fall inside one.
        let path = bam_of(20_000, "parallel_read");
        let input = path.to_str().unwrap();
        for order in [SortOrder::Coordinate, SortOrder::Queryname] {
            for max in [DEFAULT_MAX_RECORDS_IN_RAM, 4096, 777] {
                let tmp = std::env::temp_dir();
                let read = crate::bam_stream::Input::open(input).unwrap();
                let (mut header, records) = read.into_parts();
                let lengths: Vec<i32> = header.sequences.iter().map(|s| s.length).collect();
                let mut serial = SortingCollection::new(order, max, tmp.clone());
                for record in records {
                    serial
                        .add(&record, bin_zero(&lengths, record.reference_index))
                        .unwrap();
                }
                header.set_sort_order(order.name());
                let (parallel_header, parallel) = read_in_parallel(input, order, max, &tmp)
                    .unwrap()
                    .expect("a clean BAM is read in parallel");
                assert_eq!(parallel_header, header);
                assert_eq!(parallel.spills(), serial.spills(), "{order:?} {max}");
                assert_eq!(emitted(parallel), emitted(serial), "{order:?} {max}");
            }
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn a_damaged_bam_is_left_to_the_serial_read() {
        let path = bam_of(5_000, "damaged_read");
        let mut bytes = std::fs::read(&path).unwrap();
        // Inside the terminator block, then inside a block of records.
        for cut in [bytes.len() - 10, bytes.len() / 2] {
            bytes.truncate(cut);
            std::fs::write(&path, &bytes).unwrap();
            let tmp = std::env::temp_dir();
            let read = read_in_parallel(path.to_str().unwrap(), SortOrder::Coordinate, 1000, &tmp);
            assert!(read.unwrap().is_none(), "cut at {cut}");
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn more_runs_than_can_be_open_merge_in_order() {
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n");
        for i in 0..(MAX_OPEN_RUNS * 2 + 5) {
            let pos = 1 + (i * 37) % 11;
            sam.push_str(&format!(
                "r\t0\tchr1\t{pos}\t60\t4M\t*\t0\t0\tACGT\tIIII\tXI:i:{i}\n"
            ));
        }
        let (spilled, _) = sorted_names(&sam, SortOrder::Coordinate, 1);
        assert_eq!(spilled, sorted_decoded(&sam, SortOrder::Coordinate));
    }

    #[test]
    fn the_index_is_the_one_a_pointer_per_record_gives() {
        let header_text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200000000\n";
        let (header, _) = read_sam(header_text).unwrap();
        let line = |i: usize, pad: usize| {
            format!(
                "r{i}\t0\tchr1\t{}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tZZ:Z:{}\n",
                // A linear-index window each, so every start reaches the index as it is.
                1 + i * 20_000,
                "x".repeat(pad)
            )
        };
        let size = |i: usize, pad: usize| {
            let (_, records) = read_sam(&format!("{header_text}{}", line(i, pad))).unwrap();
            records[0].encode().unwrap().len()
        };
        // The header's inflated bytes come first; then records, some of them sized to end exactly
        // where a block is cut, and the last ending mid-block.
        let block = crate::output::bam_writer(Vec::new(), &header)
            .unwrap()
            .finish()
            .unwrap();
        let mut at = crate::bgzf::decompress_owned(block).unwrap().len();
        let mut sam = String::from(header_text);
        let mut on_boundary = 0;
        for i in 0..6000 {
            let base = size(i, 0);
            let to_cut = DEFAULT_UNCOMPRESSED_BLOCK_SIZE - at % DEFAULT_UNCOMPRESSED_BLOCK_SIZE;
            let pad = if i % 3 == 0 && (base..base + 400).contains(&to_cut) {
                on_boundary += 1;
                to_cut - base
            } else {
                (i * 37) % 300
            };
            at += size(i, pad);
            sam.push_str(&line(i, pad));
        }
        assert!(on_boundary >= 5, "{on_boundary}");
        assert_ne!(at % DEFAULT_UNCOMPRESSED_BLOCK_SIZE, 0);
        assert!(at > 10 * DEFAULT_UNCOMPRESSED_BLOCK_SIZE, "several blocks");
        // Blocks in flight, which a pointer per record would have waited for.
        assert!(rayon::current_num_threads() > 1);

        let (_, records) = read_sam(&sam).unwrap();
        let mut sorting =
            SortingCollection::new(SortOrder::Coordinate, 1 << 20, std::env::temp_dir());
        for record in &records {
            sorting.add(record, false).unwrap();
        }
        let mut ours = EncodedBamWriter::new(Vec::new(), &header, true).unwrap();
        let mut cut_after = 0;
        sorting
            .for_each_sorted(|encoded| {
                ours.write(encoded)?;
                cut_after += (ours.written % DEFAULT_UNCOMPRESSED_BLOCK_SIZE as u64 == 0) as usize;
                Ok::<(), io::Error>(())
            })
            .unwrap();
        assert_eq!(cut_after, on_boundary);
        assert_eq!(ours.written as usize, at);
        let (bam, index) = ours.finish().unwrap();

        let mut theirs = crate::output::bam_writer(Vec::new(), &header)
            .unwrap()
            .with_index();
        for record in &records {
            theirs.write(record).unwrap();
        }
        let (their_bam, their_index) = theirs.finish_with_index().unwrap();
        assert_eq!(bam, their_bam);
        assert_eq!(index.unwrap(), their_index);
    }
}
