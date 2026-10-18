//! `SortSam`.
//!
//! Ports `picard.sam.SortSam.doWork` at tag 3.4.0: read a SAM/BAM, set the header's sort order, and
//! rewrite the records in that order. A record transform whose byte-identity is the whole output
//! file.
//!
//! `doWork` keeps the **input header unchanged except for the `SO` field**, and adds no `@PG` and no
//! timestamp, so the output is comparable raw. The records are then sorted by the chosen order's
//! comparator: `SAMRecordCoordinateComparator` for `coordinate`, `SAMRecordQueryNameComparator` for
//! `queryname`. Both comparators treat records that agree on every field as equal, so the sort must
//! be **stable** to be byte-identical (decision 0021); Rust's `sort_by` is.
//!
//! This ports the `coordinate` and `queryname` orders, and the precondition the `duplicate` order
//! is refused on. A duplicate sort compares records by, among other things, their MATE's unclipped
//! coordinate, which it reads out of the `MC` tag; a file whose records do not carry one cannot be
//! sorted that way at all, and the reference says so before it writes a single record. The
//! ordering itself is a separate surface, and unreachable from a corpus without mate cigars.
//!
//! The sort itself is htsjdk's `SortingCollection`, which `SAMFileWriterImpl` sorts through: it
//! holds `MAX_RECORDS_IN_RAM` records, spills each full batch to `TMP_DIR` as a sorted run, and
//! merges the runs. Measured in the oracle container (3,050 records, coordinate and queryname, SAM
//! and BAM out): Picard's output with `MAX_RECORDS_IN_RAM=37` is byte-identical to its in-memory
//! sort, because the sort is stable and the merge breaks ties by run, so where the runs end
//! changes nothing that is written.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::text_parse::ParseError;
use htsjdk_bam::{coordinate, query_name};

/// The `SORT_ORDER` values this port handles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Coordinate,
    Queryname,
}

/// What `--SORT_ORDER` may name.
///
/// Two of the three are sorts this port performs. The third is the duplicate order, which is
/// refused on any file whose paired records carry no mate cigar, and that is every file this
/// corpus has.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RequestedOrder {
    Coordinate,
    Queryname,
    Duplicate,
}

impl RequestedOrder {
    pub fn parse(name: &str) -> Option<RequestedOrder> {
        match name {
            "coordinate" => Some(RequestedOrder::Coordinate),
            "queryname" => Some(RequestedOrder::Queryname),
            "duplicate" => Some(RequestedOrder::Duplicate),
            _ => None,
        }
    }

    /// The `@HD SO` the output header carries.
    pub fn name(self) -> &'static str {
        match self {
            RequestedOrder::Coordinate => "coordinate",
            RequestedOrder::Queryname => "queryname",
            RequestedOrder::Duplicate => "duplicate",
        }
    }

    /// The sort this port performs for it, where it performs one.
    pub fn sort_order(self) -> Option<SortOrder> {
        match self {
            RequestedOrder::Coordinate => Some(SortOrder::Coordinate),
            RequestedOrder::Queryname => Some(SortOrder::Queryname),
            RequestedOrder::Duplicate => None,
        }
    }
}

/// One record, as much of it as the mate-cigar precondition looks at.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MateCigarCheck {
    pub read_name: String,
    pub paired: bool,
    pub first_of_pair: bool,
    pub unmapped: bool,
    pub mate_unmapped: bool,
    pub read_length: usize,
    pub reference: Option<String>,
    pub start: i32,
    pub end: i32,
    pub mate_cigar: Option<String>,
    /// [`crate::view_sam::sam_string`] of the record.
    pub sam_string: String,
}

impl MateCigarCheck {
    /// `SAMRecord.toString()`, which is what the refusal names the record with: the whole SAM
    /// line under htsjdk 5.0.0.
    pub fn description(&self) -> String {
        self.sam_string.clone()
    }

    /// Whether sorting by duplicate order would ask this record for its mate's coordinate.
    pub fn needs_mate_cigar(&self) -> bool {
        self.paired && !self.unmapped && !self.mate_unmapped
    }
}

/// The order the comparator touches the records in, which is what decides WHICH record a refusal
/// names.
///
/// The sort compares the second record against the first before anything else, so the second is
/// the first one asked for its mate's coordinate; after that each record is touched as it is
/// inserted, in file order.
pub fn touch_order(count: usize) -> Vec<usize> {
    match count {
        0 => Vec::new(),
        1 => vec![0],
        _ => {
            let mut order = vec![1, 0];
            order.extend(2..count);
            order
        }
    }
}

/// The refusal a duplicate sort answers a file with no mate cigars with.
///
/// `SAMUtils.getMateUnclippedStart` throws on the first record that needs a mate cigar and has
/// none, and the message carries the record as `SAMRecord.toString()` writes it. The space after
/// the colon is the unclipped-coordinate methods'; the mate-alignment-end method next to them
/// writes the same sentence WITHOUT it, which is a difference a reader would not invent.
pub fn mate_cigar_refusal(records: &[MateCigarCheck]) -> Option<String> {
    mate_cigar_refusal_in_order(records, &touch_order(records.len()))
}

/// The same refusal, for a caller that knows its own touch order.
///
/// `SortSam`'s order is the `SortingCollection`'s, which compares the second record against the
/// first and so names the second. `MergeSamFiles` reaches the same throw through
/// `MergingSamRecordIterator`, which touches each stream's records in file order and therefore
/// names the first. The message is identical; which record it names is not, and only running both
/// against the reference says which is which.
pub fn mate_cigar_refusal_in_order(records: &[MateCigarCheck], order: &[usize]) -> Option<String> {
    for &index in order {
        let record = &records[index];
        if record.needs_mate_cigar() && record.mate_cigar.is_none() {
            return Some(format!(
                "Mate CIGAR (Tag MC) not found: {}",
                record.description()
            ));
        }
    }
    None
}

impl SortOrder {
    /// `SAMFileHeader.SortOrder.getSortOrder()`: the string written into the `@HD SO` field.
    pub fn name(self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
        }
    }

    pub fn comparator(self) -> fn(&BamRecord, &BamRecord) -> std::cmp::Ordering {
        match self {
            SortOrder::Coordinate => coordinate::compare,
            SortOrder::Queryname => query_name::compare,
        }
    }
}

/// `SortSam.doWork` for SAM input and output: set the sort order on the input header, sort the
/// records, and write. The header is otherwise the input header, byte for byte.
pub fn sort_sam(input_sam: &str, order: SortOrder) -> Result<String, ParseError> {
    let (mut header, mut records) = read_sam(input_sam)?;
    header.set_sort_order(order.name());
    // A stable sort so equal records keep input order (decision 0021).
    records.sort_by(order.comparator());
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

/// `SortSam.doWork` for SAM input and **BAM** output: the same sort, written as a BAM through
/// htsjdk-rs's byte-identical `BamWriter`. The returned bytes are the whole BAM file. Byte-identical
/// to Picard run with `USE_JDK_DEFLATER=true` (the BGZF blocks come from the port's zlib writer);
/// Picard's default GKL deflater is a separate surface.
pub fn sort_sam_to_bam(input_sam: &str, order: SortOrder) -> Result<Vec<u8>, ParseError> {
    let (mut header, mut records) = read_sam(input_sam)?;
    header.set_sort_order(order.name());
    records.sort_by(order.comparator());

    let mut writer =
        crate::output::bam_writer(Vec::new(), &header).expect("in-memory BAM writer never fails");
    for rec in &records {
        writer
            .write(rec)
            .expect("records that parsed re-encode as BAM");
    }
    Ok(writer
        .finish()
        .expect("in-memory BAM writer never fails to finish"))
}

/// `SortingCollection`'s default `maxRecordsInRam`, which is `MAX_RECORDS_IN_RAM`'s.
pub const DEFAULT_MAX_RECORDS_IN_RAM: usize = 500_000;

/// How many runs are merged at once. htsjdk opens every run together and fails on a file-descriptor
/// limit a whole genome reaches; merging the runs so far into one run first writes the same records
/// in the same order, since they are all earlier than the runs still to come.
const MAX_OPEN_RUNS: usize = 256;

type Compare = fn(&BamRecord, &BamRecord) -> Ordering;

/// htsjdk's `SortingCollection` of records, spilling to disk as `BAMRecordCodec` encodes them.
pub struct SortingCollection {
    compare: Compare,
    max_records_in_ram: usize,
    tmp_dir: PathBuf,
    ram: Vec<BamRecord>,
    runs: Vec<File>,
}

impl SortingCollection {
    pub fn new(compare: Compare, max_records_in_ram: usize, tmp_dir: PathBuf) -> Self {
        SortingCollection {
            compare,
            max_records_in_ram: max_records_in_ram.max(1),
            tmp_dir,
            ram: Vec::new(),
            runs: Vec::new(),
        }
    }

    /// `add`: spill first when the buffer is already full, as htsjdk does.
    pub fn add(&mut self, record: BamRecord) -> io::Result<()> {
        if self.ram.len() == self.max_records_in_ram {
            self.spill()?;
        }
        self.ram.push(record);
        Ok(())
    }

    /// How many runs have gone to disk.
    pub fn spills(&self) -> usize {
        self.runs.len()
    }

    fn spill(&mut self) -> io::Result<()> {
        // Stable, so equal records keep input order within a run (decision 0021).
        self.ram.sort_by(self.compare);
        let run = write_run(&self.tmp_dir, self.ram.drain(..).map(Ok))?;
        self.runs.push(run);
        if self.runs.len() == MAX_OPEN_RUNS {
            let merged = Merge::new(self.compare, std::mem::take(&mut self.runs))?;
            let run = write_run(&self.tmp_dir, merged)?;
            self.runs.push(run);
        }
        Ok(())
    }

    /// `doneAdding` and `iterator`: the records in order, from RAM when nothing spilled.
    pub fn into_sorted(mut self) -> io::Result<Sorted> {
        if self.runs.is_empty() {
            self.ram.sort_by(self.compare);
            return Ok(Sorted::Ram(self.ram.into_iter()));
        }
        if !self.ram.is_empty() {
            self.spill()?;
        }
        Ok(Sorted::Merge(Merge::new(self.compare, self.runs)?))
    }
}

/// Writes records to a fresh file in `dir`, unlinked at once so nothing is left behind, and
/// rewound for reading.
fn write_run(dir: &Path, records: impl Iterator<Item = io::Result<BamRecord>>) -> io::Result<File> {
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = dir.join(format!("sortingcollection.{}.{n}.tmp", std::process::id()));
    let file = File::options()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)?;
    // shortcut: unlinking an open file keeps it readable on Unix only; Windows is not a target.
    std::fs::remove_file(&path)?;
    let mut out = BufWriter::new(file);
    let mut bytes = Vec::new();
    for record in records {
        bytes.clear();
        record?
            .encode_into(&mut bytes)
            .map_err(|e| io::Error::other(format!("{e:?}")))?;
        out.write_all(&bytes)?;
    }
    let mut file = out.into_inner().map_err(|e| e.into_error())?;
    file.seek(SeekFrom::Start(0))?;
    Ok(file)
}

/// One spilled run, read back record by record.
struct Run {
    reader: BufReader<File>,
    bytes: Vec<u8>,
}

impl Iterator for Run {
    type Item = io::Result<BamRecord>;

    fn next(&mut self) -> Option<io::Result<BamRecord>> {
        match self.reader.fill_buf() {
            Ok([]) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }
        let read = (|| {
            let mut size = [0u8; 4];
            self.reader.read_exact(&mut size)?;
            self.bytes.clear();
            self.bytes.extend_from_slice(&size);
            let size = i32::from_le_bytes(size) as usize;
            self.bytes.resize(4 + size, 0);
            self.reader.read_exact(&mut self.bytes[4..])?;
            match BamRecord::decode(&self.bytes) {
                Ok(Some((record, _))) => Ok(record),
                other => Err(io::Error::other(format!("spilled record: {other:?}"))),
            }
        })();
        Some(read)
    }
}

/// The next record of one run, ordered so the heap's top is the least record, and of equal
/// records the one from the earliest run: `MergingIterator`'s tie-break.
struct Head {
    record: BamRecord,
    run: usize,
    compare: Compare,
}

impl Ord for Head {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.compare)(&other.record, &self.record).then(other.run.cmp(&self.run))
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
pub struct Merge {
    compare: Compare,
    runs: Vec<Run>,
    heap: BinaryHeap<Head>,
}

impl Merge {
    fn new(compare: Compare, files: Vec<File>) -> io::Result<Merge> {
        let mut merge = Merge {
            compare,
            runs: files
                .into_iter()
                .map(|file| Run {
                    reader: BufReader::new(file),
                    bytes: Vec::new(),
                })
                .collect(),
            heap: BinaryHeap::new(),
        };
        for run in 0..merge.runs.len() {
            merge.refill(run)?;
        }
        Ok(merge)
    }

    fn refill(&mut self, run: usize) -> io::Result<()> {
        if let Some(record) = self.runs[run].next() {
            self.heap.push(Head {
                record: record?,
                run,
                compare: self.compare,
            });
        }
        Ok(())
    }
}

impl Iterator for Merge {
    type Item = io::Result<BamRecord>;

    fn next(&mut self) -> Option<io::Result<BamRecord>> {
        let head = self.heap.pop()?;
        Some(self.refill(head.run).map(|()| head.record))
    }
}

/// The records of a [`SortingCollection`], in order.
pub enum Sorted {
    Ram(std::vec::IntoIter<BamRecord>),
    Merge(Merge),
}

impl Iterator for Sorted {
    type Item = io::Result<BamRecord>;

    fn next(&mut self) -> Option<io::Result<BamRecord>> {
        match self {
            Sorted::Ram(records) => records.next().map(Ok),
            Sorted::Merge(merge) => merge.next(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sorts through a collection holding `max` records at most.
    fn sorted_names(sam: &str, order: SortOrder, max: usize) -> (Vec<Vec<u8>>, usize) {
        let (_, records) = read_sam(sam).unwrap();
        let mut sorting = SortingCollection::new(order.comparator(), max, std::env::temp_dir());
        for record in records {
            sorting.add(record).unwrap();
        }
        let spills = sorting.spills();
        let out = sorting
            .into_sorted()
            .unwrap()
            .map(|r| r.unwrap().encode().unwrap())
            .collect();
        (out, spills)
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
            for max in [2, 7, 299, 300] {
                let (spilled, spills) = sorted_names(&sam, order, max);
                assert_eq!(spills, 300usize.div_ceil(max) - 1, "max {max}");
                assert_eq!(spilled, whole, "{order:?} with MAX_RECORDS_IN_RAM={max}");
            }
        }
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
        let (whole, _) = sorted_names(&sam, SortOrder::Coordinate, DEFAULT_MAX_RECORDS_IN_RAM);
        let (spilled, _) = sorted_names(&sam, SortOrder::Coordinate, 1);
        assert_eq!(spilled, whole);
    }

    const INPUT: &str = "@HD\tVN:1.6\tSO:unsorted\n\
        @SQ\tSN:chr1\tLN:1000\n\
        @SQ\tSN:chr2\tLN:1000\n\
        r3\t0\tchr1\t500\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        r1\t0\tchr2\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        r2\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    fn body(sam: &str) -> Vec<String> {
        sam.lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    }

    #[test]
    fn coordinate_sorts_by_reference_then_position() {
        let out = sort_sam(INPUT, SortOrder::Coordinate).unwrap();
        assert!(out.contains("@HD\tVN:1.6\tSO:coordinate"));
        // chr1:100 (r2), chr1:500 (r3), then chr2:100 (r1).
        assert_eq!(body(&out), ["r2", "r3", "r1"]);
    }

    #[test]
    fn queryname_sorts_by_name() {
        let out = sort_sam(INPUT, SortOrder::Queryname).unwrap();
        assert!(out.contains("@HD\tVN:1.6\tSO:queryname"));
        assert_eq!(body(&out), ["r1", "r2", "r3"]);
    }

    #[test]
    fn only_the_sort_order_field_of_the_header_changes() {
        let out = sort_sam(INPUT, SortOrder::Coordinate).unwrap();
        // The sequence dictionary is preserved verbatim and in order.
        assert!(out.contains("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000"));
    }
}
