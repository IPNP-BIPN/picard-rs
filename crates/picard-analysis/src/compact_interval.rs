//! Compact intervals for the interval-list tools, read line by line.
//!
//! [`htsjdk_bam::interval::Interval`] carries its contig as a `String`, and the tools cloned whole
//! lists between steps, so a million-line list peaked at ten times its size on disk. Here a contig
//! is an index into one shared [`Contigs`] table, a list is read a line at a time, and the
//! orderings, merges and renderings are `htsjdk_bam::interval`'s, ported onto that shape.
//!
//! Every tool used to read its inputs whole before parsing them, so a read failure (an unreadable
//! or non-UTF-8 file) was reported before any refusal of a line. [`read_lines`] keeps that order: a
//! refusal stops the parsing, the rest of the file is still read, and a read failure there wins.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};

use htsjdk_bam::interval::{Interval, ParseError};

/// The contig names a run has seen, each once, with its sort rank: the index of its first `@SQ` in
/// the dictionary, or `-1` for a contig the dictionary does not have (`getSequenceIndex`).
#[derive(Debug, Default)]
pub struct Contigs {
    names: Vec<Box<str>>,
    rank: Vec<i32>,
    index: HashMap<Box<str>, u32>,
}

impl Contigs {
    /// A table whose ranks follow `dictionary`.
    pub fn with_dictionary(dictionary: &[String]) -> Self {
        let mut contigs = Contigs::default();
        contigs.set_dictionary(dictionary);
        contigs
    }

    /// Ranks every contig by its first position in `dictionary`.
    pub fn set_dictionary(&mut self, dictionary: &[String]) {
        for (position, name) in dictionary.iter().enumerate() {
            let index = self.intern(name) as usize;
            if self.rank[index] == -1 {
                self.rank[index] = position as i32;
            }
        }
    }

    /// The index of `name`, added (unranked) when it is new.
    pub fn intern(&mut self, name: &str) -> u32 {
        if let Some(&index) = self.index.get(name) {
            return index;
        }
        let index = self.names.len() as u32;
        self.names.push(name.into());
        self.rank.push(-1);
        self.index.insert(name.into(), index);
        index
    }

    pub fn find(&self, name: &str) -> Option<u32> {
        self.index.get(name).copied()
    }

    pub fn name(&self, index: u32) -> &str {
        &self.names[index as usize]
    }

    pub fn rank(&self, index: u32) -> i32 {
        self.rank[index as usize]
    }
}

/// An `Interval` whose contig is an index into a [`Contigs`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Iv {
    pub contig: u32,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
    /// `null` in Java, written as `.`.
    pub name: Option<Box<str>>,
}

impl Iv {
    pub fn from_interval(interval: Interval, contigs: &mut Contigs) -> Iv {
        Iv {
            contig: contigs.intern(&interval.contig),
            start: interval.start,
            end: interval.end,
            negative_strand: interval.negative_strand,
            name: interval.name.map(String::into_boxed_str),
        }
    }

    pub fn to_interval(&self, contigs: &Contigs) -> Interval {
        Interval::with_strand_and_name(
            contigs.name(self.contig),
            self.start,
            self.end,
            self.negative_strand,
            self.name.as_deref(),
        )
    }

    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }

    /// `intersects`: same contig and overlapping coordinates.
    pub fn intersects(&self, other: &Iv) -> bool {
        self.contig == other.contig && self.start <= other.end && self.end >= other.start
    }

    /// `withinDistanceOf(other, distance)`.
    pub fn within_distance_of(&self, other: &Iv, distance: i32) -> bool {
        self.contig == other.contig
            && self.start - distance <= other.end
            && self.end + distance >= other.start
    }

    /// `intersect(that)` on two intervals known to intersect: this strand, the overlap, and the
    /// name `"{this} intersection {that}"` with a null name read as `null`.
    pub fn intersect(&self, that: &Iv) -> Iv {
        let name = format!(
            "{} intersection {}",
            self.name.as_deref().unwrap_or("null"),
            that.name.as_deref().unwrap_or("null")
        );
        Iv {
            contig: self.contig,
            start: self.start.max(that.start),
            end: self.end.min(that.end),
            negative_strand: self.negative_strand,
            name: Some(name.into_boxed_str()),
        }
    }

    pub fn strand(&self) -> char {
        if self.negative_strand {
            '-'
        } else {
            '+'
        }
    }

    /// One line as `IntervalListWriter.write` emits it.
    pub fn write_line(&self, contigs: &Contigs, out: &mut impl Write) -> io::Result<()> {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            contigs.name(self.contig),
            self.start,
            self.end,
            self.strand(),
            self.name.as_deref().unwrap_or(".")
        )
    }
}

/// `IntervalCoordinateComparator`: dictionary rank, start, end, positive strand first, then name
/// with nulls first.
pub fn compare_coordinate(contigs: &Contigs, lhs: &Iv, rhs: &Iv) -> Ordering {
    contigs
        .rank(lhs.contig)
        .cmp(&contigs.rank(rhs.contig))
        .then_with(|| lhs.start.cmp(&rhs.start))
        .then_with(|| lhs.end.cmp(&rhs.end))
        .then_with(|| lhs.negative_strand.cmp(&rhs.negative_strand))
        .then_with(|| lhs.name.cmp(&rhs.name))
}

/// `sorted()`, in place and stable.
pub fn sort(contigs: &Contigs, intervals: &mut [Iv]) {
    intervals.sort_by(|a, b| compare_coordinate(contigs, a, b));
}

/// `getUniqueIntervals(mergeAbutting, concatenateNames = true)`: sort, then fold each overlapping
/// (or, when `merge_abutting`, abutting) run into its first interval, names joined with `|`.
pub fn uniqued(contigs: &Contigs, mut intervals: Vec<Iv>, merge_abutting: bool) -> Vec<Iv> {
    sort(contigs, &mut intervals);
    // The names of the run being folded, once a second name has joined the first.
    let mut joined: Option<String> = None;
    let mut kept = 0;
    for next in 0..intervals.len() {
        if kept > 0 && {
            let (cur, next) = (&intervals[kept - 1], &intervals[next]);
            (merge_abutting && cur.within_distance_of(next, 1)) || cur.intersects(next)
        } {
            let name = intervals[next].name.take();
            let end = intervals[next].end;
            let cur = &mut intervals[kept - 1];
            cur.end = cur.end.max(end);
            if let Some(name) = name {
                match &mut joined {
                    Some(joined) => {
                        joined.push('|');
                        joined.push_str(&name);
                    }
                    None => {
                        joined = Some(match cur.name.take() {
                            Some(first) => format!("{first}|{name}"),
                            None => name.into(),
                        })
                    }
                }
            }
        } else {
            if let Some(joined) = joined.take() {
                intervals[kept - 1].name = Some(joined.into_boxed_str());
            }
            intervals.swap(kept, next);
            kept += 1;
        }
    }
    if let Some(joined) = joined {
        intervals[kept - 1].name = Some(joined.into_boxed_str());
    }
    intervals.truncate(kept);
    intervals
}

/// `IntervalList.invert`: every base of every `sequences` contig, in dictionary order, that the
/// list does not cover, named `interval-1`, `interval-2`, ... across contigs.
pub fn invert(contigs: &mut Contigs, sequences: &[(String, i32)], intervals: Vec<Iv>) -> Vec<Iv> {
    let unique = uniqued(contigs, intervals, true);
    let mut out = Vec::new();
    let mut counter = 0;
    for (name, length) in sequences {
        let contig = contigs.intern(name);
        let mut last_covered = 0;
        let mut gap = |start: i32, end: i32, out: &mut Vec<Iv>| {
            counter += 1;
            out.push(Iv {
                contig,
                start,
                end,
                negative_strand: false,
                name: Some(format!("interval-{counter}").into_boxed_str()),
            });
        };
        for interval in unique.iter().filter(|iv| iv.contig == contig) {
            if interval.start > last_covered + 1 {
                gap(last_covered + 1, interval.start - 1, &mut out);
            }
            last_covered = interval.end;
        }
        if *length > last_covered {
            gap(last_covered + 1, *length, &mut out);
        }
    }
    out
}

/// `IntervalList.padded(padding, padding)`, clamped to `[1, contig length]`; an interval left empty
/// (or on a contig `sequences` lacks) is dropped.
pub fn padded(
    contigs: &Contigs,
    sequences: &[(String, i32)],
    mut intervals: Vec<Iv>,
    padding: i32,
) -> Vec<Iv> {
    if padding == 0 {
        return intervals;
    }
    intervals.retain_mut(|iv| {
        let name = contigs.name(iv.contig);
        let length = sequences
            .iter()
            .find(|(n, _)| n == name)
            .map_or(0, |(_, l)| *l);
        iv.start = (iv.start - padding).max(1);
        iv.end = (iv.end + padding).min(length);
        iv.start <= iv.end
    });
    intervals
}

/// `IntervalList.breakIntervalsAtBandMultiples`: cut each interval at every multiple of `band`,
/// naming the pieces `"{name}.{k}"`.
pub fn break_at_band_multiples(intervals: Vec<Iv>, band: i32) -> Vec<Iv> {
    let mut out = Vec::with_capacity(intervals.len());
    for interval in intervals {
        let (first, last) = (interval.start / band, interval.end / band);
        if interval.end < interval.start || first == last {
            out.push(interval);
            continue;
        }
        let base = interval.name.as_deref().unwrap_or("null");
        let mut start = interval.start;
        for index in first..=last {
            out.push(Iv {
                start,
                end: ((index + 1) * band - 1).min(interval.end),
                name: Some(format!("{base}.{}", index - first + 1).into_boxed_str()),
                ..interval
            });
            start = (index + 1) * band;
        }
    }
    out
}

/// The `@` lines of an interval list, and how many of them came before its first other line.
#[derive(Debug, Default)]
pub struct Header {
    pub lines: Vec<String>,
    pub leading: usize,
}

impl Header {
    /// Every `@` line, each ending in a newline.
    pub fn text(&self) -> String {
        self.lines.iter().map(|line| format!("{line}\n")).collect()
    }

    /// The `@` lines before the first other line.
    pub fn leading_lines(&self) -> &[String] {
        &self.lines[..self.leading]
    }
}

/// Reads `reader` a line at a time (without its `\n` or `\r\n`, as `str::lines` splits). Once
/// `each` refuses a line the rest is still read, so that a read failure later in the file is the
/// one reported, as a whole-file read would have reported it.
pub fn read_lines<E>(
    mut reader: impl BufRead,
    mut each: impl FnMut(&str) -> Result<(), E>,
) -> io::Result<Result<(), E>> {
    let mut line = String::new();
    let mut refusal = None;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Ok(refusal.map_or(Ok(()), Err));
        }
        if refusal.is_none() {
            let text = match line.strip_suffix('\n') {
                Some(text) => text.strip_suffix('\r').unwrap_or(text),
                None => &line,
            };
            if let Err(error) = each(text) {
                refusal = Some(error);
            }
        }
    }
}

/// Reads all of `reader`, for its read failure alone.
pub fn drain(reader: impl BufRead) -> io::Result<()> {
    read_lines(reader, |_| Ok::<(), ()>(())).map(|_| ())
}

/// One interval line as `IntervalList::parse_body` reads it, or `None` for a header or blank line.
pub fn parse_line(line: &str, contigs: &mut Contigs) -> Result<Option<Iv>, ParseError> {
    if line.starts_with('@') || line.trim().is_empty() {
        return Ok(None);
    }
    let mut fields = line.split('\t');
    let mut field = || fields.next();
    let (Some(contig), Some(start), Some(end), Some(strand), Some(name)) =
        (field(), field(), field(), field(), field())
    else {
        return Err(ParseError::TooFewFields(line.to_string()));
    };
    let number = |s: &str| s.parse().map_err(|_| ParseError::BadNumber(s.into()));
    Ok(Some(Iv {
        start: number(start)?,
        end: number(end)?,
        contig: contigs.intern(contig),
        negative_strand: strand == "-",
        name: (name != ".").then(|| name.into()),
    }))
}

/// Reads an interval list: its `@` lines into a [`Header`], its intervals onto `intervals`.
pub fn read_list(
    reader: impl BufRead,
    contigs: &mut Contigs,
    intervals: &mut Vec<Iv>,
) -> io::Result<Result<Header, ParseError>> {
    let mut header = Header::default();
    let mut body = false;
    let read = read_lines(reader, |line| {
        if line.starts_with('@') {
            header.lines.push(line.to_string());
            if !body {
                header.leading += 1;
            }
            return Ok(());
        }
        body = true;
        if let Some(iv) = parse_line(line, contigs)? {
            intervals.push(iv);
        }
        Ok(())
    })?;
    Ok(read.map(|()| header))
}

/// Writes every interval as an interval-list line.
pub fn write_lines(contigs: &Contigs, intervals: &[Iv], out: &mut impl Write) -> io::Result<()> {
    for interval in intervals {
        interval.write_line(contigs, out)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::interval::IntervalList;

    const LIST: &str = "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:100\n@SQ\tSN:chr1\tLN:100\n\
        chr1\t5\t9\t+\tb\nchrU\t1\t2\t-\t.\nchr2\t40\t50\t-\tc\nchr1\t1\t4\t+\ta\n\
        chr1\t3\t20\t-\td\nchr2\t51\t60\t+\t.\nchr1\t30\t30\t+\te\n";

    fn ours(text: &str) -> (Contigs, Vec<Iv>) {
        let mut contigs = Contigs::default();
        let mut intervals = Vec::new();
        read_list(text.as_bytes(), &mut contigs, &mut intervals)
            .unwrap()
            .unwrap();
        contigs.set_dictionary(&["chr2".to_string(), "chr1".to_string()]);
        (contigs, intervals)
    }

    fn lines(contigs: &Contigs, intervals: &[Iv]) -> String {
        let mut out = Vec::new();
        write_lines(contigs, intervals, &mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// The compact sort, merge, invert and band split write what `htsjdk_bam::interval` writes.
    #[test]
    fn agrees_with_htsjdk_interval() {
        let dictionary = vec!["chr2".to_string(), "chr1".to_string()];
        let theirs = IntervalList::parse_body(dictionary.clone(), LIST).unwrap();
        let (mut contigs, intervals) = ours(LIST);
        let sorted = {
            let mut v = intervals.clone();
            sort(&contigs, &mut v);
            v
        };
        assert_eq!(lines(&contigs, &sorted), theirs.sorted().write_body());
        assert_eq!(
            lines(&contigs, &uniqued(&contigs, intervals.clone(), true)),
            theirs.uniqued(true).write_body()
        );
        let sequences = vec![("chr2".to_string(), 100), ("chr1".to_string(), 100)];
        let inverted = invert(&mut contigs, &sequences, intervals.clone());
        assert_eq!(
            lines(&contigs, &inverted),
            theirs.invert(&sequences).write_body()
        );
        let broken = htsjdk_bam::interval::break_intervals_at_band_multiples(&theirs.intervals, 7);
        let broken = IntervalList {
            dictionary,
            intervals: broken,
        };
        assert_eq!(
            lines(&contigs, &break_at_band_multiples(intervals, 7)),
            broken.write_body()
        );
    }

    #[test]
    fn a_refusal_gives_way_to_a_later_read_failure() {
        let mut text = b"chr1\t1\n".to_vec();
        text.extend_from_slice(b"chr1\t1\t2\t+\t\xff\n");
        let mut contigs = Contigs::default();
        let read = read_list(&text[..], &mut contigs, &mut Vec::new());
        assert_eq!(read.unwrap_err().kind(), io::ErrorKind::InvalidData);
        let read = read_list(&b"chr1\t1\n"[..], &mut contigs, &mut Vec::new()).unwrap();
        assert!(matches!(read, Err(ParseError::TooFewFields(_))));
    }
}
