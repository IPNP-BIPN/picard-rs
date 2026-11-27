//! An indexed BAM read over a list of regions, as `BAMFileReader.createIndexIterator` reads one:
//! the chunks the `.bai` names for every region (`getSpanOverlapping`, then `BAMFileSpan.merge`),
//! each read from its start to its end. Only those parts of the file are inflated; the caller
//! keeps its overlap filter, so it is handed the records a scan of the whole file would have kept,
//! in file order, each once.

use crate::record_validation::Validator;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::index::{read_bai, BamIndex, Chunk};
use htsjdk_bam::record::BamRecord;
use noodles_bgzf::VirtualPosition;
use std::fs::File;
use std::io::{BufReader, Read};

/// `GenomicIndexUtil.regionToBins` for a 1-based inclusive region: the bins that may hold a record
/// overlapping it.
fn region_bins(start: i32, end: i32) -> Vec<i32> {
    let max = 0x1FFF_FFFF;
    let start = if start <= 0 { 0 } else { (start - 1) & max };
    let end = if end <= 0 { max } else { (end - 1) & max };
    if start > end {
        return Vec::new();
    }
    let mut bins = vec![0];
    for (offset, shift) in [(1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)] {
        bins.extend(offset + (start >> shift)..=offset + (end >> shift));
    }
    bins
}

/// The chunks to read for `intervals` (`(sequence, start, end)`, 1-based inclusive): per interval,
/// its bins' chunks past the linear index's minimum offset; then all of them sorted and those that
/// overlap merged, so no record is read twice. htsjdk also coalesces chunks in adjacent blocks,
/// which reads more records for its filter to drop and changes none it keeps.
pub fn chunks(index: &BamIndex, intervals: &[(i32, i32, i32)]) -> Vec<Chunk> {
    let mut all = Vec::new();
    for &(sequence, start, end) in intervals {
        let Some(content) = usize::try_from(sequence)
            .ok()
            .and_then(|s| index.references.get(s))
        else {
            continue;
        };
        let bins = region_bins(start, end);
        // `LinearIndex.getMinimumOffset`.
        let window = (if start <= 0 { 0 } else { start - 1 } >> 14) as usize;
        let minimum = content.linear_index.get(window).copied().unwrap_or(0);
        for bin in &content.bins {
            if bins.contains(&bin.bin_number) {
                all.extend(bin.chunks.iter().filter(|c| c.end > minimum).copied());
            }
        }
    }
    all.sort_by_key(|c| (c.start, c.end));
    let mut merged: Vec<Chunk> = Vec::new();
    for chunk in all {
        match merged.last_mut() {
            Some(last) if chunk.start <= last.end => last.end = last.end.max(chunk.end),
            _ => merged.push(chunk),
        }
    }
    merged
}

/// `SamFiles.findIndex` for a BAM's `.bai`: `x.bai` beside `x.bam`, else `x.bam.bai`.
fn index_path(input: &str) -> Option<String> {
    let candidates = [
        input.strip_suffix(".bam").map(|stem| format!("{stem}.bai")),
        Some(format!("{input}.bai")),
    ];
    candidates
        .into_iter()
        .flatten()
        .find(|p| std::path::Path::new(p).is_file())
}

/// The records of the chunks, validated as `BAMFileIndexIterator` validates every record it
/// decodes, whether or not the caller's filter keeps it.
pub struct Query {
    reader: noodles_bgzf::io::Reader<BufReader<File>>,
    chunks: std::vec::IntoIter<Chunk>,
    end: u64,
    validator: Validator,
    bytes: Vec<u8>,
}

impl Query {
    /// The query over `intervals` of the BAM at `input`, whose header was already read; `None`
    /// when there is no `.bai` to read, for the caller to scan the file instead.
    pub fn open(input: &str, header: &SamHeader, intervals: &[(i32, i32, i32)]) -> Option<Query> {
        let index = read_bai(&std::fs::read(index_path(input)?).ok()?).ok()?;
        // The binary dictionary, which the checks fall back on when the text has no `@SQ`.
        let mut stream = crate::bam_stream::RecordStream::new(noodles_bgzf::io::Reader::new(
            File::open(input).ok()?,
        ));
        stream.header().ok()?;
        Some(Query {
            reader: noodles_bgzf::io::Reader::new(BufReader::new(File::open(input).ok()?)),
            chunks: chunks(&index, intervals).into_iter(),
            end: 0,
            validator: Validator::new(header, &stream.binary_sequences),
            bytes: Vec::new(),
        })
    }

    fn next_record(&mut self) -> std::io::Result<Option<BamRecord>> {
        while u64::from(self.reader.virtual_position()) >= self.end {
            let Some(chunk) = self.chunks.next() else {
                return Ok(None);
            };
            self.reader.seek(VirtualPosition::from(chunk.start))?;
            self.end = chunk.end;
        }
        let mut size = [0u8; 4];
        self.reader.read_exact(&mut size)?;
        let length = u32::from_le_bytes(size) as usize;
        self.bytes.clear();
        self.bytes.extend_from_slice(&size);
        self.bytes.resize(4 + length, 0);
        self.reader.read_exact(&mut self.bytes[4..])?;
        match BamRecord::decode(&self.bytes) {
            Ok(Some((record, _))) => Ok(Some(record)),
            other => Err(std::io::Error::other(format!("{other:?}"))),
        }
    }
}

impl Iterator for Query {
    type Item = BamRecord;

    fn next(&mut self) -> Option<BamRecord> {
        let record = self
            .next_record()
            .unwrap_or_else(|e| crate::metrics_cli::fail(&format!("{e}")))?;
        if let Err(message) = self.validator.check(&record) {
            crate::metrics_cli::thrown(&message);
        }
        Some(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_bins_are_htsjdks() {
        // One 16 kb window: a bin on every level.
        assert_eq!(region_bins(1, 100), vec![0, 1, 9, 73, 585, 4681]);
        // Across the first 16 kb boundary, two leaf bins.
        assert_eq!(
            region_bins(16_000, 17_000),
            vec![0, 1, 9, 73, 585, 4681, 4682]
        );
        assert_eq!(region_bins(10, 5), Vec::<i32>::new());
    }
}
