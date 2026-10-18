//! What htsjdk's writers check and leave on disk, beyond the records themselves.
//!
//! A tool that throws halfway through leaves its output file as the writer's buffers left it: the
//! writer is never closed, so nothing it still holds is written. Measured against Picard 3.5.0 on a
//! 50,000-record CleanSam run that throws at record 30,000 (decision recorded in the suite
//! `inputshapes`):
//!
//! * a SAM is its header, which `SAMTextWriter` flushes as soon as it is written, then the record
//!   text in whole chunks of the `AsciiWriter`'s 128 KiB buffer;
//! * a BAM is the BGZF blocks that filled before the throw, each written as three calls (header,
//!   deflated payload, footer) through a 128 KiB `BufferedOutputStream`, so the file ends wherever
//!   that buffer last flushed, mid-block more often than not.

use std::cmp::Ordering;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::{coordinate, query_name};
use htsjdk_bgzf::{BLOCK_FOOTER_LENGTH, BLOCK_HEADER_LENGTH, DEFAULT_UNCOMPRESSED_BLOCK_SIZE};

/// `Defaults.BUFFER_SIZE`, the size of both buffers.
pub const BUFFER_SIZE: usize = 128 * 1024;

/// The bytes of a SAM the run died writing: `header_len` bytes of header, then `text`'s records in
/// whole buffers. `text` is the whole SAM the writer had been handed, header included.
pub fn sam_on_disk(text: &[u8], header_len: usize) -> &[u8] {
    let records = text.len().saturating_sub(header_len);
    &text[..header_len + records / BUFFER_SIZE * BUFFER_SIZE]
}

/// The bytes of a BAM the run died writing, from the finished stream of everything the writer had
/// been handed: only the blocks that were full, through the buffered stream.
pub fn bam_on_disk(finished: &[u8]) -> Vec<u8> {
    // The full blocks, in order: every block but the last data block and the terminator holds
    // exactly one uncompressed block's worth.
    let mut writes = Vec::new();
    let mut p = 0;
    while p + BLOCK_HEADER_LENGTH <= finished.len() {
        let size = u16::from_le_bytes([finished[p + 16], finished[p + 17]]) as usize + 1;
        let Some(block) = finished.get(p..p + size) else {
            break;
        };
        let isize = u32::from_le_bytes(block[size - 4..].try_into().expect("four bytes")) as usize;
        if isize != DEFAULT_UNCOMPRESSED_BLOCK_SIZE {
            break;
        }
        writes.push(BLOCK_HEADER_LENGTH);
        writes.push(size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);
        writes.push(BLOCK_FOOTER_LENGTH);
        p += size;
    }
    // `BufferedOutputStream.write`: a write that does not fit flushes what is buffered first, and
    // one at least as large as the buffer goes straight through.
    let (mut on_disk, mut buffered) = (0, 0);
    for len in writes {
        if len >= BUFFER_SIZE {
            on_disk += buffered + len;
            buffered = 0;
        } else if buffered + len > BUFFER_SIZE {
            on_disk += buffered;
            buffered = len;
        } else {
            buffered += len;
        }
    }
    finished[..on_disk].to_vec()
}

/// `SAMFileWriterImpl.assertPresorted` over records in the order they are added: the index of the
/// first one `SAMSortOrderChecker.isSorted` rejects, with the exception's message.
///
/// Only `coordinate` and `queryname` have a checker here; `unsorted` and `unknown` have none in
/// htsjdk either, and `duplicate` is not checked (no tool reaches a presorted duplicate writer).
/// `filename` is `getFilename()`: `null` for a SAM, the output's `file:` URI for a BAM.
pub fn first_out_of_order(
    header: &SamHeader,
    records: &[BamRecord],
    filename: &str,
) -> Option<(usize, String)> {
    let order = header.attributes.get("SO").unwrap_or("unsorted");
    let compare: fn(&BamRecord, &BamRecord) -> Ordering = match order {
        "coordinate" => coordinate::file_order_compare,
        "queryname" => |a, b| query_name::compare_read_names(&a.read_name, &b.read_name),
        _ => return None,
    };
    let index = records
        .windows(2)
        .position(|pair| compare(&pair[0], &pair[1]) == Ordering::Greater)?
        + 1;
    // Both keys are the OFFENDING record: `isSorted` sets `prev = rec` before it returns false.
    let key = sort_key(header, &records[index], order);
    Some((
        index,
        format!(
            "java.lang.IllegalArgumentException: Alignments added out of order in \
             SAMFileWriterImpl.addAlignment for {filename}. Sort order is {order}. Offending \
             records are at [{key}] and [{key}]"
        ),
    ))
}

/// `SAMSortOrderChecker.getSortKey`.
fn sort_key(header: &SamHeader, record: &BamRecord, order: &str) -> String {
    if order == "queryname" {
        return record.read_name.clone();
    }
    let name = usize::try_from(record.reference_index)
        .ok()
        .and_then(|index| header.sequences.get(index))
        .map_or("*", |sequence| sequence.name.as_str());
    format!("{name}:{}", record.alignment_start)
}

/// `getFilename()` of a BAM writer: the output's absolute `file:` URI, as `File.toURI` spells it.
pub fn bam_filename(output: &str) -> String {
    let absolute = std::path::absolute(output).unwrap_or_else(|_| output.into());
    let mut uri = String::from("file://");
    for c in absolute.to_string_lossy().chars() {
        match c {
            ' ' => uri.push_str("%20"),
            '%' => uri.push_str("%25"),
            other => uri.push(other),
        }
    }
    uri
}

/// `IOUtil.openGzipFileForWriting`: a `CustomGzipOutputStream`, which is a `GZIPOutputStream`
/// whose deflater takes `IOUtil.compressionLevel` (the tool's `COMPRESSION_LEVEL`). It is what
/// htsjdk writes for any `.gz` or `.bfq` opened with `openFileForWriting` (a BFQ, the barcode
/// files of ExtractIlluminaBarcodes), as opposed to the BGZF of its block-compressed writers. One
/// member, whose header names no file and no time and says OS 255, which is flate2's own header;
/// the deflate stream is zlib's, the same zlib the BGZF writer's JDK deflater runs.
pub fn gzip_writer<W: std::io::Write>(inner: W, level: u32) -> flate2::write::GzEncoder<W> {
    flate2::GzBuilder::new().write(inner, flate2::Compression::new(level))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_sam_keeps_its_header_and_whole_buffers_of_records() {
        let text = vec![b'x'; 10 + 3 * BUFFER_SIZE + 7];
        assert_eq!(sam_on_disk(&text, 10).len(), 10 + 3 * BUFFER_SIZE);
        assert_eq!(sam_on_disk(&text[..10 + 5], 10).len(), 10);
    }

    #[test]
    fn a_small_bam_leaves_nothing() {
        let mut writer = htsjdk_bgzf::BgzfWriter::new(Vec::new());
        std::io::Write::write_all(&mut writer, b"BAM\x01 a header and a record").unwrap();
        assert!(bam_on_disk(&writer.into_inner().unwrap()).is_empty());
    }
}
