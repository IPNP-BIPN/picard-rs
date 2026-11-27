//! What htsjdk's writers check and leave on disk, beyond the records themselves.
//!
//! A tool that throws halfway through leaves its output file as the writer's buffers left it: the
//! writer is never closed, so nothing it still holds is written. Measured against Picard 3.5.0 on a
//! 50,000-record CleanSam run that throws at record 30,000 (decision recorded in the suite
//! `inputshapes`):
//!
//! * a SAM is its header, which `SAMTextWriter` flushes as soon as it is written, then the record
//!   text in whole chunks of the `AsciiWriter`'s 128 KiB buffer;
//! * a BAM is the BGZF blocks that filled before the throw, written by
//!   `BlockCompressedOutputStream.writeGzipBlock` through a 128 KiB `BufferedOutputStream`, so the
//!   file ends wherever that buffer last flushed, mid-block more often than not. htsjdk 5.0.0
//!   writes each block as twelve small writes for the 18-byte header, one for the deflated
//!   payload and two for the 8-byte footer (`BlockCompressedOutputStream.java:429-448`);
//!   [`bam_on_disk`] models the header and the footer as one write each.
//!
//! Neither the JDK's `BufferedOutputStream` nor any other JDK class is transcribed here (the JDK is
//! GPL2): its behaviour is named where it matters. The goldens that measure the result are
//! `tests/data/inputshapes.txt.gz`'s rows `output clean_throws_after_a_flush` (a SAM) and
//! `output_bytes revert_missing_mate_bam` (a BAM).
//!
//! Why this is computed rather than left to a Rust writer: a Rust writer's own buffering would
//! leave a file of its own length at the throw, not htsjdk's. `stream_output.rs` therefore writes
//! the whole output and then cuts it to the length these functions give.
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only view of a byte array (Java's `byte[]` that the callee may
//!   not change); `&text[..n]` is the view of its first `n` bytes, with no copy.
//! - `usize` is an unsigned index type; `/` on it truncates, as Java's `int` division does.
//! - `Option<T>` is "a value or nothing"; `?` after an `Option` returns `None` from the function
//!   at once when there is nothing.
//! - `let Some(x) = y else { ... };` unwraps `y` or runs the block, which must leave the loop or
//!   the function; Java has no single construct for it.
//! - `|a, b| ...` is a closure, Java's lambda; `fn(&A, &B) -> Ordering` is the type of a plain
//!   function pointer, which a closure that captures nothing can be.
//! - `W: std::io::Write` makes a function generic over any output stream, Java's `OutputStream`
//!   parameter.

use std::cmp::Ordering;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::{coordinate, query_name};
use htsjdk_bgzf::{BLOCK_FOOTER_LENGTH, BLOCK_HEADER_LENGTH, DEFAULT_UNCOMPRESSED_BLOCK_SIZE};

/// `Defaults.BUFFER_SIZE`, the size of both buffers.
///
/// `BUFFER_SIZE = getIntProperty("buffer_size", 1024 * 128);` (`Defaults.java:128`). The BAM's
/// `BufferedOutputStream` takes `SAMFileWriterFactory`'s `bufferSize`, which starts as this; the
/// `AsciiWriter` sizes its buffer from `Defaults.NON_ZERO_BUFFER_SIZE` (`AsciiWriter.java:38`),
/// which is the same number unless `buffer_size` is set to 0.
pub const BUFFER_SIZE: usize = 128 * 1024;

/// The bytes of a SAM the run died writing: `header_len` bytes of header, then `text`'s records in
/// whole buffers. `text` is the whole SAM the writer had been handed, header included.
///
/// The header reaches the stream because `SAMTextHeaderCodec.encode` ends with
/// `this.writer.flush();` (`SAMTextHeaderCodec.java:451`), which empties the `AsciiWriter`. After
/// that the `AsciiWriter` hands its buffer on only when it is full (`AsciiWriter.java:79-82`):
///
/// ```java
/// if (numBytes == buffer.length) {
///     os.write(buffer, 0, numBytes);
///     numBytes = 0;
/// }
/// ```
///
/// so the records reach the file in whole buffers counted from the end of the header, and the
/// partial buffer at the throw is lost. `saturating_sub` is a subtraction that stops at zero
/// instead of going below it.
pub fn sam_on_disk(text: &[u8], header_len: usize) -> &[u8] {
    let records = text.len().saturating_sub(header_len);
    &text[..header_len + records / BUFFER_SIZE * BUFFER_SIZE]
}

/// The bytes of a BAM the run died writing, from the finished stream of everything the writer had
/// been handed: only the blocks that were full, through the buffered stream.
///
/// Two steps. First the full BGZF blocks are found at the front of `finished` and turned into the
/// list of write lengths htsjdk would have made for them; the block still being filled, and the
/// terminator, were never written. Then that list is run through `BufferedOutputStream.write`'s
/// rule to find how much of it had left the buffer.
///
/// The model makes one 18-byte write per header and one 8-byte write per footer, where htsjdk
/// makes twelve and two (`BlockCompressedOutputStream.java:429-448`). The buffer flushes early
/// only for a write that does not fit, so the two can differ when the buffer has fewer free
/// bytes than a header or a footer needs at the moment one arrives: there the small writes can
/// fill it further before it flushes.
pub fn bam_on_disk(finished: &[u8]) -> Vec<u8> {
    // The full blocks, in order: every block but the last data block and the terminator holds
    // exactly one uncompressed block's worth.
    let mut sizes = Vec::new();
    let mut p = 0;
    while p + BLOCK_HEADER_LENGTH <= finished.len() {
        // Bytes 16 and 17 of a BGZF header are `BSIZE`, the block's total size minus one
        // (`codec.writeShort((short) (totalBlockSize - 1));`), little-endian.
        let size = u16::from_le_bytes([finished[p + 16], finished[p + 17]]) as usize + 1;
        // `get` returns `None` instead of panicking when the block runs past the end.
        let Some(block) = finished.get(p..p + size) else {
            break;
        };
        // The last four bytes are `ISIZE`, the block's uncompressed length
        // (`codec.writeInt(uncompressedSize);`). `try_into` turns the four-byte slice into a
        // fixed array, which cannot fail here; `expect` would panic with its message if it did.
        let isize = u32::from_le_bytes(block[size - 4..].try_into().expect("four bytes")) as usize;
        if isize != DEFAULT_UNCOMPRESSED_BLOCK_SIZE {
            break;
        }
        sizes.push(size);
        p += size;
    }
    finished[..through_the_buffer(&sizes)].to_vec()
}

/// [`bam_on_disk`] in place, for the finished BAM at `path`: the file cut to what the unclosed
/// writer leaves, reading only the blocks' headers and footers.
pub fn cut_bam_on_disk(path: &str) -> std::io::Result<()> {
    use std::os::unix::fs::FileExt;
    let file = std::fs::File::options().read(true).write(true).open(path)?;
    let length = file.metadata()?.len();
    let mut sizes = Vec::new();
    let (mut p, mut header, mut footer) = (0u64, [0u8; BLOCK_HEADER_LENGTH], [0u8; 4]);
    while p + BLOCK_HEADER_LENGTH as u64 <= length {
        file.read_exact_at(&mut header, p)?;
        let size = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
        if p + size as u64 > length {
            break;
        }
        file.read_exact_at(&mut footer, p + size as u64 - 4)?;
        if u32::from_le_bytes(footer) as usize != DEFAULT_UNCOMPRESSED_BLOCK_SIZE {
            break;
        }
        sizes.push(size);
        p += size as u64;
    }
    file.set_len(through_the_buffer(&sizes) as u64)
}

/// How much of blocks of these sizes, each written as its header, payload and footer, reaches the
/// file through the 128 KiB `BufferedOutputStream` before the writer is abandoned.
fn through_the_buffer(sizes: &[usize]) -> usize {
    let writes = sizes.iter().flat_map(|&size| {
        [
            BLOCK_HEADER_LENGTH,
            size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH,
            BLOCK_FOOTER_LENGTH,
        ]
    });
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
    on_disk
}

/// `SAMFileWriterImpl.assertPresorted` over records in the order they are added: the index of the
/// first one `SAMSortOrderChecker.isSorted` rejects, with the exception's message.
///
/// Only `coordinate` and `queryname` have a checker here; `unsorted` and `unknown` have none in
/// htsjdk either, and `duplicate` is not checked (no tool reaches a presorted duplicate writer).
/// `filename` is `getFilename()`: `null` for a SAM, the output's `file:` URI for a BAM.
///
/// htsjdk 5.0.0, `src/main/java/htsjdk/samtools/SAMSortOrderChecker.java:49-53`:
///
/// ```java
/// if (prev != null) {
///     ret = comparator.fileOrderCompare(prev, rec) <= 0;
/// }
/// prev = rec;
/// return ret;
/// ```
///
/// A record is rejected when it compares greater than the one before it, so the first pair of
/// neighbours in that relation names the record. For `queryname`, `fileOrderCompare` is
/// `return compareReadNames(samRecord1.getReadName(), samRecord2.getReadName());`
/// (`SAMRecordQueryNameComparator.java:81`), names only, which is what the closure compares.
pub fn first_out_of_order(
    header: &SamHeader,
    records: &[BamRecord],
    filename: &str,
) -> Option<(usize, String)> {
    let order = header.attributes.get("SO").unwrap_or("unsorted");
    // `match` picks the comparison for the order, or leaves with `None` (nothing to check).
    let compare: fn(&BamRecord, &BamRecord) -> Ordering = match order {
        "coordinate" => coordinate::file_order_compare,
        "queryname" => |a, b| query_name::compare_read_names(&a.read_name, &b.read_name),
        _ => return None,
    };
    // `windows(2)` walks every pair of neighbours, `position` gives the index of the first pair
    // out of order (or `None`, which `?` returns: every record was in order), and `+ 1` turns it
    // into the index of the pair's second record, the one `isSorted` rejects.
    let index = records
        .windows(2)
        .position(|pair| compare(&pair[0], &pair[1]) == Ordering::Greater)?
        + 1;
    // Both keys are the OFFENDING record: `isSorted` sets `prev = rec` before it returns false.
    // The golden shows it: `tests/data/inputshapes.txt.gz`, row `thrown clean_unsorted`, ends
    // `Offending records are at [chr1:20] and [chr1:20]`.
    let key = sort_key(header, &records[index], order);
    // A `\` at the end of a line inside a string literal drops the newline and the next line's
    // leading spaces, so the message is the one line `assertPresorted` builds
    // (`SAMFileWriterImpl.java:209-213`).
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
///
/// `SAMSortOrderChecker.java:70-73`:
///
/// ```java
/// case coordinate:
///     return rec.getReferenceName() + ":" + rec.getAlignmentStart();
/// case queryname:
///     return rec.getReadName();
/// ```
///
/// An unplaced record's reference index is -1, which `usize::try_from` refuses, giving `"*"`:
/// `getReferenceName()` returns `NO_ALIGNMENT_REFERENCE_NAME`, which is `"*"`
/// (`SAMRecord.java:131`).
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

/// `getFilename()` of a BAM writer: the output's absolute `file:` URI.
///
/// `SAMFileWriterFactory.makeBAMWriter` hands the writer `outputPath.toUri().toString()`
/// (`SAMFileWriterFactory.java:324`), and the JDK's `Path.toUri` spells an absolute Unix path as
/// `file://` followed by the path (so `file:///...`), percent-encoding the characters a URI path
/// may not hold. Only a space and `%` are encoded here.
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
///
/// `IOUtil.java:856` and `CustomGzipOutputStream.java:16-17`:
///
/// ```java
/// return new CustomGzipOutputStream(out, Defaults.BUFFER_SIZE, compressionLevel);
/// ```
///
/// ```java
/// super(outputStream, bufferSize);
/// this.def.setLevel(compressionLevel);
/// ```
///
/// The level is the one Picard's `CommandLineProgram.instanceMain` sets with
/// `IOUtil.setCompressionLevel(COMPRESSION_LEVEL);` (`CommandLineProgram.java:222`). OS 255 is
/// also what OpenJDK 17's `GZIPOutputStream` writes in that byte, Picard 3.5.0's tested JVM.
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

    #[test]
    fn a_bam_cut_in_place_is_the_one_bam_on_disk_gives() {
        let path = std::env::temp_dir().join(format!("bam-on-disk-{}.bam", std::process::id()));
        let path = path.to_str().unwrap();
        // Blocks that hardly deflate, so a few fill the 128 KiB buffer; and none.
        let mut state = 1u64;
        for blocks in [0usize, 1, 2, 5, 23] {
            let payload: Vec<u8> = (0..blocks * DEFAULT_UNCOMPRESSED_BLOCK_SIZE + 999)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    (state >> 56) as u8
                })
                .collect();
            let mut writer = htsjdk_bgzf::BgzfWriter::new(Vec::new());
            std::io::Write::write_all(&mut writer, &payload).unwrap();
            let finished = writer.into_inner().unwrap();
            std::fs::write(path, &finished).unwrap();
            cut_bam_on_disk(path).unwrap();
            let expected = bam_on_disk(&finished);
            assert_eq!(std::fs::read(path).unwrap(), expected, "{blocks}");
            // Two blocks still fit the buffer together; five do not.
            assert_eq!(expected.is_empty(), blocks < 5, "{blocks}");
        }
        std::fs::remove_file(path).unwrap();
    }
}
