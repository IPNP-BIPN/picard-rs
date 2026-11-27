//! CRAM outputs, for every tool whose writer comes from `SAMFileWriterFactory.makeWriter`.
//!
//! `makeWriter(header, presorted, output, REFERENCE_SEQUENCE)` returns a `CRAMFileWriter` when the
//! output's name ends in `.cram`, and a SAM or BAM writer otherwise. The tools here build their
//! BAM in memory and write it at the end, so a CRAM output is that BAM's header and records
//! written again as CRAM: the header the BAM writer stamped is the one the CRAM writer stamps,
//! and the records are the same records.
//!
//! The CRAM is htsjdk 5.0.0's, the version Picard 3.5.0 bundles: CRAM 3.1 under the NORMAL
//! profile, with the output's `file:` URI as the file definition's id.
//!
//! # A missing reference leaves an empty file
//!
//! `Files.newOutputStream` creates the output when the writer is made, and the header the writer
//! writes first sits in a buffer that is never flushed when the first mapped record finds no
//! reference. Picard exits with the file there and empty, and so does this.

use htsjdk_bam::reader::BamReader;

/// The memory a CRAM may take, unless `PICARD_RS_CRAM_MEMORY` (bytes) says otherwise. The CRAM
/// writer and reader take the whole file at once and hold it decoded, so past this a run would
/// exhaust the memory rather than finish.
/// shortcut: a refusal, not a CRAM; lift it once the CRAM writer and reader stream.
const CRAM_MEMORY: u64 = 16 << 30;

/// Peak memory per byte of the BAM a CRAM output is written from: 24.7 and 23.6 measured on 62
/// and 112 MiB BAMs (SortSam, chr20 downsampled).
pub const PER_BAM_BYTE_OUT: u64 = 25;

/// Peak memory per byte of a CRAM input: 76 measured on a 23 MiB CRAM (CollectQualityYieldMetrics
/// and MeanQualityByCycle), each record decoded once.
pub const PER_CRAM_BYTE_IN: u64 = 80;

/// The same for SortSam, which decodes the CRAM, re-encodes it as a BAM and decodes that again:
/// 196 and 192 measured on 12 and 23 MiB CRAMs.
pub const PER_CRAM_BYTE_SORTED: u64 = 200;

fn cram_memory() -> u64 {
    static MAX: std::sync::OnceLock<u64> = std::sync::OnceLock::new();
    *MAX.get_or_init(|| {
        std::env::var("PICARD_RS_CRAM_MEMORY")
            .ok()
            .and_then(|value| value.parse().ok())
            .unwrap_or(CRAM_MEMORY)
    })
}

/// Exits, before memory runs out, when `what`, `bytes` long at `per_byte` bytes of memory each,
/// would take more than a CRAM may. Picard streams a CRAM and has no such limit, so this is
/// picard-rs's own refusal, not one of Picard's exceptions.
pub fn refuse_past_the_cram_limit(what: &str, bytes: u64, per_byte: u64) {
    if let Some(message) = past_the_cram_limit(what, bytes, per_byte, cram_memory()) {
        eprintln!("{message}");
        std::process::exit(1);
    }
}

fn past_the_cram_limit(what: &str, bytes: u64, per_byte: u64, memory: u64) -> Option<String> {
    let needs = bytes.saturating_mul(per_byte);
    (needs > memory).then(|| {
        let gib = |b: u64| b as f64 / (1u64 << 30) as f64;
        format!(
            "picard-rs cannot hold this CRAM in memory: {what} is {bytes} bytes, which needs about \
             {:.1} GiB ({per_byte} times its size), past the {:.1} GiB a CRAM may take. picard-rs \
             reads and writes a CRAM whole until its CRAM writer and reader stream. Write a BAM \
             instead, or set PICARD_RS_CRAM_MEMORY (bytes) higher on a machine with the memory.",
            gib(needs),
            gib(memory)
        )
    })
}

/// The BAM a CRAM output is built from, in memory, refused once it outgrows the CRAM limit.
#[derive(Default)]
pub struct CramBuffer(pub Vec<u8>);

impl std::io::Write for CramBuffer {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let total = (self.0.len() + bytes.len()) as u64;
        refuse_past_the_cram_limit(
            "the BAM a CRAM output is written from",
            total,
            PER_BAM_BYTE_OUT,
        );
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// `CRAM_TYPE.hasValidFileExtension`: the name ends in `.cram`.
pub fn is_cram(path: &str) -> bool {
    path.ends_with(".cram")
}

/// The id `SAMFileWriterFactory` gives a CRAM it writes to a file: `outputFile.toUri()`.
pub fn file_id(path: &str) -> String {
    crate::fingerprint::file_uri(path)
}

/// A CRAM of `header` and `records`, against this process's reference, or the line an uncaught
/// exception prints when it fails.
pub fn encode(
    path: &str,
    header: &htsjdk_bam::header::SamHeader,
    records: &[htsjdk_bam::record::BamRecord],
) -> Result<Vec<u8>, String> {
    let mut source = writer_reference();
    htsjdk_cram::writer::write_cram(header, records, &file_id(path), source.as_mut())
        .map_err(|error| error.java_line())
}

/// The reference `makeWriter(..., REFERENCE_SEQUENCE)` builds a CRAM writer with: the FASTA, or,
/// with none, `CRAMLazyReferenceSource`, which refuses the first mapped record with its own
/// message rather than the one a reader without a reference gives.
fn writer_reference() -> Box<dyn htsjdk_cram::reader::ReferenceSource> {
    match crate::cram_input::argv_reference() {
        Some(path) => Box::new(htsjdk_cram::reader::FastaReference::open(
            std::path::Path::new(&path),
        )),
        None => Box::new(htsjdk_cram::reader::LazyReference),
    }
}

/// The finished BAM `bam` written to `path`: as it is, or, for a `.cram` path, as the CRAM of the
/// same header and records, indexed when `CREATE_INDEX` is set, which is the default every
/// `new SAMFileWriterFactory()` takes from the command line. A CRAM that cannot be written leaves
/// the file empty and exits as Picard does.
pub fn write(path: &str, bam: Vec<u8>) -> std::io::Result<()> {
    if !is_cram(path) {
        return std::fs::write(path, bam);
    }
    write_from_bam(path, &bam, create_index_argument())
}

/// `CREATE_INDEX` as this process was given it, in either spelling, the last one winning.
pub fn create_index_argument() -> bool {
    boolean_argument("CREATE_INDEX=")
}

fn boolean_argument(prefix: &str) -> bool {
    let mut value = false;
    for raw in crate::args::process().into_iter().skip(1) {
        if let Some(v) = raw.trim_start_matches('-').strip_prefix(prefix) {
            value = v.eq_ignore_ascii_case("true");
        }
    }
    value
}

/// The CRAM, and with `CREATE_MD5_FILE` its digest beside it as `Md5CalculatingOutputStream`
/// leaves it: thirty-two lower-case hex digits and no newline.
fn write_with_md5(path: &str, bytes: &[u8]) -> std::io::Result<()> {
    use md5::{Digest, Md5};
    std::fs::write(path, bytes)?;
    if boolean_argument("CREATE_MD5_FILE=") {
        let digest: String = Md5::digest(bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        std::fs::write(format!("{path}.md5"), digest)?;
    }
    Ok(())
}

/// `header` and `records` written to `path` as a CRAM.
pub fn write_records(
    path: &str,
    header: &htsjdk_bam::header::SamHeader,
    records: &[htsjdk_bam::record::BamRecord],
) -> std::io::Result<()> {
    match encode(path, header, records) {
        Ok(bytes) => write_with_md5(path, &bytes),
        Err(line) => {
            std::fs::write(path, [])?;
            crate::metrics_cli::thrown(&line)
        }
    }
}

/// [`write`] for a tool that also honours `CREATE_INDEX`, given the BAM it built: a CRAM output is
/// indexed by `CRAMBAIIndexer` into `<output>.bai`, which refuses a header that is not
/// coordinate-sorted where a BAM writer would quietly skip the index.
pub fn write_from_bam(path: &str, bam: &[u8], create_index: bool) -> std::io::Result<()> {
    let bytes = bam.len() as u64;
    refuse_past_the_cram_limit(
        "the BAM a CRAM output is written from",
        bytes,
        PER_BAM_BYTE_OUT,
    );
    let plain =
        crate::bgzf::decompress_all(bam).map_err(|e| std::io::Error::other(format!("{e:?}")))?;
    let reader = BamReader::new(&plain).map_err(|e| std::io::Error::other(format!("{e:?}")))?;
    let header = reader.header.text.clone();
    let records = reader
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| std::io::Error::other(format!("{e:?}")))?;
    write_records_indexed(path, &header, &records, create_index)
}

/// [`write_records`], with `CREATE_INDEX`.
pub fn write_records_indexed(
    path: &str,
    header: &htsjdk_bam::header::SamHeader,
    records: &[htsjdk_bam::record::BamRecord],
    create_index: bool,
) -> std::io::Result<()> {
    if !create_index {
        return write_records(path, header, records);
    }
    let index_path = format!("{path}.bai");
    // The index stream is opened, then the CRAM's, then the indexer refuses.
    if crate::metrics_cli::sort_order(header) != "coordinate" {
        std::fs::write(&index_path, [])?;
        std::fs::write(path, [])?;
        crate::metrics_cli::thrown(
            "htsjdk.samtools.SAMException: CRAM file mut be coordinate-sorted for indexing.",
        );
    }
    let mut source = writer_reference();
    let written =
        htsjdk_cram::writer::write_cram_indexed(header, records, &file_id(path), source.as_mut())
            .map_err(|error| error.java_line());
    match written {
        Ok((bytes, index)) => {
            write_with_md5(path, &bytes)?;
            std::fs::write(&index_path, index)
        }
        Err(line) => {
            std::fs::write(&index_path, [])?;
            std::fs::write(path, [])?;
            crate::metrics_cli::thrown(&line)
        }
    }
}

#[cfg(test)]
mod limit_tests {
    use super::*;

    #[test]
    fn a_cram_past_the_limit_is_refused_with_its_size() {
        assert_eq!(past_the_cram_limit("INPUT", 5 << 30, 2, 10 << 30), None);
        let message = past_the_cram_limit("INPUT", 6 << 30, 2, 10 << 30).unwrap();
        assert!(
            message.contains("INPUT is 6442450944 bytes, which needs about 12.0 GiB (2 times"),
            "{message}"
        );
        assert!(message.contains("past the 10.0 GiB"), "{message}");
        assert!(message.contains("PICARD_RS_CRAM_MEMORY"), "{message}");
    }
}
