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
