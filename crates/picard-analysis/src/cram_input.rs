//! CRAM inputs, for every tool that reads a BAM.
//!
//! htsjdk's `SamReaderFactory` tells a CRAM from a BAM by its first four bytes and decodes it
//! against the reference the tool was given, `REFERENCE_SEQUENCE` (`R`), which every Picard tool
//! accepts and hands to the factory. The tools here sniff gzip's magic for a BAM and read anything
//! else as SAM text, so a CRAM is turned into the BAM it decodes to before that sniff: the same
//! header text and the same records, BGZF-wrapped without compression. Everything after the sniff
//! then sees what htsjdk's CRAM reader hands a tool, and reads it the way it reads a BAM.
//!
//! The decoding is htsjdk 5.0.0's, the version Picard 3.5.0 bundles: `MD` and `NM` are regenerated
//! for every mapped record that does not carry its own, which a BAM read never does.
//!
//! # A missing reference fails when a record needs it
//!
//! Without a reference, htsjdk falls back to `CRAMLazyReferenceSource`, which has nothing, and the
//! first mapped slice fails with `Failure getting reference bases for sequence <name>`. The header
//! reads; unmapped records decode. The failure is raised here, before the tool has written
//! anything, with the line htsjdk's uncaught exception prints and Picard's exit status.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_cram::reader::{CramError, CramReader, FastaReference, NoReference, ReferenceSource};
use std::io::Write;

/// The four bytes a CRAM begins with.
pub const CRAM_MAGIC: &[u8; 4] = b"CRAM";

/// Does this file start like a CRAM?
pub fn is_cram(raw: &[u8]) -> bool {
    raw.starts_with(CRAM_MAGIC)
}

/// The variable through which a parent hands its reader's reference to a child it runs on the
/// same input: `CollectMultipleMetrics` reads the input once, with its own reference, for every
/// program, and the port runs each program over the input again. An empty value is no reference.
pub const READER_REFERENCE_ENV: &str = "PICARD_RS_READER_REFERENCE";

/// What a tool's reader is given instead of `REFERENCE_SEQUENCE`, when that is not what it reads
/// with.
static OVERRIDE: std::sync::OnceLock<Option<String>> = std::sync::OnceLock::new();

/// What a failure to decode is reported as, when the tool decodes inside worker threads and
/// reports the failure of the pool rather than the cause.
static WRAPPER: std::sync::OnceLock<&'static str> = std::sync::OnceLock::new();

/// The reference this tool's reader is given, set once before any input is read: `None` for a
/// reader built without one (`CalculateFingerprintMetrics`'s, `IdentifyContaminant`'s), a path
/// for a tool whose reference argument has another name (`CollectRrbsMetrics`' `REFERENCE`).
pub fn set_reader_reference(reference: Option<String>) {
    let _ = OVERRIDE.set(reference.filter(|r| r != "null"));
}

/// The exception a decoding failure is reported as instead of its own, for the fingerprint tools
/// whose `ExecutorService` turns any failure into `Failed to fingerprint`.
pub fn wrap_failures_as(line: &'static str) {
    let _ = WRAPPER.set(line);
}

/// The reference this process's reader is given: an override, else the parent's, else
/// `REFERENCE_SEQUENCE=` or `R=` in either spelling, the last one winning, with Barclay's `null`
/// meaning none.
pub fn reference_argument() -> Option<String> {
    if let Some(reference) = OVERRIDE.get() {
        return reference.clone();
    }
    if let Ok(value) = std::env::var(READER_REFERENCE_ENV) {
        return Some(value).filter(|v| !v.is_empty() && v != "null");
    }
    argv_reference()
}

/// `REFERENCE_SEQUENCE=` or `R=` as this process was given it, in either spelling, the last one
/// winning, with Barclay's `null` meaning none: what a writer is made with.
pub fn argv_reference() -> Option<String> {
    let mut found = None;
    for raw in crate::args::process().into_iter().skip(1) {
        let raw = raw.trim_start_matches('-');
        if let Some((name, value)) = raw.split_once('=') {
            if name == "REFERENCE_SEQUENCE" || name == "R" {
                found = Some(value.to_string());
            }
        }
    }
    found.filter(|v| v != "null")
}

/// The reference source `SamReaderFactory.referenceSequence(REFERENCE_SEQUENCE)` builds.
pub fn reference_source() -> Box<dyn ReferenceSource> {
    match reference_argument() {
        Some(path) => Box::new(FastaReference::open(std::path::Path::new(&path))),
        None => Box::new(NoReference),
    }
}

/// A CRAM's header and every record, decoded against this process's reference.
pub fn decode(raw: &[u8]) -> Result<(SamHeader, String, Vec<BamRecord>), CramError> {
    let mut reader = CramReader::new(raw)?;
    let mut source = reference_source();
    let records = reader.read_all(source.as_mut())?;
    Ok((reader.header.clone(), reader.header_text.clone(), records))
}

/// [`decode`], failing as an uncaught exception does.
pub fn decode_or_exit(raw: &[u8]) -> (SamHeader, Vec<BamRecord>) {
    match decode(raw) {
        Ok((header, _, records)) => (header, records),
        Err(error) => fail(&error),
    }
}

/// An uncaught exception for `error`, or the one [`wrap_failures_as`] put in its place.
fn fail(error: &CramError) -> ! {
    match WRAPPER.get() {
        Some(line) => crate::metrics_cli::thrown(line),
        None => crate::metrics_cli::thrown(&error.java_line()),
    }
}

/// The file's bytes as the BAM sniff expects them: a CRAM becomes the uncompressed-BGZF BAM it
/// decodes to, anything else is returned untouched.
pub fn uncram(raw: Vec<u8>) -> Vec<u8> {
    if !is_cram(&raw) {
        return raw;
    }
    match decode(&raw) {
        Ok((header, text, records)) => as_bgzf_bam(&header, &text, &records),
        Err(error) => fail(&error),
    }
}

/// The file at `path`, read whole and passed through [`uncram`].
pub fn read(path: &str) -> std::io::Result<Vec<u8>> {
    std::fs::read(path).map(uncram)
}

/// A BAM holding `records`, with `text` as its header text and `header`'s dictionary, in BGZF
/// blocks stored rather than deflated: it is read back at once and never written out.
pub fn as_bgzf_bam(header: &SamHeader, text: &str, records: &[BamRecord]) -> Vec<u8> {
    let mut plain = Vec::new();
    plain.extend_from_slice(b"BAM\x01");
    let text_bytes: Vec<u8> = text.chars().map(|c| c as u32 as u8).collect();
    plain.extend_from_slice(&(text_bytes.len() as i32).to_le_bytes());
    plain.extend_from_slice(&text_bytes);
    plain.extend_from_slice(&(header.sequences.len() as i32).to_le_bytes());
    for sequence in &header.sequences {
        let name: Vec<u8> = sequence.name.chars().map(|c| c as u32 as u8).collect();
        plain.extend_from_slice(&((name.len() + 1) as i32).to_le_bytes());
        plain.extend_from_slice(&name);
        plain.push(0);
        plain.extend_from_slice(&sequence.length.to_le_bytes());
    }
    for record in records {
        if let Err(error) = record.encode_into(&mut plain) {
            crate::metrics_cli::thrown(&format!(
                "htsjdk.samtools.SAMException: a CRAM record does not encode as BAM: {error:?}"
            ));
        }
    }
    let mut bgzf = htsjdk_bgzf::BgzfWriter::with_level(Vec::new(), 0);
    bgzf.write_all(&plain)
        .expect("an in-memory BGZF writer never fails");
    bgzf.into_inner()
        .expect("an in-memory BGZF writer never fails")
}
