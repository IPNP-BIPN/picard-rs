//! `SamFormatConverter`, the SAM/BAM conversions.
//!
//! Ports `picard.sam.SamFormatConverter` (which delegates to `SamFileConverter`) at tag 3.4.0, for the
//! two directions between SAM and BAM. The tool adds no `@PG` and no timestamp and applies no
//! transform: the records pass through in file order, so each output is comparable raw.
//!
//! - [`bam_to_sam`] decodes a BAM and writes it as SAM text. This is exactly the `BamReader` →
//!   `write_sam` pipeline the throughput benchmark exercises byte-for-byte over two million reads.
//! - [`sam_to_bam`] reads SAM text and writes a BAM with htsjdk-rs's `BamWriter`, which is byte-identical
//!   to htsjdk's `BAMFileWriter` (its `every_file_is_byte_identical_to_htsjdks` conformance). Because
//!   the BGZF blocks come from the port's zlib writer, the match holds against Picard run with
//!   `USE_JDK_DEFLATER=true` (java.util.zip); Picard's default GKL/igzip deflater emits different BGZF
//!   bytes and is a separate surface, as is CRAM.

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::writer::BamWriter;

/// `SamFormatConverter` for a BAM (already BGZF-decompressed) to SAM text.
pub fn bam_to_sam(decompressed_bam: &[u8]) -> Result<String, String> {
    let reader = BamReader::new(decompressed_bam).map_err(|e| format!("{e:?}"))?;
    let header = reader.header.text.clone();
    let records: Vec<BamRecord> = reader
        .map(|r| r.map_err(|e| format!("{e:?}")))
        .collect::<Result<_, _>>()?;
    write_sam(&header, &records).ok_or_else(|| "records failed to re-encode as SAM".to_string())
}

/// `SamFormatConverter` for SAM text to a BAM, byte-identical to Picard's `USE_JDK_DEFLATER=true`
/// output. The returned bytes are the whole BAM file (BGZF-framed, with the terminator block).
pub fn sam_to_bam(sam_text: &str) -> Result<Vec<u8>, String> {
    // The tool opens the input at VALIDATION_STRINGENCY.SILENT; stringency does not reach the bytes.
    let (header, records) =
        read_sam_with(sam_text, ValidationStringency::Lenient).map_err(|e| format!("{e:?}"))?;
    let mut writer = BamWriter::new(Vec::new(), &header).map_err(|e| format!("{e:?}"))?;
    for rec in &records {
        writer.write(rec).map_err(|e| format!("{e:?}"))?;
    }
    writer.finish().map_err(|e| format!("{e:?}"))
}
