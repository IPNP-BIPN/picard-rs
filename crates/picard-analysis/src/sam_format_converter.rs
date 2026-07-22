//! `SamFormatConverter`, the BAM-to-SAM path.
//!
//! Ports `picard.sam.SamFormatConverter` (which delegates to `SamFileConverter.convertSamToSam`) at
//! tag 3.4.0, for the **BAM input, SAM output** direction: decode a BAM and write it as SAM text,
//! unchanged. The tool adds no `@PG` and no timestamp, so the whole SAM is comparable raw.
//!
//! There is no transform here: the records pass through in file order. This is exactly the
//! `BamReader` → `write_sam` pipeline the throughput benchmark already exercises byte-for-byte over
//! two million reads; this exposes it as its own tool. The SAM-to-BAM and CRAM directions (which write
//! BGZF, so the deflater choice matters) are separate surfaces.

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;

/// `SamFormatConverter` for a BAM (already BGZF-decompressed) to SAM text.
pub fn bam_to_sam(decompressed_bam: &[u8]) -> Result<String, String> {
    let reader = BamReader::new(decompressed_bam).map_err(|e| format!("{e:?}"))?;
    let header = reader.header.text.clone();
    let records: Vec<BamRecord> = reader
        .map(|r| r.map_err(|e| format!("{e:?}")))
        .collect::<Result<_, _>>()?;
    write_sam(&header, &records).ok_or_else(|| "records failed to re-encode as SAM".to_string())
}
