//! `AddOATag`.
//!
//! Ports `picard.sam.AddOATag.setOATag` at tag 3.4.0, for the default path (no `INTERVAL_LIST`, so
//! every record is tagged). Stamps each record with the `OA` (original alignment) tag, a
//! semicolon-terminated `RNAME,POS,strand,CIGAR,MAPQ,NM` field, appended to any `OA` already
//! present. A per-record transform: math-free (the `NM` is copied from the record's existing tag,
//! not recomputed) and independent, so it runs on all cores and stays byte-identical.
//!
//! `doWork` copies the header unchanged and adds no `@PG` and no timestamp, so the whole output is
//! comparable raw.

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text_parse::ParseError;
use rayon::prelude::*;

const READ_UNMAPPED: u16 = 0x4;
const READ_NEGATIVE_STRAND: u16 = 0x10;

/// The existing value of a `Z` string tag, or `""`.
fn str_tag(rec: &BamRecord, name: &[u8; 2]) -> String {
    match rec.tags.get(Tag::new(name)) {
        Some(TagValue::Str(s)) => s.clone(),
        _ => String::new(),
    }
}

/// `SAMRecord.getAttribute("NM")` rendered as `AddOATag` does: the integer, or `""` when absent.
fn nm_string(rec: &BamRecord) -> String {
    match rec.tags.get(Tag::new(b"NM")) {
        Some(TagValue::Int(v)) => v.to_string(),
        _ => String::new(),
    }
}

/// `AddOATag.setOATag`. `seq_names` maps a reference index to its `@SQ` name.
fn set_oa_tag(rec: &mut BamRecord, seq_names: &[String]) {
    let strand = if rec.flags & READ_NEGATIVE_STRAND != 0 {
        '-'
    } else {
        '+'
    };

    let oa_value = if rec.flags & READ_UNMAPPED != 0 {
        format!("*,0,{strand},*,255,;")
    } else {
        let ref_name = &seq_names[rec.reference_index as usize];
        assert!(
            !ref_name.contains(','),
            "Reference name for record {} contains a comma character.",
            rec.read_name
        );
        format!(
            "{},{},{},{},{},{};",
            ref_name,
            rec.alignment_start,
            strand,
            rec.cigar.to_text(),
            rec.mapping_quality,
            nm_string(rec),
        )
    };

    let combined = str_tag(rec, b"OA") + &oa_value;
    rec.tags.insert(Tag::new(b"OA"), TagValue::Str(combined));
}

/// `AddOATag.doWork` for SAM I/O, default path: tag every record and rewrite.
pub fn add_oa_tag(input_sam: &str) -> Result<String, ParseError> {
    let (header, mut records) = read_sam(input_sam)?;
    let seq_names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    // Each record's OA is computed from its own fields and tags, so the work is embarrassingly
    // parallel and order-preserving (decision 0006).
    records
        .par_iter_mut()
        .for_each(|rec| set_oa_tag(rec, &seq_names));
    Ok(write_sam(&header, &records).expect("records that parsed re-encode as SAM text"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: &str = "@HD\tVN:1.6\tSO:coordinate\n\
        @SQ\tSN:chr1\tLN:1000\n\
        mapped\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:2\n\
        revNoNm\t16\tchr1\t200\t30\t4M\t*\t0\t0\tACGT\tIIII\n\
        unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";

    fn oa(sam: &str, name: &str) -> String {
        let row = sam.lines().find(|l| l.starts_with(name)).unwrap();
        row.split('\t')
            .find(|f| f.starts_with("OA:Z:"))
            .unwrap()
            .strip_prefix("OA:Z:")
            .unwrap()
            .to_string()
    }

    #[test]
    fn a_mapped_read_records_its_alignment_with_the_nm_tag() {
        let out = add_oa_tag(INPUT).unwrap();
        assert_eq!(oa(&out, "mapped"), "chr1,100,+,4M,60,2;");
    }

    #[test]
    fn a_reverse_read_without_nm_leaves_the_nm_field_empty() {
        let out = add_oa_tag(INPUT).unwrap();
        assert_eq!(oa(&out, "revNoNm"), "chr1,200,-,4M,30,;");
    }

    #[test]
    fn an_unmapped_read_uses_the_star_form_with_mapq_255() {
        let out = add_oa_tag(INPUT).unwrap();
        assert_eq!(oa(&out, "unmapped"), "*,0,+,*,255,;");
    }

    #[test]
    fn an_existing_oa_tag_is_appended_to() {
        let input = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
            r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tOA:Z:old,1,+,4M,10,;\n";
        let out = add_oa_tag(input).unwrap();
        assert_eq!(oa(&out, "r"), "old,1,+,4M,10,;chr1,100,+,4M,60,;");
    }
}
