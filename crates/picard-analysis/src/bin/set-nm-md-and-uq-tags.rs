//! `SetNmMdAndUqTags` as a runnable binary: the covering array's port side.
//!
//! The library entry points took SAM text and the default options, which is one of the four
//! outputs this array's rows produce. The other three come from the two arguments the tool has:
//! `SET_ONLY_UQ` leaves `MD` and `NM` as the input had them, and `IS_BISULFITE_SEQUENCE` changes
//! `NM` and `UQ` (not `MD`, whose comparison has no bisulfite branch).
//!
//! The refusal is htsjdk's rather than Picard's, and it names the order it found:
//!
//! ```text
//! htsjdk.samtools.SAMException: Input must be coordinate-sorted for this program to run. Found: queryname
//! ```
//!
//! The output is a BAM, written presorted, with the input's header unchanged and no `@PG`.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::read_sam;
use htsjdk_bam::writer::BamWriter;
use picard_analysis::set_nm_md_and_uq_tags::{fix_record, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let reference = arg(&args, "REFERENCE_SEQUENCE=")
        .or_else(|| arg(&args, "R="))
        .ok_or("REFERENCE_SEQUENCE= is required")?;
    let flag = |key: &str| arg(&args, key).map(|v| v == "true").unwrap_or(false);
    let options = Options {
        is_bisulfite_sequence: flag("IS_BISULFITE_SEQUENCE="),
        set_only_uq: flag("SET_ONLY_UQ="),
    };

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    let (header, mut records): (SamHeader, Vec<BamRecord>) = if raw.starts_with(&[0x1f, 0x8b]) {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header = reader.header.text.clone();
        let records = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        (header, records)
    } else {
        let text = String::from_utf8(raw)?;
        read_sam(&text).map_err(|e| format!("{e:?}"))?
    };

    // The check is the tool's first act after opening the reader, before the writer exists, so a
    // refused row leaves no output file behind.
    let found = header.attributes.get("SO").unwrap_or("unsorted");
    if found != "coordinate" {
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMException: Input must be \
             coordinate-sorted for this program to run. Found: {found}"
        );
        std::process::exit(1);
    }

    let contigs = htsjdk_bam::fasta::read_fasta_file(&reference).map_err(|e| format!("{e:?}"))?;
    let bases: std::collections::HashMap<&str, &[u8]> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.bases.as_slice()))
        .collect();

    for record in &mut records {
        if record.reference_index < 0 {
            continue;
        }
        let name = &header.sequences[record.reference_index as usize].name;
        let reference_bases = *bases
            .get(name.as_str())
            .ok_or_else(|| format!("no reference bases for contig {name}"))?;
        fix_record(record, reference_bases, options);
    }

    let mut writer = BamWriter::new(Vec::new(), &header).map_err(|e| format!("{e:?}"))?;
    for record in &records {
        writer.write(record).map_err(|e| format!("{e:?}"))?;
    }
    std::fs::write(&output, writer.finish().map_err(|e| format!("{e:?}"))?)?;
    Ok(())
}
