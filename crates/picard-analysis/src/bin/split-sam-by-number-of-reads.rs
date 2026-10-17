//! `SplitSamByNumberOfReads` as a runnable binary: the covering array's port side.
//!
//! The tool writes a DIRECTORY of shards, `<OUT_PREFIX>_%04d.<ext>`, so what the array compares is
//! one named file inside it rather than "the output". `SPLIT_TO_N_READS` and `SPLIT_TO_N_FILES`
//! are mutually exclusive, which no covering array can carry both of: the repository declares a
//! domain for the first and none for the second.
//!
//! The boundary is not a count of records. A shard ends at the read-count target OR at the end of
//! a queryname group, whichever comes later, so a shard can hold more reads than asked for and
//! never splits a template.

use std::io::Read;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::BamReader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam, write_sam};
use picard_analysis::args::arg;
use picard_analysis::split_sam_by_number_of_reads::{
    split_sam_by_number_of_reads, split_sam_by_number_of_reads_to_bam, SplitOptions,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SplitSamByNumberOfReads");
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let prefix = arg(&args, "OUT_PREFIX=").unwrap_or_else(|| "shard".to_string());
    let number = |key: &str| -> Result<i64, Box<dyn std::error::Error>> {
        Ok(arg(&args, key)
            .map(|v| v.parse::<i64>())
            .transpose()?
            .unwrap_or(0))
    };
    let options = SplitOptions {
        split_to_n_files: number("SPLIT_TO_N_FILES=")?,
        split_to_n_reads: number("SPLIT_TO_N_READS=")?,
        total_reads_in_input: number("TOTAL_READS_IN_INPUT=")?,
    };

    // `customCommandLineValidation`.
    if options.split_to_n_files <= 1 && options.split_to_n_reads <= 1 {
        eprintln!(
            "One of SPLIT_TO_N_FILES or SPLIT_TO_N_READS must be greater than 0. Found SPLIT_TO_N_FILES is {} and SPLIT_TO_N_READS is {}.",
            options.split_to_n_files, options.split_to_n_reads
        );
        std::process::exit(1);
    }
    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut raw = Vec::new();
    std::fs::File::open(&input)?.read_to_end(&mut raw)?;
    // `reader.type().fileExtension()`: the shards take the INPUT's format, so a SAM in gives SAM
    // shards named `.sam`, not BAM ones. A tool that always wrote BAM would put its output under a
    // name the reference never writes.
    let bam_input = raw.starts_with(&[0x1f, 0x8b]);
    // The library splits SAM text; a BAM input is decoded and re-encoded as text first, which
    // changes no record -- the shard bytes come from the writer either way.
    let text = if bam_input {
        let plain = htsjdk_bgzf::decompress_all(&raw).map_err(|e| format!("{e:?}"))?;
        let reader = BamReader::new(&plain).map_err(|e| format!("{e:?}"))?;
        let header: SamHeader = reader.header.text.clone();
        let records: Vec<BamRecord> = reader
            .map(|r| r.map_err(|e| format!("{e:?}")))
            .collect::<Result<_, _>>()?;
        write_sam(&header, &records).ok_or("records that parsed re-encode as SAM text")?
    } else {
        let text = String::from_utf8(raw)?;
        // Parse and re-encode so that a SAM input takes the same path as a BAM one.
        let (header, records) = read_sam(&text).map_err(|e| format!("{e:?}"))?;
        write_sam(&header, &records).ok_or("records that parsed re-encode as SAM text")?
    };

    std::fs::create_dir_all(&output)?;
    if bam_input {
        let shards =
            split_sam_by_number_of_reads_to_bam(&text, &options).map_err(|e| format!("{e:?}"))?;
        // `CREATE_INDEX` gives each BAM shard its `.bai`, which htsjdk writes only for a
        // coordinate-sorted file. Its virtual offsets point into this writer's blocks, so the
        // index is this BAM's and not byte for byte the reference's, as the BAM is not.
        let create_index = arg(&args, "CREATE_INDEX=")
            .is_some_and(|v| picard_analysis::args::is_true(&v))
            && text
                .lines()
                .next()
                .is_some_and(|hd| hd.starts_with("@HD") && hd.contains("\tSO:coordinate"));
        for (index, bytes) in shards.iter().enumerate() {
            let name = format!("{prefix}_{:04}.bam", index + 1);
            std::fs::write(std::path::Path::new(&output).join(name), bytes)?;
            if create_index {
                let bai = htsjdk_bam::build_index::build_bam_index(bytes)
                    .map_err(|e| format!("{e:?}"))?;
                let name = format!("{prefix}_{:04}.bai", index + 1);
                std::fs::write(std::path::Path::new(&output).join(name), bai)?;
            }
        }
    } else {
        let shards = split_sam_by_number_of_reads(&text, &options).map_err(|e| format!("{e:?}"))?;
        for (index, shard) in shards.iter().enumerate() {
            let name = format!("{prefix}_{:04}.sam", index + 1);
            std::fs::write(std::path::Path::new(&output).join(name), shard)?;
        }
    }
    Ok(())
}
