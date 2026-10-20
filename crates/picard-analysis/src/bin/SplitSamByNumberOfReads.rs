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

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use picard_analysis::args::arg;
use picard_analysis::split_sam_by_number_of_reads::{
    split_records, split_sam_by_number_of_reads, Shards, SplitOptions,
};
use picard_analysis::stream_output::{Output, OutputDir};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SplitSamByNumberOfReads");
    picard_analysis::output::install();
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

    let mut magic = Vec::new();
    std::fs::File::open(&input)?
        .take(4)
        .read_to_end(&mut magic)?;
    // `IOUtil.assertDirectoryIsWritable(OUTPUT)`, once the input is known to be readable.
    picard_analysis::stream_output::assert_directory_is_writable(&output);
    let cram_input = picard_analysis::cram_input::is_cram(&magic);
    // A CRAM's shards are named `.cram` and written by `makeSAMOrBAMWriter`, which writes BAM for
    // any name that is not `.sam`.
    let extension = if cram_input { "cram" } else { "bam" };
    // `reader.type().fileExtension()`: the shards take the INPUT's format, so a SAM in gives SAM
    // shards named `.sam`, not BAM ones. A tool that always wrote BAM would put its output under a
    // name the reference never writes.
    let bam_input = cram_input || magic.starts_with(&[0x1f, 0x8b]);
    let opened = picard_analysis::bam_stream::Input::open(&input)?;
    let streamed = opened.streams();
    let (header, records) = opened.into_parts();
    let shard_path = |index: usize| -> String {
        let extension = if bam_input { extension } else { "sam" };
        let name = format!("{prefix}_{:04}.{extension}", index + 1);
        std::path::Path::new(&output)
            .join(name)
            .display()
            .to_string()
    };
    // `CREATE_INDEX` gives each BAM shard its `.bai`, which htsjdk writes only for a
    // coordinate-sorted file. Its virtual offsets point into this writer's blocks, so the index is
    // this BAM's and not byte for byte the reference's, as the BAM is not.
    let create_index = arg(&args, "CREATE_INDEX=")
        .is_some_and(|v| picard_analysis::args::is_true(&v))
        && header.attributes.get("SO") == Some("coordinate");

    // A CRAM's shards are BAMs under a `.cram` name, and a negative TOTAL_READS_IN_INPUT is refused
    // once the input has been read: both are left to the whole-file split.
    if !cram_input && options.total_reads_in_input >= 0 {
        // Picard counts the records in a pass of its own, unless it is told how many there are.
        let (counted, records) = if options.total_reads_in_input == 0 {
            let counted = records.count() as i64;
            let again = picard_analysis::bam_stream::Input::reopen(&input, streamed)?;
            (counted, again.into_parts().1)
        } else {
            (0, records)
        };
        let mut shards = Shards::new(counted, &options).map_err(|e| format!("{e:?}"))?;
        let dir = OutputDir::create(&output)?;
        let mut out = Output::create(&shard_path(0), &header, false)?;
        let mut written = 1;
        for record in records {
            if shards.starts_new(&record) {
                out.finish(false, create_index)?;
                out = Output::create(&shard_path(written), &header, false)?;
                written += 1;
            }
            out.write(&record)?;
        }
        out.finish(false, create_index)?;
        dir.keep();
        return Ok(());
    }

    let records: Vec<BamRecord> = records.collect();
    let text = if bam_input {
        String::new()
    } else {
        write_sam(&header, &records).ok_or("records that parsed re-encode as SAM text")?
    };

    std::fs::create_dir_all(&output)?;
    if bam_input {
        let (header, shards) =
            split_records(header, records, &options).map_err(|e| format!("{e:?}"))?;
        let shards: Vec<Vec<u8>> = shards
            .iter()
            .map(|recs| -> Result<Vec<u8>, String> {
                let mut writer = picard_analysis::output::bam_writer(Vec::new(), &header)
                    .map_err(|e| format!("{e:?}"))?;
                for rec in recs {
                    writer.write(rec).map_err(|e| format!("{e:?}"))?;
                }
                writer.finish().map_err(|e| format!("{e:?}"))
            })
            .collect::<Result<_, _>>()?;
        for (index, bytes) in shards.iter().enumerate() {
            let name = format!("{prefix}_{:04}.{extension}", index + 1);
            std::fs::write(std::path::Path::new(&output).join(&name), bytes)?;
            if create_index {
                let bai = htsjdk_bam::build_index::build_bam_index(bytes)
                    .map_err(|e| format!("{e:?}"))?;
                // `createBamIndex`: `.bam` is replaced by `.bai`, anything else has it appended.
                let name = match name.strip_suffix(".bam") {
                    Some(stem) => format!("{stem}.bai"),
                    None => format!("{name}.bai"),
                };
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
