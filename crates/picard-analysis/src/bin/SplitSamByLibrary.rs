//! `SplitSamByLibrary` as a runnable binary: the covering array's port side.
//!
//! One output per `@RG`'s `LB`, named after the library, in a directory -- plus `unknown.bam` for
//! records whose read group names no library. Each output carries the INPUT's whole header, read
//! groups of other libraries included, which is Picard's own choice and not a simplification here.

use std::io::Read;

use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use picard_analysis::args::arg;
use picard_analysis::split_sam_by_library::{libraries, split_sam_by_library};
use picard_analysis::stream_output::{Output, OutputDir};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("SplitSamByLibrary");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    if let Some(stringency) = arg(&args, "VALIDATION_STRINGENCY=") {
        if !matches!(stringency.as_str(), "STRICT" | "LENIENT" | "SILENT") {
            return Err(format!("unknown VALIDATION_STRINGENCY: {stringency}").into());
        }
    }

    let mut magic = Vec::new();
    std::fs::File::open(&input)?
        .take(4)
        .read_to_end(&mut magic)?;
    let cram_input = picard_analysis::cram_input::is_cram(&magic);
    // The outputs take the INPUT's format, as `reader.type().fileExtension()` gives it.
    let bam_input = cram_input || magic.starts_with(&[0x1f, 0x8b]);
    let (header, records) = picard_analysis::bam_stream::Input::open(&input)?.into_parts();
    let extension = if bam_input { "bam" } else { "sam" };
    let path_of = |name: &str| -> String {
        let path = std::path::Path::new(&output).join(format!("{name}.{extension}"));
        path.display().to_string()
    };

    // Two libraries whose file names collide, or one named like the `unknown` output, are left to
    // the whole-file split, which writes one file after the other.
    let streams = !cram_input
        && libraries(&header).is_ok_and(|libraries| {
            let mut names: Vec<&str> = libraries.outputs.iter().map(|(n, _)| n.as_str()).collect();
            names.push("unknown");
            names.sort_unstable();
            names.windows(2).all(|pair| pair[0] != pair[1])
        });
    if streams {
        let libraries = libraries(&header).expect("checked above");
        let dir = OutputDir::create(&output)?;
        let mut outputs = libraries
            .outputs
            .iter()
            .map(|(name, header)| Output::create(&path_of(name), header, false))
            .collect::<Result<Vec<_>, _>>()?;
        let mut unknown = None;
        // A BAM input's records are written back from the bytes they were read from; a SAM one's
        // are re-encoded, as the text always was.
        for record in records {
            match libraries.route(&record) {
                Some(index) => outputs[index].write(&record)?,
                None => {
                    if unknown.is_none() {
                        let path = path_of("unknown");
                        unknown = Some(Output::create(&path, &libraries.unknown, false)?);
                    }
                    unknown.as_mut().expect("opened").write(&record)?
                }
            }
        }
        for out in outputs.into_iter().chain(unknown) {
            out.finish(false, false)?;
        }
        dir.keep();
        return Ok(());
    }

    let records: Vec<BamRecord> = records.collect();
    let text = write_sam(&header, &records).ok_or("records that parsed re-encode as SAM text")?;

    std::fs::create_dir_all(&output)?;
    // A CRAM in is CRAMs out, `makeWriter` with the reference, and `CREATE_INDEX` is honoured by
    // `CRAMBAIIndexer`, which refuses an output that is not coordinate-sorted.
    if cram_input {
        let create_index = arg(&args, "CREATE_INDEX=").as_deref() == Some("true");
        for (name, header, records) in picard_analysis::split_sam_by_library::split_outputs(&text)
            .map_err(|e| format!("{e:?}"))?
        {
            let path = std::path::Path::new(&output).join(format!("{name}.cram"));
            picard_analysis::cram_output::write_records_indexed(
                &path.display().to_string(),
                &header,
                &records,
                create_index,
            )?;
        }
        return Ok(());
    }
    if bam_input {
        let outputs = picard_analysis::split_sam_by_library::split_records(header, records)
            .map_err(|e| format!("{e:?}"))?;
        for (name, header, records) in outputs {
            let mut writer = picard_analysis::output::bam_writer(Vec::new(), &header)
                .map_err(|e| format!("{e:?}"))?;
            for record in &records {
                writer.write(record).map_err(|e| format!("{e:?}"))?;
            }
            std::fs::write(
                std::path::Path::new(&output).join(format!("{name}.bam")),
                writer.finish().map_err(|e| format!("{e:?}"))?,
            )?;
        }
    } else {
        for (name, text) in split_sam_by_library(&text).map_err(|e| format!("{e:?}"))? {
            std::fs::write(
                std::path::Path::new(&output).join(format!("{name}.sam")),
                text,
            )?;
        }
    }
    Ok(())
}
