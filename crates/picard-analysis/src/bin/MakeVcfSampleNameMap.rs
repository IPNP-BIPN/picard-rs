//! `MakeVcfSampleNameMap` as a runnable binary: the covering array's port side.
//!
//! Two lines of logic and one surprise. The tool reads each input's header, refuses any VCF that
//! does not hold exactly one sample, and writes `<path>\t<sample>` for the rest -- but the ORDER
//! of those lines is a `java.util.HashMap`'s, keyed on the path as it was WRITTEN and sized from
//! the number of inputs. The library ports that; this reads the headers.
//!
//! The path in the output is the argument, never a resolved or absolute one, so the same file
//! named twice differently is two lines.

use picard_analysis::make_vcf_sample_name_map::{build, render, wrong_sample_count_message, Entry};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn args_all(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(key).map(str::to_string))
        .collect()
}

/// The `#CHROM` line's sample columns, which is `VCFHeader.getGenotypeSamples`.
fn samples_of(text: &str) -> Vec<String> {
    for line in text.lines() {
        if line.starts_with("#CHROM") {
            return line.split('\t').skip(9).map(str::to_string).collect();
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Vec::new()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let mut inputs = args_all(&args, "INPUT=");
    inputs.extend(args_all(&args, "I="));
    if inputs.is_empty() {
        return Err("INPUT= is required".into());
    }
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;

    let mut entries = Vec::with_capacity(inputs.len());
    for path in &inputs {
        let text = std::fs::read_to_string(path)?;
        let samples = samples_of(&text);
        if samples.len() != 1 {
            eprintln!(
                "Exception in thread \"main\" picard.PicardException: {}",
                wrong_sample_count_message(path, samples.len())
            );
            std::process::exit(1);
        }
        entries.push(Entry {
            path: path.clone(),
            sample: samples[0].clone(),
        });
    }

    std::fs::write(&output, render(&build(&entries)))?;
    Ok(())
}
