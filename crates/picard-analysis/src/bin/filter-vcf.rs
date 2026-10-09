//! `FilterVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.filter.FilterVcf.doWork` at tag 3.4.0. The filters are
//! `picard_analysis::filter_vcf`, and the read and write round trip is `picard_analysis::vcf_io`.
//!
//! An OUTPUT named `.vcf` or `.bcf` needs the input's dictionary whatever `CREATE_INDEX` says,
//! because the tool sets it on the writer before it knows whether an index will be built. The
//! script filter (`JAVASCRIPT_FILE`) runs a script engine and is not ported.

use picard_analysis::filter_vcf::{add_header_lines, filter_record, Thresholds};
use picard_analysis::vcf_io::{die, header_dictionary, read_path, write_output};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn number<T: std::str::FromStr>(args: &[String], key: &str, default: T) -> T {
    match arg(args, key) {
        None => default,
        Some(value) => value.parse().unwrap_or_else(|_| {
            eprintln!("{key} is not a number: {value}");
            std::process::exit(1);
        }),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    if arg(&args, "JAVASCRIPT_FILE=").is_some() || arg(&args, "JS=").is_some() {
        return Err("JAVASCRIPT_FILE is not ported".into());
    }
    let create_index = arg(&args, "CREATE_INDEX=")
        .map(|value| value == "true")
        .unwrap_or(true);
    let thresholds = Thresholds {
        min_ab: number(&args, "MIN_AB=", 0.0),
        min_dp: number(&args, "MIN_DP=", 0),
        min_gq: number(&args, "MIN_GQ=", 0),
        max_fs: number(&args, "MAX_FS=", f64::MAX),
        min_qd: number(&args, "MIN_QD=", 0.0),
    };

    let vcf = read_path(&input).unwrap_or_else(|exception| die(&exception));
    let dictionary = header_dictionary(&vcf.file.header);
    let name = std::path::Path::new(&output)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    if (name.ends_with(".vcf") || name.ends_with(".bcf")) && dictionary.is_none() {
        die(
            "picard.PicardException: The input vcf must have a sequence dictionary in order to \
             create indexed vcf or bcfs.",
        );
    }

    let mut header = vcf.file.header.clone();
    add_header_lines(&mut header);
    let mut records = vcf.records;
    for record in &mut records {
        record.decode();
        filter_record(&mut record.variant, &thresholds)
            .unwrap_or_else(|message| die(&format!("java.lang.NumberFormatException: {message}")));
    }
    let index_dictionary = if create_index {
        dictionary.as_deref()
    } else {
        None
    };
    write_output(&output, &header, &records, index_dictionary)?;
    Ok(())
}
