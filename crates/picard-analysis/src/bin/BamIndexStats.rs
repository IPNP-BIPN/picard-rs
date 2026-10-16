//! `BamIndexStats` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.BamIndexStats.doWork` at tag 3.4.0. The counting lives in
//! `picard_analysis::bam_index_stats`; this is the argument surface, which is almost nothing on
//! purpose: the tool takes an INPUT and prints to standard output, and every other argument it
//! declares is inherited and unused.
//!
//! That is what makes it worth an array anyway. `CREATE_INDEX` and `CREATE_MD5_FILE` are declared
//! by `CommandLineProgram` and do nothing here, `VALIDATION_STRINGENCY` reaches a reader that only
//! looks at the index, and a port that refused any of them would refuse a row the reference
//! accepts. An argument that is accepted and ignored is behaviour, and the only way to know a port
//! reproduces it is to pass it.
//!
//! The index is READ from beside the input, the way `SamFiles.findIndex` looks for it: `<stem>.bai`
//! when the input ends in `.bam`, `<input>.bai` otherwise. Rebuilding it instead -- which is what
//! this binary did first -- gives the same counts for a file that has an index and the wrong
//! ANSWER for a file that does not: htsjdk throws `No index for bam file <path>` and this printed
//! statistics. The array measured it as six refused rows against three answers.
//!
//! Both refusals are the reference's own text, including the `Exception in thread "main"` prefix
//! that a Picard tool prints when it throws rather than exits: a row that refuses is a row a port
//! has to refuse the same way, message included.

use std::io::Write;

use picard_analysis::bam_index_stats::bam_index_stats_with_index;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;

    let bam = std::fs::read(&input)?;
    // `SamFiles.findIndex`: the `.bam` suffix is dropped, anything else keeps its whole name.
    let index_path = match input.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{input}.bai"),
    };
    // The stream is sniffed before the index is looked for, which is the order the two refusals
    // come out in: a SAM text file is `Invalid GZIP header` and never reaches the index check.
    if !bam.starts_with(&[0x1f, 0x8b]) {
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMFormatException: Invalid GZIP header"
        );
        std::process::exit(1);
    }
    let Ok(bai) = std::fs::read(&index_path) else {
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMException: No index for bam file {input}"
        );
        std::process::exit(1);
    };
    let out = bam_index_stats_with_index(&bam, &bai).map_err(|e| format!("{e:?}"))?;
    let mut stdout = std::io::BufWriter::new(std::io::stdout().lock());
    stdout.write_all(out.as_bytes())?;
    stdout.flush()?;
    Ok(())
}
