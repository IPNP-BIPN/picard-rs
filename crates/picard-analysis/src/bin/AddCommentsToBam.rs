//! `AddCommentsToBam` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.AddCommentsToBam.doWork` at tag 3.4.0. The reheader is
//! `picard_analysis::add_comments_to_bam`, a block copy with only the header re-encoded; this is the
//! argument surface around it. `COMMENT` is a list, each occurrence one `@CO` line in the order
//! given, and a `.sam` input is refused by its NAME before anything is read.
//!
//! The reheader asserts its INPUT is writable, so a corpus mounted read-only, as the oracle
//! contract mounts it, refuses every BAM row; the SAM rows are refused first, by name.

use picard_analysis::add_comments_to_bam::{add_comments_to_bam, AddCommentsError};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn refuse(message: &str) -> ! {
    eprintln!("Exception in thread \"main\" picard.PicardException: {message}");
    std::process::exit(1);
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let comments: Vec<String> = args
        .iter()
        .filter_map(|a| {
            a.strip_prefix("COMMENT=")
                .or_else(|| a.strip_prefix("C="))
                .map(str::to_string)
        })
        .collect();

    // `INPUT.getAbsolutePath().endsWith(".sam")`, after the readability check.
    let bam = std::fs::read(&input)?;
    if input.ends_with(".sam") {
        refuse("SAM files are not supported");
    }
    // `BamFileIoUtils.reheaderBamFile` asserts the INPUT is writable before it copies a block,
    // because the copy is written for the in-place case too, so a read-only corpus refuses every
    // BAM, naming it by its absolute path.
    if std::fs::OpenOptions::new()
        .write(true)
        .open(&input)
        .is_err()
    {
        let absolute = if std::path::Path::new(&input).is_absolute() {
            input.clone()
        } else {
            std::env::current_dir()
                .map(|cwd| cwd.join(&input).display().to_string())
                .unwrap_or_else(|_| input.clone())
        };
        eprintln!(
            "Exception in thread \"main\" htsjdk.samtools.SAMException: \
             File exists but is not writable: {absolute}"
        );
        std::process::exit(1);
    }
    let borrowed: Vec<&str> = comments.iter().map(String::as_str).collect();
    let out = match add_comments_to_bam(&bam, &borrowed) {
        Ok(out) => out,
        Err(AddCommentsError::NewlineInComment) => refuse("Comments can not contain a new line"),
        Err(other) => return Err(format!("{other:?}").into()),
    };
    std::fs::write(&output, out)?;
    Ok(())
}
