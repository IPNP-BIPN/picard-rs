//! `IntervalListToBed` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.util.IntervalListToBed.doWork` at tag 3.4.0 down to its arguments; the conversion
//! is `picard_analysis::interval_list_to_bed`.
//!
//! `SORT` (default true) is the only argument that reaches the output, and it does so only on an
//! input that is not already in coordinate order: it sorts through a `SortingCollection` under
//! `IntervalCoordinateComparator`, which keys on the sequence index taken from the header
//! dictionary, then start, then end, then positive strand first, then name. The array runs two
//! fixtures for exactly that reason, the second of which leads with `chr2` and carries both
//! strands.
//!
//! `SCORE` is written verbatim into every line's fifth column, and the array holds it at its
//! default: it is an `int` with no declared bounds, so any other value would be invented, and an
//! invented value covers nothing.
//!
//! `CREATE_INDEX` reaches nothing. The output is BED text through `IOUtil.openFileForBufferedWriting`,
//! not a `SAMFileWriter`, so there is no index to create and no refusal to raise. Neither does
//! `VALIDATION_STRINGENCY`: the input is read by `IntervalListCodec`, which does not consult it.

use std::io::Write;

use picard_analysis::interval_list_to_bed::{interval_list_to_bed, Options};

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// A Barclay `boolean` argument, which is `true`/`false` and nothing else.
fn flag(args: &[String], key: &str, default: bool) -> Result<bool, String> {
    match arg(args, key).as_deref() {
        None => Ok(default),
        Some("true") => Ok(true),
        Some("false") => Ok(false),
        Some(other) => Err(format!("unknown value for {key}{other}")),
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
    let opts = Options {
        score: match arg(&args, "SCORE=") {
            Some(value) => value.parse().map_err(|_| format!("SCORE: {value}"))?,
            None => 500,
        },
        sort: flag(&args, "SORT=", true)?,
    };

    let list = std::fs::read_to_string(&input)?;
    let bed = interval_list_to_bed(&list, &opts).map_err(|e| format!("{e:?}"))?;

    let mut out = std::io::BufWriter::new(std::fs::File::create(&output)?);
    out.write_all(bed.as_bytes())?;
    out.flush()?;
    Ok(())
}
