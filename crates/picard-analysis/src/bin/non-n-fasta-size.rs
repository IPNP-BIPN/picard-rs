//! `NonNFastaSize` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.reference.NonNFastaSize.doWork` at tag 3.4.0. The whole-genome count is
//! `picard_analysis::non_n_fasta_size`; what this adds is the argument that turns it into a
//! different question.
//!
//! `INTERVALS` is optional in the tool and always present in the array, because the corpus declares
//! an interval list for it. With one, `doWork` builds an `IntervalListReferenceSequenceMask` and
//! asks it position by position:
//!
//! ```java
//! if (referenceSequenceMask.get(sequence.getContigIndex(), i + 1)) {
//!     nonNbases += bases[i] == SequenceUtil.N ? 0 : 1;
//! }
//! ```
//!
//! So the count is over masked positions only, the mask is a union of the list's intervals (they
//! are not uniqued first, and overlapping ones simply set the same bits), and positions are
//! 1-based. Without `INTERVALS` the mask admits everything, which is the module's case.
//!
//! Two details that decide bytes: the bases are uppercased before the comparison, so a soft-masked
//! `n` counts as an `N` and is not counted; and the output is `nonNbases + "\n"`, a decimal with a
//! trailing newline.

use std::io::Write;

use htsjdk_bam::fasta::read_fasta;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// One interval of a Picard interval list: 1-based, inclusive.
struct Interval {
    contig: String,
    start: usize,
    end: usize,
}

fn read_intervals(path: &str) -> Result<Vec<Interval>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut intervals = Vec::new();
    for line in text.lines() {
        if line.starts_with('@') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        intervals.push(Interval {
            contig: fields[0].to_string(),
            start: fields[1].parse().map_err(|_| format!("start: {line}"))?,
            end: fields[2].parse().map_err(|_| format!("end: {line}"))?,
        });
    }
    Ok(intervals)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    let input = arg(&args, "INPUT=")
        .or_else(|| arg(&args, "I="))
        .ok_or("INPUT= is required")?;
    let output = arg(&args, "OUTPUT=")
        .or_else(|| arg(&args, "O="))
        .ok_or("OUTPUT= is required")?;
    let intervals = match arg(&args, "INTERVALS=") {
        Some(path) => Some(read_intervals(&path)?),
        None => None,
    };

    let fasta = std::fs::read_to_string(&input)?;
    let sequences = read_fasta(fasta.as_bytes()).map_err(|e| format!("{e:?}"))?;

    let mut non_n: u64 = 0;
    for sequence in &sequences {
        // The mask is a union of bits, so overlapping intervals count a position once. Building it
        // per contig keeps that true without sorting or uniquing the list, which the tool does not
        // do either.
        let mask: Option<Vec<bool>> = intervals.as_ref().map(|intervals| {
            let mut bits = vec![false; sequence.bases.len() + 1];
            for interval in intervals {
                if interval.contig != sequence.name {
                    continue;
                }
                let last = interval.end.min(sequence.bases.len());
                if interval.start <= last {
                    bits[interval.start..=last].fill(true);
                }
            }
            bits
        });
        for (offset, &base) in sequence.bases.iter().enumerate() {
            if let Some(mask) = &mask {
                if !mask[offset + 1] {
                    continue;
                }
            }
            // `StringUtil.toUpperCase(bases)` before the comparison against `SequenceUtil.N`.
            if !base.eq_ignore_ascii_case(&b'N') {
                non_n += 1;
            }
        }
    }

    let mut out = std::fs::File::create(&output)?;
    writeln!(out, "{non_n}")?;
    out.flush()?;
    Ok(())
}
