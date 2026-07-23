//! `SplitSamByNumberOfReads`.
//!
//! Ports `picard.sam.SplitSamByNumberOfReads.doWork` at tag 3.4.0: split one SAM/BAM into N shards
//! by read count, keeping every record verbatim under the input header (no `@PG`, no re-sort, output
//! written `presorted=true`). A shard boundary is only taken when the running count has reached the
//! per-file target **and** the current read name differs from the previous one, so a queryname group
//! is never split across shards (which is why a shard can hold up to `readsPerFile + (largest
//! same-queryname run - 1)` reads).
//!
//! Each shard's header is the input header unchanged and the tool adds nothing, so every shard SAM is
//! comparable raw; for BAM output byte-identity holds under `USE_JDK_DEFLATER=true` via `BamWriter`.
//! Scope: the split itself. `OUT_PREFIX`/extension/`OUTPUT` directory naming is a CLI concern (the
//! shards come back in order); `REFERENCE_SEQUENCE` and CRAM are out of scope.

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::{ParseError, ValidationStringency};

/// The count knobs. Exactly one of `split_to_n_files` / `split_to_n_reads` is set (`> 0`), mirroring
/// the mutex; `total_reads_in_input` overrides the counted total when `> 0`.
#[derive(Debug, Clone, Default)]
pub struct SplitOptions {
    pub split_to_n_files: i64,
    pub split_to_n_reads: i64,
    pub total_reads_in_input: i64,
}

/// Why `SplitSamByNumberOfReads` could not run.
#[derive(Debug)]
pub enum SplitError {
    Parse(ParseError),
    /// `customCommandLineValidation`: neither `SPLIT_TO_N_FILES` nor `SPLIT_TO_N_READS` exceeds 1.
    NoSplitTarget,
    /// `TOTAL_READS_IN_INPUT` was negative.
    NegativeTotalReads,
}

impl From<ParseError> for SplitError {
    fn from(e: ParseError) -> Self {
        SplitError::Parse(e)
    }
}

/// `ceil(a / b)` computed as htsjdk does it: `(int) Math.ceil(a / (double) b)`.
fn ceil_div(a: i64, b: i64) -> i64 {
    (a as f64 / b as f64).ceil() as i64
}

/// `SplitSamByNumberOfReads.doWork`: the shards in order, each a whole SAM (input header verbatim +
/// its records).
pub fn split_sam_by_number_of_reads(
    input_sam: &str,
    opts: &SplitOptions,
) -> Result<Vec<String>, SplitError> {
    if opts.total_reads_in_input < 0 {
        return Err(SplitError::NegativeTotalReads);
    }
    if opts.split_to_n_files <= 1 && opts.split_to_n_reads <= 1 {
        return Err(SplitError::NoSplitTarget);
    }

    let (header, records) = read_sam_with(input_sam, ValidationStringency::Lenient)?;

    let total_reads = if opts.total_reads_in_input == 0 {
        records.len() as i64
    } else {
        opts.total_reads_in_input
    };
    let split_to_n_files = if opts.split_to_n_files != 0 {
        opts.split_to_n_files
    } else {
        ceil_div(total_reads, opts.split_to_n_reads)
    };
    let reads_per_file = ceil_div(total_reads, split_to_n_files);

    let mut shards: Vec<Vec<&BamRecord>> = Vec::new();
    let mut current: Vec<&BamRecord> = Vec::new();
    let mut reads_written: i64 = 0;
    let mut last_read_name = String::new();

    for rec in &records {
        if reads_written >= reads_per_file && last_read_name != rec.read_name {
            shards.push(std::mem::take(&mut current));
            reads_written = 0;
        }
        current.push(rec);
        last_read_name = rec.read_name.clone();
        reads_written += 1;
    }
    shards.push(current);

    Ok(shards
        .iter()
        .map(|recs| render_shard(&header, recs))
        .collect())
}

/// One shard as SAM text: the input header verbatim, then its records.
fn render_shard(header: &SamHeader, recs: &[&BamRecord]) -> String {
    let owned: Vec<BamRecord> = recs.iter().map(|r| (*r).clone()).collect();
    write_sam(header, &owned).expect("records re-encode as SAM")
}

#[cfg(test)]
mod tests {
    use super::*;

    const H: &str = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n";

    fn rec(name: &str, start: i32) -> String {
        format!("{name}\t0\tchr1\t{start}\t60\t4M\t*\t0\t0\tACGT\tIIII\n")
    }

    #[test]
    fn splits_into_the_requested_number_of_files() {
        let input = format!(
            "{H}{}{}{}{}",
            rec("a", 1),
            rec("b", 2),
            rec("c", 3),
            rec("d", 4)
        );
        let shards = split_sam_by_number_of_reads(
            &input,
            &SplitOptions {
                split_to_n_files: 2,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(shards.len(), 2);
        // 4 reads / 2 files = 2 per shard.
        for s in &shards {
            assert_eq!(s.matches("\tchr1\t").count(), 2);
            assert!(s.starts_with("@HD\tVN:1.6\tSO:queryname"));
        }
    }

    #[test]
    fn a_queryname_group_straddling_the_boundary_stays_together() {
        // a, b, b, c with N_READS=2: reads_per_file=2, but the second "b" is not a name change, so
        // shard 1 keeps all three and shard 2 gets just "c". (Picard rejects N_READS=1, so 2 is the
        // smallest valid target that still lets a group straddle.)
        let input = format!(
            "{H}{}{}{}{}",
            rec("a", 1),
            rec("b", 2),
            rec("b", 3),
            rec("c", 4)
        );
        let shards = split_sam_by_number_of_reads(
            &input,
            &SplitOptions {
                split_to_n_reads: 2,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(shards.len(), 2);
        assert_eq!(shards[0].matches("\tchr1\t").count(), 3);
        assert_eq!(shards[1].matches("\tchr1\t").count(), 1);
    }

    #[test]
    fn no_split_target_is_rejected() {
        let input = format!("{H}{}", rec("a", 1));
        assert!(matches!(
            split_sam_by_number_of_reads(&input, &SplitOptions::default()),
            Err(SplitError::NoSplitTarget)
        ));
    }
}
