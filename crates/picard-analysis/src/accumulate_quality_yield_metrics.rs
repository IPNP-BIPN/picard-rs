//! `AccumulateQualityYieldMetrics`.
//!
//! Ports `picard.util.AccumulateQualityYieldMetrics.doWork` at tag 3.4.0: read one or more
//! `CollectQualityYieldMetrics` files (each holding a single `QualityYieldMetrics` row, typically from
//! separate shards of one read group), sum them, and write a single combined metrics file.
//!
//! `MergeableMetricBase.merge` adds every `@MergeByAdding` field, which for `QualityYieldMetrics` is
//! all ten counters (`TOTAL_READS`, `PF_READS`, `TOTAL_BASES`, `PF_BASES`, `Q20_BASES`,
//! `PF_Q20_BASES`, `Q30_BASES`, `PF_Q30_BASES`, `Q20_EQUIVALENT_YIELD`, `PF_Q20_EQUIVALENT_YIELD`);
//! `READ_LENGTH` is `@NoMergingIsDerived` and recomputed by `calculateDerivedFields` as
//! `TOTAL_READS == 0 ? 0 : (int)(TOTAL_BASES / TOTAL_READS)`. `useOriginalQualities` is
//! `@MergeByAssertEquals` but is not a written column, so it does not appear in the file.
//!
//! The input files are parsed straight from their `## METRICS CLASS` table (the column-name row then
//! the single value row); the output is written with [`htsjdk_metrics`] reusing the same
//! [`crate::quality_yield::QualityYieldMetrics`] bean. The tool writes a bare `new MetricsFile<>()`
//! (not `getMetricsFile()`), so the output carries no command-line or start-time header comments and
//! is byte-identical with no canonicalization at all.

use htsjdk_metrics::file::MetricsFile;

use crate::quality_yield::QualityYieldMetrics;

/// Why an input metrics file could not be accumulated.
#[derive(Debug, PartialEq, Eq)]
pub enum AccumulateError {
    /// No `## METRICS CLASS` table (with a column row and a value row) was found.
    NoMetricsTable,
    /// A required column was missing from the table.
    MissingColumn(String),
    /// A value cell did not parse as an integer.
    BadValue(String),
}

/// Parses the single `QualityYieldMetrics` row out of a metrics file's `## METRICS CLASS` table.
fn parse_quality_yield(metrics_file: &str) -> Result<QualityYieldMetrics, AccumulateError> {
    let mut lines = metrics_file.lines();
    // Advance to the metrics-class marker; the next two non-empty lines are the columns and values.
    while let Some(line) = lines.next() {
        if line.starts_with("## METRICS CLASS") {
            let columns = lines.next().ok_or(AccumulateError::NoMetricsTable)?;
            let values = lines.next().ok_or(AccumulateError::NoMetricsTable)?;
            let columns: Vec<&str> = columns.split('\t').collect();
            let values: Vec<&str> = values.split('\t').collect();

            let get = |name: &str| -> Result<i64, AccumulateError> {
                let idx = columns
                    .iter()
                    .position(|c| *c == name)
                    .ok_or_else(|| AccumulateError::MissingColumn(name.to_string()))?;
                let cell = values
                    .get(idx)
                    .ok_or_else(|| AccumulateError::MissingColumn(name.to_string()))?;
                cell.parse::<i64>()
                    .map_err(|_| AccumulateError::BadValue(cell.to_string()))
            };

            return Ok(QualityYieldMetrics {
                total_reads: get("TOTAL_READS")?,
                pf_reads: get("PF_READS")?,
                read_length: get("READ_LENGTH")? as i32,
                total_bases: get("TOTAL_BASES")?,
                pf_bases: get("PF_BASES")?,
                q20_bases: get("Q20_BASES")?,
                pf_q20_bases: get("PF_Q20_BASES")?,
                q30_bases: get("Q30_BASES")?,
                pf_q30_bases: get("PF_Q30_BASES")?,
                q20_equivalent_yield: get("Q20_EQUIVALENT_YIELD")?,
                pf_q20_equivalent_yield: get("PF_Q20_EQUIVALENT_YIELD")?,
            });
        }
    }
    Err(AccumulateError::NoMetricsTable)
}

/// `AccumulateQualityYieldMetrics.doWork` over the input metrics files, returning the combined file.
pub fn accumulate_quality_yield_metrics(inputs: &[&str]) -> Result<String, AccumulateError> {
    let mut total = QualityYieldMetrics::default();
    for input in inputs {
        let m = parse_quality_yield(input)?;
        // merge: add every @MergeByAdding counter.
        total.total_reads += m.total_reads;
        total.pf_reads += m.pf_reads;
        total.total_bases += m.total_bases;
        total.pf_bases += m.pf_bases;
        total.q20_bases += m.q20_bases;
        total.pf_q20_bases += m.pf_q20_bases;
        total.q30_bases += m.q30_bases;
        total.pf_q30_bases += m.pf_q30_bases;
        total.q20_equivalent_yield += m.q20_equivalent_yield;
        total.pf_q20_equivalent_yield += m.pf_q20_equivalent_yield;
    }
    // calculateDerivedFields: READ_LENGTH = TOTAL_READS == 0 ? 0 : (int)(TOTAL_BASES / TOTAL_READS).
    total.read_length = if total.total_reads == 0 {
        0
    } else {
        (total.total_bases / total.total_reads) as i32
    };

    // The tool writes a bare `new MetricsFile<>()` (not `getMetricsFile()`), so it adds no command
    // line or start-time header comments: the output is just the metrics table and is fully
    // deterministic, needing no canonicalization.
    let mut file = MetricsFile::new();
    file.add_metric(&total);
    Ok(file.write())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metrics_file(total_reads: i64, total_bases: i64) -> String {
        let read_length = if total_reads == 0 {
            0
        } else {
            total_bases / total_reads
        };
        format!(
            "## htsjdk.samtools.metrics.StringHeader\n# a command line\n\n\
             ## METRICS CLASS\tpicard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics\n\
             TOTAL_READS\tPF_READS\tREAD_LENGTH\tTOTAL_BASES\tPF_BASES\tQ20_BASES\tPF_Q20_BASES\tQ30_BASES\tPF_Q30_BASES\tQ20_EQUIVALENT_YIELD\tPF_Q20_EQUIVALENT_YIELD\n\
             {total_reads}\t{total_reads}\t{read_length}\t{total_bases}\t{total_bases}\t0\t0\t0\t0\t0\t0\n\n"
        )
    }

    #[test]
    fn it_sums_counters_and_recomputes_read_length() {
        let a = metrics_file(100, 15000);
        let b = metrics_file(50, 7600);
        let out = accumulate_quality_yield_metrics(&[&a, &b]).unwrap();
        // Summed: TOTAL_READS 150, TOTAL_BASES 22600; READ_LENGTH = 22600 / 150 = 150.
        let data = out
            .lines()
            .find(|l| l.starts_with("150\t"))
            .expect("data row");
        let cells: Vec<&str> = data.split('\t').collect();
        assert_eq!(cells[0], "150", "TOTAL_READS");
        assert_eq!(cells[2], "150", "READ_LENGTH");
        assert_eq!(cells[3], "22600", "TOTAL_BASES");
    }

    #[test]
    fn a_single_input_reproduces_its_counters() {
        let a = metrics_file(100, 15000);
        let out = accumulate_quality_yield_metrics(&[&a]).unwrap();
        assert!(out.contains("100\t100\t150\t15000\t15000\t"), "{out}");
    }

    #[test]
    fn a_file_without_a_metrics_table_is_an_error() {
        assert_eq!(
            accumulate_quality_yield_metrics(&["## htsjdk.samtools.metrics.StringHeader\n"]),
            Err(AccumulateError::NoMetricsTable)
        );
    }
}
