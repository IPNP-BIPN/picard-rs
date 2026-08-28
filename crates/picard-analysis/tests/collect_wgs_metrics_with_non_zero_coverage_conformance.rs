//! Conformance for `CollectWgsMetricsWithNonZeroCoverage` against Picard 3.4.0.
//!
//! Golden from `tools/nonzerowgs-conformance`. Each case carries the two metric rows, the one
//! histogram table under them, and the chart's version line.
//!
//! # What this suite is for
//!
//!  * **the output being two rows carrying a CATEGORY the parent has no trace of**;
//!  * **the second row's territory being the covered loci**;
//!  * **its mean coverage rising by exactly that ratio**;
//!  * **the exclusion columns NOT moving, their denominator being the bases seen**;
//!  * **a fully covered reference making the two rows identical**;
//!  * **the histogram being one table of two columns**;
//!  * **`--INCLUDE_BQ_HISTOGRAM` adding a third to it**;
//!  * **`--CHART_OUTPUT` being required**;
//!  * **the chart being written even when nothing is covered**;
//!  * **and a file whose reads are all excluded dividing by a territory of nought.**

use std::io::Read;

use picard_analysis::collect_wgs_metrics_with_non_zero_coverage::{
    drop_the_zero_bin, mean_coverage, rows, territory, NON_ZERO_COLUMN, NON_ZERO_REGIONS,
    WHOLE_GENOME, WHOLE_GENOME_COLUMN,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("collect_wgs_metrics_with_non_zero_coverage.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    s.replace("\\t", "\t").replace("\\n", "\n")
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

fn table(text: &str, kind: &str, case: &str) -> Vec<std::collections::HashMap<String, String>> {
    let body = field(text, kind, case).unwrap_or_else(|| panic!("{kind}/{case}"));
    let mut lines = body.lines().filter(|line| !line.is_empty());
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

fn number(row: &std::collections::HashMap<String, String>, name: &str) -> f64 {
    let value = row.get(name).unwrap_or_else(|| panic!("{name}"));
    if value == "?" {
        f64::NAN
    } else {
        value.parse().unwrap_or_else(|_| panic!("{name}={value}"))
    }
}

/// The depth histogram of a case, as the whole-genome column holds it.
fn histogram(text: &str, case: &str) -> Vec<u64> {
    table(text, "histogram", case)
        .iter()
        .map(|row| {
            row.get(WHOLE_GENOME_COLUMN)
                .expect("the whole-genome column")
                .parse()
                .expect("a count")
        })
        .collect()
}

/// The output is two rows, in the order the collector adds them.
#[test]
fn the_output_is_two_rows_carrying_a_category() {
    let text = corpus();
    let metrics = table(&text, "metrics", "partly-covered");
    assert_eq!(metrics.len(), 2);
    assert_eq!(metrics[0]["CATEGORY"], WHOLE_GENOME);
    assert_eq!(metrics[1]["CATEGORY"], NON_ZERO_REGIONS);
    // Which is the order the port writes them in, from one histogram.
    let ours = rows(&histogram(&text, "partly-covered"));
    assert_eq!(ours[0].0, WHOLE_GENOME);
    assert_eq!(ours[1].0, NON_ZERO_REGIONS);
}

/// The second row's territory is the covered loci, and its mean rises by exactly that ratio.
#[test]
fn the_second_row_is_the_covered_loci() {
    let text = corpus();
    for case in [
        "partly-covered",
        "unpaired-counted",
        "deep-uncapped",
        "deep-capped",
    ] {
        let metrics = table(&text, "metrics", case);
        let whole = &metrics[0];
        let non_zero = &metrics[1];
        let bins = histogram(&text, case);
        assert_eq!(
            territory(&bins),
            number(whole, "GENOME_TERRITORY") as u64,
            "{case}"
        );
        assert_eq!(
            territory(&drop_the_zero_bin(&bins)),
            number(non_zero, "GENOME_TERRITORY") as u64,
            "{case}"
        );
        // The same bases over a smaller territory.
        let ratio = number(whole, "GENOME_TERRITORY") / number(non_zero, "GENOME_TERRITORY");
        let rose = number(non_zero, "MEAN_COVERAGE") / number(whole, "MEAN_COVERAGE");
        // The written means carry six fraction digits, so the ratio of two of them is exact only
        // to about that: 2.25 against 2.25000225 on the partly covered file.
        assert!(
            (ratio - rose).abs() < 1e-4,
            "{case}: {ratio} against {rose}"
        );
        // And that is what the port's own arithmetic gives.
        assert!(
            (mean_coverage(&bins) - number(whole, "MEAN_COVERAGE")).abs() < 1e-6,
            "{case}"
        );
        assert!(
            (mean_coverage(&drop_the_zero_bin(&bins)) - number(non_zero, "MEAN_COVERAGE")).abs()
                < 1e-6,
            "{case}"
        );
    }
}

/// The exclusion columns do NOT move: their denominator is the bases seen, which the zeroed bin
/// does not change.
#[test]
fn the_exclusion_columns_do_not_move() {
    let text = corpus();
    for case in [
        "deep-capped",
        "deep-uncapped",
        "fully-covered",
        "all-duplicates",
    ] {
        let metrics = table(&text, "metrics", case);
        for column in [
            "PCT_EXC_TOTAL",
            "PCT_EXC_DUPE",
            "PCT_EXC_MAPQ",
            "PCT_EXC_UNPAIRED",
            "PCT_EXC_BASEQ",
            "PCT_EXC_OVERLAP",
            "PCT_EXC_CAPPED",
            "PCT_EXC_ADAPTER",
        ] {
            assert_eq!(metrics[0][column], metrics[1][column], "{case}/{column}");
        }
    }
    // And they are not all nought: the capped run excludes nine bases in ten on both rows.
    assert_eq!(
        table(&text, "metrics", "deep-capped")[1]["PCT_EXC_TOTAL"],
        "0.9"
    );
}

/// A fully covered reference makes the two rows identical apart from the category.
#[test]
fn a_fully_covered_reference_makes_the_rows_identical() {
    let text = corpus();
    let metrics = table(&text, "metrics", "fully-covered");
    let mut whole = metrics[0].clone();
    let mut non_zero = metrics[1].clone();
    whole.remove("CATEGORY");
    non_zero.remove("CATEGORY");
    assert_eq!(whole, non_zero);
    // Which is what says the second row is a recomputation and not a second traversal.
    let bins = histogram(&text, "fully-covered");
    assert_eq!(bins[0], 0);
    assert_eq!(drop_the_zero_bin(&bins), bins);
}

/// The histogram is one table of two columns, the second's depth-zero cell being nought.
#[test]
fn the_histogram_is_one_table_of_two_columns() {
    let text = corpus();
    let histogram = table(&text, "histogram", "partly-covered");
    assert!(histogram[0].contains_key(WHOLE_GENOME_COLUMN));
    assert!(histogram[0].contains_key(NON_ZERO_COLUMN));
    assert_eq!(histogram[0]["coverage"], "0");
    assert_eq!(histogram[0][NON_ZERO_COLUMN], "0");
    assert_ne!(histogram[0][WHOLE_GENOME_COLUMN], "0");
    // Above the zero bin the two columns agree, cell for cell.
    for row in &histogram[1..] {
        assert_eq!(
            row[WHOLE_GENOME_COLUMN], row[NON_ZERO_COLUMN],
            "at {}",
            row["coverage"]
        );
    }
    // The base-quality histogram is a third column of that same table.
    let with_bq = table(&text, "histogram", "with-bq-histogram");
    assert!(with_bq[0].contains_key("unfiltered_baseq_count"));
    assert_eq!(with_bq[0].len(), 4);
}

/// The chart argument is required, and the chart is written even when nothing is covered.
#[test]
fn the_chart_is_required_and_always_written() {
    let text = corpus();
    assert_eq!(field(&text, "error", "no-chart").as_deref(), Some("exit 1"));
    assert!(field(&text, "metrics", "no-chart").is_none());
    // The emptiness test asks the histogram whether it has bins, and a bin is created for every
    // depth up to the cap, so the "no valid bases" warning is unreachable.
    for case in ["no-reads", "all-duplicates", "all-excluded"] {
        assert_eq!(
            field(&text, "chart", case).as_deref(),
            Some("%PDF-1.4"),
            "{case}"
        );
    }
}

/// A file whose reads are all excluded divides by a territory of nought.
#[test]
fn a_territory_of_nought_writes_a_question_mark() {
    let text = corpus();
    for case in ["all-excluded", "no-reads", "all-duplicates"] {
        let metrics = table(&text, "metrics", case);
        assert_eq!(metrics[1]["GENOME_TERRITORY"], "0", "{case}");
        assert_eq!(metrics[1]["MEAN_COVERAGE"], "?", "{case}");
        assert!(number(&metrics[1], "MEAN_COVERAGE").is_nan(), "{case}");
        // The first row still has its own territory, and a mean of nought rather than a NaN.
        assert_eq!(metrics[0]["GENOME_TERRITORY"], "90", "{case}");
        assert_eq!(metrics[0]["MEAN_COVERAGE"], "0", "{case}");
        // Which is the division the port makes of the same histogram.
        let bins = histogram(&text, case);
        assert!(mean_coverage(&drop_the_zero_bin(&bins)).is_nan(), "{case}");
        assert_eq!(mean_coverage(&bins), 0.0, "{case}");
    }
}
