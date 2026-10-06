//! `BpmToNormalizationManifestCsv` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.illumina.BpmToNormalizationManifestCsv.doWork` at tag 3.4.0. The files
//! are read by `picard_analysis::infinium`, and each row is
//! `picard_analysis::illumina_arrays::normalization_line`. The cluster file is read FIRST, so a
//! run with two bad files names the cluster file; and a locus whose index is past the cluster
//! file's last code throws mid-write, after the header and the rows before it.

use picard_analysis::fingerprint::reference_path;
use picard_analysis::illumina_arrays::{
    normalization_line, NormalizationRow, NORMALIZATION_HEADER,
};
use picard_analysis::infinium::{Bpm, Egt, ReadError};
use picard_analysis::metrics_cli::{thrown, Args};

/// An exception as the JVM prints it: a `PicardException` unless the message names its class.
fn raise(message: &str) -> ! {
    if message.starts_with("java.") {
        thrown(message)
    } else {
        thrown(&format!("picard.PicardException: {message}"))
    }
}

fn main() {
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("CF", "CLUSTER_FILE")]);
    let input = args.required("INPUT");
    let cluster = args.required("CLUSTER_FILE");
    let output = args.required("OUTPUT");
    let read = |path: &str| {
        std::fs::read(path)
            .unwrap_or_else(|e| raise(&format!("java.io.FileNotFoundException: {e}")))
    };
    let egt = match Egt::parse(&read(&cluster)) {
        Ok(egt) => egt,
        Err(ReadError::Io(_)) => raise(&format!(
            "Error reading cluster file '{}'",
            reference_path(&cluster)
        )),
        Err(ReadError::Picard(m)) => raise(&m),
    };
    let bpm = match Bpm::parse(&read(&input)) {
        Ok(bpm) => bpm,
        Err(ReadError::Io(_)) => raise(&format!(
            "Error reading bpm file '{}'",
            reference_path(&input)
        )),
        Err(ReadError::Picard(m)) => raise(&m),
    };
    let mut text = format!("{NORMALIZATION_HEADER}\n");
    for locus in &bpm.loci {
        let Some(score) = egt.total_score.get(locus.index) else {
            let _ = std::fs::write(&output, &text);
            raise(&format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {} out of bounds for length {}",
                locus.index,
                egt.total_score.len()
            ));
        };
        text.push_str(&normalization_line(&NormalizationRow {
            index: locus.index + 1,
            name: locus.name.clone(),
            chromosome: locus.chrom.clone(),
            position: locus.map_info,
            gentrain_score: *score,
            snp: locus.snp.clone(),
            illumina_strand: locus.ilmn_strand.clone(),
            customer_strand: locus.customer_strand.clone(),
            normalization_id: locus.normalization_id,
        }));
        text.push('\n');
    }
    if let Err(e) = std::fs::write(&output, text) {
        raise(&format!(
            "Error writing bpm.csv file '{}': {e}",
            reference_path(&output)
        ));
    }
}
