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
use picard_analysis::infinium::{Bpm, Egt, ReadError, Reader};
use picard_analysis::metrics_cli::{thrown, Args};
use std::io::Write;

/// An exception as the JVM prints it: a `PicardException` unless the message names its class.
fn raise(message: &str) -> ! {
    if message.starts_with("java.") {
        thrown(message)
    } else {
        thrown(&format!("picard.PicardException: {message}"))
    }
}

fn main() {
    picard_analysis::argument_check::check("BpmToNormalizationManifestCsv");
    picard_analysis::output::install();
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("CF", "CLUSTER_FILE")]);
    let input = args.required("INPUT");
    let cluster = args.required("CLUSTER_FILE");
    let output = args.required("OUTPUT");
    let open = |path: &str| {
        Reader::open(path).unwrap_or_else(|e| raise(&format!("java.io.FileNotFoundException: {e}")))
    };
    let egt = match Egt::parse(&mut open(&cluster)) {
        Ok(egt) => egt,
        Err(ReadError::Io(_)) => raise(&format!(
            "Error reading cluster file '{}'",
            reference_path(&cluster)
        )),
        Err(ReadError::Picard(m)) => raise(&m),
    };
    let mut reader = open(&input);
    let bpm_error = |e: ReadError| -> ! {
        match e {
            ReadError::Io(_) => raise(&format!(
                "Error reading bpm file '{}'",
                reference_path(&input)
            )),
            ReadError::Picard(m) => raise(&m),
        }
    };
    let bpm = Bpm::parse(&mut reader).unwrap_or_else(|e| bpm_error(e));
    // A file that cannot be made is refused once every row is checked, as the whole text was
    // written at the end.
    let mut out =
        std::fs::File::create(&output).map(|file| std::io::BufWriter::with_capacity(1 << 16, file));
    let mut put = |text: &str| {
        if let Ok(w) = &mut out {
            if let Err(e) = w.write_all(text.as_bytes()) {
                out = Err(e);
            }
        }
    };
    put(&format!("{NORMALIZATION_HEADER}\n"));
    for i in 0..bpm.len() {
        let locus = bpm.locus(&mut reader, i).unwrap_or_else(|e| bpm_error(e));
        let Some(score) = egt.total_score.get(locus.index) else {
            if let Ok(w) = &mut out {
                let _ = w.flush();
            }
            raise(&format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {} out of bounds for length {}",
                locus.index,
                egt.total_score.len()
            ));
        };
        let mut line = normalization_line(&NormalizationRow {
            index: locus.index + 1,
            name: locus.name,
            chromosome: locus.chrom,
            position: locus.map_info,
            gentrain_score: *score,
            snp: locus.snp,
            illumina_strand: locus.ilmn_strand,
            customer_strand: locus.customer_strand,
            normalization_id: locus.normalization_id,
        });
        line.push('\n');
        put(&line);
    }
    if let Err(e) = out.and_then(|w| w.into_inner().map_err(|e| e.into_error())) {
        raise(&format!(
            "Error writing bpm.csv file '{}': {e}",
            reference_path(&output)
        ));
    }
}
