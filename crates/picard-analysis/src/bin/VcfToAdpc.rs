//! `VcfToAdpc` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.VcfToAdpc.doWork` at tag 3.4.0. The record layout is
//! `picard_analysis::vcf_to_adpc`; this reads each sample of each VCF in turn.
//!
//! Every failure inside the run is caught, logged as the exception's `toString` and answered with
//! exit code 1, and the `.adpc.bin` keeps whatever records were written before it: the writer is
//! closed on the way out. The samples and markers files are written only at the end.

use htsjdk_vcf::variant::Value;
use picard_analysis::fingerprint::reference_path;
use picard_analysis::metrics_cli::Args;
use picard_analysis::vcf_io::{log_error, read_path, unroll_paths};
use picard_analysis::vcf_to_adpc::{
    illumina_genotype, markers_file, samples_file, write_record, IlluminaGenotype, Record, HEADER,
};

fn text(value: &Value) -> String {
    match value {
        Value::Str(s) => s.clone(),
        Value::Int(i) => i.to_string(),
        Value::Double(d) => d.to_string(),
        Value::Bool(b) => b.to_string(),
        Value::List(l) => l.iter().map(text).collect::<Vec<_>>().join(","),
        Value::Missing => ".".to_string(),
    }
}

/// `Float.parseFloat`, with the exception it throws.
fn float(s: &str) -> Result<f32, String> {
    s.trim()
        .parse()
        .map_err(|_| format!("java.lang.NumberFormatException: For input string: \"{s}\""))
}

fn run(vcfs: &[String], out: &mut Vec<u8>) -> Result<(Vec<String>, usize), String> {
    let mut samples: Vec<String> = Vec::new();
    let mut loci: Option<usize> = None;
    for path in vcfs {
        let vcf = read_path(path).map_err(|e| e.to_string())?;
        for (n, sample) in vcf.file.header.samples.iter().enumerate() {
            samples.push(sample.clone());
            let mut count = 0;
            for record in &vcf.records {
                let ctx = &record.variant;
                let info = |key: &str| -> Result<String, String> {
                    ctx.attributes
                        .iter()
                        .find(|(k, _)| k == key)
                        .map(|(_, v)| text(v))
                        .ok_or_else(|| {
                            format!("picard.PicardException: Unable to find attribute {key} in VCF.  Is this an Arrays VCF file?")
                        })
                };
                let gc = float(&info("GC_SCORE")?)?;
                let g = &ctx.genotypes[n];
                let genotype = if g.is_called() {
                    let a = info("ALLELE_A")?;
                    let b = info("ALLELE_B")?;
                    if g.alleles.len() != 2 {
                        return Err(format!("picard.PicardException: Unexpected number of called alleles in variant context {} found alleles", ctx.contig));
                    }
                    let first = g.alleles[0].base_string();
                    let second = g.alleles[1].base_string();
                    illumina_genotype(Some((&first, &second)), &a, &b).ok_or_else(|| {
                        "picard.PicardException: Error matching called alleles to Illumina alleles.".to_string()
                    })?
                } else {
                    IlluminaGenotype::Nn
                };
                let format = |key: &str| g.get(key).map(text);
                let unsigned = |key: &str| -> Result<u16, String> {
                    let v = format(key).ok_or_else(|| format!("picard.PicardException: Unable to find attribute {key} in VCF Genotype field.  Is this an Arrays VCF file?"))?;
                    let n: i32 = v.trim().parse().map_err(|_| {
                        format!("java.lang.NumberFormatException: For input string: \"{v}\"")
                    })?;
                    if n < 0 {
                        return Err(format!("picard.PicardException: Value for key {key} ({n}) is <= 0!  Invalid value for unsigned int"));
                    }
                    Ok(n.min(65535) as u16)
                };
                let optional = |key: &str| -> Result<f32, String> {
                    match format(key) {
                        Some(v) if v != "?" => float(&v),
                        _ => Ok(f32::NAN),
                    }
                };
                let rec = Record {
                    a_intensity: unsigned("X")?,
                    b_intensity: unsigned("Y")?,
                    a_normalized: optional("NORMX")?,
                    b_normalized: optional("NORMY")?,
                    gc_score: gc,
                    genotype,
                };
                out.extend(write_record(&rec));
                count += 1;
            }
            if count == 0 {
                return Err(format!(
                    "picard.PicardException: Found no records in VCF' {}'",
                    reference_path(path)
                ));
            }
            match loci {
                None => loci = Some(count),
                Some(l) if l != count => {
                    return Err(
                        "picard.PicardException: VCFs have differing number of loci".to_string()
                    )
                }
                _ => {}
            }
        }
    }
    Ok((samples, loci.unwrap_or(0)))
}

fn main() {
    let args = Args::from_env(&[("O", "OUTPUT")]);
    let vcfs = unroll_paths(&args.collection("VCF", &[])).unwrap_or_default();
    let output = args.required("OUTPUT");
    let samples_path = args.required("SAMPLES_FILE");
    let markers_path = args.required("NUM_MARKERS_FILE");
    let mut out = HEADER.to_vec();
    let result = run(&vcfs, &mut out);
    let _ = std::fs::write(&output, &out);
    match result {
        Ok((samples, loci)) => {
            let names: Vec<&str> = samples.iter().map(String::as_str).collect();
            let _ = std::fs::write(&samples_path, samples_file(&names));
            let _ = std::fs::write(&markers_path, markers_file(loci));
        }
        Err(e) => {
            log_error("VcfToAdpc", &e);
            std::process::exit(1);
        }
    }
}
