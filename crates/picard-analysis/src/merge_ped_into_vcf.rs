//! `MergePedIntoVcf`: a zCall PED merged into a genotyping-array VCF.
//!
//! Reading and writing the VCF are not ported. What is ported is the merge: the PED and MAP read
//! in step, the thresholds table, the allele translation, and the two genotype fields.
//!
//! Ported from `picard.arrays.MergePedIntoVcf` and `picard.arrays.illumina.ZCallPedFile` in
//! Picard 3.4.0.

use std::collections::BTreeMap;

/// `ZCallPedFile.OFFSET`: the six PED fields before the alleles, which are ignored whatever they
/// hold.
pub const PED_OFFSET: usize = 6;
/// `MergePedIntoVcf.MISSING_VALUE_v4`, written where a threshold pair is `NA`.
pub const MISSING_VALUE: &str = ".";
/// `parseZCallThresholds`, where one of a pair is `NA` and the other is not.
pub const HALF_NA_MESSAGE: &str = "Thresholds should either both exist or both not exist.";
/// `ZCallPedFile.fromFile`, on a PED of more than one sample.
pub const MULTI_SAMPLE_PED_MESSAGE: &str = "Only single-sample .ped files are supported.";
/// The same, on an allele of more than one character.
pub const LONG_ALLELE_MESSAGE: &str = "Malformed file: each allele should be a single character.";
/// `doWork`, on a VCF of more than one sample.
pub const MULTI_SAMPLE_VCF_MESSAGE: &str = "MergePedIntoVCF only works with single-sample VCFs.";

/// Why a merge failed.
///
/// Every one of these is caught by `doWork`, printed as a stack trace and followed by a return of
/// ZERO, so the tool reports success and writes no file. The port answers the error; the caller is
/// what decides to report it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MergeError {
    HalfNa,
    MultiSamplePed,
    LongAllele,
    MultiSampleVcf,
    /// The genotype carried no extended attributes, so the map the tool writes into is immutable.
    NoExtendedAttributes,
    /// The looked-up allele is built as non-reference and the context does not hold it.
    AlleleNotInContext(String),
}

impl MergeError {
    pub fn message(&self) -> String {
        match self {
            MergeError::HalfNa => HALF_NA_MESSAGE.to_string(),
            MergeError::MultiSamplePed => MULTI_SAMPLE_PED_MESSAGE.to_string(),
            MergeError::LongAllele => LONG_ALLELE_MESSAGE.to_string(),
            MergeError::MultiSampleVcf => MULTI_SAMPLE_VCF_MESSAGE.to_string(),
            MergeError::NoExtendedAttributes => {
                "java.lang.UnsupportedOperationException".to_string()
            }
            MergeError::AlleleNotInContext(allele) => {
                format!("Allele in genotype {allele} not in the variant context")
            }
        }
    }
}

/// `parseZCallThresholds`: one line per SNP, three tab-separated fields.
///
/// A pair of `NA` becomes the missing value rather than being dropped; ONE `NA` of a pair is
/// refused. The map is keyed by the SNP's ID, and in the reference it is STATIC, so a second run
/// in one process still holds the first's.
pub fn parse_thresholds(text: &str) -> Result<BTreeMap<String, [String; 2]>, MergeError> {
    let mut thresholds = BTreeMap::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let (x, y) = (fields[1], fields[2]);
        if x == "NA" || y == "NA" {
            if x != "NA" || y != "NA" {
                return Err(MergeError::HalfNa);
            }
            thresholds.insert(
                fields[0].to_string(),
                [MISSING_VALUE.to_string(), MISSING_VALUE.to_string()],
            );
        } else {
            thresholds.insert(fields[0].to_string(), [x.to_string(), y.to_string()]);
        }
    }
    Ok(thresholds)
}

/// `ZCallPedFile.fromFile`: the PED's alleles keyed by the MAP's SNP names, read in step.
///
/// The first six PED fields are ignored whatever they hold, and the pairing is by INDEX: the
/// map's nth line names the PED's nth allele pair. Each allele must be one character.
pub fn parse_ped(ped: &str, map: &str) -> Result<BTreeMap<String, String>, MergeError> {
    if ped.lines().filter(|line| !line.is_empty()).count() > 1 {
        return Err(MergeError::MultiSamplePed);
    }
    let fields: Vec<&str> = ped.split_whitespace().collect();
    let names: Vec<&str> = map
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| line.split_whitespace().nth(1).unwrap_or(""))
        .collect();
    let mut alleles = BTreeMap::new();
    for (i, name) in names.iter().enumerate() {
        let index = i * 2 + PED_OFFSET;
        let (first, second) = (fields[index], fields[index + 1]);
        if first.len() != 1 || second.len() != 1 {
            return Err(MergeError::LongAllele);
        }
        alleles.insert((*name).to_string(), format!("{first}{second}"));
    }
    Ok(alleles)
}

/// `translateAllele`: `A` and `B` are looked up in the record's own two, and `0` is a no-call.
///
/// The looked-up allele is always built as NON-reference, so a PED calling `A` where the record's
/// `ALLELE_A` is its `REF` makes an allele the context does not hold.
pub fn translate_allele(
    allele: char,
    allele_a: &str,
    allele_b: &str,
) -> Result<Option<String>, String> {
    match allele {
        'A' => Ok(Some(allele_a.to_string())),
        'B' => Ok(Some(allele_b.to_string())),
        '0' => Ok(None),
        other => Err(format!("Illegal allele: {other}")),
    }
}

/// The two alleles a PED pair translates to, and whether the context can hold them.
pub fn zcall_alleles(
    pair: &str,
    allele_a: &str,
    allele_b: &str,
    reference: &str,
) -> Result<[Option<String>; 2], MergeError> {
    let mut out: [Option<String>; 2] = [None, None];
    for (i, allele) in pair.chars().take(2).enumerate() {
        let translated =
            translate_allele(allele, allele_a, allele_b).map_err(MergeError::AlleleNotInContext)?;
        if let Some(base) = &translated {
            // Built as non-reference: naming the reference base is what the context refuses.
            if base == reference {
                return Err(MergeError::AlleleNotInContext(base.clone()));
            }
        }
        out[i] = translated;
    }
    Ok(out)
}

/// The `GT` string a pair of alleles writes, a missing one being `.`.
pub fn genotype_string(alleles: &[Option<String>; 2], alternate: &str) -> String {
    let code = |allele: &Option<String>| match allele {
        None => ".".to_string(),
        Some(base) if base == alternate => "1".to_string(),
        Some(_) => "0".to_string(),
    };
    format!("{}/{}", code(&alleles[0]), code(&alleles[1]))
}
