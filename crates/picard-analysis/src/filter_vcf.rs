//! `FilterVcf`.
//!
//! Ports `picard.vcf.filter.FilterVcf.doWork` and `FilterApplyingVariantIterator` at tag 3.4.0,
//! with the three site filters (`AlleleBalanceFilter`, `FisherStrandFilter`, `QdFilter`) and the
//! two genotype filters (`GenotypeQualityFilter`, `DepthFilter`) it builds.
//!
//! # Every record is rebuilt
//!
//! The iterator replaces each record's filters (`passFilters()` or `filters(set)`, so an input
//! filter such as `q10` is dropped and `.` becomes `PASS`) and rebuilds every genotype with a
//! filter of its own, `PASS` when none fired. A genotype filter set to `PASS` is no filter at all
//! (`GenotypeBuilder.filter` stores null), so a record whose genotypes all pass has no FT key and
//! the column is not written, while one filtered genotype brings FT into every column.
//!
//! # A missing GQ or DP is -1, and -1 is below a threshold of 0
//!
//! `Genotype.getGQ()` and `getDP()` return -1 when the field is absent, and the filters compare
//! with `<`, so at the default thresholds of 0 a genotype without GQ is `LowGQ` and one without DP
//! is `LowDP`. A no-call usually has neither and is filtered twice.

use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};

/// The thresholds `FilterVcf` passes to its filters.
#[derive(Debug, Clone, Copy)]
pub struct Thresholds {
    pub min_ab: f64,
    pub min_dp: i32,
    pub min_gq: i32,
    pub max_fs: f64,
    pub min_qd: f64,
}

impl Default for Thresholds {
    fn default() -> Self {
        Self {
            min_ab: 0.0,
            min_dp: 0,
            min_gq: 0,
            max_fs: f64::MAX,
            min_qd: 0.0,
        }
    }
}

/// `CommonInfo.getAttributeAsDouble`: the default when absent, the number as it is, a string
/// parsed. `None` stands for the `NumberFormatException` a string that is not a number throws.
fn attribute_as_double(vc: &VariantContext, key: &str, default: f64) -> Option<f64> {
    match vc.attributes.iter().find(|(k, _)| k == key).map(|(_, v)| v) {
        None => Some(default),
        Some(Value::Double(d)) => Some(*d),
        Some(Value::Int(i)) => Some(*i as f64),
        Some(Value::Str(s)) => s.trim().parse::<f64>().ok(),
        Some(_) => None,
    }
}

/// `AlleleBalanceFilter.filter`: per het allele pair, the smaller AD share over the pair's total.
fn allele_balance(vc: &VariantContext, min_ab: f64) -> Option<&'static str> {
    if !vc.genotypes.iter().any(Genotype::is_het) {
        return None;
    }
    // `HashMap<List<Allele>, Counts>`, visited in an order that cannot matter: any pair below the
    // limit filters the site.
    let mut counts: Vec<(Vec<_>, i64, i64)> = Vec::new();
    for gt in vc.genotypes.iter() {
        if gt.is_no_call() || !gt.is_het() {
            continue;
        }
        let Some(ad) = &gt.ad else { continue };
        let index = |allele| vc.alleles.iter().position(|a| a == allele);
        let (Some(first), Some(second)) = (index(&gt.alleles[0]), index(&gt.alleles[1])) else {
            continue;
        };
        let key = gt.alleles.clone();
        let at = match counts.iter().position(|(k, _, _)| *k == key) {
            Some(at) => at,
            None => {
                counts.push((key, 0, 0));
                counts.len() - 1
            }
        };
        counts[at].1 += i64::from(ad.get(first).copied().unwrap_or(0));
        counts[at].2 += i64::from(ad.get(second).copied().unwrap_or(0));
    }
    counts
        .iter()
        .any(|&(_, one, two)| {
            let total = one + two;
            total > 0 && (one.min(two) as f64 / total as f64) < min_ab
        })
        .then_some("AlleleBalance")
}

/// `FilterApplyingVariantIterator.next`, with the filters `FilterVcf` builds.
pub fn filter_record(vc: &mut VariantContext, t: &Thresholds) -> Result<(), String> {
    let mut site: Vec<String> = Vec::new();
    if let Some(name) = allele_balance(vc, t.min_ab) {
        site.push(name.to_string());
    }
    let fs = attribute_as_double(vc, "FS", 0.0).ok_or("FS is not a number")?;
    if fs > t.max_fs {
        site.push("StrandBias".to_string());
    }
    let qd = attribute_as_double(vc, "QD", -1.0).ok_or("QD is not a number")?;
    if qd >= 0.0 && qd < t.min_qd {
        site.push("LowQD".to_string());
    }

    let mut variant_samples = 0usize;
    let mut all_filtered = true;
    for gt in vc.genotypes.iter_mut() {
        let mut filters: Vec<&str> = Vec::new();
        if gt.gq.unwrap_or(-1) < t.min_gq {
            filters.push("LowGQ");
        }
        if gt.dp.unwrap_or(-1) < t.min_dp {
            filters.push("LowDP");
        }
        if gt.is_called() && !gt.is_hom_ref() {
            variant_samples += 1;
            all_filtered &= !filters.is_empty();
        }
        // `GenotypeBuilder.filters(List)`: one filter as it is, several sorted and joined.
        filters.sort_unstable();
        gt.filters = (!filters.is_empty()).then(|| filters.join(";"));
    }
    if variant_samples > 0 && all_filtered {
        site.push("AllGtsFiltered".to_string());
    }
    site.sort_unstable();
    site.dedup();
    vc.filters = Some(site);
    Ok(())
}

/// `VCFHeader.addMetaDataLine` for each line `FilterVcf` adds: a line whose key and ID the header
/// already has is not added again, so an input's own `FT` description is the one written.
pub fn add_header_lines(header: &mut VcfHeader) {
    let filters = [
        (
            "AllGtsFiltered",
            "Site filtered out because all genotypes are filtered out.",
        ),
        (
            "AlleleBalance",
            "Heterozygote allele balance below required threshold.",
        ),
        (
            "StrandBias",
            "Site exhibits excessive allele/strand correlation.",
        ),
        ("LowQD", "Site exhibits QD value below a hard limit."),
    ];
    let has_format_ft = header.lines.iter().any(
        |line| matches!(line, HeaderLine::Compound { key, id, .. } if key == "FORMAT" && id == "FT"),
    );
    if !has_format_ft {
        header.lines.push(HeaderLine::format(
            "FT",
            Cardinality::Unbounded,
            LineType::String,
            "Genotype filters.",
        ));
    }
    for (id, description) in filters {
        let present = header
            .lines
            .iter()
            .any(|line| matches!(line, HeaderLine::Filter { id: have, .. } if have == id));
        if !present {
            header.lines.push(HeaderLine::filter(id, description));
        }
    }
}
