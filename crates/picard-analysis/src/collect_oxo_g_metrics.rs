//! `CollectOxoGMetrics`: the eight counters behind the 8-oxo-G argument, and the rates they give.
//!
//! Walking the loci is not ported. What is ported is which read reaches which counter, which
//! context a site is filed under, the order the rows come out in, and the arithmetic of the
//! derived columns.
//!
//! Ported from `picard.analysis.CollectOxoGMetrics` in Picard 3.4.0.

use crate::make_vcf_sample_name_map::hash_map_order;

/// `MINIMUM_QUALITY_SCORE`, which drops single bases.
pub const DEFAULT_MINIMUM_QUALITY_SCORE: i32 = 20;
/// `MINIMUM_MAPPING_QUALITY`, which drops whole reads.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 30;
/// The insert-size window, whose ends both drop whole pairs.
pub const DEFAULT_MINIMUM_INSERT_SIZE: i32 = 60;
pub const DEFAULT_MAXIMUM_INSERT_SIZE: i32 = 600;
/// `CONTEXT_SIZE`, the bases on EACH side of the assayed one.
pub const DEFAULT_CONTEXT_SIZE: usize = 1;

/// The floor under the two reference-bias rates, which caps their Q at a hundred.
pub const MINIMUM_REFERENCE_BIAS_RATE: f64 = 1e-10;

/// `makeContextStrings`: every kmer of the right length whose middle base is a `C`.
///
/// A `G` site does not get contexts of its own: it is folded into the reverse complement of its
/// own, which is why sixteen rows cover both halves of the genome at a context size of one.
pub fn contexts(context_size: usize) -> Vec<String> {
    let width = 1 + 2 * context_size;
    let mut out = Vec::new();
    let bases = *b"ACGT";
    let mut kmer = vec![b'A'; width];
    loop {
        if kmer[context_size] == b'C' {
            out.push(String::from_utf8(kmer.clone()).expect("ascii"));
        }
        let mut position = width;
        loop {
            if position == 0 {
                out.sort();
                return out;
            }
            position -= 1;
            let index = bases
                .iter()
                .position(|b| *b == kmer[position])
                .expect("a base");
            if index + 1 < bases.len() {
                kmer[position] = bases[index + 1];
                break;
            }
            kmer[position] = bases[0];
        }
    }
}

/// The order the rows reach the file in: a `HashMap`'s over the context strings, one entry per
/// library per context, the libraries of one context together.
///
/// The contexts go into a `HashSet` and then into a `ListMap`, and both are default-sized maps
/// over the same keys, so the second's order is the first's and the printed order is the table's:
/// sixteen contexts sit in a table of thirty-two, and `CCA` comes out first.
pub fn row_order(contexts: &[String], libraries: &[String]) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for context in hash_map_order(contexts, 16) {
        for library in libraries {
            out.push((context.clone(), library.clone()));
        }
    }
    out
}

/// `SequenceUtil.reverseComplement`.
pub fn reverse_complement(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            other => other as char,
        })
        .collect()
}

/// The context a site is filed under, or `None` if it is not assayed at all.
///
/// A site is skipped when it is within the context of either end of the contig, and when its base
/// is neither a `C` nor a `G`. A `G` is filed under the reverse complement of its own context,
/// which is what puts the two strands of one context in one row.
pub fn context_at(reference: &[u8], position: usize, context_size: usize) -> Option<String> {
    if position <= context_size || position > reference.len() - context_size {
        return None;
    }
    let index = position - 1;
    let base = reference[index].to_ascii_uppercase();
    if base != b'C' && base != b'G' {
        return None;
    }
    let window = String::from_utf8(
        reference[index - context_size..=index + context_size]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect(),
    )
    .expect("ascii");
    Some(if base == b'C' {
        window
    } else {
        reverse_complement(&window)
    })
}

/// `customCommandLineValidation`, whose two messages name the context they refuse.
pub fn validate_context(context: &str, context_size: usize) -> Result<(), String> {
    let width = 1 + 2 * context_size;
    if context.len() != width {
        return Err(format!(
            "Context {context} is not {width} long as implied by CONTEXT_SIZE={context_size}"
        ));
    }
    if context.as_bytes()[context.len() / 2] != b'C' {
        return Err(format!(
            "Middle base of context sequence {context} must be C"
        ));
    }
    Ok(())
}

/// The eight counters one row keeps.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Counts {
    pub ref_c_control_a: i64,
    pub ref_c_oxidated_a: i64,
    pub ref_c_control_c: i64,
    pub ref_c_oxidated_c: i64,
    pub ref_g_control_a: i64,
    pub ref_g_oxidated_a: i64,
    pub ref_g_control_c: i64,
    pub ref_g_oxidated_c: i64,
}

/// `computeAlleleFraction`, on one read at one site.
///
/// The read's own orientation is undone first, and then the question is which END of the pair
/// carries the base: read one carrying a `G` as read, or read two carrying a `C`, is the oxidised
/// state on the reference side, and the same pairing with the alternate base is the oxidised state
/// on the alternate side. Everything else is the control. A base that is neither the reference nor
/// the one alternate the tool looks for is counted nowhere at all.
pub fn accept(
    counts: &mut Counts,
    reference_base: u8,
    read_base: u8,
    read_number: u8,
    negative_strand: bool,
) {
    let reference_base = reference_base.to_ascii_uppercase();
    let read_base = read_base.to_ascii_uppercase();
    let alternate = if reference_base == b'C' { b'A' } else { b'T' };
    let as_read = if negative_strand {
        match read_base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        }
    } else {
        read_base
    };
    let on_c = reference_base == b'C';
    if read_base == reference_base {
        match (as_read, read_number) {
            (b'G', 1) | (b'C', 2) => {
                if on_c {
                    counts.ref_c_oxidated_c += 1;
                } else {
                    counts.ref_g_oxidated_c += 1;
                }
            }
            (b'G', 2) | (b'C', 1) => {
                if on_c {
                    counts.ref_c_control_c += 1;
                } else {
                    counts.ref_g_control_c += 1;
                }
            }
            _ => {}
        }
    } else if read_base == alternate {
        match (as_read, read_number) {
            (b'T', 1) | (b'A', 2) => {
                if on_c {
                    counts.ref_c_oxidated_a += 1;
                } else {
                    counts.ref_g_oxidated_a += 1;
                }
            }
            (b'T', 2) | (b'A', 1) => {
                if on_c {
                    counts.ref_c_control_a += 1;
                } else {
                    counts.ref_g_control_a += 1;
                }
            }
            _ => {}
        }
    }
}

/// The columns `finish` derives from the counters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Metrics {
    pub total_bases: i64,
    pub ref_oxo_bases: i64,
    pub ref_nonoxo_bases: i64,
    pub ref_total_bases: i64,
    pub alt_nonoxo_bases: i64,
    pub alt_oxo_bases: i64,
    pub oxidation_error_rate: f64,
    pub oxidation_q: f64,
    pub c_ref_ref_bases: i64,
    pub g_ref_ref_bases: i64,
    pub c_ref_alt_bases: i64,
    pub g_ref_alt_bases: i64,
    pub c_ref_oxo_error_rate: f64,
    pub c_ref_oxo_q: f64,
    pub g_ref_oxo_error_rate: f64,
    pub g_ref_oxo_q: f64,
}

/// `Calculator.finish`.
///
/// The oxidation rate takes the oxidised count LESS the control one, on the argument that damage
/// from other causes falls evenly on the two ends, and then floors the difference at ONE BASE
/// rather than at nought: a library with no alternate at all therefore reports one over its total
/// and a finite Q instead of nothing. The two reference-bias rates are floored at 1e-10, which
/// caps their Q at a hundred, and a context no read covered divides nought by nought.
pub fn finish(counts: &Counts) -> Metrics {
    let total_bases = counts.ref_c_control_c
        + counts.ref_c_oxidated_c
        + counts.ref_c_control_a
        + counts.ref_c_oxidated_a
        + counts.ref_g_control_c
        + counts.ref_g_oxidated_c
        + counts.ref_g_control_a
        + counts.ref_g_oxidated_a;
    let ref_oxo_bases = counts.ref_c_oxidated_c + counts.ref_g_oxidated_c;
    let ref_nonoxo_bases = counts.ref_c_control_c + counts.ref_g_control_c;
    let alt_nonoxo_bases = counts.ref_c_control_a + counts.ref_g_control_a;
    let alt_oxo_bases = counts.ref_c_oxidated_a + counts.ref_g_oxidated_a;
    let oxidation_error_rate =
        std::cmp::max(alt_oxo_bases - alt_nonoxo_bases, 1) as f64 / total_bases as f64;
    let c_ref_ref_bases = counts.ref_c_control_c + counts.ref_c_oxidated_c;
    let g_ref_ref_bases = counts.ref_g_control_c + counts.ref_g_oxidated_c;
    let c_ref_alt_bases = counts.ref_c_control_a + counts.ref_c_oxidated_a;
    let g_ref_alt_bases = counts.ref_g_control_a + counts.ref_g_oxidated_a;
    let c_rate = c_ref_alt_bases as f64 / (c_ref_alt_bases + c_ref_ref_bases) as f64;
    let g_rate = g_ref_alt_bases as f64 / (g_ref_alt_bases + g_ref_ref_bases) as f64;
    let c_ref_oxo_error_rate = (c_rate - g_rate).max(MINIMUM_REFERENCE_BIAS_RATE);
    let g_ref_oxo_error_rate = (g_rate - c_rate).max(MINIMUM_REFERENCE_BIAS_RATE);
    Metrics {
        total_bases,
        ref_oxo_bases,
        ref_nonoxo_bases,
        ref_total_bases: ref_oxo_bases + ref_nonoxo_bases,
        alt_nonoxo_bases,
        alt_oxo_bases,
        oxidation_error_rate,
        oxidation_q: -10.0 * oxidation_error_rate.log10(),
        c_ref_ref_bases,
        g_ref_ref_bases,
        c_ref_alt_bases,
        g_ref_alt_bases,
        c_ref_oxo_error_rate,
        c_ref_oxo_q: -10.0 * c_ref_oxo_error_rate.log10(),
        g_ref_oxo_error_rate,
        g_ref_oxo_q: -10.0 * g_ref_oxo_error_rate.log10(),
    }
}
