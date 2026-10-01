//! Port of `picard.vcf.SortVcf` (Picard 3.4.0): the merged header and the sort.
//!
//! # The header goes through `VCFUtils.smartMergeHeaders` even for one input
//!
//! The output header is `new VCFHeader(smartMergeHeaders(inputHeaders, false), sampleList)`, and
//! the merge is not the identity on a single header. It walks the header's lines in SORTED order
//! into a map keyed by the line's key, plus `-<ID>` for a line that has an ID, and keeps the first
//! line under each key. So two unstructured lines with the same key (`##source=a` and
//! `##source=b`) are one line in the output -- the one that sorts first -- and two compound lines
//! with the same ID but different counts come out as the first, its count promoted to `.`.
//!
//! # The samples are written in sorted order
//!
//! `sampleList` is `getSampleNamesInOrder()`, which is the SORTED names, so a file whose columns
//! were not in sorted order comes out with them rearranged. The codec decoded that file's genotypes
//! when it read them (its names were not "already sorted"), so the writer encodes each column from
//! the parsed genotypes, in the header's new order.
//!
//! # The sort is by contig index and start, and stable
//!
//! `VariantContextComparator` compares the contig's index among the header's contig lines, then
//! the start, and nothing else. With fewer records than `MAX_RECORDS_IN_RAM` the collection never
//! spills, so the records are sorted in memory by `Arrays.parallelSort`, which is stable: two
//! records at one position keep their input order, and each keeps its still-lazy genotype block.

use htsjdk_vcf::header::{Cardinality, HeaderLine};

/// A line's key in `smartMergeHeaders`' map: `getKey()`, plus `-<ID>` for a `VCFIDHeaderLine`.
fn merge_key(line: &HeaderLine) -> String {
    let id = match line {
        HeaderLine::Unstructured { .. } => None,
        HeaderLine::Compound { id, .. } | HeaderLine::Filter { id, .. } => Some(id.clone()),
        HeaderLine::Contig { fields, .. } | HeaderLine::Structured { fields, .. } => fields
            .iter()
            .find(|(key, _)| key == "ID")
            .map(|(_, value)| value.clone()),
    };
    match id {
        Some(id) => format!("{}-{id}", line.key()),
        None => line.key().to_string(),
    }
}

/// The order `getMetaDataInSortedOrder` walks: contigs by index among themselves, everything by
/// its rendered string otherwise. Only the order WITHIN one merge key decides anything, and two
/// lines under one key are never a contig and a non-contig.
fn sorted_position(line: &HeaderLine) -> (String, i32) {
    match line {
        HeaderLine::Contig { index, .. } => (String::from("contig="), *index),
        other => (other.render(), 0),
    }
}

/// `smartMergeHeaders` over one header, returning its lines in their original order minus the
/// ones the merge drops, or the `IllegalStateException` message it throws.
pub fn smart_merge_header_lines(lines: &[HeaderLine]) -> Result<Vec<HeaderLine>, String> {
    let mut order: Vec<usize> = (0..lines.len()).collect();
    order.sort_by(|&a, &b| sorted_position(&lines[a]).cmp(&sorted_position(&lines[b])));

    // key -> the index of the line kept for it, in the map's (sorted) insertion order.
    let mut kept: Vec<(String, usize)> = Vec::new();
    let mut merged: Vec<HeaderLine> = lines.to_vec();
    for index in order {
        let line = &lines[index];
        let key = merge_key(line);
        let Some(&(_, other_index)) = kept.iter().find(|(k, _)| *k == key) else {
            kept.push((key, index));
            continue;
        };
        let other = &merged[other_index];
        if std::mem::discriminant(line) != std::mem::discriminant(other) {
            return Err(format!(
                "Incompatible header types: {} {}",
                line.render(),
                other.render()
            ));
        }
        if let (
            HeaderLine::Compound {
                number: line_number,
                line_type: line_type_,
                ..
            },
            HeaderLine::Compound {
                number: other_number,
                line_type: other_type,
                ..
            },
        ) = (line, other)
        {
            // `equalsExcludingDescription`: the ID is the key's, so what is left is count and type.
            if line_number != other_number || line_type_ != other_type {
                if line_type_ == other_type {
                    // "Promoting header field Number to ." -- on the line already kept.
                    if let HeaderLine::Compound { number, .. } = &mut merged[other_index] {
                        *number = Cardinality::Unbounded;
                    }
                } else {
                    use htsjdk_vcf::header::LineType::{Float, Integer};
                    let promotable = matches!(
                        (line_type_, other_type),
                        (Integer, Float) | (Float, Integer)
                    );
                    if !promotable {
                        return Err(format!(
                            "Incompatible header types, collision between these two types: {} {}",
                            line.render(),
                            other.render()
                        ));
                    }
                }
            }
        }
    }
    let mut survivors: Vec<usize> = kept.into_iter().map(|(_, index)| index).collect();
    survivors.sort_unstable();
    Ok(survivors
        .into_iter()
        .map(|index| merged[index].clone())
        .collect())
}

/// Sort by `VariantContextComparator`: the contig's index, then the start; stable.
///
/// A record on a contig the header does not have is the comparator's `NullPointerException`, and
/// is returned as its message.
pub fn sort_records<T>(
    records: &mut [T],
    contigs: &[String],
    key: impl Fn(&T) -> (&str, i64),
) -> Result<(), String> {
    let index_of = |contig: &str| contigs.iter().position(|c| c == contig);
    if records.len() > 1 {
        for record in records.iter() {
            if index_of(key(record).0).is_none() {
                return Err("java.lang.NullPointerException: Cannot invoke \
                     \"java.lang.Integer.intValue()\" because the return value of \
                     \"java.util.Map.get(Object)\" is null"
                    .to_string());
            }
        }
    }
    records.sort_by_key(|record| {
        let (contig, start) = key(record);
        (index_of(contig).unwrap_or(0), start)
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_vcf::header::LineType;

    fn unstructured(key: &str, value: &str) -> HeaderLine {
        HeaderLine::Unstructured {
            key: key.to_string(),
            value: value.to_string(),
        }
    }

    #[test]
    fn two_lines_under_one_key_are_the_one_that_sorts_first() {
        let lines = [
            unstructured("source", "handwritten"),
            HeaderLine::contig("chr1", 10, 0),
            unstructured("source", "a second"),
        ];
        let merged = smart_merge_header_lines(&lines).unwrap();
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[1].render(), "source=a second");
    }

    #[test]
    fn a_count_disagreement_promotes_the_kept_line_to_unbounded() {
        let lines = [
            HeaderLine::info("X", Cardinality::Fixed(2), LineType::Integer, "b"),
            HeaderLine::info("X", Cardinality::Fixed(1), LineType::Integer, "a"),
        ];
        let merged = smart_merge_header_lines(&lines).unwrap();
        assert_eq!(merged.len(), 1);
        assert_eq!(
            merged[0].render(),
            "INFO=<ID=X,Number=.,Type=Integer,Description=\"a\">"
        );
    }

    #[test]
    fn the_sort_is_stable_within_a_position() {
        let contigs = vec!["chr1".to_string(), "chr2".to_string()];
        let mut records = vec![("chr2", 5, 'a'), ("chr1", 9, 'b'), ("chr1", 9, 'c')];
        sort_records(&mut records, &contigs, |r| (r.0, r.1)).unwrap();
        let order: String = records.iter().map(|r| r.2).collect();
        assert_eq!(order, "bca");
    }
}
