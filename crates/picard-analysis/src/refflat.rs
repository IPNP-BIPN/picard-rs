//! Reading a refFlat gene-annotation file into an [`OverlapDetector`] of genes.
//!
//! Ported from `picard.annotation.RefFlatReader`, tag 3.4.0.
//!
//! The refFlat columns are, in order: gene name, transcript name, chromosome, strand, txStart,
//! txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds. **The coordinates are 0-based
//! half-open, and this reader is where they become 1-based inclusive** — `+1` on every start,
//! nothing on the ends. That single conversion is the reader's whole reason to exist as a
//! distinct layer; everything downstream assumes it has already happened.
//!
//! Two grouping rules that decide the gene set:
//!
//!   * lines are grouped by **gene name**, and every transcript line for a gene must agree on
//!     strand and chromosome or the whole gene is rejected;
//!   * a transcript whose chromosome is not in the sequence dictionary is dropped **before**
//!     grouping, so a gene all of whose transcripts are on unknown contigs never appears. The
//!     order of `getOverlaps` is unobservable (htsjdk-rs decision 0020), but *which* genes exist
//!     is not, and this filter is part of that.
//!
//! And two validation rules that reject a malformed transcript rather than silently mangling it:
//! an exon with `start > end`, and two exons that overlap or abut (`exons[i-1].end >= exons[i].start`).

use htsjdk_bam::overlap::OverlapDetector;

use crate::annotation::Gene;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RefFlatError {
    WrongColumnCount { line: usize, found: usize },
    StrandDisagreement(String),
    ChromosomeDisagreement(String),
    ExonCountMismatch(String),
    ExonHasNoExtent(String),
    ExonsOverlap(String),
    BadNumber(String),
}

/// One parsed refFlat row, coordinates still 0-based half-open as read.
struct Row {
    gene: String,
    transcript: String,
    chromosome: String,
    strand: String,
    tx_start: i32,
    tx_end: i32,
    cds_start: i32,
    cds_end: i32,
    exon_count: usize,
    exon_starts: Vec<i32>,
    exon_ends: Vec<i32>,
}

fn parse_row(line: &str, line_number: usize) -> Result<Row, RefFlatError> {
    let f: Vec<&str> = line.split('\t').collect();
    // The refFlat has 11 columns.
    if f.len() != 11 {
        return Err(RefFlatError::WrongColumnCount {
            line: line_number,
            found: f.len(),
        });
    }
    let num = |s: &str| {
        s.parse::<i32>()
            .map_err(|_| RefFlatError::BadNumber(s.to_string()))
    };
    // exonStarts and exonEnds are comma-terminated lists. Java's String.split(",") drops the
    // trailing empty field, so "1,2," splits to ["1","2"]; a plain Rust split would keep the
    // trailing "" and has to filter it out to match.
    let list = |s: &str| -> Result<Vec<i32>, RefFlatError> {
        s.split(',').filter(|p| !p.is_empty()).map(num).collect()
    };
    Ok(Row {
        gene: f[0].to_string(),
        transcript: f[1].to_string(),
        chromosome: f[2].to_string(),
        strand: f[3].to_string(),
        tx_start: num(f[4])?,
        tx_end: num(f[5])?,
        cds_start: num(f[6])?,
        cds_end: num(f[7])?,
        exon_count: num(f[8])? as usize,
        exon_starts: list(f[9])?,
        exon_ends: list(f[10])?,
    })
}

/// `RefFlatReader.load`.
///
/// `recognized` decides whether a chromosome is in the sequence dictionary; a transcript on an
/// unrecognized chromosome is skipped. Pass `|_| true` to keep everything.
pub fn load(
    text: &str,
    recognized: impl Fn(&str) -> bool,
) -> Result<OverlapDetector<Gene>, RefFlatError> {
    // Group rows by gene, preserving first-seen order. htsjdk groups into a HashMap, whose
    // iteration order is unobservable downstream (decision 0020); insertion order here is a
    // deterministic stand-in that changes nothing about the resulting gene set.
    let mut order: Vec<String> = Vec::new();
    let mut by_gene: std::collections::HashMap<String, Vec<Row>> = std::collections::HashMap::new();

    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let row = parse_row(line, i + 1)?;
        if !recognized(&row.chromosome) {
            continue;
        }
        let gene = row.gene.clone();
        by_gene.entry(gene.clone()).or_insert_with(|| {
            order.push(gene);
            Vec::new()
        });
        by_gene.get_mut(&row.gene).unwrap().push(row);
    }

    let mut detector = OverlapDetector::create();
    for gene_name in &order {
        let rows = &by_gene[gene_name];
        // A gene that fails validation is skipped, matching htsjdk's `catch (AnnotationException)`
        // around a single gene rather than aborting the whole load.
        if let Ok(gene) = make_gene(rows) {
            let (contig, start, end) = (gene.contig.clone(), gene.start, gene.end);
            detector.add(&contig, start, end, gene);
        }
    }
    Ok(detector)
}

/// `makeGeneFromRefFlatLines` plus `makeTranscriptFromRefFlatLine`, where the 0-based to 1-based
/// conversion lives.
fn make_gene(rows: &[Row]) -> Result<Gene, RefFlatError> {
    let first = &rows[0];
    let negative = first.strand == "-";

    // The gene extent is the min 1-based start and the max inclusive end over its transcripts.
    let mut start = i32::MAX;
    let mut end = i32::MIN;
    for row in rows {
        start = start.min(row.tx_start + 1);
        end = end.max(row.tx_end);
    }

    let mut gene = Gene::new(&first.chromosome, start, end, negative, &first.gene);

    for row in rows {
        if row.strand != first.strand {
            return Err(RefFlatError::StrandDisagreement(first.gene.clone()));
        }
        if row.chromosome != first.chromosome {
            return Err(RefFlatError::ChromosomeDisagreement(first.gene.clone()));
        }
        let description = format!("{}:{}", row.gene, row.transcript);
        if row.exon_count != row.exon_starts.len() || row.exon_count != row.exon_ends.len() {
            return Err(RefFlatError::ExonCountMismatch(description));
        }

        let mut tx = crate::annotation::Transcript::new(
            &row.transcript,
            row.tx_start + 1,
            row.tx_end,
            row.cds_start + 1,
            row.cds_end,
        );
        for i in 0..row.exon_count {
            let e_start = row.exon_starts[i] + 1;
            let e_end = row.exon_ends[i];
            if e_start > e_end {
                return Err(RefFlatError::ExonHasNoExtent(description));
            }
            if i > 0 && tx.exons[i - 1].end >= e_start {
                return Err(RefFlatError::ExonsOverlap(description));
            }
            tx.add_exon(e_start, e_end);
        }
        gene.transcripts.push(tx);
    }

    Ok(gene)
}

#[cfg(test)]
mod tests {
    use super::*;

    const LINE: &str = "GENEA\ttxA\tchr1\t+\t999\t2000\t1199\t1800\t2\t999,1499,\t1200,2000,";

    #[test]
    fn coordinates_become_one_based_inclusive() {
        let d = load(LINE, |_| true).unwrap();
        let genes = d.get_overlaps("chr1", 1000, 1000);
        assert_eq!(genes.len(), 1);
        let g = genes[0];
        // txStart 999 (0-based) -> 1000 (1-based); txEnd 2000 stays 2000.
        assert_eq!((g.start, g.end), (1000, 2000));
        let tx = &g.transcripts[0];
        assert_eq!((tx.transcription_start, tx.transcription_end), (1000, 2000));
        assert_eq!((tx.coding_start, tx.coding_end), (1200, 1800));
        // exonStarts 999,1499 -> 1000,1500; ends unchanged.
        assert_eq!(
            tx.exons[0],
            crate::annotation::Exon {
                start: 1000,
                end: 1200
            }
        );
        assert_eq!(
            tx.exons[1],
            crate::annotation::Exon {
                start: 1500,
                end: 2000
            }
        );
    }

    #[test]
    fn an_unrecognized_chromosome_is_skipped() {
        let d = load(LINE, |c| c == "chr2").unwrap();
        assert!(d.get_all().is_empty(), "chr1 is not recognized here");
    }

    /// Two transcript lines for one gene are grouped, and the extent spans both.
    #[test]
    fn transcripts_of_one_gene_are_grouped() {
        let text = "GENEA\ttx1\tchr1\t+\t100\t500\t100\t500\t1\t100,\t500,\n\
                    GENEA\ttx2\tchr1\t+\t400\t900\t400\t900\t1\t400,\t900,";
        let d = load(text, |_| true).unwrap();
        let g = d.get_overlaps("chr1", 200, 200);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].transcripts.len(), 2);
        assert_eq!((g[0].start, g[0].end), (101, 900), "min start, max end");
    }

    #[test]
    fn a_gene_whose_transcripts_disagree_on_strand_is_skipped() {
        let text = "G\ttx1\tchr1\t+\t100\t500\t100\t500\t1\t100,\t500,\n\
                    G\ttx2\tchr1\t-\t100\t500\t100\t500\t1\t100,\t500,";
        // The invalid gene is dropped, not fatal, so the load succeeds with no genes.
        let d = load(text, |_| true).unwrap();
        assert!(d.get_all().is_empty());
    }

    #[test]
    fn overlapping_exons_reject_the_gene() {
        let text = "G\ttx\tchr1\t+\t100\t500\t100\t500\t2\t100,200,\t250,500,";
        let d = load(text, |_| true).unwrap();
        assert!(
            d.get_all().is_empty(),
            "exon 1 ends at 250, exon 2 starts at 201: overlap"
        );
    }

    #[test]
    fn a_wrong_column_count_is_a_hard_error() {
        assert!(matches!(
            load("too\tfew\tcolumns", |_| true),
            Err(RefFlatError::WrongColumnCount { found: 3, .. })
        ));
    }

    #[test]
    fn the_trailing_comma_in_exon_lists_is_dropped() {
        // "100," must parse to one exon, not one exon and an empty field.
        let d = load(LINE, |_| true).unwrap();
        assert_eq!(
            d.get_overlaps("chr1", 1000, 1000)[0].transcripts[0]
                .exons
                .len(),
            2
        );
    }
}
