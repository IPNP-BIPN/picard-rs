//! Conformance for `CreateExtendedIlluminaManifest` against Picard 3.4.0.
//!
//! Golden from `tools/arrays-conformance/`: two manifests written out, and two loci the reference
//! cannot place at all.
//!
//! # What this suite is for
//!
//!  * **the output being the input plus seven columns**;
//!  * **the alleles being placed on the build's strand**, so a negative-strand locus is
//!    complemented;
//!  * **an ambiguous SNP being decided by the probes**, because complementing tells its two
//!    alleles apart from nothing;
//!  * **an rsID coming from dbSNP by position**, and `null` standing for a position it does not
//!    carry;
//!  * **the report counting what happened**, ambiguity included;
//!  * **and a locus the reference cannot place stopping the run**, which is what the reference
//!    does rather than flagging it.

use std::io::Read;

use picard_analysis::create_extended_illumina_manifest::{
    process_snp, render_row, report, rs_id, Flag, Record, Statistics, Strand, COLUMNS,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/create_extended_manifest.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

/// The reference base at a position, the contigs being `ACGT` repeating.
fn base(position: i32) -> String {
    ["A", "C", "G", "T"][((position - 1) % 4) as usize].to_string()
}

/// The fixture's four loci, as Illumina's own manifest holds them.
fn manifest() -> Vec<Record> {
    [
        ("rs1", "[A/G]", "1", 1001, 11, 0, "+"),
        ("rs2", "[T/C]", "1", 2001, 12, 13, "-"),
        ("rs3", "[A/C]", "2", 3001, 14, 0, "+"),
        ("rs4", "[A/T]", "2", 4001, 15, 16, "+"),
    ]
    .iter()
    .map(
        |(name, snp, chr, position, address_a, address_b, strand)| Record {
            ilmn_id: format!("{name}_ilmn"),
            name: name.to_string(),
            ilmn_strand: "TOP".to_string(),
            snp: snp.to_string(),
            address_a: address_a.to_string(),
            allele_a_probe_seq: "ACGTACGT".to_string(),
            address_b: if *address_b == 0 {
                String::new()
            } else {
                address_b.to_string()
            },
            allele_b_probe_seq: if *address_b == 0 {
                String::new()
            } else {
                "ACGTACGA".to_string()
            },
            genome_build: "37".to_string(),
            chr: chr.to_string(),
            map_info: *position,
            ref_strand: Strand::parse(strand),
        },
    )
    .collect()
}

/// The input's own columns of one row, as the fixture writes them.
fn input_columns(record: &Record) -> Vec<String> {
    vec![
        record.ilmn_id.clone(),
        record.name.clone(),
        "TOP".to_string(),
        record.snp.clone(),
        record.address_a.clone(),
        record.allele_a_probe_seq.clone(),
        record.address_b.clone(),
        record.allele_b_probe_seq.clone(),
        "37".to_string(),
        record.chr.clone(),
        record.map_info.to_string(),
        "diploid".to_string(),
        "Homo sapiens".to_string(),
        "source".to_string(),
        "1".to_string(),
        "TOP".to_string(),
        "ACGTACGTACGT".to_string(),
        "ACGT".to_string(),
        "1".to_string(),
        "3".to_string(),
        match record.ref_strand {
            Strand::Positive => "+".to_string(),
            Strand::Negative => "-".to_string(),
            Strand::None => String::new(),
        },
        "0".to_string(),
    ]
}

/// The whole written manifest for one run, header and rows.
fn written(known: &[(String, i32, String)]) -> String {
    let mut lines = vec![COLUMNS.join(",")];
    for record in manifest() {
        let mut extension = process_snp(&record, &base(record.map_info));
        extension.rs_id = rs_id(known, &record.chr, record.map_info);
        lines.push(render_row(&input_columns(&record), &extension));
    }
    lines.join("\n")
}

/// The output is the input plus seven columns, and the alleles in them are the build's.
#[test]
fn the_output_is_the_input_plus_seven_columns() {
    let text = corpus();
    let known = vec![
        ("1".to_string(), 1001, "rs1001".to_string()),
        ("2".to_string(), 3001, "rs3001".to_string()),
    ];
    assert_eq!(
        written(&known),
        field(&text, "manifest", "four-loci").expect("the golden")
    );
    assert_eq!(COLUMNS.len(), 29);
    assert_eq!(
        COLUMNS[22..],
        [
            "build37Chr",
            "build37Pos",
            "build37RefAllele",
            "build37AlleleA",
            "build37AlleleB",
            "build37Rsid",
            "build37Flag"
        ]
    );
}

/// A locus read on the negative strand has its alleles complemented.
#[test]
fn the_alleles_are_placed_on_the_builds_strand() {
    let records = manifest();
    // `[T/C]` on the negative strand is `A` and `G` on the build.
    let negative = process_snp(&records[1], &base(2001));
    assert_eq!(
        (negative.allele_a.as_str(), negative.allele_b.as_str()),
        ("A", "G")
    );
    // The same pair on the positive strand would have been left alone.
    let positive = Record {
        ref_strand: Strand::Positive,
        ..records[1].clone()
    };
    let placed = process_snp(&positive, &base(2001));
    assert_eq!(
        (placed.allele_a.as_str(), placed.allele_b.as_str()),
        ("T", "C")
    );
}

/// An ambiguous SNP is decided by the probes, because the strand cannot decide it.
#[test]
fn an_ambiguous_snp_is_decided_by_the_probes() {
    let records = manifest();
    // `[A/T]` is its own complement, so the pair survives the strand unchanged and the last base
    // of each probe sequence stands in: `ACGTACGT` reads `T` and `ACGTACGA` reads `A`, and neither
    // matches the pair as written, so the pair becomes the probes'.
    assert!(records[3].is_ambiguous());
    let ambiguous = process_snp(&records[3], &base(4001));
    assert_eq!(
        (ambiguous.allele_a.as_str(), ambiguous.allele_b.as_str()),
        ("T", "A")
    );
    assert_eq!(ambiguous.flag, Flag::Pass);

    // A manifest with no B probe cannot do that, and the locus is flagged rather than guessed at.
    let without = Record {
        allele_b_probe_seq: String::new(),
        ..records[3].clone()
    };
    assert_eq!(
        process_snp(&without, &base(4001)).flag,
        Flag::MissingAlleleBProbeseq
    );

    // A SNP that is not ambiguous keeps its own pair whatever the probes read.
    assert!(!records[0].is_ambiguous());
    let plain = process_snp(&records[0], &base(1001));
    assert_eq!(
        (plain.allele_a.as_str(), plain.allele_b.as_str()),
        ("A", "G")
    );
}

/// The rsID comes from dbSNP by position, and a position it does not carry is `null`.
#[test]
fn the_rs_id_comes_from_dbsnp() {
    let text = corpus();
    assert_eq!(
        written(&[]),
        field(&text, "manifest", "no-known-sites").expect("the golden")
    );

    let known = vec![("1".to_string(), 1001, "rs1001".to_string())];
    assert_eq!(rs_id(&known, "1", 1001), "rs1001");
    // The same position on another contig is not the same position.
    assert_eq!(rs_id(&known, "2", 1001), "null");
    assert_eq!(rs_id(&known, "1", 2001), "null");
}

/// The report counts what happened, ambiguity included.
#[test]
fn the_report_counts_the_flags() {
    let text = corpus();
    let mut statistics = Statistics::default();
    for record in manifest() {
        let extension = process_snp(&record, &base(record.map_info));
        statistics.update(&record, &extension, "37");
    }
    assert_eq!(
        report(
            "extended.csv",
            "<dir>/fixture.csv",
            "<dir>/fixture.egt",
            true,
            "37",
            &statistics
        ),
        field(&text, "report", "four-loci").expect("the golden")
    );

    // Four assays, all of them SNPs, all passing, and one of them ambiguous on the positive
    // strand: `rs4`, whose alleles the probes decided.
    assert_eq!(statistics.assays, 4);
    assert_eq!(statistics.snps, 4);
    assert_eq!(statistics.assays_flagged, 0);
    assert_eq!(statistics.ambiguous_snps_on_positive_strand, 1);
    assert_eq!(statistics.ambiguous_snps_on_negative_strand, 0);
    assert_eq!(statistics.indels, 0);
}

/// A locus the reference cannot place stops the run, rather than being flagged.
#[test]
fn a_locus_off_the_reference_stops_the_run() {
    let text = corpus();
    // The tool asks the reference for the base at the locus before it decides anything, and a
    // contig the dictionary does not carry comes back as nothing at all.
    assert_eq!(
        field(&text, "error", "a-contig-that-is-not-there").as_deref(),
        Some(
            "java.lang.NullPointerException:Cannot invoke \
             \"htsjdk.samtools.SAMSequenceRecord.getSequenceLength()\" because the return value \
             of \"htsjdk.samtools.SAMSequenceDictionary.getSequence(String)\" is null"
        )
    );
    // A position past the end of a contig it does carry is refused by the reference reader.
    assert_eq!(
        field(&text, "error", "a-position-past-the-end").as_deref(),
        Some("htsjdk.samtools.SAMException:Malformed query; start point 99000 lies after end point 5000")
    );
    // Neither is a flag: the manifest never gets written, so there is nothing to flag it in.
    assert_eq!(field(&text, "manifest", "a-contig-that-is-not-there"), None);
    assert!(Flag::MissingAlleleBProbeseq.is_fail());
    assert!(!Flag::Dupe.is_fail());
}
