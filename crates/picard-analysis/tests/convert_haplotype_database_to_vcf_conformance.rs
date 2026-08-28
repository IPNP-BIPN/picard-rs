//! Conformance for `ConvertHaplotypeDatabaseToVcf` against Picard 3.4.0.
//!
//! Each case carries the haplotype database the tool read and the VCF it wrote. The reference is a
//! repeating `ACGT` on chr1 and a repeating `TTTTGGGGCCCCAAAA` on chr2, both sixty bases, so the
//! base at any position is known by arithmetic and the port needs no FASTA.
//!
//! # What this suite is for
//!
//!  * **REF being the allele that disagrees with the reference, in both directions**;
//!  * **AF being the frequency of the ALT the tool wrote**;
//!  * **a block of one being unphased and carrying no phase set**;
//!  * **a block of more carrying one, and that phase set being the first SNP's position rather
//!    than the anchor's**;
//!  * **a matching major allele reversing the genotype inside a phased block only**;
//!  * **the records coming out in dictionary order**;
//!  * **a row neither of whose alleles matches being refused**;
//!  * **and the three shapes of malformed table each being refused differently.**

use std::io::Read;

use picard_analysis::haplotype_map::{
    as_vcf, block_as_vcf, format_frequency, parse_haplotype_database, HaplotypeBlock, Snp,
    ALLELE_DISAGREEMENT_PREFIX, CHROMOSOME_MISMATCH_PREFIX, HET_GENOTYPE_FOR_PHASING,
    INVALID_RECORD_PREFIX, MISSING_HEADER_PREFIX, NO_HAPLOTYPE_PREFIX, VCF_SOURCE,
};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join("convert_haplotype_database_to_vcf.txt.gz");
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

fn field(text: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| unescape(&line[prefix.len()..]))
}

/// The dictionary's contig order, which is the FASTA's.
fn sequence_order() -> Vec<String> {
    vec!["chr1".to_string(), "chr2".to_string()]
}

/// The reference the dump wrote, which is periodic and therefore known by arithmetic.
fn reference_base(contig: &str, position: i32) -> Option<u8> {
    let index = (position - 1) as usize;
    match contig {
        "chr1" => Some(b"ACGT"[index % 4]),
        "chr2" => Some(b"TTTTGGGGCCCCAAAA"[index % 16]),
        _ => None,
    }
}

/// One record of a written VCF, as its fields.
#[derive(Debug, PartialEq)]
struct Written {
    chromosome: String,
    position: i32,
    id: String,
    reference: String,
    alternate: String,
    info: String,
    format: String,
    sample: String,
}

fn written(text: &str, case: &str) -> Vec<Written> {
    field(text, "out", case)
        .unwrap_or_else(|| panic!("{case} wrote a VCF"))
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let columns: Vec<&str> = line.split('\t').collect();
            Written {
                chromosome: columns[0].to_string(),
                position: columns[1].parse().expect("a position"),
                id: columns[2].to_string(),
                reference: columns[3].to_string(),
                alternate: columns[4].to_string(),
                info: columns[7].to_string(),
                format: columns[8].to_string(),
                sample: columns[9].to_string(),
            }
        })
        .collect()
}

fn blocks(text: &str, name: &str) -> Vec<HaplotypeBlock> {
    parse_haplotype_database(&field(text, "db", name).unwrap_or_else(|| panic!("{name}")))
        .unwrap_or_else(|error| panic!("{name}: {error}"))
}

/// The port's records rendered the way the tool writes them.
fn rendered(text: &str, database: &str) -> Vec<Written> {
    as_vcf(&blocks(text, database), &sequence_order(), reference_base)
        .expect("the records")
        .into_iter()
        .map(|record| {
            let phased = record.phase_set.is_some();
            Written {
                chromosome: record.chromosome,
                position: record.position,
                id: record.id,
                reference: record.reference.to_string(),
                alternate: record.alternate.to_string(),
                info: format!("AF={:.3}", record.allele_frequency),
                format: if phased {
                    "GT:PS".to_string()
                } else {
                    "GT".to_string()
                },
                sample: match record.phase_set {
                    Some(set) => format!("{}:{set}", record.genotype),
                    None => record.genotype,
                },
            }
        })
        .collect()
}

/// Every case that writes a VCF writes the records the port produces.
#[test]
fn every_case_writes_the_same_records() {
    let text = corpus();
    for (case, database) in [
        ("major-is-reference", "plain"),
        ("minor-is-reference", "swapped"),
        ("phased-block", "block"),
        ("sorted-output", "unsorted"),
    ] {
        assert_eq!(rendered(&text, database), written(&text, case), "{case}");
    }
    // The table with a header and no rows writes a VCF holding nothing.
    assert!(written(&text, "no-rows").is_empty());
}

/// REF is the allele that DISAGREES with the reference, in both directions: the two one-row cases
/// come out with the same alleles and only the frequency tells them apart.
#[test]
fn ref_is_the_allele_that_disagrees() {
    let text = corpus();
    let major = &written(&text, "major-is-reference")[0];
    let minor = &written(&text, "minor-is-reference")[0];
    assert_eq!(reference_base("chr1", 1), Some(b'A'));
    assert_eq!(major.reference, "C");
    assert_eq!(major.alternate, "A");
    assert_eq!(minor.reference, major.reference);
    assert_eq!(minor.alternate, major.alternate);
    // ALT is the allele that matches, and AF is ALT's frequency.
    assert_eq!(major.info, "AF=0.750");
    assert_eq!(minor.info, "AF=0.250");
    // Neither is phased, both blocks holding one SNP.
    assert_eq!(major.format, "GT");
    assert_eq!(major.sample, "0/1");
    assert_eq!(minor.sample, "0/1");
}

/// A block of more than one SNP is phased on its FIRST SNP by position, not on the row named as
/// the anchor, and the row whose major allele matched has its genotype reversed.
#[test]
fn the_phase_set_is_the_first_snps_position() {
    let text = corpus();
    let records = written(&text, "phased-block");
    assert_eq!(records.len(), 3);
    // The fixture named rs2, at position 2, as the anchor of all three.
    let database = field(&text, "db", "block").expect("the database");
    assert!(database.contains("\trs2\n") || database.contains("\trs2\t"));
    for record in &records {
        assert_eq!(record.format, "GT:PS");
        assert!(record.sample.ends_with(":1"), "{record:?}");
    }
    assert_eq!(records[0].sample, "1|0:1");
    assert_eq!(records[1].sample, "0|1:1");
    assert_eq!(records[2].sample, "1|0:1");
}

/// The same block of one SNP whose major allele matches keeps `0/1`, so the reversal belongs to
/// the phasing and not to the swap.
#[test]
fn the_reversal_belongs_to_the_phasing() {
    let snp = Snp {
        name: "rs1".to_string(),
        chromosome: "chr1".to_string(),
        position: 1,
        major_allele: b'A',
        minor_allele: b'C',
        minor_allele_frequency: 0.25,
        panels: None,
    };
    let mut alone = HaplotypeBlock::default();
    alone.add_snp(snp.clone()).expect("added");
    let records = block_as_vcf(&alone, reference_base).expect("records");
    assert_eq!(records[0].genotype, "0/1");
    assert_eq!(records[0].phase_set, None);
    // The same SNP beside another is reversed and phased.
    let mut pair = alone.clone();
    pair.add_snp(Snp {
        name: "rs2".to_string(),
        position: 2,
        major_allele: b'G',
        minor_allele: b'C',
        ..snp
    })
    .expect("added");
    let records = block_as_vcf(&pair, reference_base).expect("records");
    assert_eq!(records[0].genotype, "1|0");
    assert_eq!(records[1].genotype, "0|1");
    assert_eq!(records[0].phase_set, Some(1));
    assert_eq!(records[1].phase_set, Some(1));
}

/// The records come out in dictionary order whatever order the table listed them in.
#[test]
fn the_records_come_out_in_dictionary_order() {
    let text = corpus();
    let records = written(&text, "sorted-output");
    let positions: Vec<(String, i32)> = records
        .iter()
        .map(|record| (record.chromosome.clone(), record.position))
        .collect();
    assert_eq!(
        positions,
        vec![
            ("chr1".to_string(), 1),
            ("chr1".to_string(), 3),
            ("chr2".to_string(), 5)
        ]
    );
    // The table listed chr2 first.
    let database = field(&text, "db", "unsorted").expect("the database");
    let first = database
        .lines()
        .find(|line| !line.starts_with('@') && !line.starts_with('#'))
        .expect("a row");
    assert!(first.starts_with("chr2\t"), "{first}");
}

/// A row neither of whose alleles matches the reference is refused, by a message naming the SNP
/// as contig and position and neither allele.
#[test]
fn a_row_that_matches_neither_allele_is_refused() {
    let text = corpus();
    let error = field(&text, "error", "neither-allele-matches").expect("the refusal");
    assert_eq!(
        error,
        format!("java.lang.RuntimeException:{ALLELE_DISAGREEMENT_PREFIX}chr1:2")
    );
    // chr1 position 2 is C, and the row named A and T.
    assert_eq!(reference_base("chr1", 2), Some(b'C'));
    let blocks = blocks(&text, "disagreeing");
    assert_eq!(
        as_vcf(&blocks, &sequence_order(), reference_base),
        Err(format!("{ALLELE_DISAGREEMENT_PREFIX}chr1:2"))
    );
}

/// The three shapes of malformed table are each refused differently, and by the READER rather
/// than by the conversion.
#[test]
fn a_malformed_table_is_refused_by_the_reader() {
    let text = corpus();
    let no_header = field(&text, "error", "no-header").expect("the refusal");
    assert!(no_header.contains(MISSING_HEADER_PREFIX), "{no_header}");
    assert!(parse_haplotype_database("#CHROMOSOME\tPOSITION\n")
        .unwrap_err()
        .starts_with(MISSING_HEADER_PREFIX));

    let short = field(&text, "error", "short-row").expect("the refusal");
    assert_eq!(
        short,
        format!("picard.PicardException:{INVALID_RECORD_PREFIX}4 fields: chr1\t1\trs1\tA")
    );
    assert_eq!(
        parse_haplotype_database("@HD\tVN:1.6\nchr1\t1\trs1\tA\n").unwrap_err(),
        format!("{INVALID_RECORD_PREFIX}4 fields: chr1\t1\trs1\tA")
    );

    let dangling = field(&text, "error", "dangling-anchor").expect("the refusal");
    assert_eq!(
        dangling,
        format!("picard.PicardException:{NO_HAPLOTYPE_PREFIX}rsX")
    );
    assert_eq!(
        parse_haplotype_database("@HD\tVN:1.6\nchr1\t1\trs1\tA\tC\t0.10\trsX\n").unwrap_err(),
        format!("{NO_HAPLOTYPE_PREFIX}rsX")
    );

    // And a block spanning two contigs is refused by the block itself.
    let across = parse_haplotype_database(
        "@HD\tVN:1.6\nchr1\t1\trs1\tA\tC\t0.10\trs1\nchr2\t5\trs2\tG\tC\t0.20\trs1\n",
    )
    .unwrap_err();
    assert!(across.starts_with(CHROMOSOME_MISMATCH_PREFIX), "{across}");
}

/// The header the tool writes names its own source and a sample that exists only to carry a
/// phase, while the reference line is the FASTA's URI rather than the literal the code sets.
#[test]
fn the_header_names_the_source_and_the_phasing_sample() {
    let text = corpus();
    let vcf = field(&text, "out", "major-is-reference").expect("the VCF");
    assert!(vcf.contains(&format!("##source={VCF_SOURCE}")), "{vcf}");
    assert!(vcf.contains(HET_GENOTYPE_FOR_PHASING), "{vcf}");
    assert!(vcf.contains("##reference=file://"), "{vcf}");
    assert!(!vcf.contains(&format!("##reference={VCF_SOURCE}")), "{vcf}");
}

/// The frequency formatter drops trailing zeros, which is what the table's own rows show.
#[test]
fn the_frequency_drops_its_trailing_zeros() {
    assert_eq!(format_frequency(0.10), "0.1");
    assert_eq!(format_frequency(0.25), "0.25");
    assert_eq!(format_frequency(1.0), "1");
    assert_eq!(format_frequency(0.0), "0");
}
