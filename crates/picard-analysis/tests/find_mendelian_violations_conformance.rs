//! Conformance for `FindMendelianViolations` against Picard 3.4.0.
//!
//! Golden from `tools/mendelian-conformance/`: twelve runs over one trio, one site at a time.
//!
//! # What this suite is for
//!
//!  * **a violation being named as well as counted**, per parent;
//!  * **the quality floor reading the likelihood rather than the quality** where two homozygous
//!    reference parents have a child that is not;
//!  * **the depth being a separate floor on a separate field**;
//!  * **the allele balance being asked of the child and of nobody else**;
//!  * **a male child's X being haploid**, where a heterozygous call is not judged at all;
//!  * **the mitochondrion being skipped by default**;
//!  * **and `--VCF_DIR` refusing the run in the reference**, which is a bug the golden records.

use std::io::Read;

use picard_analysis::find_mendelian_violations::{
    accumulate, collect, is_mendelian_violation, render, Genotype, Metrics, Options, Outcome, Sex,
    Site, Trio, Violation,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/mendelian_violations.txt.gz");
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

/// A genotype with the defaults the fixture fills in, and the likelihoods its call implies.
fn genotype(call: &str) -> Genotype {
    let alleles: Vec<usize> = match call {
        "./." => {
            return Genotype {
                alleles: None,
                gq: 60,
                dp: 30,
                ad: Some(vec![15, 15]),
                pl: None,
            }
        }
        _ => call
            .split('/')
            .map(|allele| allele.parse().expect("an index"))
            .collect(),
    };
    let likelihoods = match call {
        "0/0" => vec![0, 50, 100],
        "1/1" => vec![100, 50, 0],
        _ => vec![50, 0, 50],
    };
    Genotype {
        alleles: Some(alleles),
        gq: 60,
        dp: 30,
        ad: Some(vec![15, 15]),
        pl: Some(likelihoods),
    }
}

fn site(contig: &str, father: Genotype, mother: Genotype, child: Genotype) -> Site {
    Site {
        contig: contig.to_string(),
        position: 100,
        filtered: false,
        alleles: vec!["A".to_string(), "C".to_string()],
        father,
        mother,
        child,
    }
}

fn trio(sex: Sex) -> Trio {
    Trio {
        family_id: "fam".to_string(),
        mother: "mother".to_string(),
        father: "father".to_string(),
        offspring: "child".to_string(),
        offspring_sex: sex,
    }
}

/// Run one case and render the table it writes.
fn run(sites: &[Site], sex: Sex, options: &Options) -> String {
    let (metrics, _) = collect(sites, &trio(sex), options);
    render(&[metrics])
}

/// A child that could have come from its parents, and one that could not.
#[test]
fn a_violation_is_named_as_well_as_counted() {
    let text = corpus();
    let options = Options::default();

    let possible = [site(
        "chr1",
        genotype("0/0"),
        genotype("0/1"),
        genotype("0/1"),
    )];
    assert_eq!(
        run(&possible, Sex::Male, &options),
        field(&text, "metrics", "a-possible-child").expect("the golden")
    );

    let denovo = [site(
        "chr1",
        genotype("0/0"),
        genotype("0/0"),
        genotype("1/1"),
    )];
    assert_eq!(
        run(&denovo, Sex::Male, &options),
        field(&text, "metrics", "a-violation").expect("the golden")
    );
    let (_, violations) = collect(&denovo, &trio(Sex::Male), &options);
    assert_eq!(violations, vec![Violation::DiploidDenovo]);
    assert_eq!(violations[0].name(), "Diploid_Denovo");

    // The other parent's, which is a different column of the same table.
    let from_the_father = [site(
        "chr1",
        genotype("0/0"),
        genotype("1/1"),
        genotype("0/0"),
    )];
    assert_eq!(
        run(&from_the_father, Sex::Male, &options),
        field(&text, "metrics", "a-violation-from-the-father").expect("the golden")
    );
    let (_, violations) = collect(&from_the_father, &trio(Sex::Male), &options);
    assert_eq!(violations, vec![Violation::HomRefHomVarHom]);
}

/// The floor that a de-novo call is judged by is the likelihood, not the quality.
#[test]
fn the_quality_floor_reads_the_likelihood() {
    let text = corpus();
    // A confident call of the variant: GQ 10, and a likelihood of ten for being reference.
    let low_quality = Genotype {
        alleles: Some(vec![1, 1]),
        gq: 10,
        dp: 30,
        ad: Some(vec![0, 30]),
        pl: Some(vec![10, 5, 0]),
    };
    let sites = [site("chr1", genotype("0/0"), genotype("0/0"), low_quality)];

    // At the default floor the site is not looked at, so it is not even a variant site.
    let mut metrics = Metrics::default();
    assert_eq!(
        accumulate(
            &sites[0],
            &trio(Sex::Male),
            &Options::default(),
            &mut metrics
        ),
        Outcome::Skipped
    );
    assert_eq!(metrics.num_variant_sites, 0);
    assert_eq!(
        run(&sites, Sex::Male, &Options::default()),
        field(&text, "metrics", "a-low-quality-child").expect("the golden")
    );

    // Lower the floor under the likelihood and the same call is a violation. The GQ of ten never
    // enters into it: this branch reads PL[0] and nothing else.
    let lower = Options {
        min_gq: 5,
        ..Options::default()
    };
    assert_eq!(
        run(&sites, Sex::Male, &lower),
        field(&text, "metrics", "a-low-quality-child-with-a-lower-floor").expect("the golden")
    );
}

/// The depth is its own floor on its own field.
#[test]
fn the_depth_is_a_separate_floor() {
    let text = corpus();
    let shallow = Genotype {
        alleles: Some(vec![1, 1]),
        gq: 60,
        dp: 3,
        ad: Some(vec![0, 3]),
        pl: Some(vec![100, 50, 0]),
    };
    let sites = [site("chr1", genotype("0/0"), genotype("0/0"), shallow)];
    let deep_enough = Options {
        min_dp: 10,
        ..Options::default()
    };
    assert_eq!(
        run(&sites, Sex::Male, &deep_enough),
        field(&text, "metrics", "a-shallow-child").expect("the golden")
    );
    // The same call at the default floor of zero is a violation, so it is the depth that dropped
    // it and not the call.
    let (metrics, violations) = collect(&sites, &trio(Sex::Male), &Options::default());
    assert_eq!(metrics.num_variant_sites, 1);
    assert_eq!(violations, vec![Violation::DiploidDenovo]);
}

/// The allele balance is asked of the child, and of nobody else.
#[test]
fn a_lopsided_parent_is_still_a_parent() {
    let text = corpus();
    let lopsided = Genotype {
        alleles: Some(vec![0, 1]),
        gq: 60,
        dp: 30,
        ad: Some(vec![29, 1]),
        pl: Some(vec![50, 0, 50]),
    };
    let sites = [site("chr1", lopsided, genotype("0/0"), genotype("0/0"))];
    // One in thirty is well under the default third, and the site is still counted: the father's
    // het is nobody's business but his.
    assert_eq!(
        run(&sites, Sex::Male, &Options::default()),
        field(&text, "metrics", "a-lopsided-het").expect("the golden")
    );

    // The same depths on the CHILD drop the site altogether.
    let child_is_lopsided = [site(
        "chr1",
        genotype("0/0"),
        genotype("0/0"),
        Genotype {
            alleles: Some(vec![0, 1]),
            gq: 60,
            dp: 30,
            ad: Some(vec![29, 1]),
            pl: Some(vec![50, 0, 50]),
        },
    )];
    let mut metrics = Metrics::default();
    assert_eq!(
        accumulate(
            &child_is_lopsided[0],
            &trio(Sex::Male),
            &Options::default(),
            &mut metrics
        ),
        Outcome::Skipped
    );
}

/// A male child's X is haploid, and a heterozygous call there is not judged.
#[test]
fn the_sex_chromosomes_are_counted_differently() {
    let text = corpus();
    let options = Options::default();
    let on_the_x = [site(
        "chrX",
        genotype("0/0"),
        genotype("0/1"),
        genotype("0/1"),
    )];

    assert_eq!(
        run(&on_the_x, Sex::Male, &options),
        field(&text, "metrics", "a-male-child-on-the-x").expect("the golden")
    );
    assert_eq!(
        run(&on_the_x, Sex::Female, &options),
        field(&text, "metrics", "a-female-child-on-the-x").expect("the golden")
    );

    // Both count the site; neither is a violation, but for different reasons. The male child's
    // call is heterozygous where it should have one allele, so the tool declines to judge it; the
    // female child's is a call her mother could have donated.
    let mut male = Metrics::default();
    assert_eq!(
        accumulate(&on_the_x[0], &trio(Sex::Male), &options, &mut male),
        Outcome::Counted
    );
    assert_eq!(male.num_variant_sites, 1);

    // A male child whose only allele is one neither parent carries is a haploid de-novo call, and
    // the donor asked is the MOTHER, the father having no X to give.
    let haploid_denovo = [site(
        "chrX",
        genotype("0/0"),
        genotype("0/0"),
        genotype("1/1"),
    )];
    let (_, violations) = collect(&haploid_denovo, &trio(Sex::Male), &options);
    assert_eq!(violations, vec![Violation::HaploidDenovo]);
    assert_eq!(violations[0].name(), "Haploid_Denovo");
}

/// The mitochondrion is left out unless it is asked for.
#[test]
fn the_mitochondrion_is_skipped() {
    let text = corpus();
    let sites = [site(
        "chrM",
        genotype("0/0"),
        genotype("0/0"),
        genotype("1/1"),
    )];
    assert_eq!(
        run(&sites, Sex::Male, &Options::default()),
        field(&text, "metrics", "the-mitochondrion").expect("the golden")
    );
    // The same site on an autosome is the violation the contig hid.
    let elsewhere = [site(
        "chr1",
        genotype("0/0"),
        genotype("0/0"),
        genotype("1/1"),
    )];
    let (metrics, _) = collect(&elsewhere, &trio(Sex::Male), &Options::default());
    assert_eq!(metrics.total_mendelian_violations, 1);
}

/// Writing the offending records out refuses the run in the reference.
#[test]
fn the_records_cannot_be_written_out() {
    let text = corpus();
    // `--VCF_DIR` subsets each record to the trio's three samples, and the header it writes for
    // that subset carries a sample twice, so the file it has just written cannot be read back.
    let recorded = field(&text, "error", "with-the-offending-records").expect("the golden");
    assert_eq!(
        recorded,
        "htsjdk.tribble.TribbleException$InvalidHeader:Your input file has a malformed header: \
         BUG: VCF header has duplicate sample names"
    );
    // The same run without the directory writes the table, so it is the writing that fails and
    // not the counting.
    assert_eq!(
        field(
            &text,
            "metrics",
            "with-the-offending-records-and-a-tab-report"
        )
        .expect("the golden"),
        field(&text, "metrics", "a-violation").expect("the golden")
    );
}

/// The inheritance test itself, in both directions.
#[test]
fn an_allele_may_come_from_either_parent() {
    let mother = genotype("0/0");
    let father = genotype("1/1");
    // One from each, whichever way round they are taken.
    assert!(!is_mendelian_violation(&mother, &father, &genotype("0/1")));
    assert!(!is_mendelian_violation(&father, &mother, &genotype("0/1")));
    // And a child homozygous for either allele has one that came from nowhere.
    assert!(is_mendelian_violation(&mother, &father, &genotype("0/0")));
    assert!(is_mendelian_violation(&mother, &father, &genotype("1/1")));
}
