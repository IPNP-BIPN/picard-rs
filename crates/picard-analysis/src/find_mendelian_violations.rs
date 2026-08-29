//! `FindMendelianViolations`: the calls in a trio that cannot have been inherited.
//!
//! A child's two alleles come one from each parent, so a genotype that cannot be built that way is
//! a violation. What makes the tool more than that arithmetic is everything it declines to look
//! at: a call below `--MIN_GQ`, a call below `--MIN_DP`, a heterozygous child whose allele depths
//! do not look heterozygous, a contig in `--SKIP_CHROMS`, and a site where the trio is not variant
//! at all.
//!
//! The sex chromosomes are counted differently. A male child's X is haploid outside the
//! pseudo-autosomal regions, so the only parent that can have donated is the mother, and a
//! heterozygous call there is not a violation but a call the tool refuses to judge.
//!
//! Ported from `picard.vcf.MendelianViolations.MendelianViolationDetector`,
//! `picard.vcf.MendelianViolations.MendelianViolationMetrics` and
//! `picard.vcf.MendelianViolations.FindMendelianViolations` in Picard 3.4.0.

/// The sex of a child, as the pedigree spells it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

impl Sex {
    /// The pedigree's fifth column: `1` is male, `2` is female, anything else unknown.
    pub fn from_pedigree(code: &str) -> Sex {
        match code {
            "1" => Sex::Male,
            "2" => Sex::Female,
            _ => Sex::Unknown,
        }
    }

    /// How the metrics file spells it.
    pub fn name(&self) -> &'static str {
        match self {
            Sex::Male => "Male",
            Sex::Female => "Female",
            Sex::Unknown => "Unknown",
        }
    }
}

/// One call: the alleles by index into the site's, and the fields the filters read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genotype {
    /// `None` for a no-call; otherwise the allele indices, `0` being the reference.
    pub alleles: Option<Vec<usize>>,
    pub gq: i32,
    pub dp: i32,
    /// One depth per allele of the site, reference first.
    pub ad: Option<Vec<i32>>,
    /// One likelihood per possible genotype, the homozygous reference first.
    pub pl: Option<Vec<i32>>,
}

impl Genotype {
    pub fn is_no_call(&self) -> bool {
        self.alleles.is_none()
    }
    pub fn is_called(&self) -> bool {
        self.alleles.is_some()
    }
    pub fn alleles(&self) -> &[usize] {
        self.alleles.as_deref().unwrap_or(&[])
    }
    pub fn is_hom(&self) -> bool {
        let alleles = self.alleles();
        !alleles.is_empty() && alleles.iter().all(|allele| *allele == alleles[0])
    }
    pub fn is_het(&self) -> bool {
        self.is_called() && !self.is_hom()
    }
    pub fn is_hom_ref(&self) -> bool {
        self.is_hom() && self.alleles()[0] == 0
    }
    pub fn is_hom_var(&self) -> bool {
        self.is_hom() && self.alleles()[0] != 0
    }
    /// A heterozygous call with no reference allele, which the tool refuses to judge.
    pub fn is_het_non_ref(&self) -> bool {
        self.is_het() && !self.alleles().contains(&0)
    }
}

/// One site of one trio.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub contig: String,
    pub position: i32,
    pub filtered: bool,
    /// The site's alleles, reference first, each as written.
    pub alleles: Vec<String>,
    pub father: Genotype,
    pub mother: Genotype,
    pub child: Genotype,
}

/// A trio, as the pedigree names it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Trio {
    pub family_id: String,
    pub mother: String,
    pub father: String,
    pub offspring: String,
    pub offspring_sex: Sex,
}

/// A region of the female chromosome the male child is diploid over after all.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

/// What a run was asked for.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    pub min_gq: i32,
    pub min_dp: i32,
    pub min_het_fraction: f64,
    pub skip_chroms: Vec<String>,
    pub male_chroms: Vec<String>,
    pub female_chroms: Vec<String>,
    pub par_intervals: Vec<Interval>,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            min_gq: 30,
            min_dp: 0,
            min_het_fraction: 0.3,
            skip_chroms: vec!["MT".to_string(), "chrM".to_string()],
            male_chroms: vec!["Y".to_string(), "chrY".to_string()],
            female_chroms: vec!["X".to_string(), "chrX".to_string()],
            par_intervals: Vec::new(),
        }
    }
}

/// The kinds of violation, named the way the record written out names them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Violation {
    DiploidDenovo,
    HomVarHomVarHet,
    HomRefHomVarHom,
    HomHetHom,
    HaploidDenovo,
    HaploidOther,
    Other,
}

impl Violation {
    pub fn name(&self) -> &'static str {
        match self {
            Violation::DiploidDenovo => "Diploid_Denovo",
            Violation::HomVarHomVarHet => "HomVar_HomVar_Het",
            Violation::HomRefHomVarHom => "HomRef_HomVar_Hom",
            Violation::HomHetHom => "Hom_Het_Hom",
            Violation::HaploidDenovo => "Haploid_Denovo",
            Violation::HaploidOther => "Haploid_Other",
            Violation::Other => "Other",
        }
    }
}

/// What one site did to a trio's counters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    /// The site was not looked at, whichever of the reasons it was.
    Skipped,
    /// The site counted as a variant site and was possible.
    Counted,
    /// The site counted as a variant site and was not.
    Violated(Violation),
}

/// One row of the metrics file.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Metrics {
    pub family_id: String,
    pub mother: String,
    pub father: String,
    pub offspring: String,
    pub offspring_sex: String,
    pub num_variant_sites: u64,
    pub num_diploid_denovo: u64,
    pub num_homvar_homvar_het: u64,
    pub num_homref_homvar_hom: u64,
    pub num_hom_het_hom: u64,
    pub num_haploid_denovo: u64,
    pub num_haploid_other: u64,
    pub num_other: u64,
    pub total_mendelian_violations: u64,
}

impl Metrics {
    /// The total is derived rather than counted, and it is the sum of the seven kinds.
    pub fn calculate_derived_fields(&mut self) {
        self.total_mendelian_violations = self.num_diploid_denovo
            + self.num_homvar_homvar_het
            + self.num_homref_homvar_hom
            + self.num_hom_het_hom
            + self.num_haploid_denovo
            + self.num_haploid_other
            + self.num_other;
    }
}

/// Whether the trio is variant at all: one called genotype that is not homozygous reference.
pub fn is_variant(genotypes: &[&Genotype]) -> bool {
    genotypes
        .iter()
        .any(|genotype| genotype.is_called() && !genotype.is_hom_ref())
}

/// Whether the child's alleles can be built one from each parent.
pub fn is_mendelian_violation(mother: &Genotype, father: &Genotype, child: &Genotype) -> bool {
    let alleles = child.alleles();
    if alleles.len() < 2 {
        return false;
    }
    let (first, second) = (alleles[0], alleles[1]);
    if mother.alleles().contains(&first) && father.alleles().contains(&second) {
        return false;
    }
    if father.alleles().contains(&first) && mother.alleles().contains(&second) {
        return false;
    }
    true
}

/// Whether a position falls in one of the pseudo-autosomal regions.
pub fn is_in_pseudo_autosomal_region(intervals: &[Interval], contig: &str, position: i32) -> bool {
    intervals
        .iter()
        .any(|par| par.contig == contig && position >= par.start && position <= par.end)
}

/// One site, judged for one trio, with the trio's counters advanced.
///
/// The order of the tests is the tool's own, and it is what decides the difference between a clean
/// site and a site that was never looked at: a call the quality floor drops is not a variant site,
/// so it neither counts as a violation nor as one of the sites a violation could have come from.
pub fn accumulate(site: &Site, trio: &Trio, options: &Options, metrics: &mut Metrics) -> Outcome {
    if site.filtered {
        return Outcome::Skipped;
    }
    if site.alleles.len() < 2 {
        return Outcome::Skipped;
    }
    if options.skip_chroms.contains(&site.contig) {
        return Outcome::Skipped;
    }

    let mother = &site.mother;
    let father = &site.father;
    let child = &site.child;

    // A genotype with a non-SNP allele, or one with no reference allele of its own, takes the
    // whole trio out of the count.
    let reference_is_a_base = site.alleles[0].len() == 1;
    for genotype in [mother, father, child] {
        if genotype.is_het_non_ref() {
            return Outcome::Skipped;
        }
        let all_bases = genotype
            .alleles()
            .iter()
            .all(|allele| site.alleles[*allele].len() == 1);
        if !reference_is_a_base || !all_bases {
            return Outcome::Skipped;
        }
    }

    // More than two alleles between the three of them, the reference included, and the site is not
    // bi-allelic in this trio.
    let mut seen: Vec<usize> = vec![0];
    for genotype in [mother, father, child] {
        for allele in genotype.alleles() {
            if !seen.contains(allele) {
                seen.push(*allele);
            }
        }
    }
    if seen.len() > 2 {
        return Outcome::Skipped;
    }

    if !is_variant(&[mother, father, child]) {
        return Outcome::Skipped;
    }

    // The allele balance is asked of the CHILD and of nobody else, so a parent whose het call is
    // as lopsided as it likes changes nothing.
    if child.is_het() {
        let Some(depths) = &child.ad else {
            return Outcome::Skipped;
        };
        let called: Vec<i32> = child
            .alleles()
            .iter()
            .map(|allele| depths[*allele])
            .collect();
        let total = called[0] + called[1];
        let minimum = called[0].min(called[1]);
        if f64::from(minimum) / f64::from(total) < options.min_het_fraction {
            return Outcome::Skipped;
        }
    }

    // Whether the child is haploid here, and if it is, which parent can have donated.
    let mut haploid = false;
    let mut haploid_parent: Option<&Genotype> = None;
    // A male child's X is haploid, and only outside the pseudo-autosomal regions: inside one he is
    // diploid like his sister, and so is a child whose sex the pedigree does not give.
    if options.female_chroms.contains(&site.contig)
        && trio.offspring_sex == Sex::Male
        && !is_in_pseudo_autosomal_region(&options.par_intervals, &site.contig, site.position)
    {
        haploid = true;
        haploid_parent = Some(mother);
    }
    if options.male_chroms.contains(&site.contig) {
        if trio.offspring_sex == Sex::Male {
            haploid = true;
            haploid_parent = Some(father);
        } else {
            return Outcome::Skipped;
        }
    }

    // The quality floor. The parents are always asked for their GQ, and so is the child, EXCEPT
    // where two homozygous reference parents have a child that is not: there the number read is
    // the likelihood of being homozygous reference rather than the quality of the call made.
    if haploid {
        let parent = haploid_parent.expect("a haploid parent");
        if parent.is_no_call() || parent.gq < options.min_gq {
            return Outcome::Skipped;
        }
    } else if mother.is_no_call()
        || mother.gq < options.min_gq
        || father.is_no_call()
        || father.gq < options.min_gq
    {
        return Outcome::Skipped;
    }
    if child.is_no_call() {
        return Outcome::Skipped;
    }
    if mother.is_hom_ref() && father.is_hom_ref() && !child.is_hom_ref() {
        let Some(likelihoods) = &child.pl else {
            return Outcome::Skipped;
        };
        if likelihoods[0] < options.min_gq {
            return Outcome::Skipped;
        }
    } else if child.gq < options.min_gq {
        return Outcome::Skipped;
    }

    // And the depth, which is a different field from the quality and asked of the same genotypes.
    if haploid {
        let parent = haploid_parent.expect("a haploid parent");
        if child.dp < options.min_dp || parent.dp < options.min_dp {
            return Outcome::Skipped;
        }
    } else if child.dp < options.min_dp || mother.dp < options.min_dp || father.dp < options.min_dp
    {
        return Outcome::Skipped;
    }

    metrics.num_variant_sites += 1;

    if haploid {
        // A heterozygous call where the child has one allele is not judged at all.
        if child.is_het() {
            return Outcome::Counted;
        }
        let parent = haploid_parent.expect("a haploid parent");
        if !parent.alleles().contains(&child.alleles()[0]) {
            return if child.is_hom_ref() {
                metrics.num_haploid_other += 1;
                Outcome::Violated(Violation::HaploidOther)
            } else {
                metrics.num_haploid_denovo += 1;
                Outcome::Violated(Violation::HaploidDenovo)
            };
        }
        return Outcome::Counted;
    }

    if !is_mendelian_violation(mother, father, child) {
        return Outcome::Counted;
    }
    let violation = if mother.is_hom_ref() && father.is_hom_ref() && !child.is_hom_ref() {
        metrics.num_diploid_denovo += 1;
        Violation::DiploidDenovo
    } else if mother.is_hom_var() && father.is_hom_var() && child.is_het() {
        metrics.num_homvar_homvar_het += 1;
        Violation::HomVarHomVarHet
    } else if child.is_hom()
        && ((mother.is_hom_ref() && father.is_hom_var())
            || (mother.is_hom_var() && father.is_hom_ref()))
    {
        metrics.num_homref_homvar_hom += 1;
        Violation::HomRefHomVarHom
    } else if child.is_hom()
        && ((mother.is_hom() && father.is_het()) || (mother.is_het() && father.is_hom()))
    {
        metrics.num_hom_het_hom += 1;
        Violation::HomHetHom
    } else {
        metrics.num_other += 1;
        Violation::Other
    };
    Outcome::Violated(violation)
}

/// Every site of one file, for one trio.
pub fn collect(sites: &[Site], trio: &Trio, options: &Options) -> (Metrics, Vec<Violation>) {
    let mut metrics = Metrics {
        family_id: trio.family_id.clone(),
        mother: trio.mother.clone(),
        father: trio.father.clone(),
        offspring: trio.offspring.clone(),
        offspring_sex: trio.offspring_sex.name().to_string(),
        ..Metrics::default()
    };
    let mut violations = Vec::new();
    for site in sites {
        if let Outcome::Violated(violation) = accumulate(site, trio, options, &mut metrics) {
            violations.push(violation);
        }
    }
    metrics.calculate_derived_fields();
    (metrics, violations)
}

/// The metrics file's header.
pub const HEADER: &str = "FAMILY_ID\tMOTHER\tFATHER\tOFFSPRING\tOFFSPRING_SEX\tNUM_VARIANT_SITES\t\
NUM_DIPLOID_DENOVO\tNUM_HOMVAR_HOMVAR_HET\tNUM_HOMREF_HOMVAR_HOM\tNUM_HOM_HET_HOM\t\
NUM_HAPLOID_DENOVO\tNUM_HAPLOID_OTHER\tNUM_OTHER\tTOTAL_MENDELIAN_VIOLATIONS";

/// The table as the metrics file writes it.
pub fn render(rows: &[Metrics]) -> String {
    let mut text = String::from(HEADER);
    for row in rows {
        text.push('\n');
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.family_id,
            row.mother,
            row.father,
            row.offspring,
            row.offspring_sex,
            row.num_variant_sites,
            row.num_diploid_denovo,
            row.num_homvar_homvar_het,
            row.num_homref_homvar_hom,
            row.num_hom_het_hom,
            row.num_haploid_denovo,
            row.num_haploid_other,
            row.num_other,
            row.total_mendelian_violations
        ));
    }
    text
}
