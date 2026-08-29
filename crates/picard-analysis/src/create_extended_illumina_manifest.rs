//! `CreateExtendedIlluminaManifest`: Illumina's manifest, plus where each locus sits on the build.
//!
//! The output is the input with seven columns added, and the columns are the answer to one
//! question per locus: which contig and position it has on the target build, which base the build
//! carries there, which two alleles the assay reads there, what dbSNP calls it, and whether any of
//! that worked. A locus that fails is FLAGGED and written out anyway: the flag is what the
//! downstream tool acts on, and dropping the row would leave nothing to act on.
//!
//! The alleles are not a copy of the manifest's. The SNP column is written on the assay's own
//! strand, so on a locus whose reference strand is negative the pair is complemented. And where
//! the two alleles are each other's complement, the pair cannot be told apart that way at all: the
//! last base of each probe sequence decides instead.
//!
//! Ported from `picard.arrays.illumina.CreateExtendedIlluminaManifest`,
//! `picard.arrays.illumina.Build37ExtendedIlluminaManifestRecordCreator`,
//! `picard.arrays.illumina.Build37ExtendedIlluminaManifestRecord` and
//! `picard.arrays.illumina.IlluminaManifestRecord` in Picard 3.4.0.

/// The tool's own version, which the report names.
pub const VERSION: &str = "2.0";

/// The two SNP spellings that stand for an insertion and a deletion.
pub const INDELS: [&str; 2] = ["[D/I]", "[I/D]"];

/// The SNPs whose two alleles are each other's complement, and so cannot be placed on a strand.
pub const AMBIGUOUS_SNPS: [&str; 4] = ["[A/T]", "[T/A]", "[C/G]", "[G/C]"];

/// Which strand of the build the assay reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Positive,
    Negative,
    None,
}

impl Strand {
    pub fn parse(text: &str) -> Strand {
        match text {
            "+" => Strand::Positive,
            "-" => Strand::Negative,
            _ => Strand::None,
        }
    }
}

/// Why a locus was flagged, or that it was not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flag {
    Pass,
    Dupe,
    IlluminaFlagged,
    LiftoverFailed,
    UnsupportedGenomeBuild,
    ProbeSequenceMismatch,
    MissingAlleleBProbeseq,
}

impl Flag {
    pub fn name(&self) -> &'static str {
        match self {
            Flag::Pass => "PASS",
            Flag::Dupe => "DUPE",
            Flag::IlluminaFlagged => "ILLUMINA_FLAGGED",
            Flag::LiftoverFailed => "LIFTOVER_FAILED",
            Flag::UnsupportedGenomeBuild => "UNSUPPORTED_GENOME_BUILD",
            Flag::ProbeSequenceMismatch => "PROBE_SEQUENCE_MISMATCH",
            Flag::MissingAlleleBProbeseq => "MISSING_ALLELE_B_PROBESEQ",
        }
    }

    pub fn is_fail(&self) -> bool {
        !matches!(self, Flag::Pass | Flag::Dupe)
    }
}

/// One row of Illumina's own manifest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub ilmn_id: String,
    pub name: String,
    pub ilmn_strand: String,
    /// `[A/G]` and the like, upper case.
    pub snp: String,
    pub address_a: String,
    pub allele_a_probe_seq: String,
    pub address_b: String,
    pub allele_b_probe_seq: String,
    pub genome_build: String,
    pub chr: String,
    pub map_info: i32,
    pub ref_strand: Strand,
}

impl Record {
    pub fn is_indel(&self) -> bool {
        INDELS.contains(&self.snp.as_str())
    }
    pub fn is_snp(&self) -> bool {
        !self.is_indel()
    }
    /// Whether the two alleles are each other's complement.
    pub fn is_ambiguous(&self) -> bool {
        AMBIGUOUS_SNPS.contains(&self.snp.as_str())
    }
    /// The first allele of the SNP column.
    pub fn allele_a(&self) -> String {
        self.snp[1..2].to_string()
    }
    /// The second.
    pub fn allele_b(&self) -> String {
        self.snp[3..4].to_string()
    }
}

/// The seven columns the extension adds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Extension {
    pub b37_chr: String,
    pub b37_pos: i32,
    pub ref_allele: String,
    pub allele_a: String,
    pub allele_b: String,
    /// `null` where dbSNP has nothing at the position.
    pub rs_id: String,
    pub flag: Flag,
}

/// The complement of one base.
fn complement(base: &str) -> String {
    base.chars()
        .rev()
        .map(|character| match character {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

/// The two alleles the assay reads on the build, and the flag that comes with them.
///
/// The pair is taken from the SNP column and complemented where the assay reads the negative
/// strand. An ambiguous SNP survives that unchanged, so the probes decide: the last base of each
/// probe sequence is the allele it reads, and where NEITHER matches the pair as written, on the
/// positive strand, the pair is replaced by the probes' own. A manifest with no B probe cannot do
/// that at all, and the locus is flagged for it.
pub fn process_snp(record: &Record, reference_base: &str) -> Extension {
    let mut allele_a = record.allele_a();
    let mut allele_b = record.allele_b();
    if record.ref_strand == Strand::Negative {
        allele_a = complement(&allele_a);
        allele_b = complement(&allele_b);
    }
    let mut flag = Flag::Pass;
    if record.is_ambiguous() {
        if record.allele_b_probe_seq.is_empty() {
            flag = Flag::MissingAlleleBProbeseq;
        } else {
            let probe_a = last_base(&record.allele_a_probe_seq);
            let probe_b = last_base(&record.allele_b_probe_seq);
            if probe_a != allele_a && probe_b != allele_b && record.ref_strand == Strand::Positive {
                allele_a = probe_a;
                allele_b = probe_b;
            }
        }
    }
    Extension {
        b37_chr: record.chr.clone(),
        b37_pos: record.map_info,
        ref_allele: reference_base.to_string(),
        allele_a,
        allele_b,
        rs_id: String::new(),
        flag,
    }
}

/// The base a probe reads, which is the last of its sequence.
fn last_base(probe: &str) -> String {
    probe[probe.len() - 1..].to_string()
}

/// The name dbSNP gives a position, or `null` where it gives none.
///
/// The lookup is by position alone: the alleles are not compared, so a site dbSNP carries under a
/// different variant still lends its name.
pub fn rs_id(known_sites: &[(String, i32, String)], chr: &str, position: i32) -> String {
    known_sites
        .iter()
        .find(|(contig, at, _)| contig == chr && *at == position)
        .map(|(_, _, name)| name.clone())
        .unwrap_or_else(|| "null".to_string())
}

/// The counters the report prints.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Statistics {
    pub assays: i32,
    pub assays_flagged: i32,
    pub assays_duplicated: i32,
    pub snps: i32,
    pub snps_flagged: i32,
    pub snps_duplicated: i32,
    pub snps_illumina_flagged: i32,
    pub snp_probe_sequence_mismatch: i32,
    pub snp_missing_allele_b_probe_sequence: i32,
    pub ambiguous_snps_on_positive_strand: i32,
    pub ambiguous_snps_on_negative_strand: i32,
    pub indels: i32,
    pub indels_flagged: i32,
    pub indels_duplicated: i32,
    pub indels_illumina_flagged: i32,
    pub indel_probe_sequence_mismatch: i32,
    pub indel_source_sequence_invalid: i32,
    pub indels_not_found: i32,
    pub indel_conflict: i32,
    pub on_target_build: i32,
    pub on_unsupported_genome_build: i32,
    pub liftover_failed: i32,
    pub ref_strand_mismatch: i32,
}

impl Statistics {
    /// One record's contribution, by the same tests the reference applies in the same order.
    pub fn update(&mut self, record: &Record, extension: &Extension, target_build: &str) {
        self.assays += 1;
        if record.is_snp() {
            self.snps += 1;
        } else {
            self.indels += 1;
        }
        if record.genome_build == target_build {
            self.on_target_build += 1;
        }
        match extension.flag {
            Flag::UnsupportedGenomeBuild => self.on_unsupported_genome_build += 1,
            Flag::LiftoverFailed => self.liftover_failed += 1,
            _ => {}
        }
        if !extension.flag.is_fail() {
            if extension.flag == Flag::Dupe {
                self.assays_duplicated += 1;
                if record.is_snp() {
                    self.snps_duplicated += 1;
                } else {
                    self.indels_duplicated += 1;
                }
            }
            // An ambiguous SNP is counted by the strand the MANIFEST declares, whether or not the
            // probes ended up deciding the alleles.
            if record.is_ambiguous() {
                match record.ref_strand {
                    Strand::Negative => self.ambiguous_snps_on_negative_strand += 1,
                    Strand::Positive => self.ambiguous_snps_on_positive_strand += 1,
                    Strand::None => {}
                }
            }
            return;
        }
        self.assays_flagged += 1;
        if record.is_snp() {
            self.snps_flagged += 1;
            match extension.flag {
                Flag::IlluminaFlagged => self.snps_illumina_flagged += 1,
                Flag::ProbeSequenceMismatch => self.snp_probe_sequence_mismatch += 1,
                Flag::MissingAlleleBProbeseq => self.snp_missing_allele_b_probe_sequence += 1,
                _ => {}
            }
        } else {
            self.indels_flagged += 1;
            if extension.flag == Flag::IlluminaFlagged {
                self.indels_illumina_flagged += 1;
            }
        }
    }
}

/// The report, with the blank lines the file separates its sections with.
pub fn report(
    output_name: &str,
    input_path: &str,
    cluster_path: &str,
    flag_duplicates: bool,
    target_build: &str,
    statistics: &Statistics,
) -> String {
    let mut lines = vec![format!(
        "CreateExtendedIlluminaManifest (version: {VERSION}) Report For: {output_name}"
    )];
    lines.push(format!("Using Illumina Manifest: {input_path}"));
    if flag_duplicates {
        lines.push("Duplicates were flagged".to_string());
    }
    lines.push(format!("Using Illumina EGT: {cluster_path}"));
    lines.push(format!("Total Number of Assays: {}", statistics.assays));
    lines.push(format!(
        "Number of Assays on Build {target_build}: {}",
        statistics.on_target_build
    ));
    lines.push(format!(
        "Number of Assays on unsupported genome build: {}",
        statistics.on_unsupported_genome_build
    ));
    lines.push(format!(
        "Number of Assays failing liftover: {}",
        statistics.liftover_failed
    ));
    lines.push(format!(
        "Number of Assays on Build {target_build} or successfully lifted over: {}",
        statistics.on_target_build
    ));
    lines.push(format!(
        "Number of Passing Assays: {}",
        statistics.assays - statistics.assays_flagged
    ));
    lines.push(format!(
        "Number of Duplicated Assays: {}",
        statistics.assays_duplicated
    ));
    lines.push(format!(
        "Number of Failing Assays: {}",
        statistics.assays_flagged
    ));
    lines.push(format!("Number of SNPs: {}", statistics.snps));
    lines.push(format!(
        "Number of Passing SNPs: {}",
        statistics.snps - statistics.snps_flagged
    ));
    lines.push(format!(
        "Number of Duplicated SNPs: {}",
        statistics.snps_duplicated
    ));
    lines.push(format!(
        "Number of Failing SNPs: {}",
        statistics.snps_flagged
    ));
    lines.push(format!(
        "Number of SNPs failed by Illumina: {}",
        statistics.snps_illumina_flagged
    ));
    lines.push(format!(
        "Number of SNPs failed for refStrand mismatch: {}",
        statistics.ref_strand_mismatch
    ));
    lines.push(format!(
        "Number of SNPs failed for missing AlleleB ProbeSeq: {}",
        statistics.snp_missing_allele_b_probe_sequence
    ));
    lines.push(format!(
        "Number of SNPs failed for alleleA probe sequence mismatch: {}",
        statistics.snp_probe_sequence_mismatch
    ));
    lines.push(format!(
        "Number of ambiguous SNPs on Positive Strand: {}",
        statistics.ambiguous_snps_on_positive_strand
    ));
    lines.push(format!(
        "Number of ambiguous SNPs on Negative Strand: {}",
        statistics.ambiguous_snps_on_negative_strand
    ));
    lines.push(format!("Number of Indels: {}", statistics.indels));
    lines.push(format!(
        "Number of Passing Indels: {}",
        statistics.indels - statistics.indels_flagged
    ));
    lines.push(format!(
        "Number of Duplicated Indels: {}",
        statistics.indels_duplicated
    ));
    lines.push(format!(
        "Number of Failing Indels: {}",
        statistics.indels_flagged
    ));
    lines.push(format!(
        "Number of Indels failed by Illumina: {}",
        statistics.indels_illumina_flagged
    ));
    lines.push(format!(
        "Number of Indels failed for probe sequence mismatch: {}",
        statistics.indel_probe_sequence_mismatch
    ));
    lines.push(format!(
        "Number of Indels failed for source sequence invalid: {}",
        statistics.indel_source_sequence_invalid
    ));
    lines.push(format!(
        "Number of Indels not found: {}",
        statistics.indels_not_found
    ));
    lines.push(format!(
        "Number of Indels failed for conflict: {}",
        statistics.indel_conflict
    ));
    lines.join("\n")
}

/// The columns of the extended manifest, the input's followed by the seven that are added.
pub const COLUMNS: [&str; 29] = [
    "IlmnID",
    "Name",
    "IlmnStrand",
    "SNP",
    "AddressA_ID",
    "AlleleA_ProbeSeq",
    "AddressB_ID",
    "AlleleB_ProbeSeq",
    "GenomeBuild",
    "Chr",
    "MapInfo",
    "Ploidy",
    "Species",
    "Source",
    "SourceVersion",
    "SourceStrand",
    "SourceSeq",
    "TopGenomicSeq",
    "BeadSetID",
    "Exp_Clusters",
    "RefStrand",
    "Intensity_Only",
    "build37Chr",
    "build37Pos",
    "build37RefAllele",
    "build37AlleleA",
    "build37AlleleB",
    "build37Rsid",
    "build37Flag",
];

/// One row of the written manifest: the input's own columns, then the seven computed ones.
pub fn render_row(input_columns: &[String], extension: &Extension) -> String {
    let mut columns = input_columns.to_vec();
    columns.push(extension.b37_chr.clone());
    columns.push(extension.b37_pos.to_string());
    columns.push(extension.ref_allele.clone());
    columns.push(extension.allele_a.clone());
    columns.push(extension.allele_b.clone());
    columns.push(extension.rs_id.clone());
    columns.push(extension.flag.name().to_string());
    columns.join(",")
}
