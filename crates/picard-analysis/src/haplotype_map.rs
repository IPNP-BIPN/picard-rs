//! `HaplotypeMap`: a table of major and minor alleles grouped into blocks by an anchor SNP, and
//! the two things Picard turns it into.
//!
//! Ported from `picard.fingerprint.HaplotypeMap`, `picard.fingerprint.HaplotypeBlock` and
//! `picard.fingerprint.Snp` in Picard 3.4.0.

use std::collections::BTreeMap;

/// `HaplotypeMap.fromHaplotypeDatabase`, on a file whose first line is not a `@` header.
pub const MISSING_HEADER_PREFIX: &str = "Haplotype map file must contain header: ";
/// The same, on a row with too few or too many fields.
pub const INVALID_RECORD_PREFIX: &str = "Invalid haplotype map record contains ";
/// The same, on a row whose anchor names no row of its own.
pub const NO_HAPLOTYPE_PREFIX: &str = "No haplotype found for anchor snp ";
/// `HaplotypeBlock.addSnp`, on a SNP whose contig is not the block's.
pub const CHROMOSOME_MISMATCH_PREFIX: &str = "Snp chromosome ";
/// `HaplotypeMap.asVcf`, on a SNP neither of whose alleles is the reference base.
pub const ALLELE_DISAGREEMENT_PREFIX: &str =
    "One of the two alleles should agree with the reference: ";
/// `HaplotypeMap.writeAsVcf`'s sample, which exists only so a genotype can carry a phase.
pub const HET_GENOTYPE_FOR_PHASING: &str = "HetGenotypeForPhasing";
/// And the `##source` value it writes, which the `##reference` line does NOT keep.
pub const VCF_SOURCE: &str = "HaplotypeMap::writeAsVcf";

/// One row of the table.
#[derive(Debug, Clone, PartialEq)]
pub struct Snp {
    pub name: String,
    pub chromosome: String,
    pub position: i32,
    pub major_allele: u8,
    pub minor_allele: u8,
    pub minor_allele_frequency: f64,
    pub panels: Option<Vec<String>>,
}

impl Snp {
    /// `Snp.getAlleleString`, which lower-cases the SECOND allele and leaves the first alone.
    pub fn allele_string(&self) -> String {
        format!(
            "{}{}",
            self.major_allele as char,
            (self.minor_allele as char).to_ascii_lowercase()
        )
    }
}

/// One block, whose SNPs share a contig and are held in position order.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct HaplotypeBlock {
    pub snps: Vec<Snp>,
}

impl HaplotypeBlock {
    /// `HaplotypeBlock.addSnp`, which refuses a SNP from another contig.
    pub fn add_snp(&mut self, snp: Snp) -> Result<(), String> {
        if let Some(first) = self.snps.first() {
            if first.chromosome != snp.chromosome {
                return Err(format!(
                    "{CHROMOSOME_MISMATCH_PREFIX}{} does not agree with chromosome of existing snp(s): {}",
                    snp.chromosome, first.chromosome
                ));
            }
        }
        self.snps.push(snp);
        Ok(())
    }

    /// `new TreeSet<>(block.getSnps())`, which orders by contig and then by position.
    pub fn sorted_snps(&self) -> Vec<Snp> {
        let mut sorted = self.snps.clone();
        sorted.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.position.cmp(&b.position))
        });
        sorted
    }
}

/// `fromHaplotypeDatabase`: the `@` header, then one row per SNP, grouped by anchor.
///
/// A row whose anchor is empty or names itself STARTS a block; every other row is held back and
/// added to its anchor's block afterwards, which is why an anchor that names no row of its own is
/// only caught at the end. The blocks come out in the anchors' insertion order.
pub fn parse_haplotype_database(text: &str) -> Result<Vec<HaplotypeBlock>, String> {
    let mut lines = text.lines().peekable();
    let mut header = String::new();
    while let Some(line) = lines.peek() {
        if !line.starts_with('@') {
            break;
        }
        header.push_str(line);
        header.push('\n');
        lines.next();
    }
    if header.is_empty() {
        return Err(MISSING_HEADER_PREFIX.to_string());
    }
    let mut anchors: Vec<String> = Vec::new();
    let mut blocks: BTreeMap<String, HaplotypeBlock> = BTreeMap::new();
    let mut held: Vec<(String, Snp)> = Vec::new();
    for line in lines {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 || fields.len() > 8 {
            return Err(format!(
                "{INVALID_RECORD_PREFIX}{} fields: {line}",
                fields.len()
            ));
        }
        let snp = Snp {
            name: fields[2].to_string(),
            chromosome: fields[0].to_string(),
            position: fields[1].parse().map_err(|_| "a position")?,
            major_allele: fields[3].as_bytes()[0],
            minor_allele: fields[4].as_bytes()[0],
            minor_allele_frequency: fields[5].parse().map_err(|_| "a frequency")?,
            panels: fields
                .get(7)
                .filter(|panels| !panels.is_empty())
                .map(|panels| panels.split(',').map(str::to_string).collect()),
        };
        let anchor = fields.get(6).copied().unwrap_or("");
        if anchor.trim().is_empty() || anchor == snp.name {
            anchors.push(snp.name.clone());
            let mut block = HaplotypeBlock::default();
            block.add_snp(snp.clone())?;
            blocks.insert(snp.name, block);
        } else {
            held.push((anchor.to_string(), snp));
        }
    }
    for (anchor, snp) in held {
        let block = blocks
            .get_mut(&anchor)
            .ok_or_else(|| format!("{NO_HAPLOTYPE_PREFIX}{anchor}"))?;
        block.add_snp(snp)?;
    }
    Ok(anchors
        .into_iter()
        .filter_map(|anchor| blocks.remove(&anchor))
        .collect())
}

/// One row of the file `writeToFile` writes.
#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub chromosome: String,
    pub position: i32,
    pub name: String,
    pub major_allele: u8,
    pub minor_allele: u8,
    pub minor_allele_frequency: f64,
    /// The block's FIRST SNP by position, or nothing for that first SNP itself.
    pub anchor: Option<String>,
    pub panels: Option<String>,
}

/// `writeToFile`'s rows: the anchor column is REWRITTEN rather than carried.
///
/// Each block's first SNP by position gets an empty anchor and every later row gets that first
/// row's name, whatever the input named. The rows are then sorted by the dictionary's contig
/// order and the position, so a block's SNPs need not stay together.
pub fn rows(blocks: &[HaplotypeBlock], sequence_order: &[String]) -> Vec<Row> {
    let mut entries: Vec<Row> = Vec::new();
    for block in blocks {
        let mut anchor: Option<String> = None;
        for snp in block.sorted_snps() {
            entries.push(Row {
                chromosome: snp.chromosome.clone(),
                position: snp.position,
                name: snp.name.clone(),
                major_allele: snp.major_allele,
                minor_allele: snp.minor_allele,
                minor_allele_frequency: snp.minor_allele_frequency,
                anchor: anchor.clone(),
                panels: snp.panels.as_ref().map(|panels| panels.join(",")),
            });
            if anchor.is_none() {
                anchor = Some(snp.name);
            }
        }
    }
    let index = |contig: &str| {
        sequence_order
            .iter()
            .position(|name| name == contig)
            .map_or(-1, |position| position as i32)
    };
    entries.sort_by(|a, b| {
        index(&a.chromosome)
            .cmp(&index(&b.chromosome))
            .then(a.position.cmp(&b.position))
            .then(a.name.cmp(&b.name))
    });
    entries
}

/// `FormatUtil.format(double)`: at most six fraction digits, no grouping, trailing zeros dropped.
pub fn format_frequency(value: f64) -> String {
    let rendered = format!("{:.6}", value);
    let trimmed = rendered.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-" {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

/// One VCF record `asVcf` builds for one SNP.
#[derive(Debug, Clone, PartialEq)]
pub struct VcfRecord {
    pub chromosome: String,
    pub position: i32,
    pub id: String,
    pub reference: char,
    pub alternate: char,
    pub allele_frequency: f64,
    pub genotype: String,
    /// The block's first SNP's position, for a block of more than one SNP.
    pub phase_set: Option<i32>,
}

/// `asVcf` for one block, against a reference that answers the base at a position.
///
/// REF is the allele that DISAGREES with the reference. The code asks which allele matches the
/// reference base and then writes THAT one as ALT, so a row whose major allele is the reference
/// base and a row whose minor allele is come out with the same REF and ALT, and only the
/// frequency tells them apart: AF is the frequency of the ALT the tool wrote, which is `1 - MAF`
/// when the major allele matched and `MAF` when the minor one did.
///
/// A block of one SNP is unphased and carries no phase set. A block of more carries a phase set
/// on every row, and that phase set is the position of the block's FIRST SNP by position and not
/// of the row named as the anchor. A row inside such a block whose major allele matched has its
/// genotype reversed to `1|0`.
pub fn block_as_vcf(
    block: &HaplotypeBlock,
    reference_base: impl Fn(&str, i32) -> Option<u8>,
) -> Result<Vec<VcfRecord>, String> {
    let sorted = block.sorted_snps();
    let phased = sorted.len() > 1;
    let anchor_position = sorted.first().map(|snp| snp.position);
    let mut records = Vec::with_capacity(sorted.len());
    for snp in &sorted {
        let base = reference_base(&snp.chromosome, snp.position)
            .ok_or_else(|| format!("no reference base at {}:{}", snp.chromosome, snp.position))?;
        let swap = if base.eq_ignore_ascii_case(&snp.major_allele) {
            true
        } else if base.eq_ignore_ascii_case(&snp.minor_allele) {
            false
        } else {
            return Err(format!(
                "{ALLELE_DISAGREEMENT_PREFIX}{}:{}",
                snp.chromosome, snp.position
            ));
        };
        let alleles = snp.allele_string().into_bytes();
        let reference = if swap { alleles[1] } else { alleles[0] };
        let alternate = if swap { alleles[0] } else { alleles[1] };
        let frequency = if swap {
            1.0 - snp.minor_allele_frequency
        } else {
            snp.minor_allele_frequency
        };
        let genotype = if phased && swap {
            "1|0"
        } else if phased {
            "0|1"
        } else {
            "0/1"
        };
        records.push(VcfRecord {
            chromosome: snp.chromosome.clone(),
            position: snp.position,
            id: snp.name.clone(),
            reference: reference.to_ascii_uppercase() as char,
            alternate: alternate.to_ascii_uppercase() as char,
            allele_frequency: frequency,
            genotype: genotype.to_string(),
            phase_set: if phased { anchor_position } else { None },
        });
    }
    Ok(records)
}

/// `writeAsVcf`: every block's records, sorted by the dictionary's contig order and the position.
pub fn as_vcf(
    blocks: &[HaplotypeBlock],
    sequence_order: &[String],
    reference_base: impl Fn(&str, i32) -> Option<u8> + Copy,
) -> Result<Vec<VcfRecord>, String> {
    let mut records = Vec::new();
    for block in blocks {
        records.extend(block_as_vcf(block, reference_base)?);
    }
    let index = |contig: &str| {
        sequence_order
            .iter()
            .position(|name| name == contig)
            .map_or(-1, |position| position as i32)
    };
    records.sort_by(|a, b| {
        index(&a.chromosome)
            .cmp(&index(&b.chromosome))
            .then(a.position.cmp(&b.position))
    });
    Ok(records)
}
