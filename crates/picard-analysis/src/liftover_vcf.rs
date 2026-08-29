//! `LiftoverVcf`: moving variants from one reference to another through a chain file.
//!
//! The arithmetic is the easy half. What the tool is really about is what happens when the new
//! reference disagrees with the old one: a variant whose reference allele the target does not
//! carry cannot simply be renumbered, and the answers are a rejection, or, where the alleles turn
//! out to be the other way round, a lift with the alleles swapped and the genotypes with them.
//!
//! A rejected variant goes to the reject file rather than being dropped, carrying a filter that
//! names the reason and, where the reason is the reference allele, the locus and the alleles that
//! were attempted.
//!
//! Ported from `picard.vcf.LiftoverVcf` and `picard.util.LiftoverUtils` in Picard 3.4.0, and from
//! `htsjdk.samtools.liftover.Chain`, `htsjdk.samtools.liftover.LiftOver` and
//! `htsjdk.variant.vcf.VCFEncoder` in htsjdk 4.2.0.

/// One continuous block of a chain, in zero-based half-open coordinates on both sides.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Block {
    pub from_start: i32,
    pub to_start: i32,
    pub length: i32,
}

impl Block {
    pub fn from_end(&self) -> i32 {
        self.from_start + self.length
    }
    pub fn to_end(&self) -> i32 {
        self.to_start + self.length
    }
}

/// One chain: the source is the `t` side and the target is the `q` side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chain {
    pub from_name: String,
    pub to_name: String,
    pub to_size: i32,
    /// A chain whose target side is written on the reverse strand.
    pub to_opposite_strand: bool,
    pub id: i32,
    pub blocks: Vec<Block>,
}

/// An interval, one-based and inclusive, with the strand it ended up on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
}

impl Interval {
    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }
}

/// Read a chain file.
///
/// The header names both sides and the block lines walk them forward together: each line is the
/// block's length and then the gap on each side before the next block, and the last line of a
/// chain is a length on its own.
pub fn parse_chains(text: &str) -> Vec<Chain> {
    let mut chains = Vec::new();
    let mut current: Option<Chain> = None;
    let mut from_block_start = 0;
    let mut to_block_start = 0;
    let mut saw_last_line = false;
    for line in text.lines() {
        if line.trim().is_empty() {
            if let Some(chain) = current.take() {
                chains.push(chain);
            }
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields[0] == "chain" {
            from_block_start = fields[5].parse().expect("a start");
            to_block_start = fields[10].parse().expect("a start");
            saw_last_line = false;
            current = Some(Chain {
                from_name: fields[2].to_string(),
                to_name: fields[7].to_string(),
                to_size: fields[8].parse().expect("a size"),
                to_opposite_strand: fields[9] == "-",
                id: fields[12].parse().expect("an id"),
                blocks: Vec::new(),
            });
            continue;
        }
        let Some(chain) = current.as_mut() else {
            continue;
        };
        if saw_last_line {
            continue;
        }
        let size: i32 = fields[0].parse().expect("a block length");
        chain.blocks.push(Block {
            from_start: from_block_start,
            to_start: to_block_start,
            length: size,
        });
        if fields.len() == 1 {
            saw_last_line = true;
        } else {
            from_block_start += fields[1].parse::<i32>().expect("a gap") + size;
            to_block_start += fields[2].parse::<i32>().expect("a gap") + size;
        }
    }
    if let Some(chain) = current.take() {
        chains.push(chain);
    }
    chains
}

/// How much of an interval one chain covers, and where in its blocks.
#[derive(Debug, Clone, PartialEq, Eq)]
struct TargetIntersection<'a> {
    chain: &'a Chain,
    intersection_length: i32,
    start_offset: i32,
    offset_from_end: i32,
    first_block: usize,
    last_block: usize,
}

fn target_intersection<'a>(
    chain: &'a Chain,
    interval: &Interval,
) -> Option<TargetIntersection<'a>> {
    if chain.from_name != interval.contig {
        return None;
    }
    let start = interval.start - 1;
    let end = interval.end;
    let mut intersection_length = 0;
    let mut first_block = None;
    let mut last_block = 0;
    let mut start_offset = 0;
    let mut offset_from_end = 0;
    for (index, block) in chain.blocks.iter().enumerate() {
        if block.from_start >= end {
            break;
        }
        if block.from_end() <= start {
            continue;
        }
        if first_block.is_none() {
            first_block = Some(index);
            start_offset = if start > block.from_start {
                start - block.from_start
            } else {
                0
            };
        }
        last_block = index;
        offset_from_end = if block.from_end() > end {
            block.from_end() - end
        } else {
            0
        };
        intersection_length += end.min(block.from_end()) - start.max(block.from_start);
    }
    let first_block = first_block?;
    if intersection_length == 0 {
        return None;
    }
    Some(TargetIntersection {
        chain,
        intersection_length,
        start_offset,
        offset_from_end,
        first_block,
        last_block,
    })
}

/// Where an interval of the source lands on the target, if it lands at all.
///
/// A reversed chain writes its target coordinates on the other strand, so the position is counted
/// from the far end of the contig and the strand comes back flipped.
pub fn lift_over(chains: &[Chain], interval: &Interval, min_match: f64) -> Option<Interval> {
    let minimum = min_match * f64::from(interval.length());
    let mut best: Option<TargetIntersection> = None;
    for chain in chains {
        let Some(candidate) = target_intersection(chain, interval) else {
            continue;
        };
        if f64::from(candidate.intersection_length) < minimum {
            continue;
        }
        if best.is_some() {
            // In basic liftover, more than one hit is no hit.
            return None;
        }
        best = Some(candidate);
    }
    let found = best?;
    let chain = found.chain;
    let mut to_start = chain.blocks[found.first_block].to_start + found.start_offset;
    let mut to_end = chain.blocks[found.last_block].to_end() - found.offset_from_end;
    if to_end <= to_start || to_start < 0 {
        return None;
    }
    if chain.to_opposite_strand {
        let negative_start = chain.to_size - to_end;
        let negative_end = chain.to_size - to_start;
        to_start = negative_start;
        to_end = negative_end;
    }
    Some(Interval {
        contig: chain.to_name.clone(),
        start: to_start + 1,
        end: to_end,
        negative_strand: if chain.to_opposite_strand {
            !interval.negative_strand
        } else {
            interval.negative_strand
        },
    })
}

/// One INFO entry, which is either a key with a value or a flag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Attribute {
    pub key: String,
    /// `None` is a flag, written as the key alone.
    pub value: Option<String>,
}

/// One variant, as much of it as the tool moves.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    pub position: i32,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
    pub quality: String,
    pub filters: Vec<String>,
    pub attributes: Vec<Attribute>,
    pub format: Vec<String>,
    /// One per sample, each a list of fields; the first is the genotype.
    pub samples: Vec<Vec<String>>,
}

impl Variant {
    pub fn end(&self) -> i32 {
        self.position + self.reference.len() as i32 - 1
    }
    pub fn is_snp(&self) -> bool {
        self.reference.len() == 1 && self.alternates.iter().all(|allele| allele.len() == 1)
    }
    pub fn is_biallelic(&self) -> bool {
        self.alternates.len() == 1
    }
    pub fn attribute(&self, key: &str) -> Option<&Attribute> {
        self.attributes.iter().find(|entry| entry.key == key)
    }
    fn set_attribute(&mut self, key: &str, value: Option<String>) {
        self.remove_attribute(key);
        self.attributes.push(Attribute {
            key: key.to_string(),
            value,
        });
    }
    fn remove_attribute(&mut self, key: &str) {
        self.attributes.retain(|entry| entry.key != key);
    }
}

/// The filters a rejected variant carries.
pub const FILTER_NO_TARGET: &str = "NoTarget";
pub const FILTER_MISMATCHING_REF_ALLELE: &str = "MismatchedRefAllele";
pub const FILTER_INDEL_STRADDLES_TWO_INTERVALS: &str = "IndelStraddlesTwoIntervals";
pub const FILTER_CANNOT_LIFTOVER_REV_COMP: &str = "CannotLiftOver";

/// The INFO keys the tool adds.
pub const SWAPPED_ALLELES: &str = "SwappedAlleles";
pub const REV_COMPED_ALLELES: &str = "ReverseComplementedAlleles";
pub const ORIGINAL_CONTIG: &str = "OriginalContig";
pub const ORIGINAL_START: &str = "OriginalStart";
pub const ORIGINAL_ALLELES: &str = "OriginalAlleles";
pub const ATTEMPTED_LOCUS: &str = "AttemptedLocus";
pub const ATTEMPTED_ALLELES: &str = "AttemptedAlleles";

/// The INFO fields a swap makes false, and the ones it turns into their complement.
pub const DEFAULT_TAGS_TO_REVERSE: [&str; 1] = ["AF"];
pub const DEFAULT_TAGS_TO_DROP: [&str; 1] = ["MAX_AF"];

/// The reverse complement of one allele.
pub fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            other => other,
        })
        .collect()
}

/// Move one variant onto the target, without asking whether the reference agrees.
///
/// A reversed target complements every allele and says so in the INFO; a forward one takes the
/// attribute back off, in case the variant was lifted once before.
pub fn lift_variant(
    source: &Variant,
    target: &Interval,
    write_original_position: bool,
    write_original_alleles: bool,
) -> Variant {
    let mut lifted = source.clone();
    lifted.contig = target.contig.clone();
    lifted.position = target.start;
    if target.negative_strand {
        lifted.reference = reverse_complement(&source.reference);
        lifted.alternates = source
            .alternates
            .iter()
            .map(|allele| reverse_complement(allele))
            .collect();
    }

    // A variant is never carried over as already swapped, and it is marked as complemented only
    // where this lift complemented it.
    lifted.remove_attribute(SWAPPED_ALLELES);
    if target.negative_strand {
        lifted.set_attribute(REV_COMPED_ALLELES, None);
    } else {
        lifted.remove_attribute(REV_COMPED_ALLELES);
    }

    if write_original_position {
        lifted.set_attribute(ORIGINAL_CONTIG, Some(source.contig.clone()));
        lifted.set_attribute(ORIGINAL_START, Some(source.position.to_string()));
    }
    // The original alleles are recorded only where they CHANGED, so a plain lift records nothing
    // however loudly it is asked to.
    if write_original_alleles
        && (lifted.reference != source.reference || lifted.alternates != source.alternates)
    {
        let mut alleles = vec![source.reference.clone()];
        alleles.extend(source.alternates.clone());
        lifted.set_attribute(ORIGINAL_ALLELES, Some(alleles.join(",")));
    }
    lifted
}

/// `VCFEncoder.formatVCFDouble`, which is how a computed number reaches the file.
pub fn format_vcf_double(value: f64) -> String {
    if value < 1.0 {
        if value < 0.01 {
            if value.abs() >= 1e-20 {
                return format!("{value:.3e}");
            }
            return "0.00".to_string();
        }
        return format!("{value:.3}");
    }
    format!("{value:.2}")
}

/// Turn a bi-allelic SNP the other way round.
///
/// The genotypes follow the alleles rather than being re-sorted, so a `0/1` becomes `1/0`, and the
/// depths and likelihoods are reversed with them. An INFO field that is a fraction of the alternate
/// becomes its complement, which is where a computed number enters a file whose other numbers came
/// through as text.
pub fn swap_ref_alt(variant: &Variant, tags_to_reverse: &[&str], tags_to_drop: &[&str]) -> Variant {
    let mut swapped = variant.clone();
    swapped.set_attribute(SWAPPED_ALLELES, None);
    swapped.reference = variant.alternates[0].clone();
    swapped.alternates = vec![variant.reference.clone()];

    let genotype_index = variant.format.iter().position(|field| field == "GT");
    let ad_index = variant.format.iter().position(|field| field == "AD");
    let pl_index = variant.format.iter().position(|field| field == "PL");
    for sample in &mut swapped.samples {
        if let Some(index) = genotype_index {
            if let Some(genotype) = sample.get(index) {
                let separator = if genotype.contains('|') { '|' } else { '/' };
                let swapped_alleles: Vec<String> = genotype
                    .split(['/', '|'])
                    .map(|allele| match allele {
                        "0" => "1".to_string(),
                        "1" => "0".to_string(),
                        other => other.to_string(),
                    })
                    .collect();
                sample[index] = swapped_alleles.join(&separator.to_string());
            }
        }
        for index in [ad_index, pl_index].into_iter().flatten() {
            if let Some(field) = sample.get(index) {
                let mut values: Vec<&str> = field.split(',').collect();
                values.reverse();
                sample[index] = values.join(",");
            }
        }
    }

    for key in tags_to_drop {
        swapped.remove_attribute(key);
    }
    for key in tags_to_reverse {
        let Some(attribute) = swapped.attribute(key) else {
            continue;
        };
        let Some(text) = attribute.value.clone() else {
            continue;
        };
        if text == "." {
            continue;
        }
        let value: f64 = text.parse().unwrap_or(-1.0);
        swapped.set_attribute(key, Some(format_vcf_double(1.0 - value)));
    }
    swapped
}

/// What one variant became.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    /// Written to the output, at its new coordinates.
    Lifted(Variant),
    /// Written to the reject file, filtered and carrying whatever the reason recorded.
    Rejected(Variant),
}

/// What a run was asked for.
#[derive(Debug, Clone, PartialEq)]
pub struct Options {
    pub recover_swapped_ref_alt: bool,
    pub write_original_position: bool,
    pub write_original_alleles: bool,
    pub min_match: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            recover_swapped_ref_alt: false,
            write_original_position: false,
            write_original_alleles: false,
            min_match: 1.0,
        }
    }
}

/// Reject one variant, which keeps the SOURCE record and adds a filter to it.
fn reject(source: &Variant, filter: &str) -> Variant {
    let mut rejected = source.clone();
    rejected.filters = vec![filter.to_string()];
    rejected
}

/// Whether the target reference carries what the lifted variant says it does, and what to do about
/// it if it does not.
pub fn try_to_add(
    lifted: &Variant,
    reference: &[u8],
    source: &Variant,
    options: &Options,
) -> Outcome {
    let start = (lifted.position - 1) as usize;
    let end = lifted.end() as usize;
    let carried = String::from_utf8_lossy(&reference[start..end]).to_string();
    if carried.eq_ignore_ascii_case(&lifted.reference) {
        return Outcome::Lifted(lifted.clone());
    }

    // The reference disagrees. On a bi-allelic SNP whose ALT is what the target carries, the
    // alleles are simply the other way round, and recovering that is a lift rather than a
    // rejection.
    if lifted.is_biallelic()
        && lifted.is_snp()
        && carried.eq_ignore_ascii_case(&lifted.alternates[0])
        && options.recover_swapped_ref_alt
    {
        return Outcome::Lifted(swap_ref_alt(
            lifted,
            &DEFAULT_TAGS_TO_REVERSE,
            &DEFAULT_TAGS_TO_DROP,
        ));
    }

    // Otherwise it is rejected, and what is written is the SOURCE record with the locus and the
    // alleles that were attempted recorded on it.
    let mut rejected = reject(source, FILTER_MISMATCHING_REF_ALLELE);
    rejected.set_attribute(
        ATTEMPTED_LOCUS,
        Some(format!(
            "{}:{}-{}",
            lifted.contig,
            lifted.position,
            lifted.end()
        )),
    );
    rejected.set_attribute(
        ATTEMPTED_ALLELES,
        // The reference allele is written with the star that marks it as one.
        Some(format!(
            "{}*->{}",
            lifted.reference,
            lifted.alternates.join(",")
        )),
    );
    Outcome::Rejected(rejected)
}

/// Lift one variant, from the chain to the answer.
pub fn liftover(
    source: &Variant,
    chains: &[Chain],
    reference: &[u8],
    options: &Options,
) -> Outcome {
    let interval = Interval {
        contig: source.contig.clone(),
        start: source.position,
        end: source.end(),
        negative_strand: false,
    };
    let Some(target) = lift_over(chains, &interval, options.min_match) else {
        return Outcome::Rejected(reject(source, FILTER_NO_TARGET));
    };
    // An interval that grew or shrank straddled two blocks, and there is no telling what the
    // alleles should be afterwards.
    if source.reference.len() as i32 != target.length() {
        return Outcome::Rejected(reject(source, FILTER_INDEL_STRADDLES_TWO_INTERVALS));
    }
    let lifted = lift_variant(
        source,
        &target,
        options.write_original_position,
        options.write_original_alleles,
    );
    try_to_add(&lifted, reference, source, options)
}

/// A whole file: what was lifted, sorted by the TARGET's coordinates, and what was rejected, in
/// the order it was read.
pub fn run(
    sources: &[Variant],
    chains: &[Chain],
    reference: &[u8],
    options: &Options,
) -> (Vec<Variant>, Vec<Variant>) {
    let mut lifted = Vec::new();
    let mut rejected = Vec::new();
    for source in sources {
        match liftover(source, chains, reference, options) {
            Outcome::Lifted(variant) => lifted.push(variant),
            Outcome::Rejected(variant) => rejected.push(variant),
        }
    }
    lifted.sort_by(|left, right| {
        left.contig
            .cmp(&right.contig)
            .then(left.position.cmp(&right.position))
    });
    (lifted, rejected)
}

/// One variant as the file writes it, INFO keys in the order the encoder sorts them.
pub fn render(variant: &Variant) -> String {
    let mut attributes: Vec<&Attribute> = variant.attributes.iter().collect();
    attributes.sort_by(|left, right| left.key.cmp(&right.key));
    let info = if attributes.is_empty() {
        ".".to_string()
    } else {
        attributes
            .iter()
            .map(|entry| match &entry.value {
                Some(value) => format!("{}={}", entry.key, value),
                None => entry.key.clone(),
            })
            .collect::<Vec<_>>()
            .join(";")
    };
    let filter = if variant.filters.is_empty() {
        "PASS".to_string()
    } else {
        variant.filters.join(";")
    };
    let mut line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        variant.contig,
        variant.position,
        variant.id,
        variant.reference,
        variant.alternates.join(","),
        variant.quality,
        filter,
        info
    );
    if !variant.format.is_empty() {
        line.push('\t');
        line.push_str(&variant.format.join(":"));
        for sample in &variant.samples {
            line.push('\t');
            line.push_str(&sample.join(":"));
        }
    }
    line
}
