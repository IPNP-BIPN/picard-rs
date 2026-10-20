//! `LiftoverVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.vcf.LiftoverVcf.doWork` and `picard.util.LiftoverUtils` at tag 3.4.0 over
//! htsjdk-rs's VCF reader, encoder and `LiftOver`:
//!
//! * the `CREATE_INDEX`/`DISABLE_SORT` refusal, which is a logged error and an exit code of one;
//! * the two headers: the output's is the input's lines without `##reference`, the target's
//!   contig lines, the `Swapped`/`ReverseComplemented` flags (and the three original-position
//!   lines under `WRITE_ORIGINAL_POSITION`), a new `##reference`, and the samples sorted; the
//!   reject file's is the input's with the four filters and the two attempted-locus lines;
//! * per record: no target, a target of another length (an indel across two blocks), a contig the
//!   reference lacks (an error unless `WARN_ON_MISSING_CONTIG`), then `liftVariant` (the reverse
//!   strand complements every allele, re-checks and left-aligns an indel, and maps the genotypes
//!   through the new alleles), then `tryToAddVariant`'s reference check with the swapped-SNP
//!   recovery;
//! * rejects go out as they come, carrying the source record with the filter ADDED to its own;
//!   the lifted records are sorted by the target dictionary and start unless `DISABLE_SORT`, and
//!   every write refuses an INFO key its header lacks unless `ALLOW_MISSING_FIELDS_IN_HEADER`.
//!
//! The input is streamed. Both outputs are held as their encoded lines in a `SortingCollection`
//! that spills past `MAX_RECORDS_IN_RAM`, as Picard's `VCFRecordCodec` does, and each line is
//! decoded again on its way to the writer. Neither file is opened before the whole input is read,
//! as with the whole-file read, and the rejects are still written first.
//!
//! A record untouched by the lift keeps its genotype block as the file had it, when the file's
//! samples were already sorted: htsjdk's lazy genotypes, which nothing here decodes.

use std::collections::HashMap;
use std::sync::OnceLock;

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::interval::Interval;
use htsjdk_bam::liftover::LiftOver;
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::attributes::{as_double, as_string};
use htsjdk_vcf::encoder::{MissingFields, VcfEncoder};
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::header_parse::VcfVersion;
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};
use picard_analysis::liftover_vcf::LineSpool;
use picard_analysis::metrics_cli::{absolute, Args};
use picard_analysis::sort_sam::DEFAULT_MAX_RECORDS_IN_RAM;
use picard_analysis::vcf_io::{
    add_other_meta_data_line, die, is_bcf, log_error, parse_sam_dictionary,
    set_sequence_dictionary, Record, Sequence, VcfReader, VcfWriter, WriteError,
};

const TOOL: &str = "LiftoverVcf";

const FILTER_CANNOT_LIFTOVER_REV_COMP: &str = "CannotLiftOver";
const FILTER_NO_TARGET: &str = "NoTarget";
const FILTER_MISMATCHING_REF_ALLELE: &str = "MismatchedRefAllele";
const FILTER_INDEL_STRADDLES_TWO_INTERVALS: &str = "IndelStraddlesMultipleIntervals";
const ORIGINAL_CONTIG: &str = "OriginalContig";
const ORIGINAL_START: &str = "OriginalStart";
const ORIGINAL_ALLELES: &str = "OriginalAlleles";
const ATTEMPTED_LOCUS: &str = "AttemptedLocus";
const ATTEMPTED_ALLELES: &str = "AttemptedAlleles";
const SWAPPED_ALLELES: &str = "SwappedAlleles";
const REV_COMPED_ALLELES: &str = "ReverseComplementedAlleles";
const END_KEY: &str = "END";

/// The input, re-read for a failure further on once something has been refused.
static INPUT: OnceLock<String> = OnceLock::new();

/// The whole-file read's failure, when the input has one: it came before anything the records
/// could be refused for.
fn read_failure_first() {
    if let Some(Ok(mut vcf)) = INPUT.get().map(|path| VcfReader::open(path)) {
        if let Err(exception) = vcf.check() {
            die(&exception);
        }
    }
}

/// [`die`] for a refusal once the records are being read.
fn refuse(exception: &str) -> ! {
    read_failure_first();
    die(exception)
}

fn has_attribute(vc: &VariantContext, key: &str) -> bool {
    vc.attributes.iter().any(|(k, _)| k == key)
}

fn set_attribute(vc: &mut VariantContext, key: &str, value: Value) {
    match vc.attributes.iter_mut().find(|(k, _)| k == key) {
        Some(entry) => entry.1 = value,
        None => vc.attributes.push((key.to_string(), value)),
    }
}

fn remove_attribute(vc: &mut VariantContext, key: &str) {
    vc.attributes.retain(|(k, _)| k != key);
}

/// `VariantContextBuilder.filter(reason)`: added to whatever the record already carries.
fn add_filter(vc: &mut VariantContext, reason: &str) {
    let filters = vc.filters.get_or_insert_with(Vec::new);
    if !filters.iter().any(|f| f == reason) {
        filters.push(reason.to_string());
    }
}

/// `SequenceUtil.complement`, case kept.
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'T' => b'A',
        b't' => b'a',
        other => other,
    }
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|b| complement(*b)).collect()
}

fn is_span_del(allele: &Allele) -> bool {
    !allele.is_symbolic() && !allele.is_no_call() && allele.display_string() == "*"
}

/// `LiftoverUtils.reverseComplement(Allele)`.
fn reverse_complement_allele(allele: &Allele) -> Allele {
    if allele.is_symbolic() || allele.is_no_call() || is_span_del(allele) {
        return allele.clone();
    }
    let bases = reverse_complement(allele.base_string().as_bytes());
    Allele::create(&bases, allele.is_reference())
        .unwrap_or_else(|e| refuse(&format!("java.lang.IllegalArgumentException: {e:?}")))
}

/// `VariantContext.isSNP()`: every allele one base, at least one alternate.
fn is_snp(vc: &VariantContext) -> bool {
    vc.alleles.len() > 1
        && vc
            .alleles
            .iter()
            .all(|a| !a.is_symbolic() && !is_span_del(a) && a.len() == 1)
}

fn is_biallelic(vc: &VariantContext) -> bool {
    vc.alleles.len() == 2
}

/// `isIndelForLiftover`.
fn is_indel_for_liftover(vc: &VariantContext) -> bool {
    if vc.reference().len() != 1 {
        return vc
            .alternate_alleles()
            .iter()
            .any(|a| !a.is_symbolic() && !is_span_del(a));
    }
    vc.alleles
        .iter()
        .filter(|a| !a.is_symbolic() && !is_span_del(a))
        .any(|a| a.len() != 1)
}

/// `referenceAlleleDiffersFromReferenceForIndel`.
fn reference_differs(alleles: &[Allele], reference: &[u8], start: i64, end: i64) -> bool {
    let ref_string = slice(reference, start, end);
    let ref_allele = alleles
        .iter()
        .find(|a| a.is_reference())
        .unwrap_or_else(|| {
            refuse("java.lang.IllegalStateException: Error: no reference allele was present")
        });
    !ref_string.eq_ignore_ascii_case(&ref_allele.base_string())
}

/// `StringUtil.bytesToString(bases, start - 1, end - start + 1)`.
fn slice(reference: &[u8], start: i64, end: i64) -> String {
    let from = (start - 1).max(0) as usize;
    let to = (end.max(start - 1) as usize).min(reference.len());
    String::from_utf8_lossy(&reference[from.min(to)..to]).into_owned()
}

/// `leftAlignVariant`: returns the new start, stop and alleles.
fn left_align(
    start: i64,
    end: i64,
    alleles: &[Allele],
    reference: &[u8],
) -> (i64, i64, Vec<Allele>) {
    // The map is keyed by allele; the order it is walked in changes nothing, every step being
    // applied to every allele.
    let mut bases: Vec<(Allele, Vec<u8>)> = alleles
        .iter()
        .filter(|a| !is_span_del(a) && !a.is_symbolic())
        .map(|a| (a.clone(), a.base_string().into_bytes()))
        .collect();
    let mut the_start = start;
    let mut the_end = end;
    loop {
        let (old_start, old_end) = (the_start, the_end);
        let mut last: Vec<u8> = bases.iter().map(|(_, b)| b[b.len() - 1]).collect();
        last.sort_unstable();
        last.dedup();
        if last.len() == 1 && the_end > 1 {
            for (_, b) in bases.iter_mut() {
                b.pop();
            }
            the_end -= 1;
        }
        if bases.iter().any(|(_, b)| b.is_empty()) {
            let (extra, left) = if the_start > 1 {
                let base = reference[(the_start - 2) as usize];
                the_start -= 1;
                (base, true)
            } else {
                let base = reference[the_end as usize];
                the_end += 1;
                (base, false)
            };
            for (_, b) in bases.iter_mut() {
                if left {
                    b.insert(0, extra);
                } else {
                    b.push(extra);
                }
            }
        }
        if the_start == old_start && the_end == old_end {
            break;
        }
    }
    loop {
        let all_long = bases.iter().all(|(_, b)| b.len() >= 2);
        let mut first: Vec<u8> = bases.iter().map(|(_, b)| b[0]).collect();
        first.sort_unstable();
        first.dedup();
        if !(all_long && first.len() == 1) {
            break;
        }
        for (_, b) in bases.iter_mut() {
            b.remove(0);
        }
        the_start += 1;
    }
    let fixed: Vec<Allele> = alleles
        .iter()
        .map(|a| {
            if is_span_del(a) || a.is_symbolic() {
                return a.clone();
            }
            let new_bases = &bases.iter().find(|(k, _)| k == a).expect("allele").1;
            Allele::create(new_bases, a.is_reference())
                .unwrap_or_else(|e| refuse(&format!("java.lang.IllegalArgumentException: {e:?}")))
        })
        .collect();
    (the_start, the_end, fixed)
}

/// `fixGenotypes`: every called allele mapped to its counterpart in the new list.
fn fix_genotypes(
    originals: &[Genotype],
    original_alleles: &[Allele],
    new_alleles: &[Allele],
) -> Option<Vec<Genotype>> {
    if original_alleles == new_alleles {
        return None;
    }
    let map: HashMap<&Allele, &Allele> = original_alleles.iter().zip(new_alleles.iter()).collect();
    let fixed: Vec<Genotype> = originals
        .iter()
        .map(|g| {
            let mut g = g.clone();
            g.alleles = g
                .alleles
                .iter()
                .map(|a| {
                    if a.is_no_call() {
                        a.clone()
                    } else {
                        (*map.get(a).unwrap_or_else(|| {
                            refuse(&format!(
                                "java.lang.IllegalStateException: Allele not found: {}",
                                a.display_string()
                            ))
                        }))
                        .clone()
                    }
                })
                .collect();
            g
        })
        .collect();
    Some(fixed)
}

/// `LiftoverUtils.liftVariant`, or `None` where an indel on the reverse strand does not fit.
///
/// The flag says whether the genotypes were rebuilt, which is what decodes a lazy record.
fn lift_variant(
    source: &VariantContext,
    target: &Interval,
    reference: &[u8],
    write_original_position: bool,
    write_original_alleles: bool,
) -> Option<(VariantContext, bool)> {
    let mut vc = source.clone();
    let mut touched = false;
    vc.contig = target.contig.clone();
    vc.start = i64::from(target.start);
    vc.stop = i64::from(target.end);
    if target.negative_strand {
        remove_attribute(&mut vc, END_KEY);
        let original = source.alleles.clone();
        vc.alleles = original.iter().map(reverse_complement_allele).collect();
        if is_indel_for_liftover(source) {
            if reference_differs(&vc.alleles, reference, vc.start, vc.stop) {
                return None;
            }
            let (start, stop, alleles) = left_align(vc.start, vc.stop, &vc.alleles, reference);
            vc.start = start;
            vc.stop = stop;
            vc.alleles = alleles;
        }
        if let Some(fixed) = fix_genotypes(&source.genotypes, &original, &vc.alleles) {
            vc.genotypes = fixed.into();
            touched = true;
        }
    }
    // `builder.filters(source.getFilters())`: an unfiltered source comes out as PASS.
    vc.filters = Some(source.filters.clone().unwrap_or_default());
    vc.log10_p_error = source.log10_p_error;
    if has_attribute(source, END_KEY) && vc.alleles.iter().all(|a| !a.is_symbolic()) {
        let stop = vc.stop;
        set_attribute(&mut vc, END_KEY, Value::Int(stop));
    } else {
        remove_attribute(&mut vc, END_KEY);
    }
    remove_attribute(&mut vc, SWAPPED_ALLELES);
    if target.negative_strand {
        set_attribute(&mut vc, REV_COMPED_ALLELES, Value::Bool(true));
    } else {
        remove_attribute(&mut vc, REV_COMPED_ALLELES);
    }
    vc.id = source.id.clone();
    if write_original_position {
        set_attribute(&mut vc, ORIGINAL_CONTIG, Value::Str(source.contig.clone()));
        set_attribute(&mut vc, ORIGINAL_START, Value::Int(source.start));
    }
    if write_original_alleles && source.alleles != vc.alleles {
        let alleles = source
            .alleles
            .iter()
            .map(|a| {
                Value::Str(if a.is_no_call() {
                    ".".to_string()
                } else {
                    a.display_string()
                })
            })
            .collect();
        set_attribute(&mut vc, ORIGINAL_ALLELES, Value::List(alleles));
    }
    Some((vc, touched))
}

/// `LiftoverUtils.swapRefAlt`.
fn swap_ref_alt(vc: &VariantContext, to_reverse: &[String], to_drop: &[String]) -> VariantContext {
    let mut swapped = vc.clone();
    set_attribute(&mut swapped, SWAPPED_ALLELES, Value::Bool(true));
    let new_ref = Allele::create(vc.alleles[1].base_string().as_bytes(), true)
        .unwrap_or_else(|e| refuse(&format!("java.lang.IllegalArgumentException: {e:?}")));
    let new_alt = Allele::create(vc.alleles[0].base_string().as_bytes(), false)
        .unwrap_or_else(|e| refuse(&format!("java.lang.IllegalArgumentException: {e:?}")));
    swapped.alleles = vec![new_ref.clone(), new_alt.clone()];
    let genotypes: Vec<Genotype> = vc
        .genotypes
        .iter()
        .map(|g| {
            let mut g = g.clone();
            g.alleles = g
                .alleles
                .iter()
                .map(|a| {
                    if a.is_no_call() {
                        a.clone()
                    } else if *a == vc.alleles[0] {
                        new_alt.clone()
                    } else {
                        new_ref.clone()
                    }
                })
                .collect();
            g.ad = match &g.ad {
                Some(ad) if ad.len() == 2 => Some(vec![ad[1], ad[0]]),
                _ => None,
            };
            g.pl = match &g.pl {
                Some(pl) if pl.len() == 3 => Some(vec![pl[2], pl[1], pl[0]]),
                _ => None,
            };
            g
        })
        .collect();
    swapped.genotypes = genotypes.into();
    for (key, value) in &vc.attributes {
        if to_drop.contains(key) {
            remove_attribute(&mut swapped, key);
        } else if to_reverse.contains(key) && as_string(Some(value), "") != "." {
            let x = as_double(Some(value), -1.0)
                .unwrap_or_else(|e| refuse(&format!("{}: {}", e.class(), e.message())));
            set_attribute(&mut swapped, key, Value::Double(1.0 - x));
        }
    }
    swapped
}

/// `VCFHeader.addMetaDataLine` for an INFO or FILTER line: kept out when one of the same kind and
/// ID is already there.
fn add_line(header: &mut VcfHeader, line: HeaderLine) {
    let same = |existing: &HeaderLine| match (existing, &line) {
        (
            HeaderLine::Compound { key: a, id: x, .. },
            HeaderLine::Compound { key: b, id: y, .. },
        ) => a == b && x == y,
        (HeaderLine::Filter { id: x, .. }, HeaderLine::Filter { id: y, .. }) => x == y,
        _ => false,
    };
    if !header.lines.iter().any(same) {
        header.lines.push(line);
    }
}

/// `VCFEncoder.encode` of one record as `VCFWriter.add` makes it: a record still lazy has its site
/// columns encoded over the header without samples and its genotype block copied, as
/// `vcf_io::VcfWriter` does.
struct LineEncoder {
    header: VcfHeader,
    sites: VcfHeader,
}

impl LineEncoder {
    fn new(header: &VcfHeader) -> Self {
        LineEncoder {
            header: header.clone(),
            sites: VcfHeader {
                lines: header.lines.clone(),
                samples: Vec::new(),
            },
        }
    }

    fn encode(
        &self,
        record: &Record,
        missing: MissingFields,
        out: &mut String,
    ) -> Result<(), WriteError> {
        out.clear();
        match &record.lazy_genotypes {
            Some(block) => {
                let mut site = record.variant.clone();
                site.genotypes.clear();
                VcfEncoder::new(&self.sites)
                    .with_missing_fields(missing)
                    .encode_into(&site, out)?;
                out.push('\t');
                out.push_str(block);
            }
            None => {
                record.variant.genotypes.decode();
                VcfEncoder::new(&self.header)
                    .with_missing_fields(missing)
                    .encode_into(&record.variant, out)?
            }
        }
        Ok(())
    }

    /// The line decoded again: its site columns parsed, its genotype block kept to be copied.
    fn decode(&self, line: &[u8], version: VcfVersion) -> Record {
        let line = std::str::from_utf8(line).expect("a line this tool encoded");
        let decoded = htsjdk_vcf::record_parse::decode_line(line, &self.header, 0, version)
            .ok()
            .flatten()
            .expect("a line this tool encoded decodes");
        Record {
            variant: decoded.variant,
            lazy_genotypes: decoded.genotype_block,
            fully_decoded: false,
        }
    }
}

/// One output's records until the input is read: encoded lines, or for a `.bcf` the records
/// themselves.
enum Held {
    Lines(LineSpool),
    // shortcut: a BCF output keeps its records in memory, since a line decoded again carries
    // strings where the BCF2 writer was given typed values; spill them as lines once the two are
    // shown to write the same bytes.
    Records(Vec<(u32, i64, Record)>),
}

impl Held {
    fn new(path: &str, sort: bool, max_records_in_ram: usize, tmp_dir: &std::path::Path) -> Self {
        if is_bcf(path) {
            Held::Records(Vec::new())
        } else {
            Held::Lines(LineSpool::new(
                sort,
                max_records_in_ram,
                tmp_dir.to_path_buf(),
            ))
        }
    }

    fn add(&mut self, key: (u32, i64), record: Record, line: &str) {
        match self {
            Held::Lines(spool) => spool
                .add(key.0, key.1, line.as_bytes())
                .unwrap_or_else(|e| {
                    refuse(&format!("htsjdk.samtools.util.RuntimeIOException: {e}"))
                }),
            Held::Records(records) => records.push((key.0, key.1, record)),
        }
    }
}

/// `VCFWriter.writeHeader` then `add` per record, the records in order; every encode refuses an
/// INFO key the header lacks unless missing fields are allowed.
fn write_all(
    path: &str,
    header: &VcfHeader,
    held: Held,
    sort: bool,
    missing: MissingFields,
    index_dictionary: Option<&[Sequence]>,
    version: VcfVersion,
) {
    let fail = |error: WriteError| -> ! {
        match error {
            e @ (WriteError::Encode(_) | WriteError::Exception(_)) => die(&e.exception()),
            WriteError::Io(e) => die(&format!("htsjdk.tribble.TribbleException: {e}")),
            // The BCF2 writer's own failure, which it reports as an I/O one.
            WriteError::Other(e) => die(&format!("htsjdk.samtools.util.RuntimeIOException: {e}")),
        }
    };
    let mut writer =
        VcfWriter::create(path, header, index_dictionary, missing).unwrap_or_else(|e| fail(e));
    match held {
        Held::Lines(spool) => {
            let encoder = LineEncoder::new(header);
            let lines = spool
                .into_lines()
                .unwrap_or_else(|e| die(&format!("htsjdk.samtools.util.RuntimeIOException: {e}")));
            for line in lines {
                let line = line.unwrap_or_else(|e| {
                    die(&format!("htsjdk.samtools.util.RuntimeIOException: {e}"))
                });
                writer
                    .add(&encoder.decode(&line, version))
                    .unwrap_or_else(|e| fail(e));
            }
        }
        Held::Records(mut records) => {
            if sort {
                records.sort_by_key(|(contig, start, _)| (*contig, *start));
            }
            for (_, _, record) in &records {
                writer.add(record).unwrap_or_else(|e| fail(e));
            }
        }
    }
    writer.finish().unwrap_or_else(|e| fail(e));
}

fn main() {
    picard_analysis::argument_check::check("LiftoverVcf");
    picard_analysis::output::install();
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("C", "CHAIN"),
        ("R", "REFERENCE_SEQUENCE"),
        ("WMC", "WARN_ON_MISSING_CONTIG"),
        ("LFI", "LOG_FAILED_INTERVALS"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let chain = args.required("CHAIN");
    let reject = args.required("REJECT");
    let reference = args.required("REFERENCE_SEQUENCE");
    let warn_on_missing_contig = args.bool("WARN_ON_MISSING_CONTIG", false);
    let _ = args.bool("LOG_FAILED_INTERVALS", true);
    let write_original_position = args.bool("WRITE_ORIGINAL_POSITION", false);
    let write_original_alleles = args.bool("WRITE_ORIGINAL_ALLELES", false);
    let min_match = args.double("LIFTOVER_MIN_MATCH", 1.0);
    let allow_missing = args.bool("ALLOW_MISSING_FIELDS_IN_HEADER", false);
    let recover_swapped = args.bool("RECOVER_SWAPPED_REF_ALT", false);
    let tags_to_reverse = args.collection("TAGS_TO_REVERSE", &["AF"]);
    let tags_to_drop = args.collection("TAGS_TO_DROP", &["MAX_AF"]);
    let disable_sort = args.bool("DISABLE_SORT", false);
    let create_index = args.bool("CREATE_INDEX", false);
    let max_records_in_ram = args
        .int("MAX_RECORDS_IN_RAM", DEFAULT_MAX_RECORDS_IN_RAM as i64)
        .max(1) as usize;
    // shortcut: the first TMP_DIR only, where htsjdk picks the one with the most free space.
    let tmp_dir = args
        .all("TMP_DIR")
        .first()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(std::env::temp_dir);

    if create_index && disable_sort {
        log_error(
            TOOL,
            "CREATE_INDEX=true and DISABLE_SORT=true are mutually exclusive.",
        );
        std::process::exit(1);
    }

    let lift_over =
        LiftOver::load(&std::fs::read_to_string(&chain).unwrap_or_else(|e| die(&format!("{e}"))))
            .unwrap_or_else(|e| die(&format!("htsjdk.samtools.SAMException: {e:?}")));
    // A file whose samples are not in sorted order has its genotypes decoded by the codec, and
    // the reader leaves those records without a block to copy.
    let mut vcf = VcfReader::open(&input).unwrap_or_else(|exception| die(&exception));
    INPUT.set(input.clone()).expect("set once");
    let file = vcf.file.clone();
    let file = &file;

    // The target reference and its dictionary, which must sit beside it.
    let dict_path = std::path::Path::new(&reference).with_extension("dict");
    let dictionary = match std::fs::read_to_string(&dict_path) {
        Ok(text) => parse_sam_dictionary(&text),
        Err(_) => {
            read_failure_first();
            log_error(
                TOOL,
                &format!(
                    "Reference {} must have an associated Dictionary .dict file in the same \
                     directory.",
                    absolute(&reference)
                ),
            );
            std::process::exit(1);
        }
    };
    let contigs = read_fasta_file(&reference).unwrap_or_else(|e| refuse(&format!("{e:?}")));
    // The first contig of each name the dictionary has, moved rather than copied.
    let mut ref_seqs: HashMap<String, Vec<u8>> = HashMap::new();
    for contig in contigs {
        if dictionary.iter().any(|s| s.name == contig.name) {
            ref_seqs.entry(contig.name).or_insert(contig.bases);
        }
    }

    // The output header.
    let mut out_header = VcfHeader {
        lines: file
            .header
            .lines
            .iter()
            .filter(|l| l.key() != "reference")
            .cloned()
            .collect(),
        samples: {
            let mut s = file.header.samples.clone();
            s.sort_by(|a, b| picard_analysis::vcf_io::java_compare(a, b));
            s
        },
    };
    set_sequence_dictionary(&mut out_header, &dictionary);
    if write_original_position {
        add_line(
            &mut out_header,
            HeaderLine::info(
                ORIGINAL_CONTIG,
                Cardinality::Fixed(1),
                LineType::String,
                "The name of the source contig/chromosome prior to liftover.",
            ),
        );
        add_line(
            &mut out_header,
            HeaderLine::info(
                ORIGINAL_START,
                Cardinality::Fixed(1),
                LineType::String,
                "The position of the variant on the source contig prior to liftover.",
            ),
        );
        add_line(
            &mut out_header,
            HeaderLine::info(
                ORIGINAL_ALLELES,
                Cardinality::R,
                LineType::String,
                "A list of the original alleles (including REF) of the variant prior to \
                 liftover.  If the alleles were not changed during liftover, this attribute will \
                 be omitted.",
            ),
        );
    }
    add_line(
        &mut out_header,
        HeaderLine::info(
            SWAPPED_ALLELES,
            Cardinality::Fixed(0),
            LineType::Flag,
            "The REF and the ALT alleles have been swapped in liftover due to changes in the \
             reference. It is possible that not all INFO annotations reflect this swap, and in \
             the genotypes, only the GT, PL, and AD fields have been modified. You should check \
             the TAGS_TO_REVERSE parameter that was used during the LiftOver to be sure.",
        ),
    );
    add_line(
        &mut out_header,
        HeaderLine::info(
            REV_COMPED_ALLELES,
            Cardinality::Fixed(0),
            LineType::Flag,
            "The REF and the ALT alleles have been reverse complemented in liftover since the \
             mapping from the previous reference to the current one was on the negative strand.",
        ),
    );
    add_other_meta_data_line(
        &mut out_header,
        "reference",
        // Mapped back to the path the reference was given: a .vcf.gz carries it inside its blocks,
        // where the harness cannot rewrite it.
        &format!(
            "file:{}",
            picard_analysis::fingerprint::reference_path(&absolute(&reference))
        ),
    );

    // The reject header.
    let mut reject_header = file.header.clone();
    for (id, description) in [
        (
            FILTER_CANNOT_LIFTOVER_REV_COMP,
            "Liftover of a variant that needed reverse-complementing failed for unknown reasons.",
        ),
        (
            FILTER_NO_TARGET,
            "Variant could not be lifted between genome builds.",
        ),
        (
            FILTER_MISMATCHING_REF_ALLELE,
            "Reference allele does not match reference genome sequence after liftover.",
        ),
        (
            FILTER_INDEL_STRADDLES_TWO_INTERVALS,
            "Reference allele in Indel is straddling multiple intervals in the chain, and so the \
             results are not well defined.",
        ),
    ] {
        add_line(&mut reject_header, HeaderLine::filter(id, description));
    }
    add_line(
        &mut reject_header,
        HeaderLine::info(
            ATTEMPTED_LOCUS,
            Cardinality::Fixed(1),
            LineType::String,
            "The locus of the variant in the TARGET prior to failing due to reference allele \
             mismatching to the target reference.",
        ),
    );
    add_line(
        &mut reject_header,
        HeaderLine::info(
            ATTEMPTED_ALLELES,
            Cardinality::Fixed(1),
            LineType::String,
            "The alleles of the variant in the TARGET prior to failing due to reference allele \
             mismatching to the target reference.",
        ),
    );

    let missing = if allow_missing {
        MissingFields::Allow
    } else {
        MissingFields::Refuse
    };
    let out_encoder = LineEncoder::new(&out_header);
    let reject_encoder = LineEncoder::new(&reject_header);
    let index_of = |name: &str| dictionary.iter().position(|s| s.name == name).unwrap_or(0) as u32;
    let mut lifted = Held::new(&output, !disable_sort, max_records_in_ram, &tmp_dir);
    let mut rejected = Held::new(&reject, false, max_records_in_ram, &tmp_dir);
    let mut line = String::new();
    // Held as encoded with missing fields allowed, so that the writer refuses them in its turn;
    // with DISABLE_SORT the output is written on the fly, so a record the header cannot carry
    // fails at the record rather than after the loop.
    let lifted_missing = if disable_sort {
        missing
    } else {
        MissingFields::Allow
    };
    let mut add_lifted = |record: Record, line: &mut String| {
        if let Err(e) = out_encoder.encode(&record, lifted_missing, line) {
            refuse(&e.exception());
        }
        let key = (index_of(&record.variant.contig), record.variant.start);
        lifted.add(key, record, line);
    };
    let mut add_rejected = |record: Record, line: &mut String| {
        if let Err(e) = reject_encoder.encode(&record, MissingFields::Allow, line) {
            refuse(&e.exception());
        }
        rejected.add((0, 0), record, line);
    };
    let reject_record = |record: &Record, reason: &str| {
        let mut r = record.clone();
        add_filter(&mut r.variant, reason);
        r
    };

    while let Some(record) = vcf.try_next() {
        let record = &record.unwrap_or_else(|exception| die(&exception));
        let ctx = &record.variant;
        let source = Interval {
            contig: ctx.contig.clone(),
            start: ctx.start as i32,
            end: ctx.stop as i32,
            negative_strand: false,
            name: Some(format!("{}:{}-{}", ctx.contig, ctx.start, ctx.stop)),
        };
        let Some(target) = lift_over.lift_over_with_min_match(&source, min_match) else {
            add_rejected(reject_record(record, FILTER_NO_TARGET), &mut line);
            continue;
        };
        if ctx.reference().len() as i64 != i64::from(target.end - target.start + 1) {
            add_rejected(
                reject_record(record, FILTER_INDEL_STRADDLES_TWO_INTERVALS),
                &mut line,
            );
            continue;
        }
        let Some(reference_bases) = ref_seqs.get(&target.contig) else {
            add_rejected(reject_record(record, FILTER_NO_TARGET), &mut line);
            if !warn_on_missing_contig {
                read_failure_first();
                log_error(
                    TOOL,
                    &format!(
                        "Encountered a contig, {} that is not part of the target reference.",
                        target.contig
                    ),
                );
                std::process::exit(1);
            }
            continue;
        };
        let Some((vc, touched)) = lift_variant(
            ctx,
            &target,
            reference_bases,
            write_original_position,
            write_original_alleles,
        ) else {
            add_rejected(
                reject_record(record, FILTER_CANNOT_LIFTOVER_REV_COMP),
                &mut line,
            );
            continue;
        };

        // tryToAddVariant.
        let ref_string = slice(reference_bases, vc.start, vc.stop);
        if !ref_string.eq_ignore_ascii_case(&vc.reference().base_string()) {
            if is_biallelic(&vc)
                && is_snp(&vc)
                && ref_string.eq_ignore_ascii_case(&vc.alleles[1].base_string())
                && recover_swapped
            {
                add_lifted(
                    Record {
                        variant: swap_ref_alt(&vc, &tags_to_reverse, &tags_to_drop),
                        lazy_genotypes: None,
                        fully_decoded: false,
                    },
                    &mut line,
                );
                continue;
            }
            let mut r = reject_record(record, FILTER_MISMATCHING_REF_ALLELE);
            set_attribute(
                &mut r.variant,
                ATTEMPTED_LOCUS,
                Value::Str(format!("{}:{}-{}", vc.contig, vc.start, vc.stop)),
            );
            let alleles = |a: &Allele| {
                format!(
                    "{}{}",
                    a.display_string(),
                    if a.is_reference() { "*" } else { "" }
                )
            };
            set_attribute(
                &mut r.variant,
                ATTEMPTED_ALLELES,
                Value::Str(format!(
                    "{}->{}",
                    alleles(vc.reference()),
                    vc.alternate_alleles()
                        .iter()
                        .map(alleles)
                        .collect::<Vec<_>>()
                        .join(",")
                )),
            );
            add_rejected(r, &mut line);
        } else {
            let lazy_genotypes = if touched {
                None
            } else {
                record.lazy_genotypes.clone()
            };
            add_lifted(
                Record {
                    variant: vc,
                    lazy_genotypes,
                    fully_decoded: false,
                },
                &mut line,
            );
        }
    }

    let version = file.codec_version;
    // `unsetOption(INDEX_ON_THE_FLY)`: the rejects are never indexed.
    write_all(
        &reject,
        &reject_header,
        rejected,
        false,
        missing,
        None,
        version,
    );
    // `modifyOption(INDEX_ON_THE_FLY, !DISABLE_SORT)` over the target reference's dictionary.
    let index_dictionary = (!disable_sort).then_some(dictionary.as_slice());
    write_all(
        &output,
        &out_header,
        lifted,
        !disable_sort,
        missing,
        index_dictionary,
        version,
    );
}
