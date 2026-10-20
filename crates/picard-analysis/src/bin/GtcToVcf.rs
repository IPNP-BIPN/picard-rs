//! `GtcToVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.GtcToVcf.doWork` at tag 3.4.0. The files are read by
//! `picard_analysis::infinium`, and each record's text is `picard_analysis::gtc_to_vcf`'s; this
//! is the argument surface, the header and the order.
//!
//! The records are sorted by the reference's dictionary (`SortingCollection` with a
//! `VariantContextComparator`), not by contig name, and a site whose assay alleles are three is
//! filtered `TRIALLELIC` as it is written. The call rate in the header is the GTC's call count over
//! its SNPs less the assays the cluster file zeroed out, as a Java `double`.

use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use picard_analysis::extract_run::java_double_to_string;
use picard_analysis::fingerprint::reference_path;
use picard_analysis::gtc_to_vcf::{variant_line, Call, Cluster, Transformation};
use picard_analysis::infinium::{Bpm, Egt, ExtendedManifest, Gtc, ReadError, Reader};
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};
use picard_analysis::vcf_io::parse_sam_dictionary;
use std::io::{Read, Seek, Write};

const TOOL: &str = "GtcToVcf";

/// The twenty-three controls of `ArraysControlInfo.CONTROL_INFO`, in order.
const CONTROLS: [(&str, &str); 23] = [
    ("DNP(High)", "Staining"),
    ("DNP(Bgnd)", "Staining"),
    ("Biotin(High)", "Staining"),
    ("Biotin(Bgnd)", "Staining"),
    ("Extension(A)", "Extension"),
    ("Extension(T)", "Extension"),
    ("Extension(C)", "Extension"),
    ("Extension(G)", "Extension"),
    ("TargetRemoval", "TargetRemoval"),
    ("Hyb(High)", "Hybridization"),
    ("Hyb(Medium)", "Hybridization"),
    ("Hyb(Low)", "Hybridization"),
    ("String(PM)", "Stringency"),
    ("String(MM)", "Stringency"),
    ("NSB(Bgnd)Red", "Non-SpecificBinding"),
    ("NSB(Bgnd)Purple", "Non-SpecificBinding"),
    ("NSB(Bgnd)Blue", "Non-SpecificBinding"),
    ("NSB(Bgnd)Green", "Non-SpecificBinding"),
    ("NP(A)", "Non-Polymorphic"),
    ("NP(T)", "Non-Polymorphic"),
    ("NP(C)", "Non-Polymorphic"),
    ("NP(G)", "Non-Polymorphic"),
    ("Restore", "Restoration"),
];

fn raise(message: &str) -> ! {
    if message.starts_with("java.") {
        thrown(message)
    } else {
        thrown(&format!("picard.PicardException: {message}"))
    }
}

/// `FilenameUtils.removeExtension` on a name.
fn without_extension(name: &str) -> &str {
    match name.rfind('.') {
        Some(dot) if !name[dot..].contains('/') => &name[..dot],
        _ => name,
    }
}

fn file_name(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

fn unstructured(key: &str, value: &str) -> HeaderLine {
    HeaderLine::Unstructured {
        key: key.to_string(),
        value: value.to_string(),
    }
}

fn read_gtc(path: &str, bpm: &Bpm) -> Gtc {
    let mut reader = Reader::open(path)
        .unwrap_or_else(|e| raise(&format!("java.io.FileNotFoundException: {e}")));
    match Gtc::parse(&mut reader, bpm) {
        Ok(g) => g,
        Err(ReadError::Io(_)) => raise(&format!(
            "Error processing GTC File: {}",
            reference_path(path)
        )),
        Err(ReadError::Picard(m)) => raise(&m),
    }
}

fn main() {
    picard_analysis::argument_check::check("GtcToVcf");
    picard_analysis::output::install();
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("R", "REFERENCE_SEQUENCE")]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let manifest_path = args.required("EXTENDED_ILLUMINA_MANIFEST");
    let cluster_path = args.required("CLUSTER_FILE");
    let bpm_path = args.required("ILLUMINA_BEAD_POOL_MANIFEST_FILE");
    let reference = args.required("REFERENCE_SEQUENCE");
    let expected_gender = args.get("EXPECTED_GENDER").map(str::to_string);
    let sample_alias = args.get("SAMPLE_ALIAS").map(str::to_string);
    let pipeline_version = args.get("PIPELINE_VERSION").map(str::to_string);
    let analysis_version = args.get("ANALYSIS_VERSION_NUMBER").map(str::to_string);
    let gender_gtc = args.get("GENDER_GTC").map(str::to_string);
    let fingerprint_vcf = args
        .get("FINGERPRINT_GENOTYPES_VCF_FILE")
        .map(str::to_string);
    let strict_zeroed = args.bool("DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS", false);

    // `customCommandLineValidation`: the dictionary must say GRCh37.
    let dict_path = match reference.rsplit_once('.') {
        Some((stem, _)) => format!("{stem}.dict"),
        None => format!("{reference}.dict"),
    };
    let dictionary = std::fs::read_to_string(&dict_path)
        .map(|t| parse_sam_dictionary(&t))
        .unwrap_or_default();
    match dictionary.first().and_then(|s| s.assembly.clone()) {
        None => refuse_validation(
            TOOL,
            &["Assembly tag ('AS') is required in the sequence dictionary).".to_string()],
        ),
        Some(a) if a != "GRCh37" => refuse_validation(TOOL, &[format!("The selected reference sequence ('{a}') is not supported.  This tool is currently only implemented to support NCBI Build 37 / HG19 Reference Sequence.")]),
        _ => {}
    }

    let fingerprint_sex = match &fingerprint_vcf {
        Some(path) => {
            let vcf = picard_analysis::vcf_io::read_path(path).unwrap_or_else(|e| raise(&e));
            let gender = vcf.file.header.lines.iter().find_map(|l| match l {
                HeaderLine::Unstructured { key, value } if key == "gender" => Some(value.clone()),
                _ => None,
            });
            match gender.as_deref() {
                None => "Unknown".to_string(),
                Some(g @ ("Male" | "Female" | "Unknown" | "NotReported")) => g.to_string(),
                Some(other) => raise(&format!(
                    "java.lang.IllegalArgumentException: No enum constant picard.pedigree.Sex.{other}"
                )),
            }
        }
        None => "Unknown".to_string(),
    };
    let mut bpm_reader = Reader::open(&bpm_path)
        .unwrap_or_else(|e| raise(&format!("java.io.FileNotFoundException: {e}")));
    let bpm = match Bpm::parse(&mut bpm_reader) {
        Ok(b) => b,
        Err(ReadError::Io(_)) => raise(&format!(
            "Error reading bpm file '{}'",
            reference_path(&bpm_path)
        )),
        Err(ReadError::Picard(m)) => raise(&m),
    };
    let gtc_gender = gender_gtc.as_ref().map(|p| read_gtc(p, &bpm).gender);
    let gtc = read_gtc(&input, &bpm);
    drop(bpm_reader);
    let egt = match Reader::open(&cluster_path)
        .map_err(drop)
        .and_then(|mut r| Egt::parse(&mut r).map_err(drop))
    {
        Ok(e) => e,
        Err(()) => raise(&format!(
            "Error processing GTC File: {}",
            reference_path(&input)
        )),
    };
    let mut manifest = ExtendedManifest::open(&manifest_path).unwrap_or_else(|e| raise(&e));
    let descriptor_name = manifest.descriptor_file_name.clone();
    let manifest_version = manifest.extended_manifest_version.clone();
    // The rows are read as they are used, where the reference read them all first: a refusal
    // gives way to a row the reader would have refused.
    let gtc_manifest = without_extension(gtc.snp_manifest.as_deref().unwrap_or(""));
    let descriptor = without_extension(&descriptor_name);
    if !gtc_manifest.eq_ignore_ascii_case(descriptor) {
        let message = format!(
            "The GTC's manifest name {gtc_manifest} does not match the Illumina manifest name {descriptor}"
        );
        raise(&manifest.first_refusal().unwrap_or(message));
    }

    // `fillContexts`, each line spilled as it is made: the `SortingCollection` the reference sorts
    // them through, as compact keys over a file that is gone once the run ends.
    // `CommandLineProgram` makes each TMP_DIR and goes on whatever happens; a sort that stays in
    // memory never writes there, so a TMP_DIR that cannot take the spill falls back to the
    // system's.
    let spill_in = |dir: std::path::PathBuf| {
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join(format!("GtcToVcf-{}.lines", std::process::id()));
        std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .map(|file| (file, path))
    };
    let (spill, spill_path) = args
        .all("TMP_DIR")
        .first()
        .map(std::path::PathBuf::from)
        .and_then(|dir| spill_in(dir).ok())
        .map_or_else(|| spill_in(std::env::temp_dir()), Ok)
        .unwrap_or_else(|e| raise(&format!("htsjdk.samtools.SAMException: {e}")));
    let _ = std::fs::remove_file(&spill_path);
    let mut spilled = std::io::BufWriter::with_capacity(1 << 16, &spill);
    let mut spilled_at = 0u64;
    let mut zeroed = 0;
    // (contig, position, where the line starts in the spill, its length)
    let mut rows: Vec<(usize, i32, u64, u32)> = Vec::new();
    let mut next_index = 0usize;
    while let Some(record) = manifest.next() {
        let record = record.unwrap_or_else(|e| raise(&e));
        let gtc_index = next_index;
        next_index += 1;
        let mut refuse = |message: &str| -> ! {
            raise(
                &manifest
                    .first_refusal()
                    .unwrap_or_else(|| message.to_string()),
            )
        };
        let Some(egt_index) = egt.index_of(&record.name) else {
            refuse(&format!(
                "Found no record in cluster file for manifest entry '{}'",
                record.name
            ));
        };
        let total = egt.total_score[egt_index];
        if total == 0.0 {
            zeroed += 1;
        }
        if record.flag.is_fail() {
            continue;
        }
        let at = |v: &Option<Vec<i32>>| v.as_ref().and_then(|v| v.get(gtc_index).copied());
        let atf = |v: &Option<Vec<f32>>| v.as_ref().and_then(|v| v.get(gtc_index).copied());
        let oob_message = format!(
            "java.lang.ArrayIndexOutOfBoundsException: Index {gtc_index} out of bounds for length {}",
            gtc.number_of_snps
        );
        let Some(genotype) = gtc
            .genotypes
            .as_ref()
            .and_then(|g| g.get(gtc_index).copied())
        else {
            refuse(&oob_message)
        };
        if !(0..=3).contains(&genotype) {
            refuse(&format!(
                "Unexpected genotype call [{genotype}] for SNP: {}",
                record.name
            ));
        }
        if record.allele_a == record.allele_b {
            refuse(&format!(
                "Found same allele ({}) for A and B ",
                if record.allele_a.is_empty() {
                    "."
                } else {
                    &record.allele_a
                }
            ));
        }
        let (Some(score), Some(raw_x), Some(raw_y), Some(b_allele_freq), Some(log_r_ratio)) = (
            atf(&gtc.genotype_scores),
            at(&gtc.raw_x),
            at(&gtc.raw_y),
            atf(&gtc.b_allele_freqs),
            atf(&gtc.log_r_ratios),
        ) else {
            refuse(&oob_message)
        };
        let call = Call {
            genotype: genotype as u8,
            score,
            raw_x,
            raw_y,
            b_allele_freq,
            log_r_ratio,
        };
        if total == 0.0 && genotype != 0 && strict_zeroed {
            refuse("Found a call on a zeroed out Assay!!");
        }
        let cluster = Cluster {
            total_score: total,
            n: egt.n[egt_index],
            dev_r: egt.dev_r[egt_index],
            mean_r: egt.mean_r[egt_index],
            dev_theta: egt.dev_theta[egt_index],
            mean_theta: egt.mean_theta[egt_index],
        };
        // The normalization the GTC reader applied, given back as the transformation the line
        // re-applies: the locus's own, or none (the identity) when the manifest gave it no id.
        let transformation = bpm
            .all_normalization_ids
            .get(gtc_index)
            .and_then(|id| bpm.unique_normalization_ids.iter().position(|u| u == id))
            .and_then(|i| gtc.transformations.get(i).copied())
            .unwrap_or(Transformation {
                offset_x: 0.0,
                offset_y: 0.0,
                scale_x: 1.0,
                scale_y: 1.0,
                shear: 0.0,
                theta: 0.0,
            });
        let line = variant_line(&record, &call, &cluster, &transformation);
        let contig = dictionary
            .iter()
            .position(|s| s.name == record.b37_chr)
            .unwrap_or(usize::MAX);
        if let Err(e) = spilled.write_all(line.as_bytes()) {
            raise(&format!("htsjdk.samtools.SAMException: {e}"));
        }
        rows.push((contig, record.b37_pos, spilled_at, line.len() as u32));
        spilled_at += line.len() as u64;
    }
    if let Err(e) = spilled.flush() {
        raise(&format!("htsjdk.samtools.SAMException: {e}"));
    }
    drop(spilled);
    rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let call_rate = f64::from(gtc.num_calls) / f64::from(gtc.number_of_snps - zeroed);
    let mut header = VcfHeader::new();
    let mut add = |key: &str, value: &str| header.lines.push(unstructured(key, value));
    let input_name = file_name(&input);
    let chip_well_barcode =
        input_name[..input_name.rfind('.').unwrap_or(input_name.len())].to_string();
    add("fileDate", &picard_analysis::command_line::java_date_now());
    add("source", "GtcToVcf");
    let descriptor_name = &descriptor_name;
    add(
        "arrayType",
        &descriptor_name[..descriptor_name.rfind('.').unwrap_or(descriptor_name.len())],
    );
    add("extendedManifestFile", file_name(&manifest_path));
    add("extendedIlluminaManifestVersion", &manifest_version);
    add("chipWellBarcode", &chip_well_barcode);
    if let Some(v) = &analysis_version {
        add("analysisVersionNumber", v);
    }
    add("sampleAlias", sample_alias.as_deref().unwrap_or("null"));
    if let Some(g) = &expected_gender {
        add("expectedGender", g);
    }
    if let Some(p) = &pipeline_version {
        add("pipelineVersion", p);
    }
    let control_x = gtc.raw_control_x.clone().unwrap_or_default();
    let control_y = gtc.raw_control_y.clone().unwrap_or_default();
    let measurements = control_x.len() / CONTROLS.len();
    for (i, (control, category)) in CONTROLS.iter().enumerate() {
        let offset = i * measurements;
        let (Some(red), Some(green)) = (control_x.get(offset), control_y.get(offset)) else {
            raise(&format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {offset} out of bounds for length {}",
                control_x.len()
            ));
        };
        add(control, &format!("{control}|{category}|{red}|{green}"));
    }
    add("fingerprintGender", &fingerprint_sex);
    let autocall_gender = match &gtc_gender {
        Some(g) => g.clone(),
        None => gtc.gender.clone(),
    };
    add(
        "autocallGender",
        autocall_gender.as_deref().unwrap_or("null"),
    );
    add(
        "autocallDate",
        gtc.auto_call_date.as_deref().unwrap_or("null"),
    );
    add("imagingDate", gtc.imaging_date.as_deref().unwrap_or("null"));
    add("clusterFile", file_name(&cluster_path));
    add("manifestFile", descriptor_name);
    add("content", file_name(&manifest_path));
    add(
        "autocallVersion",
        gtc.auto_call_version.as_deref().unwrap_or("null"),
    );
    add("reference", &reference_path(&reference));
    // `CommandLineProgram.getVersion()`, which is the jar's implementation version with its prefix.
    add("picardVersion", "Version:3.5.0");
    let percentile = |p: &Option<[i32; 3]>| -> String {
        p.map(|v| v[2].to_string())
            .unwrap_or_else(|| raise("java.lang.NullPointerException"))
    };
    add("p95Red", &percentile(&gtc.red_percentiles));
    add("p95Green", &percentile(&gtc.green_percentiles));
    add("scannerName", gtc.scanner_name.as_deref().unwrap_or("null"));
    add("gtcCallRate", &java_double_to_string(call_rate));
    add("gtcCallRateDetail", "The gtcCallRate is the Call Rate reported in the Illumina GTC file, corrected for the presence of Zeroed-Out SNPs");
    let fmt = |id: &str, t: LineType, d: &str| HeaderLine::format(id, Cardinality::Fixed(1), t, d);
    let info = |id: &str, t: LineType, d: &str| HeaderLine::info(id, Cardinality::Fixed(1), t, d);
    let mut compound = vec![
        fmt("GT", LineType::String, "Genotype"),
        HeaderLine::info(
            "AC",
            Cardinality::A,
            LineType::Integer,
            "Allele count in genotypes, for each ALT allele, in the same order as listed",
        ),
        HeaderLine::info(
            "AF",
            Cardinality::A,
            LineType::Float,
            "Allele Frequency, for each ALT allele, in the same order as listed",
        ),
        info(
            "AN",
            LineType::Integer,
            "Total number of alleles in called genotypes",
        ),
        fmt("IGC", LineType::Float, "Illumina GenCall Confidence Score"),
        fmt("X", LineType::Integer, "Raw X intensity"),
        fmt("Y", LineType::Integer, "Raw Y intensity"),
        fmt("NORMX", LineType::Float, "Normalized X intensity"),
        fmt("NORMY", LineType::Float, "Normalized Y intensity"),
        fmt("R", LineType::Float, "Normalized R value"),
        fmt("THETA", LineType::Float, "Normalized Theta value"),
        fmt("BAF", LineType::Float, "B Allele Frequency"),
        fmt("LRR", LineType::Float, "Log R Ratio"),
        info("ALLELE_A", LineType::String, "A allele"),
        info("ALLELE_B", LineType::String, "B allele"),
        info("ILLUMINA_STRAND", LineType::String, "Probe strand"),
        info("PROBE_A", LineType::String, "Probe base pair sequence"),
        info(
            "PROBE_B",
            LineType::String,
            "Probe base pair sequence; not missing for strand-ambiguous SNPs",
        ),
        info(
            "BEADSET_ID",
            LineType::Integer,
            "Bead set ID for normalization",
        ),
        info(
            "ILLUMINA_CHR",
            LineType::String,
            "Chromosome in Illumina manifest",
        ),
        info(
            "ILLUMINA_POS",
            LineType::Integer,
            "Position in Illumina manifest",
        ),
        info(
            "ILLUMINA_BUILD",
            LineType::String,
            "Genome Build in Illumina manifest",
        ),
        info("SOURCE", LineType::String, "Probe source"),
        info("GC_SCORE", LineType::Float, "Gentrain Score"),
    ];
    for gt in ["AA", "AB", "BB"] {
        compound.push(info(
            &format!("N_{gt}"),
            LineType::Integer,
            &format!("Number of {gt} calls in training set"),
        ));
        compound.push(info(
            &format!("devR_{gt}"),
            LineType::Float,
            &format!("Standard deviation of normalized R for {gt} cluster"),
        ));
        compound.push(info(
            &format!("devTHETA_{gt}"),
            LineType::Float,
            &format!("Standard deviation of normalized THETA for {gt} cluster"),
        ));
        compound.push(info(
            &format!("devX_{gt}"),
            LineType::Float,
            &format!("Standard deviation of normalized X for {gt} cluster"),
        ));
        compound.push(info(
            &format!("devY_{gt}"),
            LineType::Float,
            &format!("Standard deviation of normalized Y for {gt} cluster"),
        ));
        compound.push(info(
            &format!("meanR_{gt}"),
            LineType::Float,
            &format!("Mean of normalized R for {gt} cluster"),
        ));
        compound.push(info(
            &format!("meanTHETA_{gt}"),
            LineType::Float,
            &format!("Mean of normalized THETA for {gt} cluster"),
        ));
        compound.push(info(
            &format!("meanX_{gt}"),
            LineType::Float,
            &format!("Mean of normalized X for {gt} cluster"),
        ));
        compound.push(info(
            &format!("meanY_{gt}"),
            LineType::Float,
            &format!("Mean of normalized Y for {gt} cluster"),
        ));
    }
    compound.push(info("refSNP", LineType::String, "dbSNP rsID"));
    compound.push(HeaderLine::filter("DUPE", "Duplicate assays position."));
    compound.push(HeaderLine::filter("TRIALLELIC", "Tri-allelic assay."));
    compound.push(HeaderLine::filter(
        "FAIL_REF",
        "Assay failed to map to reference.",
    ));
    compound.push(HeaderLine::filter(
        "ZEROED_OUT_ASSAY",
        "Assay Zeroed out (marked as uncallable) in the Illumina Cluster File",
    ));
    header.lines.extend(compound);
    header.samples = vec![chip_well_barcode.clone()];
    picard_analysis::vcf_io::set_sequence_dictionary(&mut header, &dictionary);

    if let Err(e) = write_vcf(&output, &header.write(), &spill, &rows, &dictionary) {
        raise(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}

/// `vcf_io::write_text_indexed` over lines read back from the spill in `rows`' order: the file,
/// then its `.idx` unless the name is block-compressed or BCF.
fn write_vcf(
    path: &str,
    header: &str,
    spill: &std::fs::File,
    rows: &[(usize, i32, u64, u32)],
    dictionary: &[picard_analysis::vcf_io::Sequence],
) -> std::io::Result<()> {
    use htsjdk_tribble::index_write::{BalanceApproach, DynamicIndexCreator, Feature};
    let lower = path.to_ascii_lowercase();
    let indexed = ![".gz", ".bgz", ".bcf"]
        .iter()
        .any(|ext| lower.ends_with(ext));
    let mut creator = DynamicIndexCreator::new(BalanceApproach::ForSeekTime);
    let mut out = std::io::BufWriter::with_capacity(1 << 16, std::fs::File::create(path)?);
    out.write_all(header.as_bytes())?;
    let mut at = header.len() as u64;
    let mut source = std::io::BufReader::with_capacity(1 << 16, spill);
    source.rewind()?;
    let mut source_at = 0u64;
    let mut line = Vec::new();
    for &(_, _, start, len) in rows {
        // In the order the manifest gave, every seek is to where the reader already is.
        source.seek_relative(start as i64 - source_at as i64)?;
        line.resize(len as usize, 0);
        source.read_exact(&mut line)?;
        source_at = start + u64::from(len);
        if indexed {
            let text = std::str::from_utf8(&line).unwrap_or("");
            let fields: Vec<&str> = text.splitn(9, '\t').collect();
            let start: i64 = fields.get(1).and_then(|p| p.parse().ok()).unwrap_or(0);
            let reference = fields.get(3).map_or(1, |r| r.len() as i64);
            let end = fields
                .get(7)
                .and_then(|info| {
                    info.split(';')
                        .find_map(|entry| entry.strip_prefix("END="))
                        .and_then(|end| end.parse().ok())
                })
                .unwrap_or(start + reference - 1);
            creator.add_feature(
                &Feature {
                    contig: fields[0].to_string(),
                    start: start as i32,
                    end: end as i32,
                },
                at as i64,
            );
        }
        out.write_all(&line)?;
        out.write_all(b"\n")?;
        at += u64::from(len) + 1;
    }
    out.into_inner().map_err(|e| e.into_error())?;
    if !indexed {
        return Ok(());
    }
    let output = std::path::Path::new(path);
    let timestamp = std::fs::metadata(output)?
        .modified()?
        .duration_since(std::time::UNIX_EPOCH)
        .map(|elapsed| elapsed.as_millis() as i64)
        .unwrap_or(0);
    let indexed = picard_analysis::vcf_io::Indexed {
        dictionary,
        path: output,
        timestamp,
    };
    if let Some(index) = finish_tribble(creator, at as i64, &indexed) {
        std::fs::write(format!("{path}.idx"), index)?;
    }
    Ok(())
}

/// `vcf_io::finish_tribble`, which that module keeps to itself: the `.idx` of the features
/// `creator` was given, the file being `file_length` bytes long.
fn finish_tribble(
    creator: htsjdk_tribble::index_write::DynamicIndexCreator,
    file_length: i64,
    indexed: &picard_analysis::vcf_io::Indexed,
) -> Option<Vec<u8>> {
    use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR, VERSION};
    use htsjdk_tribble::index_write::BuiltIndex;

    let mut properties: Vec<(String, String)> = indexed
        .dictionary
        .iter()
        .map(|sequence| {
            (
                format!("DICT:{}", sequence.name),
                sequence.length.to_string(),
            )
        })
        .collect();
    properties.extend(creator.properties());
    let (index_type, contigs, interval_contigs) = match creator.finalize(file_length) {
        Ok(BuiltIndex::Linear(contigs)) => (LINEAR, by_name(contigs, |c| &c.name), Vec::new()),
        Ok(BuiltIndex::IntervalTree(intervals)) => {
            (INTERVAL_TREE, Vec::new(), by_name(intervals, |c| &c.name))
        }
        Err(_) => return None,
    };
    let absolute = std::path::absolute(indexed.path).unwrap_or_else(|_| indexed.path.into());
    let index = TribbleIndex {
        index_type,
        version: VERSION,
        indexed_path: format!("file://{}", absolute.display()),
        indexed_file_size: file_length,
        indexed_file_timestamp: indexed.timestamp,
        indexed_file_md5: String::new(),
        flags: 0,
        properties,
        contigs,
        interval_contigs,
    };
    index.write().ok()
}

/// `LinkedHashMap.put` over a list: a repeated key keeps its first position and its last value.
fn by_name<T>(items: Vec<T>, name: impl Fn(&T) -> &String) -> Vec<T> {
    let mut kept: Vec<T> = Vec::with_capacity(items.len());
    for item in items {
        match kept.iter().position(|k| name(k) == name(&item)) {
            Some(at) => kept[at] = item,
            None => kept.push(item),
        }
    }
    kept
}
