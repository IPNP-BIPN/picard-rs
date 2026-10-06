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
use picard_analysis::infinium::{Bpm, Egt, ExtendedManifest, Gtc, ReadError};
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};
use picard_analysis::vcf_io::parse_sam_dictionary;

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
    let bytes = std::fs::read(path)
        .unwrap_or_else(|e| raise(&format!("java.io.FileNotFoundException: {e}")));
    match Gtc::parse(&bytes, bpm) {
        Ok(g) => g,
        Err(ReadError::Io(_)) => raise(&format!(
            "Error processing GTC File: {}",
            reference_path(path)
        )),
        Err(ReadError::Picard(m)) => raise(&m),
    }
}

fn main() {
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
    let bpm_bytes = std::fs::read(&bpm_path)
        .unwrap_or_else(|e| raise(&format!("java.io.FileNotFoundException: {e}")));
    let bpm = match Bpm::parse(&bpm_bytes) {
        Ok(b) => b,
        Err(ReadError::Io(_)) => raise(&format!(
            "Error reading bpm file '{}'",
            reference_path(&bpm_path)
        )),
        Err(ReadError::Picard(m)) => raise(&m),
    };
    let gtc_gender = gender_gtc.as_ref().map(|p| read_gtc(p, &bpm).gender);
    let gtc = read_gtc(&input, &bpm);
    let egt = match Egt::parse(&std::fs::read(&cluster_path).unwrap_or_default()) {
        Ok(e) => e,
        Err(_) => raise(&format!(
            "Error processing GTC File: {}",
            reference_path(&input)
        )),
    };
    let manifest_text = std::fs::read_to_string(&manifest_path).unwrap_or_default();
    let manifest =
        ExtendedManifest::parse(&manifest_text, &manifest_path).unwrap_or_else(|e| raise(&e));
    let gtc_manifest = without_extension(gtc.snp_manifest.as_deref().unwrap_or(""));
    let descriptor = without_extension(&manifest.descriptor_file_name);
    if !gtc_manifest.eq_ignore_ascii_case(descriptor) {
        raise(&format!(
            "The GTC's manifest name {gtc_manifest} does not match the Illumina manifest name {descriptor}"
        ));
    }

    // `fillContexts`.
    let mut zeroed = 0;
    let mut rows: Vec<(usize, i32, String)> = Vec::new();
    for (gtc_index, record) in manifest.records.iter().enumerate() {
        let Some(egt_index) = egt.rs_names.iter().position(|n| *n == record.name) else {
            raise(&format!(
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
        let oob = || -> ! {
            raise(&format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {gtc_index} out of bounds for length {}",
                gtc.number_of_snps
            ))
        };
        let genotype = gtc
            .genotypes
            .as_ref()
            .and_then(|g| g.get(gtc_index).copied())
            .unwrap_or_else(|| oob());
        if !(0..=3).contains(&genotype) {
            raise(&format!(
                "Unexpected genotype call [{genotype}] for SNP: {}",
                record.name
            ));
        }
        if record.allele_a == record.allele_b {
            raise(&format!(
                "Found same allele ({}) for A and B ",
                if record.allele_a.is_empty() {
                    "."
                } else {
                    &record.allele_a
                }
            ));
        }
        let call = Call {
            genotype: genotype as u8,
            score: atf(&gtc.genotype_scores).unwrap_or_else(|| oob()),
            raw_x: at(&gtc.raw_x).unwrap_or_else(|| oob()),
            raw_y: at(&gtc.raw_y).unwrap_or_else(|| oob()),
            b_allele_freq: atf(&gtc.b_allele_freqs).unwrap_or_else(|| oob()),
            log_r_ratio: atf(&gtc.log_r_ratios).unwrap_or_else(|| oob()),
        };
        if total == 0.0 && genotype != 0 && strict_zeroed {
            raise("Found a call on a zeroed out Assay!!");
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
        let line = variant_line(record, &call, &cluster, &transformation);
        let contig = dictionary
            .iter()
            .position(|s| s.name == record.b37_chr)
            .unwrap_or(usize::MAX);
        rows.push((contig, record.b37_pos, line));
    }
    rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let call_rate = f64::from(gtc.num_calls) / f64::from(gtc.number_of_snps - zeroed);
    let mut header = VcfHeader::new();
    let mut add = |key: &str, value: &str| header.lines.push(unstructured(key, value));
    let input_name = file_name(&input);
    let chip_well_barcode =
        input_name[..input_name.rfind('.').unwrap_or(input_name.len())].to_string();
    add("fileDate", "<now>");
    add("source", "GtcToVcf");
    let descriptor_name = &manifest.descriptor_file_name;
    add(
        "arrayType",
        &descriptor_name[..descriptor_name.rfind('.').unwrap_or(descriptor_name.len())],
    );
    add("extendedManifestFile", file_name(&manifest_path));
    add(
        "extendedIlluminaManifestVersion",
        &manifest.extended_manifest_version,
    );
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
    add("picardVersion", "Version:3.4.0");
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

    let mut text = header.write();
    for (_, _, line) in &rows {
        text.push_str(line);
        text.push('\n');
    }
    if let Err(e) = std::fs::write(&output, text) {
        raise(&format!("htsjdk.samtools.SAMException: {e}"));
    }
}
