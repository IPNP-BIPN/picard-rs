//! `CreateExtendedIlluminaManifest` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.illumina.CreateExtendedIlluminaManifest.doWork` at tag 3.4.0 for a
//! manifest on build 37, the only target the tool supports. Each locus is
//! `picard_analysis::create_extended_illumina_manifest::process_snp`, and the counts are its
//! `Statistics`; this reads the manifest, looks the loci up in dbSNP, flags duplicates and writes
//! the three files.
//!
//! A duplicate is an assay sharing another's build-37 position and alleles; of each group, the one
//! the cluster file scores highest stays `PASS` (the first, when scores tie) and the rest are
//! `DUPE`. The heading rows are written padded to the width of the first one, which is why
//! `[Heading]` comes out as `[Heading],`.

use picard_analysis::create_extended_illumina_manifest::{
    process_snp, render_row, Flag, Record, Statistics, Strand, VERSION,
};
use picard_analysis::infinium::{check_utf8, text_lines, Egt, Reader};
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};
use picard_analysis::vcf_io::VcfReader;
use std::io::Write;

const TOOL: &str = "CreateExtendedIlluminaManifest";

const EXTENDED: [&str; 7] = [
    "build37Chr",
    "build37Pos",
    "build37RefAllele",
    "build37AlleleA",
    "build37AlleleB",
    "build37Rsid",
    "build37Flag",
];

fn absolute(path: &str) -> String {
    std::path::absolute(path)
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| path.to_string())
}

fn file_name(path: &str) -> &str {
    path.rsplit('/').next().unwrap_or(path)
}

/// `addHeaderLine`: the fields, padded with empty ones to the heading's width.
fn header_line(width: usize, fields: &[&str]) -> String {
    let mut row: Vec<&str> = fields.to_vec();
    while row.len() < width {
        row.push("");
    }
    row.join(",")
}

/// The lines of a file [`check_utf8`] passed, without their closing `\r`s.
fn lines(path: &str) -> impl Iterator<Item = String> {
    text_lines(path)
        .unwrap_or_else(|e| thrown(&format!("{e}")))
        .map(|l| l.trim_end_matches('\r').to_string())
}

/// The assay rows: up to `count` of the lines after the `[Assay]` header, until one is not as
/// wide as the header, each split and as a record.
fn assays<'a>(
    path: &str,
    columns: &'a [String],
    count: usize,
) -> impl Iterator<Item = (Vec<String>, Record)> + 'a {
    let mut lines = lines(path);
    lines.by_ref().find(|l| l.trim() == "[Assay]");
    lines.next();
    lines
        .take(count)
        .map(|line| line.split(',').map(str::to_string).collect::<Vec<String>>())
        .take_while(|row| row.len() == columns.len())
        .map(|row| {
            let column = |name: &str| columns.iter().position(|c| c == name);
            let get = |name: &str| column(name).map(|i| row[i].clone()).unwrap_or_default();
            let record = Record {
                ilmn_id: get("IlmnID"),
                name: get("Name"),
                ilmn_strand: get("IlmnStrand"),
                snp: get("SNP").to_uppercase(),
                address_a: get("AddressA_ID"),
                allele_a_probe_seq: get("AlleleA_ProbeSeq"),
                address_b: get("AddressB_ID"),
                allele_b_probe_seq: get("AlleleB_ProbeSeq"),
                genome_build: get("GenomeBuild"),
                chr: get("Chr"),
                map_info: get("MapInfo").trim().parse().unwrap_or(0),
                ref_strand: Strand::parse(&get("RefStrand")),
            };
            (row, record)
        })
}

fn main() {
    picard_analysis::argument_check::check("CreateExtendedIlluminaManifest");
    picard_analysis::output::install();
    let args = Args::from_env(&[("I", "INPUT"), ("O", "OUTPUT"), ("R", "REFERENCE_SEQUENCE")]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let report_path = args.required("REPORT_FILE");
    let bad_path = args.get("BAD_ASSAYS_FILE").map(str::to_string);
    let cluster = args.get("CLUSTER_FILE").map(str::to_string);
    let dbsnp = args.get("DBSNP_FILE").map(str::to_string);
    let flag_duplicates = args.bool("FLAG_DUPLICATES", true);
    let target_build = args.get("TARGET_BUILD").unwrap_or("37").to_string();
    let reference = args.required("REFERENCE_SEQUENCE");
    let supported_builds = args.all("SUPPORTED_BUILD");
    let supported_references = args.all("SUPPORTED_REFERENCE_FILE");
    let supported_chains = args.all("SUPPORTED_CHAIN_FILE");

    let mut errors = Vec::new();
    if target_build != "37" {
        errors.push("Currently this tool only supports Build 37".to_string());
    }
    if flag_duplicates && cluster.is_none() {
        errors.push("In order to flag duplicates, a CLUSTER_FILE must be supplied".to_string());
    }
    // The three supported-build lists go together, one entry each per build.
    if !supported_builds.is_empty()
        || !supported_references.is_empty()
        || !supported_chains.is_empty()
    {
        if supported_builds.is_empty()
            || supported_references.is_empty()
            || supported_chains.is_empty()
        {
            errors.push("Parameters for 'SUPPORTED_BUILD', 'SUPPORTED_REFERENCE_FILE', and 'SUPPORTED_CHAIN_FILE' must ALL be specified or not at all.".to_string());
        } else {
            if supported_builds.len() != supported_references.len() {
                errors.push("The number of inputs for 'SUPPORTED_BUILD' does not match the number of inputs for 'SUPPORTED_REFERENCE_FILE'".to_string());
            }
            if supported_builds.len() != supported_chains.len() {
                errors.push("The number of inputs for 'SUPPORTED_BUILD' does not match the number of inputs for 'SUPPORTED_CHAIN_FILE'".to_string());
            }
        }
    }
    if !errors.is_empty() {
        refuse_validation(TOOL, &errors);
    }

    let fasta = std::fs::read(&reference).unwrap_or_else(|e| thrown(&format!("{e}")));
    let sequences: Vec<(String, Vec<u8>)> = htsjdk_bam::fasta::read_fasta(&fasta[..])
        .unwrap_or_else(|e| thrown(&format!("{e:?}")))
        .into_iter()
        .map(|s| (s.name, s.bases))
        .collect();

    // `IlluminaManifest`: the heading rows up to "Loci Count", then the assay rows, read twice: to
    // find the duplicates and the dbSNP loci, then to write each row.
    check_utf8(&input).unwrap_or_else(|e| thrown(&format!("{e}")));
    let mut heading: Vec<Vec<String>> = Vec::new();
    let mut count = 0usize;
    for line in lines(&input) {
        let row: Vec<String> = line.split(',').map(str::to_string).collect();
        let tag = row[0].trim().to_string();
        heading.push(row.clone());
        if tag == "Loci Count" {
            count = row.get(1).and_then(|v| v.trim().parse().ok()).unwrap_or(0);
            break;
        }
    }
    let mut assay_lines = lines(&input);
    assay_lines
        .by_ref()
        .find(|l| l.trim() == "[Assay]")
        .unwrap_or_else(|| thrown("picard.PicardException: no [Assay] section"));
    let columns: Vec<String> = assay_lines
        .next()
        .expect("a line after [Assay]")
        .trim()
        .split(',')
        .map(str::to_string)
        .collect();
    drop(assay_lines);
    let reference_base = |record: &Record| -> String {
        let bases = sequences
            .iter()
            .find(|(n, _)| *n == record.chr)
            .map(|(_, b)| b)
            .unwrap_or_else(|| thrown("java.lang.NullPointerException"));
        bases
            .get((record.map_info - 1).max(0) as usize)
            .map(|b| (b.to_ascii_uppercase() as char).to_string())
            .unwrap_or_else(|| {
                thrown(&format!(
                    "htsjdk.samtools.SAMException: Malformed query; start point {} lies after end point {}",
                    record.map_info,
                    bases.len()
                ))
            })
    };

    // The first pass keeps, of each passing locus, its build-37 span (by contig number) and its
    // duplicate key.
    let mut contigs: Vec<String> = Vec::new();
    let mut spans: Vec<(usize, i64, i64)> = Vec::new();
    let mut keys: Vec<(String, usize, String)> = Vec::new();
    let mut rows = 0usize;
    for (i, (_, record)) in assays(&input, &columns, count).enumerate() {
        rows += 1;
        let e = process_snp(&record, &reference_base(&record));
        if e.flag.is_fail() {
            continue;
        }
        let contig = match contigs.iter().position(|c| *c == e.b37_chr) {
            Some(at) => at,
            None => {
                contigs.push(e.b37_chr.clone());
                contigs.len() - 1
            }
        };
        let start = i64::from(e.b37_pos);
        spans.push((
            contig,
            start,
            start + e.allele_a.len().max(e.allele_b.len()) as i64,
        ));
        if flag_duplicates {
            let mut key = format!("{}:{}.{}", e.b37_chr, e.b37_pos, e.ref_allele);
            if e.allele_a != e.ref_allele {
                key.push_str(&format!(".{}", e.allele_a));
            }
            if e.allele_b != e.allele_a && e.allele_b != e.ref_allele {
                key.push_str(&format!(".{}", e.allele_b));
            }
            keys.push((key, i, record.name));
        }
    }

    // dbSNP: every position a record spans, mapped to its ID, for the loci the manifest has.
    let mut known: foldhash::HashMap<String, String> = foldhash::HashMap::default();
    if let Some(path) = &dbsnp {
        let mut vcf = VcfReader::open(path).unwrap_or_else(|e| thrown(&e));
        // The loci a record may overlap, found by an overlap index: a scan of every locus for
        // every record was quadratic on a chip.
        let loci = picard_analysis::overlap::Spans::new(spans.drain(..));
        while let Some(record) = vcf.try_next() {
            let v = record.unwrap_or_else(|e| thrown(&e)).variant;
            let Some(contig) = contigs.iter().position(|c| *c == v.contig) else {
                continue;
            };
            if !loci.overlaps(&contig, v.start, v.stop) {
                continue;
            }
            for pos in v.start..=v.stop {
                known.insert(format!("{}.{pos}", v.contig), v.id.clone());
            }
        }
    }
    drop(spans);

    // `flagDuplicates`: of each group of loci sharing a key, all but the best scored.
    let mut dupe = vec![false; rows];
    if flag_duplicates {
        let path = cluster.clone().unwrap_or_default();
        let egt = Reader::open(&path)
            .map_err(drop)
            .and_then(|mut r| Egt::parse(&mut r).map_err(drop))
            .unwrap_or_else(|()| {
                thrown(&format!(
                    "picard.PicardException: Error reading cluster file '{}'",
                    absolute(&path)
                ))
            });
        let score = |name: &str| -> f32 {
            egt.index_of(name)
                .map(|at| egt.total_score[at])
                .unwrap_or_else(|| thrown("java.lang.NullPointerException"))
        };
        // Each group's loci together, in the order the manifest has them.
        keys.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        for group in keys.chunk_by(|a, b| a.0 == b.0).filter(|g| g.len() > 1) {
            let mut best = 0;
            for at in 1..group.len() {
                if score(&group[at].2) > score(&group[best].2) {
                    best = at;
                }
            }
            for (at, (_, i, _)) in group.iter().enumerate() {
                if at != best {
                    dupe[*i] = true;
                }
            }
        }
    }
    drop(keys);

    // The manifest.
    let width = heading.first().map_or(1, Vec::len);
    let mut out = String::new();
    for row in &heading[..heading.len().saturating_sub(1)] {
        let fields: Vec<&str> = row.iter().map(String::as_str).collect();
        out.push_str(&header_line(width, &fields));
        out.push('\n');
    }
    let reference_abs = absolute(&reference);
    for (key, value) in [
        (
            "CreateExtendedIlluminaManifest.version",
            Some(VERSION.to_string()),
        ),
        ("Target Build", Some(target_build.clone())),
        ("Target Reference File", Some(reference_abs.clone())),
        ("Cluster File", cluster.as_deref().map(absolute)),
        ("dbSNP File", dbsnp.as_deref().map(absolute)),
    ] {
        if let Some(v) = value {
            out.push_str(&header_line(width, &[key, &v]));
            out.push('\n');
        }
    }
    // A row each for the supported builds, their references and their chain files.
    if !supported_builds.is_empty() {
        for (label, values) in [
            ("Supported Build", supported_builds.clone()),
            (
                "Supported Reference File",
                supported_references
                    .iter()
                    .map(|path| absolute(path))
                    .collect(),
            ),
            (
                "Supported Chain File",
                supported_chains.iter().map(|path| absolute(path)).collect(),
            ),
        ] {
            let mut fields: Vec<&str> = vec![label];
            fields.extend(values.iter().map(String::as_str));
            out.push_str(&header_line(width, &fields));
            out.push('\n');
        }
    }
    if let Some(last) = heading.last() {
        let fields: Vec<&str> = last.iter().map(String::as_str).collect();
        out.push_str(&header_line(width, &fields));
        out.push('\n');
    }
    out.push_str(&header_line(width, &["[Assay]"]));
    out.push('\n');
    let mut names: Vec<&str> = columns.iter().map(String::as_str).collect();
    names.extend(EXTENDED);
    out.push_str(&names.join(","));
    out.push('\n');
    let mut statistics = Statistics::default();
    let mut bad: Vec<String> = Vec::new();
    let mut file = std::fs::File::create(&output)
        .map(|f| std::io::BufWriter::with_capacity(1 << 16, f))
        .unwrap_or_else(|e| thrown(&format!("picard.PicardException: {e}")));
    let mut put = |text: &str| {
        file.write_all(text.as_bytes())
            .unwrap_or_else(|e| thrown(&format!("picard.PicardException: {e}")))
    };
    put(&out);
    for (i, (row, record)) in assays(&input, &columns, count).enumerate() {
        let mut extension = process_snp(&record, &reference_base(&record));
        extension.rs_id = known
            .get(&format!("{}.{}", record.chr, record.map_info))
            .cloned()
            .unwrap_or_else(|| "null".to_string());
        if dupe[i] {
            extension.flag = Flag::Dupe;
        }
        let extension = &extension;
        let record = &record;
        statistics.update(record, extension, &target_build);
        if extension.flag.is_fail() {
            bad.push(
                [
                    record.ilmn_id.as_str(),
                    &record.name,
                    &record.genome_build,
                    &record.chr,
                    &record.map_info.to_string(),
                    extension.flag.name(),
                ]
                .join(","),
            );
        }
        let mut line = render_row(&row, extension);
        line.push('\n');
        put(&line);
    }
    if let Err(e) = file.into_inner().map_err(|e| e.into_error()) {
        thrown(&format!("picard.PicardException: {e}"));
    }
    if let Some(path) = &bad_path {
        let mut text = format!(
            "## The following assays were marked by CreateExtendedIlluminaManifest as Unparseable (input file: {})\n#IlmnId,Name,GenomeBuild,Chr,MapInfo,FailureFlag\n",
            absolute(&input)
        );
        for line in &bad {
            text.push_str(line);
            text.push('\n');
        }
        let _ = std::fs::write(path, text);
    }

    // The report.
    let s = &statistics;
    let mut r = format!(
        "CreateExtendedIlluminaManifest (version: {VERSION}) Report For: {}\nGenerated on: {}\nUsing Illumina Manifest: {}\n",
        file_name(&output),
        picard_analysis::command_line::java_date_now(),
        absolute(&input)
    );
    if flag_duplicates {
        r.push_str("Duplicates were flagged\n");
    }
    if let Some(c) = &cluster {
        r.push_str(&format!("Using Illumina EGT: {}\n", absolute(c)));
    }
    r.push('\n');
    let lines = [
        format!("Total Number of Assays: {}", s.assays),
        format!(
            "Number of Assays on Build {target_build}: {}",
            s.on_target_build
        ),
        format!(
            "Number of Assays on unsupported genome build: {}",
            s.on_unsupported_genome_build
        ),
        format!("Number of Assays failing liftover: {}", s.liftover_failed),
        String::new(),
        format!(
            "Number of Assays on Build {target_build} or successfully lifted over: {}",
            s.on_target_build - s.on_unsupported_genome_build - s.liftover_failed
        ),
        format!("Number of Passing Assays: {}", s.assays - s.assays_flagged),
        format!("Number of Duplicated Assays: {}", s.assays_duplicated),
        format!("Number of Failing Assays: {}", s.assays_flagged),
        String::new(),
        format!("Number of SNPs: {}", s.snps),
        format!("Number of Passing SNPs: {}", s.snps - s.snps_flagged),
        format!("Number of Duplicated SNPs: {}", s.snps_duplicated),
        String::new(),
        format!("Number of Failing SNPs: {}", s.snps_flagged),
        format!(
            "Number of SNPs failed by Illumina: {}",
            s.snps_illumina_flagged
        ),
        format!(
            "Number of SNPs failed for refStrand mismatch: {}",
            s.ref_strand_mismatch
        ),
        format!(
            "Number of SNPs failed for missing AlleleB ProbeSeq: {}",
            s.snp_missing_allele_b_probe_sequence
        ),
        format!(
            "Number of SNPs failed for alleleA probe sequence mismatch: {}",
            s.snp_probe_sequence_mismatch
        ),
        format!(
            "Number of ambiguous SNPs on Positive Strand: {}",
            s.ambiguous_snps_on_positive_strand
        ),
        format!(
            "Number of ambiguous SNPs on Negative Strand: {}",
            s.ambiguous_snps_on_negative_strand
        ),
        String::new(),
        format!("Number of Indels: {}", s.indels),
        format!("Number of Passing Indels: {}", s.indels - s.indels_flagged),
        format!("Number of Duplicated Indels: {}", s.indels_duplicated),
        String::new(),
        format!("Number of Failing Indels: {}", s.indels_flagged),
        format!(
            "Number of Indels failed by Illumina: {}",
            s.indels_illumina_flagged
        ),
        format!(
            "Number of Indels failed for probe sequence mismatch: {}",
            s.indel_probe_sequence_mismatch
        ),
        format!(
            "Number of Indels failed for source sequence invalid: {}",
            s.indel_source_sequence_invalid
        ),
        format!("Number of Indels not found: {}", s.indels_not_found),
        format!("Number of Indels failed for conflict: {}", s.indel_conflict),
    ];
    for line in lines {
        r.push_str(&line);
        r.push('\n');
    }
    if let Err(e) = std::fs::write(&report_path, r) {
        thrown(&format!("picard.PicardException: {e}"));
    }
}
