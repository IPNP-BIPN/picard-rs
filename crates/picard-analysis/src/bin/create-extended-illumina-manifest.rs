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
    process_snp, render_row, Extension, Flag, Record, Statistics, Strand, VERSION,
};
use picard_analysis::infinium::Egt;
use picard_analysis::metrics_cli::{refuse_validation, thrown, Args};
use picard_analysis::vcf_io::read_path;

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

fn main() {
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

    let mut errors = Vec::new();
    if target_build != "37" {
        errors.push("Currently this tool only supports Build 37".to_string());
    }
    if flag_duplicates && cluster.is_none() {
        errors.push("In order to flag duplicates, a CLUSTER_FILE must be supplied".to_string());
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

    // `IlluminaManifest`: the heading rows up to "Loci Count", then the assay rows.
    let text = std::fs::read_to_string(&input).unwrap_or_else(|e| thrown(&format!("{e}")));
    let lines: Vec<&str> = text.lines().map(|l| l.trim_end_matches('\r')).collect();
    let mut heading: Vec<Vec<String>> = Vec::new();
    let mut count = 0usize;
    for line in &lines {
        let row: Vec<String> = line.split(',').map(str::to_string).collect();
        let tag = row[0].trim().to_string();
        heading.push(row.clone());
        if tag == "Loci Count" {
            count = row.get(1).and_then(|v| v.trim().parse().ok()).unwrap_or(0);
            break;
        }
    }
    let assay = lines
        .iter()
        .position(|l| l.trim() == "[Assay]")
        .unwrap_or_else(|| thrown("picard.PicardException: no [Assay] section"));
    let columns: Vec<String> = lines[assay + 1]
        .trim()
        .split(',')
        .map(str::to_string)
        .collect();
    let column = |name: &str| columns.iter().position(|c| c == name);
    let mut rows: Vec<(Vec<String>, Record, Extension)> = Vec::new();
    for line in lines.iter().skip(assay + 2).take(count) {
        let row: Vec<String> = line.split(',').map(str::to_string).collect();
        if row.len() != columns.len() {
            break;
        }
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
        let bases = sequences
            .iter()
            .find(|(n, _)| *n == record.chr)
            .map(|(_, b)| b)
            .unwrap_or_else(|| thrown("java.lang.NullPointerException"));
        let base = bases
            .get((record.map_info - 1).max(0) as usize)
            .map(|b| (b.to_ascii_uppercase() as char).to_string())
            .unwrap_or_else(|| {
                thrown(&format!(
                    "htsjdk.samtools.SAMException: Malformed query; start point {} lies after end point {}",
                    record.map_info,
                    bases.len()
                ))
            });
        let extension = process_snp(&record, &base);
        rows.push((row, record, extension));
    }

    // dbSNP: every position a record spans, mapped to its ID, for the loci the manifest has.
    let mut known: Vec<(String, String)> = Vec::new();
    if let Some(path) = &dbsnp {
        let vcf = read_path(path).unwrap_or_else(|e| thrown(&e));
        for v in vcf.records.iter().map(|r| &r.variant) {
            let wanted = rows.iter().any(|(_, _, e)| {
                !e.flag.is_fail() && e.b37_chr == v.contig && {
                    let start = i64::from(e.b37_pos);
                    let end = start + e.allele_a.len().max(e.allele_b.len()) as i64;
                    v.start <= end && v.stop >= start
                }
            });
            if !wanted {
                continue;
            }
            for pos in v.start..=v.stop {
                let key = format!("{}.{pos}", v.contig);
                match known.iter_mut().find(|(k, _)| *k == key) {
                    Some(slot) => slot.1 = v.id.clone(),
                    None => known.push((key, v.id.clone())),
                }
            }
        }
    }
    for (_, record, extension) in &mut rows {
        let key = format!("{}.{}", record.chr, record.map_info);
        extension.rs_id = known
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, id)| id.clone())
            .unwrap_or_else(|| "null".to_string());
    }

    // `flagDuplicates`.
    if flag_duplicates {
        let path = cluster.clone().unwrap_or_default();
        let egt = Egt::parse(&std::fs::read(&path).unwrap_or_default()).unwrap_or_else(|_| {
            thrown(&format!(
                "picard.PicardException: Error reading cluster file '{}'",
                absolute(&path)
            ))
        });
        let mut groups: Vec<(String, Vec<usize>)> = Vec::new();
        for (i, (_, _, e)) in rows.iter().enumerate() {
            if e.flag.is_fail() {
                continue;
            }
            let mut key = format!("{}:{}.{}", e.b37_chr, e.b37_pos, e.ref_allele);
            if e.allele_a != e.ref_allele {
                key.push_str(&format!(".{}", e.allele_a));
            }
            if e.allele_b != e.allele_a && e.allele_b != e.ref_allele {
                key.push_str(&format!(".{}", e.allele_b));
            }
            match groups.iter_mut().find(|(k, _)| *k == key) {
                Some((_, list)) => list.push(i),
                None => groups.push((key, vec![i])),
            }
        }
        let score = |i: usize| -> f32 {
            let name = &rows[i].1.name;
            egt.rs_names
                .iter()
                .position(|n| n == name)
                .map(|at| egt.total_score[at])
                .unwrap_or_else(|| thrown("java.lang.NullPointerException"))
        };
        let mut dupes: Vec<usize> = Vec::new();
        for (_, list) in groups.into_iter().filter(|(_, l)| l.len() > 1) {
            let mut best = list[0];
            for &i in &list[1..] {
                if score(i) > score(best) {
                    best = i;
                }
            }
            dupes.extend(list.into_iter().filter(|&i| i != best));
        }
        for i in dupes {
            rows[i].2.flag = Flag::Dupe;
        }
    }

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
    for (row, record, extension) in &rows {
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
        out.push_str(&render_row(row, extension));
        out.push('\n');
    }
    if let Err(e) = std::fs::write(&output, out) {
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
        "CreateExtendedIlluminaManifest (version: {VERSION}) Report For: {}\nGenerated on: <now>\nUsing Illumina Manifest: {}\n",
        file_name(&output),
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
