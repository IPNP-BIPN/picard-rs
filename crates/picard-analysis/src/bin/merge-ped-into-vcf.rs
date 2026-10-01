//! `MergePedIntoVcf` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.arrays.MergePedIntoVcf.doWork` and `picard.arrays.ZCallPedFile` at tag 3.4.0.
//! The VCF round trip is `picard_analysis::vcf_io`, and AC, AF and AN are htsjdk-rs's
//! `calculateChromosomeCounts`.
//!
//! # A failure is an exit code of zero
//!
//! `doWork` wraps everything after the readability checks in `catch (Exception e) {
//! e.printStackTrace(); }` and returns 0. So every refusal this tool has is a success to its
//! caller, and what tells them apart is the file:
//!
//! * a failure before the writer is built -- the thresholds (one `NA` of a pair, a short line),
//!   the PED (two lines, a field that is not one character, too few fields for the MAP), a VCF of
//!   more than one sample -- leaves no file at all;
//! * the writer is built in a try-with-resources, so a failure INSIDE the record loop closes it:
//!   the file holds the header and every record before the one that threw. A record throws when
//!   the MAP names no call for its ID, when a PED letter is not `A`, `B` or `0`, when an array
//!   allele is not a valid allele, when its genotype has fewer than two alleles or no extended
//!   attribute (`GenotypeBuilder` hands back an immutable empty map, and `put` throws), when a
//!   translated allele is not one of the context's (an allele built non-reference on the
//!   reference's bases is a different allele), and when the encoder meets a key the header does
//!   not declare.
//!
//! # The merged genotype
//!
//! `GenotypeBuilder.create(name, alleles, attributes)` takes the EXTENDED attributes only, so GQ,
//! DP, AD, PL and FT of the original genotype are gone; GTA (the original call, by allele index)
//! and GTZ (zCall's) join them. The thresholds are INFO strings, written as they were read, and
//! `NA NA` becomes `.`.

use std::collections::HashMap;

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::chromosome_counts::{
    calculate_chromosome_counts, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY,
};
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Genotype, Value, VariantContext};
use picard_analysis::metrics_cli::{absolute, thrown, Args};
use picard_analysis::vcf_io::{
    add_other_meta_data_line, header_dictionary, read_path, write_output, Record,
};

/// Java's `String.split(regex)` for a one-character pattern: trailing empty strings removed.
fn java_split(text: &str, separator: impl Fn(char) -> bool) -> Vec<&str> {
    let mut parts: Vec<&str> = text.split(separator).collect();
    while parts.len() > 1 && parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    if parts.len() == 1 && parts[0].is_empty() && !text.is_empty() {
        parts.clear();
    }
    parts
}

/// `\s` in a Java regular expression.
fn is_java_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')
}

/// `BufferedReader.readLine` over a whole file: `\n`, `\r` and `\r\n` end a line.
fn lines(text: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let bytes = text.as_bytes();
    let (mut start, mut i) = (0, 0);
    while i < bytes.len() {
        if bytes[i] == b'\n' || bytes[i] == b'\r' {
            out.push(&text[start..i]);
            if bytes[i] == b'\r' && i + 1 < bytes.len() && bytes[i + 1] == b'\n' {
                i += 1;
            }
            start = i + 1;
        }
        i += 1;
    }
    if start < bytes.len() {
        out.push(&text[start..]);
    }
    out
}

fn read(path: &str) -> Result<String, String> {
    std::fs::read(path)
        .map(|bytes| String::from_utf8_lossy(&bytes).into_owned())
        .map_err(|e| format!("java.io.IOException: {e}"))
}

/// `parseZCallThresholds`.
fn parse_thresholds(path: &str) -> Result<HashMap<String, [String; 2]>, String> {
    let text = read(path)?;
    let mut thresholds = HashMap::new();
    for line in lines(&text) {
        let tokens = java_split(line, |c| c == '\t');
        if tokens.len() < 3 {
            return Err(format!(
                "java.lang.ArrayIndexOutOfBoundsException: Index {} out of bounds for length {}",
                tokens.len().max(1),
                tokens.len()
            ));
        }
        let (x, y) = (tokens[1], tokens[2]);
        if x == "NA" || y == "NA" {
            if x != "NA" || y != "NA" {
                return Err(
                    "picard.PicardException: Thresholds should either both exist or \
                            both not exist."
                        .to_string(),
                );
            }
            thresholds.insert(tokens[0].to_string(), [".".to_string(), ".".to_string()]);
        } else {
            thresholds.insert(tokens[0].to_string(), [x.to_string(), y.to_string()]);
        }
    }
    Ok(thresholds)
}

/// `ZCallPedFile.fromFile`: the PED's allele pairs keyed by the MAP's names, paired by index.
fn parse_ped(ped_path: &str, map_path: &str) -> Result<HashMap<String, String>, String> {
    let ped = read(ped_path)?;
    if lines(&ped).len() > 1 {
        return Err(
            "picard.PicardException: Only single-sample .ped files are supported.".to_string(),
        );
    }
    let fields = java_split(&ped, is_java_space);
    let map = read(map_path)?;
    let out_of_bounds = |index: usize, length: usize| {
        format!(
            "java.lang.ArrayIndexOutOfBoundsException: Index {index} out of bounds for length \
             {length}"
        )
    };
    let mut alleles = HashMap::new();
    for (i, line) in lines(&map).into_iter().enumerate() {
        let index = i * 2 + 6;
        let first = fields
            .get(index)
            .ok_or_else(|| out_of_bounds(index, fields.len()))?;
        let second = fields
            .get(index + 1)
            .ok_or_else(|| out_of_bounds(index + 1, fields.len()))?;
        if first.chars().count() != 1 || second.chars().count() != 1 {
            return Err(
                "picard.PicardException: Malformed file: each allele should be a \
                        single character."
                    .to_string(),
            );
        }
        let name_fields = java_split(line, is_java_space);
        let name = name_fields
            .get(1)
            .ok_or_else(|| out_of_bounds(1, name_fields.len()))?;
        alleles.insert(name.to_string(), format!("{first}{second}"));
    }
    Ok(alleles)
}

/// `formatAllele`: `*` is the spanning deletion, a starred allele is the reference, and anything
/// else is built as non-reference.
fn format_allele(text: &str) -> Result<Allele, String> {
    let made = if text == "*" {
        Allele::from_str("*", false)
    } else if text.contains('*') {
        Allele::from_str(text.replace('*', " ").trim(), true)
    } else {
        Allele::from_str(text, false)
    };
    made.map_err(|e| format!("java.lang.IllegalArgumentException: {e:?}"))
}

/// `String.valueOf(attributes.get(key))`.
fn attribute_string(vc: &VariantContext, key: &str) -> String {
    match vc.attributes.iter().find(|(k, _)| k == key).map(|(_, v)| v) {
        None => "null".to_string(),
        Some(Value::Str(s)) => s.clone(),
        Some(Value::List(items)) => {
            let parts: Vec<String> = items
                .iter()
                .map(|v| v.format().unwrap_or_default())
                .collect();
            format!("[{}]", parts.join(", "))
        }
        Some(other) => other.format().unwrap_or_default(),
    }
}

/// `translateAllele`.
fn translate_allele(allele_a: &str, allele_b: &str, allele: char) -> Result<Allele, String> {
    match allele {
        'A' => format_allele(allele_a),
        'B' => format_allele(allele_b),
        '0' => Ok(Allele::no_call()),
        other => Err(format!("picard.PicardException: Illegal allele: {other}")),
    }
}

/// `Map.put` over an ordered list of pairs.
fn put(pairs: &mut Vec<(String, Value)>, key: &str, value: Value) {
    match pairs.iter_mut().find(|(k, _)| k == key) {
        Some(entry) => entry.1 = value,
        None => pairs.push((key.to_string(), value)),
    }
}

/// One iteration of `writeVcf`'s loop: the merged record, or why it threw.
fn merge_record(
    vc: &VariantContext,
    thresholds: &HashMap<String, [String; 2]>,
    calls: &HashMap<String, String>,
) -> Result<VariantContext, String> {
    let mut merged = vc.clone();
    if let Some([x, y]) = thresholds.get(&vc.id) {
        put(&mut merged.attributes, "zthresh_X", Value::Str(x.clone()));
        put(&mut merged.attributes, "zthresh_Y", Value::Str(y.clone()));
    }
    let original = vc
        .genotypes
        .first()
        .ok_or("java.lang.IndexOutOfBoundsException: Index 0 out of bounds for length 0")?;

    // `VCFEncoder.buildAlleleStrings`: an index per allele of the context, `.` for a no-call.
    let allele_string = |allele: &Allele| -> String {
        if allele.is_no_call() {
            return ".".to_string();
        }
        vc.alleles
            .iter()
            .position(|a| a == allele)
            .map_or("null".to_string(), |i| i.to_string())
    };

    let pair = calls.get(&vc.id).ok_or_else(|| {
        format!(
            "picard.PicardException: No zCall alleles found for snp {}",
            vc.id
        )
    })?;
    let allele_a = attribute_string(vc, "ALLELE_A");
    let allele_b = attribute_string(vc, "ALLELE_B");
    let mut letters = pair.chars();
    let first = translate_allele(&allele_a, &allele_b, letters.next().unwrap_or(' '))?;
    let second = translate_allele(&allele_a, &allele_b, letters.next().unwrap_or(' '))?;
    let zcall = vec![first, second];

    if original.alleles.len() < 2 {
        return Err(format!(
            "java.lang.IndexOutOfBoundsException: Index {} out of bounds for length {}",
            original.alleles.len(),
            original.alleles.len()
        ));
    }
    // `getExtendedAttributes()` is the codec's map, or an immutable empty one when the genotype
    // had none, and `put` on that throws.
    if original.extended.is_empty() {
        return Err("java.lang.UnsupportedOperationException".to_string());
    }
    let mut extended = original.extended.clone();
    let gta = format!(
        "{}/{}",
        allele_string(&original.alleles[0]),
        allele_string(&original.alleles[1])
    );
    let gtz = format!("{}/{}", allele_string(&zcall[0]), allele_string(&zcall[1]));
    put(&mut extended, "GTA", Value::Str(gta));
    put(&mut extended, "GTZ", Value::Str(gtz));

    let mut genotype = Genotype::new(&original.sample_name, zcall);
    genotype.extended = extended;

    // `VariantContext.validateGenotypes`, run by the `make()` inside calculateChromosomeCounts.
    for allele in &genotype.alleles {
        if !allele.is_no_call() && !vc.alleles.contains(allele) {
            return Err(format!(
                "java.lang.IllegalStateException: Allele in genotype {} not in the variant \
                 context",
                allele.display_string()
            ));
        }
    }
    merged.genotypes = vec![genotype];

    // `calculateChromosomeCounts(builder, false)`.
    let counts = calculate_chromosome_counts(&merged, false, &[]);
    for (key, value) in counts.attributes() {
        put(&mut merged.attributes, &key, value);
    }
    if counts.allele_number.is_some() && merged.alternate_alleles().is_empty() {
        merged
            .attributes
            .retain(|(k, _)| k != ALLELE_COUNT_KEY && k != ALLELE_FREQUENCY_KEY);
    }
    Ok(merged)
}

/// `VCFHeader.addMetaDataLine` for an INFO or FORMAT line: kept only when no line of that kind
/// already has its ID.
fn add_compound(header: &mut VcfHeader, line: HeaderLine) {
    let HeaderLine::Compound { key, id, .. } = &line else {
        return;
    };
    let taken = header.lines.iter().any(|existing| {
        matches!(existing, HeaderLine::Compound { key: k, id: i, .. } if k == key && i == id)
    });
    if !taken {
        header.lines.push(line);
    }
}

/// Everything inside `doWork`'s try block.
fn run(args: &Args, output: &str) -> Result<(), String> {
    let vcf_path = args.required("ORIGINAL_VCF");
    let thresholds_path = args.required("ZCALL_THRESHOLDS_FILE");
    let version = args.required("ZCALL_VERSION");

    let thresholds = parse_thresholds(&thresholds_path)?;
    let calls = parse_ped(&args.required("PED_FILE"), &args.required("MAP_FILE"))?;
    let vcf = read_path(&vcf_path)?;
    if vcf.file.header.samples.len() > 1 {
        return Err(
            "picard.PicardException: MergePedIntoVCF only works with single-sample VCFs."
                .to_string(),
        );
    }

    // `addAdditionalHeaderFields`.
    let mut header = vcf.file.header.clone();
    add_other_meta_data_line(&mut header, "zcallVersion", &version);
    let name = std::path::Path::new(&thresholds_path)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    add_other_meta_data_line(&mut header, "zcallThresholds", &name);
    for (id, description) in [
        ("zthresh_X", "zCall X threshold"),
        ("zthresh_Y", "zCall Y threshold"),
    ] {
        add_compound(
            &mut header,
            HeaderLine::info(id, Cardinality::Fixed(1), LineType::Float, description),
        );
    }
    for (id, description) in [
        ("GTA", "Illumina Autocall Genotype"),
        ("GTZ", "zCall Genotype"),
    ] {
        add_compound(
            &mut header,
            HeaderLine::format(id, Cardinality::Fixed(1), LineType::String, description),
        );
    }

    // `VariantContextWriterBuilder.build()` opens the file, and only then refuses to index on
    // the fly without a dictionary.
    let Some(dictionary) = header_dictionary(&vcf.file.header) else {
        std::fs::write(output, b"").map_err(|e| format!("{e}"))?;
        return Err(
            "java.lang.IllegalArgumentException: A reference dictionary is required for \
                    creating Tribble indices on the fly"
                .to_string(),
        );
    };

    let encoder = VcfEncoder::new(&header);
    let mut written: Vec<Record> = Vec::new();
    let mut failure = None;
    for record in &vcf.records {
        let merged = match merge_record(&record.variant, &thresholds, &calls) {
            Ok(merged) => merged,
            Err(e) => {
                failure = Some(e);
                break;
            }
        };
        if let Err(e) = encoder.encode(&merged) {
            failure = Some(format!("java.lang.IllegalStateException: {e:?}"));
            break;
        }
        written.push(Record {
            variant: merged,
            lazy_genotypes: None,
        });
    }
    write_output(output, &header, &written, Some(&dictionary)).map_err(|e| format!("{e}"))?;
    match failure {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

fn main() {
    let args = Args::from_env(&[
        ("VCF", "ORIGINAL_VCF"),
        ("PED", "PED_FILE"),
        ("MAP", "MAP_FILE"),
        ("ZCALL_T_FILE", "ZCALL_THRESHOLDS_FILE"),
        ("O", "OUTPUT"),
    ]);
    let output = args.required("OUTPUT");
    // `IOUtil.assertFileIsReadable` on the three side files, before the try block.
    for name in ["PED_FILE", "MAP_FILE", "ZCALL_THRESHOLDS_FILE"] {
        let path = args.required(name);
        if !std::path::Path::new(&path).is_file() {
            thrown(&format!(
                "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
                absolute(&path)
            ));
        }
    }
    if let Err(e) = run(&args, &output) {
        // `e.printStackTrace()`, and then `return 0`.
        eprintln!("{e}");
    }
}
