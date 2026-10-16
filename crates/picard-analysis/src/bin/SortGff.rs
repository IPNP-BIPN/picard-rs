//! `SortGff` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.annotation.SortGff.doWork` at tag 3.4.0, with the parts of htsjdk 4.2.0 it
//! drives: `Gff3Codec` at `DecodeDepth.SHALLOW` (its `canDecode`, its directives, its comments and
//! sequence regions, `parseLine` and `validateFeature`), `TribbleIndexedFeatureReader`'s
//! whole-file iterator (which stamps a `TribbleException` with its input source), and
//! `Gff3Writer`. The comparator and the stable sort are `picard_analysis::sort_gff`.
//!
//! What the order of the checks decides:
//!
//! * `canDecode` comes before anything is read. It wants a GFF3 extension, a version directive
//!   on the first line, and a feature line after the comments that follow it; a file with
//!   directives and no feature is refused exactly like a file that is not GFF3.
//! * `SortingCollection` refuses a non-positive `nRecordsInMemory` once the dictionary is read.
//!   A positive one only decides where records wait: a spill writes them with `Gff3Writer` and
//!   reads them back with a fresh codec, which round-trips every field, and the merge breaks ties
//!   by file, so the order is the one stable sort.
//! * The writer is opened LAST, after the whole input has been read, and it is the writer that
//!   refuses an output whose name is not `.gff3`, `.gff`, `.gff3.gz` or `.gff.gz`.
//!
//! The output is the writer's own version directive (3.1.25, whatever the input said), every
//! comment of the input in order wherever it stood, every sequence region, and then the features,
//! with a flush directive before each one that no earlier feature's family can still reach.

use std::collections::HashMap;
use std::io::Write;

use picard_analysis::metrics_cli::{absolute, fail, thrown, Args};
use picard_analysis::sort_gff::{sort, Feature};
use picard_analysis::vcf_io::extract_dictionary;

const GFF3_EXTENSIONS: [&str; 4] = [".gff3", ".gff", ".gff3.gz", ".gff.gz"];

/// One decoded line: `Gff3BaseData`.
struct Gff3Feature {
    contig: String,
    source: String,
    kind: String,
    start: i32,
    end: i32,
    score: f64,
    strand: char,
    phase: i32,
    /// A `LinkedHashMap`: a repeated key keeps its first position and its last value.
    attributes: Vec<(String, Vec<String>)>,
    id: Option<String>,
}

impl Gff3Feature {
    fn attribute(&self, key: &str) -> Option<&Vec<String>> {
        self.attributes
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v)
    }

    /// `getAttribute`, which answers an empty list for a key the feature does not carry.
    fn attribute_or_empty(&self, key: &str) -> Vec<String> {
        self.attribute(key).cloned().unwrap_or_default()
    }
}

/// `SequenceRegion`.
struct SequenceRegion {
    contig: String,
    start: i32,
    end: i32,
    circular: bool,
}

/// A `TribbleException` thrown while decoding: the iterator appends its source to the message.
enum DecodeError {
    Tribble(String),
    /// Anything else, reported as the throwable it is.
    Other(String),
}

/// `Gff3Codec` at `DecodeDepth.SHALLOW`, with the state `SortGff` reads back from it.
#[derive(Default)]
struct Codec {
    current_line: i32,
    reached_fasta: bool,
    comments: Vec<String>,
    regions: Vec<SequenceRegion>,
}

impl Codec {
    /// `decode` for one line: a feature, or `None` for a directive, a comment or anything after
    /// the FASTA section begins.
    fn decode(&mut self, line: &str) -> Result<Option<Gff3Feature>, DecodeError> {
        self.current_line += 1;
        if self.reached_fasta {
            return Ok(None);
        }
        if line.starts_with('>') {
            self.reached_fasta = true;
            return Ok(None);
        }
        if line.starts_with('#') && !line.starts_with("##") {
            self.comments.push(line[1..].to_string());
            return Ok(None);
        }
        if line.starts_with("##") {
            self.parse_directive(line)?;
            return Ok(None);
        }
        let feature = parse_line(line, self.current_line)?;
        self.validate(&feature)?;
        Ok(Some(feature))
    }

    fn parse_directive(&mut self, line: &str) -> Result<(), DecodeError> {
        match to_directive(line) {
            Some(Directive::Version) | Some(Directive::Flush) => {}
            Some(Directive::Fasta) => self.reached_fasta = true,
            Some(Directive::SequenceRegion) => {
                let split = java_split_whitespace(line);
                let contig = url_decode(split.get(1).copied().unwrap_or_default())?;
                let number = |index: usize| -> Result<i32, DecodeError> {
                    let text = split.get(index).copied().unwrap_or_default();
                    java_parse_int(text).ok_or_else(|| {
                        DecodeError::Other(format!(
                            "java.lang.NumberFormatException: For input string: \"{text}\""
                        ))
                    })
                };
                let start = number(2)?;
                let end = number(3)?;
                if self.regions.iter().any(|r| r.contig == contig) {
                    return Err(DecodeError::Tribble(format!(
                        "directive for sequence-region {contig} included more than once."
                    )));
                }
                self.regions.push(SequenceRegion {
                    contig,
                    start,
                    end,
                    circular: false,
                });
            }
            // `logger.warn("ignoring directive " + line)`, on stderr and nowhere else.
            None => {}
        }
        Ok(())
    }

    /// `validateFeature`: a feature on a contig with a sequence region must lie inside it, or,
    /// once the region's landmark feature has said it is circular, overlap it.
    fn validate(&mut self, feature: &Gff3Feature) -> Result<(), DecodeError> {
        let Some(region) = self.regions.iter_mut().find(|r| r.contig == feature.contig) else {
            return Ok(());
        };
        if feature.start == region.start && feature.end == region.end {
            let value = extract_single_attribute(feature.attribute("Is_circular"))?;
            region.circular = value.is_some_and(|v| v.eq_ignore_ascii_case("true"));
        }
        let inside = if region.circular {
            feature.start <= region.end && region.start <= feature.end
        } else {
            region.start <= feature.start && feature.end <= region.end
        };
        if !inside {
            return Err(DecodeError::Tribble(format!(
                "feature at {}:{}-{} not contained in specified sequence region ({}:{}-{}",
                feature.contig, feature.start, feature.end, region.contig, region.start, region.end
            )));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Directive {
    Version,
    SequenceRegion,
    Flush,
    Fasta,
}

/// `\s` in a Java regular expression.
fn is_java_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\x0B' | '\x0C' | '\r')
}

/// `Gff3Directive.toDirective`: the first directive whose pattern matches the WHOLE line.
fn to_directive(line: &str) -> Option<Directive> {
    // `##gff-version\s+3(?:\.\d*)*$`
    if let Some(rest) = line.strip_prefix("##gff-version") {
        let trimmed = rest.trim_start_matches(is_java_space);
        if trimmed.len() < rest.len() {
            if let Some(tail) = trimmed.strip_prefix('3') {
                if tail.is_empty()
                    || (tail.starts_with('.')
                        && tail.chars().all(|c| c == '.' || c.is_ascii_digit()))
                {
                    return Some(Directive::Version);
                }
            }
        }
    }
    // `##sequence-region\s+.+ \d+ \d+$`
    if let Some(rest) = line.strip_prefix("##sequence-region") {
        let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
        if let Some((before, last)) = rest.rsplit_once(' ') {
            if let Some((head, middle)) = before.rsplit_once(' ') {
                let mut chars = head.chars();
                if digits(last)
                    && digits(middle)
                    && chars.next().is_some_and(is_java_space)
                    && chars.next().is_some()
                {
                    return Some(Directive::SequenceRegion);
                }
            }
        }
    }
    if line == "###" {
        return Some(Directive::Flush);
    }
    if line == "##FASTA" {
        return Some(Directive::Fasta);
    }
    None
}

/// `String.split("\\s+")`: a leading separator makes a leading empty string, trailing empties go.
fn java_split_whitespace(line: &str) -> Vec<&str> {
    let mut out: Vec<&str> = line.split(is_java_space).collect();
    // Runs of separators: `\s+` consumes them all, so the empties between them never appear.
    let leading_empty = out.first().is_some_and(|s| s.is_empty()) && !line.is_empty();
    out = out
        .into_iter()
        .enumerate()
        .filter(|(i, s)| !s.is_empty() || (*i == 0 && leading_empty))
        .map(|(_, s)| s)
        .collect();
    out
}

/// `ParsingUtils.split(input, delim)`: every field, empty ones and a trailing one included.
fn parsing_utils_split(input: &str, delimiter: char) -> Vec<&str> {
    input.split(delimiter).collect()
}

/// `Integer.parseInt`.
fn java_parse_int(text: &str) -> Option<i32> {
    text.parse::<i32>().ok()
}

/// `Double.parseDouble`: leading and trailing whitespace is ignored, a `d` or `f` suffix is
/// allowed, and the special values are spelled `NaN` and `Infinity`.
fn java_parse_double(text: &str) -> Option<f64> {
    let trimmed = text.trim_matches(|c: char| c <= ' ');
    let body = trimmed.trim_start_matches(['+', '-']);
    match body {
        "NaN" => return Some(f64::NAN),
        "Infinity" => {
            return Some(if trimmed.starts_with('-') {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            })
        }
        _ => {}
    }
    let numeric = trimmed
        .strip_suffix(['d', 'D', 'f', 'F'])
        .unwrap_or(trimmed);
    if numeric.is_empty()
        || !numeric
            .chars()
            .all(|c| c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '+' | '-'))
    {
        return None;
    }
    numeric.parse::<f64>().ok()
}

/// `Double.toString`: plain between 10^-3 and 10^7, computerized scientific notation outside.
fn java_double_to_string(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    let magnitude = value.abs();
    if magnitude == 0.0 {
        return format!("{sign}0.0");
    }
    let scientific = format!("{magnitude:e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if (1e-3..1e7).contains(&magnitude) {
        if exponent >= 0 {
            let whole_len = exponent as usize + 1;
            let mut padded = digits.clone();
            while padded.len() < whole_len {
                padded.push('0');
            }
            let (whole, fraction) = padded.split_at(whole_len);
            let fraction = if fraction.is_empty() { "0" } else { fraction };
            format!("{sign}{whole}.{fraction}")
        } else {
            let zeros = "0".repeat((-exponent - 1) as usize);
            format!("{sign}0.{zeros}{digits}")
        }
    } else {
        let (first, rest) = digits.split_at(1);
        let rest = if rest.is_empty() { "0" } else { rest };
        format!("{sign}{first}.{rest}E{exponent}")
    }
}

/// `URLDecoder.decode(s, "UTF-8")`.
fn url_decode(text: &str) -> Result<String, DecodeError> {
    let bytes = text.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'+' => {
                out.push(b' ');
                i += 1;
            }
            b'%' => {
                let mut decoded = Vec::new();
                while i + 2 < bytes.len() && bytes[i] == b'%' {
                    let hex = std::str::from_utf8(&bytes[i + 1..i + 3]).unwrap_or("");
                    let value = u8::from_str_radix(hex, 16).map_err(|_| {
                        DecodeError::Other(format!(
                            "java.lang.IllegalArgumentException: URLDecoder: Illegal hex \
                             characters in escape (%) pattern - Error at index 0 in: \"{hex}\""
                        ))
                    })?;
                    decoded.push(value);
                    i += 3;
                }
                if i < bytes.len() && bytes[i] == b'%' {
                    return Err(DecodeError::Other(
                        "java.lang.IllegalArgumentException: URLDecoder: Incomplete trailing \
                         escape (%) pattern"
                            .to_string(),
                    ));
                }
                out.extend_from_slice(String::from_utf8_lossy(&decoded).as_bytes());
            }
            b => {
                out.push(b);
                i += 1;
            }
        }
    }
    Ok(String::from_utf8_lossy(&out).into_owned())
}

/// `Gff3Writer.encodeString`: `URLEncoder.encode(s, "UTF-8")` with its `+` turned back into the
/// space it stood for.
fn url_encode(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for byte in text.bytes() {
        match byte {
            b'a'..=b'z' | b'A'..=b'Z' | b'0'..=b'9' | b'.' | b'-' | b'*' | b'_' | b' ' => {
                out.push(byte as char)
            }
            other => out.push_str(&format!("%{other:02X}")),
        }
    }
    out
}

/// `Gff3Codec.extractSingleAttribute`.
fn extract_single_attribute(values: Option<&Vec<String>>) -> Result<Option<String>, DecodeError> {
    match values {
        None => Ok(None),
        Some(v) if v.is_empty() => Ok(None),
        Some(v) if v.len() != 1 => Err(DecodeError::Tribble(
            "Attribute has multiple values when only one expected".to_string(),
        )),
        Some(v) => Ok(Some(v[0].clone())),
    }
}

/// `Gff3Codec.parseAttributes`.
fn parse_attributes(text: &str) -> Result<Vec<(String, Vec<String>)>, DecodeError> {
    let mut attributes: Vec<(String, Vec<String>)> = Vec::new();
    if text == "." {
        return Ok(attributes);
    }
    for attribute in parsing_utils_split(text, ';') {
        let key_value = parsing_utils_split(attribute, '=');
        if key_value.len() != 2 {
            return Err(DecodeError::Tribble(format!(
                "Attribute string {text} is invalid"
            )));
        }
        let key = url_decode(java_trim(key_value[0]))?;
        let mut values = Vec::new();
        for value in parsing_utils_split(java_trim(key_value[1]), ',') {
            values.push(url_decode(java_trim(value))?);
        }
        match attributes.iter_mut().find(|(k, _)| *k == key) {
            Some(entry) => entry.1 = values,
            None => attributes.push((key, values)),
        }
    }
    Ok(attributes)
}

/// `String.trim`: every character at or below the space.
fn java_trim(text: &str) -> &str {
    text.trim_matches(|c: char| c <= ' ')
}

/// `Gff3Codec.parseLine`, and the `Gff3BaseData` it builds.
fn parse_line(line: &str, current_line: i32) -> Result<Gff3Feature, DecodeError> {
    let split = parsing_utils_split(line, '\t');
    if split.len() != 9 {
        return Err(DecodeError::Tribble(format!(
            "Found an invalid number of columns in the given Gff3 file at line + {current_line} - \
             Given: {} Expected: 9 : {line}",
            split.len()
        )));
    }
    let not_a_number = || {
        DecodeError::Tribble(format!(
            "Cannot read integer value for start/end position from line {current_line}.  Line is: \
             {line}"
        ))
    };
    let contig = url_decode(split[0])?;
    let source = url_decode(split[1])?;
    let kind = url_decode(split[2])?;
    let start = java_parse_int(split[3]).ok_or_else(not_a_number)?;
    let end = java_parse_int(split[4]).ok_or_else(not_a_number)?;
    let score = if split[5] == "." {
        -1.0
    } else {
        java_parse_double(split[5]).ok_or_else(not_a_number)?
    };
    let phase = if split[7] == "." {
        -1
    } else {
        java_parse_int(split[7]).ok_or_else(not_a_number)?
    };
    // `Strand.decode`: one of `+`, `-`, `.`, and anything else is NONE.
    let strand = match split[6] {
        "+" => '+',
        "-" => '-',
        _ => '.',
    };
    let attributes = parse_attributes(split[8])?;
    let mut feature = Gff3Feature {
        contig,
        source,
        kind,
        start,
        end,
        score,
        strand,
        phase,
        attributes,
        id: None,
    };
    feature.id = extract_single_attribute(feature.attribute("ID"))?;
    extract_single_attribute(feature.attribute("Name"))?;
    Ok(feature)
}

/// The lines `LongLineBufferedReader.readLine` yields: split on `\n`, `\r` or `\r\n`.
fn read_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let bytes = text.as_bytes();
    let mut start = 0;
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'\n' => {
                lines.push(&text[start..i]);
                i += 1;
                start = i;
            }
            b'\r' => {
                lines.push(&text[start..i]);
                i += 1;
                if i < bytes.len() && bytes[i] == b'\n' {
                    i += 1;
                }
                start = i;
            }
            _ => i += 1,
        }
    }
    if start < bytes.len() {
        lines.push(&text[start..]);
    }
    lines
}

/// The file's text, gunzipped when its name says it is compressed.
fn read_text(path: &str) -> String {
    let raw = std::fs::read(path).unwrap_or_else(|e| fail(&format!("{e}")));
    let raw = if path.ends_with(".gz") || path.ends_with(".bgz") {
        htsjdk_bgzf::decompress_all(&raw).unwrap_or_else(|e| fail(&format!("{e:?}")))
    } else {
        raw
    };
    String::from_utf8_lossy(&raw).into_owned()
}

/// `Gff3Codec.canDecode`.
fn can_decode(path: &str) -> bool {
    if !GFF3_EXTENSIONS.iter().any(|e| path.ends_with(e)) {
        return false;
    }
    let text = read_text(path);
    let lines = read_lines(&text);
    let mut iter = lines.into_iter();
    let Some(mut line) = iter.next() else {
        thrown(
            "java.lang.NullPointerException: Cannot invoke \"java.lang.CharSequence.length()\" \
             because \"this.text\" is null",
        );
    };
    if to_directive(line) != Some(Directive::Version) {
        return false;
    }
    while line.starts_with('#') {
        match iter.next() {
            Some(next) => line = next,
            None => return false,
        }
    }
    let fields = parsing_utils_split(line, '\t');
    if fields.len() != 9 {
        return false;
    }
    if java_parse_int(fields[3]).is_none() || java_parse_int(fields[4]).is_none() {
        return false;
    }
    matches!(fields[6], "+" | "-" | "." | "?")
}

/// `Gff3Writer.addFeature`.
fn write_feature(out: &mut String, feature: &Gff3Feature) {
    let score = if feature.score < 0.0 {
        ".".to_string()
    } else {
        java_double_to_string(feature.score)
    };
    let phase = if feature.phase < 0 {
        ".".to_string()
    } else {
        feature.phase.to_string()
    };
    out.push_str(&format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
        url_encode(&feature.contig),
        url_encode(&feature.source),
        url_encode(&feature.kind),
        feature.start,
        feature.end,
        score,
        feature.strand,
        phase
    ));
    if feature.attributes.is_empty() {
        out.push('.');
    }
    let attributes: Vec<String> = feature
        .attributes
        .iter()
        .map(|(key, values)| {
            // The key is written as it is; only the values are escaped.
            let values: Vec<String> = values.iter().map(|v| url_encode(v)).collect();
            format!("{key}={}", values.join(","))
        })
        .collect();
    out.push_str(&attributes.join(";"));
    out.push('\n');
}

fn main() {
    let args = Args::from_env(&[
        ("I", "INPUT"),
        ("O", "OUTPUT"),
        ("SD", "SEQUENCE_DICTIONARY"),
    ]);
    let input = args.required("INPUT");
    let output = args.required("OUTPUT");
    let dictionary_path = args.get("SEQUENCE_DICTIONARY").map(str::to_string);
    let records_in_memory = args.int("nRecordsInMemory", 50000);

    // `IOUtil.assertFileIsReadable(INPUT)`, `assertFileIsWritable(OUTPUT)`.
    if !std::path::Path::new(&input).is_file() {
        thrown(&format!(
            "htsjdk.samtools.SAMException: Cannot read non-existent file: file://{}",
            absolute(&input)
        ));
    }

    if !can_decode(&input) {
        thrown(&format!(
            "java.lang.IllegalArgumentException: Input file {input} cannot be read by Gff3Codec"
        ));
    }

    let dictionary: Option<Vec<String>> = dictionary_path.as_ref().and_then(|path| {
        extract_dictionary(path)
            .unwrap_or_else(|e| fail(&format!("{e}")))
            .map(|sequences| sequences.into_iter().map(|s| s.name).collect())
    });

    if records_in_memory <= 0 {
        thrown("java.lang.IllegalArgumentException: maxRecordsInRam must be > 0");
    }

    // The whole-file iterator, and the latest start of every ID's family as it is read.
    let source = absolute(&input);
    let text = read_text(&input);
    let mut codec = Codec::default();
    let mut features: Vec<Gff3Feature> = Vec::new();
    let mut latest_start_map: HashMap<String, i32> = HashMap::new();
    for line in read_lines(&text) {
        let feature = match codec.decode(line) {
            Ok(Some(feature)) => feature,
            Ok(None) => continue,
            Err(DecodeError::Tribble(message)) => thrown(&format!(
                "htsjdk.tribble.TribbleException: {message}, for input source: {source}"
            )),
            Err(DecodeError::Other(message)) => thrown(&message),
        };
        if let Some(id) = &feature.id {
            let entry = latest_start_map.entry(id.clone()).or_insert(feature.start);
            *entry = (*entry).max(feature.start);
        }
        let current_latest_start = match &feature.id {
            None => feature.start,
            Some(id) => latest_start_map[id],
        };
        for parent in feature.attribute_or_empty("Parent") {
            let entry = latest_start_map
                .entry(parent)
                .or_insert(current_latest_start);
            *entry = (*entry).max(current_latest_start);
        }
        features.push(feature);
    }

    // `FeatureComparator` over a stable sort, which a spill and its merge leave unchanged.
    let keys: Vec<Feature> = features
        .iter()
        .enumerate()
        .map(|(index, f)| Feature {
            contig: f.contig.clone(),
            start: f.start,
            end: f.end,
            index,
        })
        .collect();
    let order = sort(&keys, dictionary.as_deref());

    // `new Gff3Writer(OUTPUT.toPath())`.
    if !GFF3_EXTENSIONS.iter().any(|e| output.ends_with(e)) {
        thrown(&format!(
            "htsjdk.tribble.TribbleException: File {output} does not have extension consistent \
             with gff3"
        ));
    }
    let mut out = String::from("##gff-version 3.1.25\n");
    for comment in &codec.comments {
        out.push('#');
        out.push_str(comment);
        out.push('\n');
    }
    for region in &codec.regions {
        out.push_str(&format!(
            "##sequence-region {} {} {}\n",
            url_encode(&region.contig),
            region.start,
            region.end
        ));
    }

    let mut latest_start = 0;
    let mut latest_chrom: Option<String> = None;
    for (n, key) in order.iter().enumerate() {
        let feature = &features[key.index];
        if n > 0
            && (feature.start > latest_start
                || latest_chrom.as_deref() != Some(feature.contig.as_str()))
        {
            out.push_str("###\n");
        }
        write_feature(&mut out, feature);
        // `updateLatestStart`.
        if latest_chrom.as_deref() != Some(feature.contig.as_str()) {
            latest_chrom = Some(feature.contig.clone());
            latest_start = 0;
        }
        if let Some(id) = &feature.id {
            latest_start = latest_start.max(latest_start_map[id]);
        }
        let parents = feature.attribute_or_empty("Parent");
        let from_parents = parents
            .iter()
            .map(|p| latest_start_map[p])
            .max()
            .unwrap_or(feature.start);
        latest_start = latest_start.max(from_parents);
    }

    let written = if output.ends_with(".gz") {
        std::fs::File::create(&output).and_then(|file| {
            let mut writer = htsjdk_bgzf::BgzfWriter::new(file);
            writer.write_all(out.as_bytes())?;
            writer.finish()
        })
    } else {
        std::fs::write(&output, out.as_bytes())
    };
    if let Err(e) = written {
        thrown(&format!(
            "picard.PicardException: Error opening  {output} to write to: {e}"
        ));
    }
}
