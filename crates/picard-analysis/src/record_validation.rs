//! The validation htsjdk 5.0.0 runs on every record a BAM reader returns.
//!
//! `BAMFileReader.BAMFileIterator.advance` numbers each record it decodes from 1 and, unless the
//! reader is `SILENT`, calls `SAMRecord.isValid(stringency == STRICT)` on it and hands the errors
//! to `SAMUtils.processValidationErrors`: `STRICT` throws the first as a `SAMFormatException`,
//! `LENIENT` prints every one to stderr. Picard sets that stringency for every reader from its
//! `VALIDATION_STRINGENCY` argument, `STRICT` by default.
//!
//! One check is lazy. `BAMRecord.getCigar` decodes the CIGAR on first use and, for a mapped read,
//! validates it then (`SAMRecord.validateCigar`) with no record number. `isValid` first asks for
//! the CIGAR when it compares its length with a non-empty sequence, so for most reads that check
//! runs in the middle of `isValid`, before the errors `isValid` itself has gathered are processed.

use std::sync::atomic::{AtomicBool, AtomicU8, Ordering};

use htsjdk_bam::cigar::{Cigar, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::reader::{BamReader, ReadError};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// htsjdk's `ValidationStringency`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stringency {
    Strict,
    Lenient,
    Silent,
}

/// 0 until read from the command line or set by a tool.
static STRINGENCY: AtomicU8 = AtomicU8::new(0);

/// The stringency every reader of this process validates with: what a tool set, or else this
/// process's `VALIDATION_STRINGENCY` argument, `STRICT` when it is absent.
pub fn stringency() -> Stringency {
    match STRINGENCY.load(Ordering::Relaxed) {
        1 => Stringency::Strict,
        2 => Stringency::Lenient,
        3 => Stringency::Silent,
        _ => {
            let args: Vec<String> = crate::args::process().into_iter().skip(1).collect();
            let found = from_args(&args);
            set_stringency(found);
            found
        }
    }
}

/// Set by a tool that never asks a record for its CIGAR.
static CIGAR_NEVER_DECODED: AtomicBool = AtomicBool::new(false);

/// For a tool that never calls `getCigar` (nor `getAlignmentEnd` or anything built on it): a read
/// with no sequence then never has its CIGAR decoded, so never validated.
pub fn cigar_never_decoded() {
    CIGAR_NEVER_DECODED.store(true, Ordering::Relaxed);
}

/// `SamReaderFactory.validationStringency` for a tool that overrides its argument.
pub fn set_stringency(stringency: Stringency) {
    let code = match stringency {
        Stringency::Strict => 1,
        Stringency::Lenient => 2,
        Stringency::Silent => 3,
    };
    STRINGENCY.store(code, Ordering::Relaxed);
}

/// The last `VALIDATION_STRINGENCY` the arguments give, in the legacy `NAME=value` form or the
/// `--NAME value` one. A value that is not a level is left to the tool's own parser to refuse.
pub fn from_args(args: &[String]) -> Stringency {
    let mut found = Stringency::Strict;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        let value = match arg.trim_start_matches('-').split_once('=') {
            Some(("VALIDATION_STRINGENCY", value)) => Some(value),
            _ if arg.starts_with('-') && arg.trim_start_matches('-') == "VALIDATION_STRINGENCY" => {
                iter.next().map(String::as_str)
            }
            _ => None,
        };
        match value {
            Some("STRICT") => found = Stringency::Strict,
            Some("LENIENT") => found = Stringency::Lenient,
            Some("SILENT") => found = Stringency::Silent,
            _ => {}
        }
    }
    found
}

/// A `SAMValidationError`, as `toString` prints it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationError {
    pub kind: &'static str,
    pub read_name: String,
    pub message: String,
    /// Set by `processValidationErrors`; printed only when positive.
    pub record_number: i64,
}

impl std::fmt::Display for ValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let severity = match self.kind {
            "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND" | "ADJACENT_INDEL_IN_CIGAR" => "WARNING",
            _ => "ERROR",
        };
        write!(f, "{severity}::{}:", self.kind)?;
        if self.record_number > 0 {
            write!(f, "Record {}, ", self.record_number)?;
        }
        write!(f, "Read name {}, {}", self.read_name, self.message)
    }
}

/// What a reader's header tells the checks: the dictionary htsjdk resolves indexes against (the
/// text's, or the binary one when the text has none) and the read group IDs.
pub struct Validator {
    sequences: Vec<(String, i32)>,
    read_groups: Vec<String>,
    record_number: i64,
    stringency: Stringency,
}

impl Validator {
    pub fn new(header: &SamHeader, binary_sequences: &[(String, i32)]) -> Self {
        let sequences = if header.sequences.is_empty() {
            binary_sequences.to_vec()
        } else {
            header
                .sequences
                .iter()
                .map(|s| (s.name.clone(), s.length))
                .collect()
        };
        Validator {
            sequences,
            read_groups: header.read_groups.iter().map(|rg| rg.id.clone()).collect(),
            record_number: 0,
            stringency: stringency(),
        }
    }

    /// `advance`'s checks on the next record: the exception it throws, if any, as
    /// `Exception in thread "main"` prints it. `LENIENT` errors are printed here.
    pub fn check(&mut self, rec: &BamRecord) -> Result<(), String> {
        self.record_number += 1;
        // The BAMRecord constructor resolves both indexes, at any stringency.
        for index in [rec.reference_index, rec.mate_reference_index] {
            if index != -1 && self.sequence(index).is_none() {
                return Err(format!(
                    "java.lang.IllegalArgumentException: Reference name for '{index}' not found in sequence dictionary."
                ));
            }
        }
        if self.stringency == Stringency::Silent {
            return Ok(());
        }
        let errors = self.is_valid(rec, self.stringency == Stringency::Strict)?;
        self.process(errors, self.record_number)?;
        // A read whose sequence is empty never reaches the CIGAR in `isValid`; it is decoded when
        // the tool first asks for it. shortcut: checked right after `isValid`, as a tool that asks
        // before reading on does; a tool that asks later would print or throw later.
        if rec.read_bases.is_empty() && !CIGAR_NEVER_DECODED.load(Ordering::Relaxed) {
            self.process(self.validate_cigar(rec), -1)?;
        }
        Ok(())
    }

    /// [`Self::check`] on a record [`Self::passes`] passed: it is counted, and nothing else.
    pub fn skip(&mut self) {
        self.record_number += 1;
    }

    /// Whether [`Self::check`] passes `rec` in silence: nothing thrown and nothing printed.
    pub fn passes(&self, rec: &BamRecord) -> bool {
        if [rec.reference_index, rec.mate_reference_index]
            .into_iter()
            .any(|index| index != -1 && self.sequence(index).is_none())
        {
            return false;
        }
        if self.stringency == Stringency::Silent {
            return true;
        }
        matches!(self.is_valid(rec, self.stringency == Stringency::Strict), Ok(errors) if errors.is_empty())
            && !(rec.read_bases.is_empty()
                && !CIGAR_NEVER_DECODED.load(Ordering::Relaxed)
                && !self.validate_cigar(rec).is_empty())
    }

    /// `SAMUtils.processValidationErrors`.
    fn process(&self, mut errors: Vec<ValidationError>, record_number: i64) -> Result<(), String> {
        if errors.is_empty() {
            return Ok(());
        }
        for e in &mut errors {
            e.record_number = record_number;
        }
        match self.stringency {
            Stringency::Strict => Err(format!(
                "htsjdk.samtools.SAMFormatException: SAM validation error: {}",
                errors[0]
            )),
            Stringency::Lenient => {
                for e in &errors {
                    eprintln!("Ignoring SAM validation error: {e}");
                }
                Ok(())
            }
            Stringency::Silent => Ok(()),
        }
    }

    fn sequence(&self, index: i32) -> Option<&(String, i32)> {
        usize::try_from(index)
            .ok()
            .and_then(|i| self.sequences.get(i))
    }

    /// `SAMRecord.isValid(firstOnly)`. `Err` is the CIGAR check thrown from inside it.
    fn is_valid(&self, rec: &BamRecord, first_only: bool) -> Result<Vec<ValidationError>, String> {
        let mut ret = Vec::new();
        let flags = rec.flags;
        let add = |ret: &mut Vec<ValidationError>, kind: &'static str, message: String| {
            ret.push(ValidationError {
                kind,
                read_name: rec.read_name.clone(),
                message,
                record_number: 0,
            });
            first_only
        };
        let mate_name = if rec.mate_reference_index == -1 {
            "*"
        } else {
            self.sequences[rec.mate_reference_index as usize].0.as_str()
        };
        if flags & 0x1 == 0 {
            for (bit, kind, what) in [
                (0x2, "INVALID_FLAG_PROPER_PAIR", "Proper pair"),
                (0x8, "INVALID_FLAG_MATE_UNMAPPED", "Mate unmapped"),
                (0x20, "INVALID_FLAG_MATE_NEG_STRAND", "Mate negative strand"),
                (0x40, "INVALID_FLAG_FIRST_OF_PAIR", "First of pair"),
                (0x80, "INVALID_FLAG_SECOND_OF_PAIR", "Second of pair"),
            ] {
                if flags & bit != 0
                    && add(
                        &mut ret,
                        kind,
                        format!("{what} flag should not be set for unpaired read."),
                    )
                {
                    return Ok(ret);
                }
            }
            if rec.mate_reference_index != -1 {
                if add(
                    &mut ret,
                    "INVALID_MATE_REF_INDEX",
                    "MRNM should not be set for unpaired read.".into(),
                ) {
                    return Ok(ret);
                }
                if add(
                    &mut ret,
                    "INVALID_UNPAIRED_MATE_REFERENCE",
                    format!("Unpaired read mate reference is {mate_name} not * for unpaired read"),
                ) {
                    return Ok(ret);
                }
            }
        } else {
            let errors = self.reference_and_position(
                rec,
                rec.mate_reference_index,
                rec.mate_alignment_start,
                true,
                first_only,
            );
            if !errors.is_empty() {
                if first_only {
                    return Ok(errors);
                }
                ret.extend(errors);
            }
            if rec.mate_reference_index == -1
                && flags & 0x8 == 0
                && add(
                    &mut ret,
                    "INVALID_FLAG_MATE_UNMAPPED",
                    "Mapped mate should have mate reference name".into(),
                )
            {
                return Ok(ret);
            }
            if flags & 0xC0 == 0
                && add(
                    &mut ret,
                    "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND",
                    "Paired read should be marked as first of pair or second of pair.".into(),
                )
            {
                return Ok(ret);
            }
        }
        // MAX_INSERT_SIZE is Integer.MAX_VALUE, so only its negation's neighbour is out of range.
        if rec.inferred_insert_size == i32::MIN
            && add(
                &mut ret,
                "INVALID_INSERT_SIZE",
                "Insert size out of range".into(),
            )
        {
            return Ok(ret);
        }
        if flags & 0x4 != 0 {
            for (hit, kind, message) in [
                (
                    flags & 0x100 != 0,
                    "INVALID_FLAG_NOT_PRIM_ALIGNMENT",
                    "Secondary alignment flag should not be set for unmapped read.",
                ),
                (
                    flags & 0x800 != 0,
                    "INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT",
                    "Supplementary alignment flag should not be set for unmapped read.",
                ),
                (
                    rec.mapping_quality != 0,
                    "INVALID_MAPPING_QUALITY",
                    "MAPQ should be 0 for unmapped read.",
                ),
            ] {
                if hit && add(&mut ret, kind, message.into()) {
                    return Ok(ret);
                }
            }
        } else {
            // MAPQ >= 256 cannot be stored in a BAM's single byte.
            for (hit, kind, message) in [
                (
                    rec.cigar.elements.is_empty(),
                    "INVALID_CIGAR",
                    "CIGAR should have > zero elements for mapped read.",
                ),
                (
                    self.sequences.is_empty(),
                    "MISSING_SEQUENCE_DICTIONARY",
                    "Empty sequence dictionary.",
                ),
                (
                    rec.reference_index == -1,
                    "INVALID_FLAG_READ_UNMAPPED",
                    "Mapped read should have valid reference name",
                ),
            ] {
                if hit && add(&mut ret, kind, message.into()) {
                    return Ok(ret);
                }
            }
        }
        if let Some(TagValue::Str(id)) = rec.tags.get(Tag::new(b"RG")) {
            if !self.read_groups.iter().any(|rg| rg == id)
                && add(
                    &mut ret,
                    "READ_GROUP_NOT_FOUND",
                    format!("RG ID on SAMRecord not found in header: {id}"),
                )
            {
                return Ok(ret);
            }
        }
        // htsjdk passes no firstOnly here: every error of the read's own position is gathered.
        let errors = self.reference_and_position(
            rec,
            rec.reference_index,
            rec.alignment_start,
            false,
            false,
        );
        if !errors.is_empty() {
            ret.extend(errors);
            if first_only {
                return Ok(ret);
            }
        }
        // MISMATCH_READ_LENGTH_AND_QUALS_LENGTH and MISMATCH_SEQ_QUAL_LENGTH cannot happen in a
        // BAM: the qualities are stored at the sequence's length, or absent.
        if mate_name == "*"
            && rec.mate_alignment_start != 0
            && add(
                &mut ret,
                "INVALID_UNALIGNED_MATE_START",
                // htsjdk quotes the read's own start here.
                format!(
                    "The unaligned mate start position is {}, should be 0",
                    rec.alignment_start
                ),
            )
        {
            return Ok(ret);
        }
        if !rec.read_bases.is_empty() {
            self.process(self.validate_cigar(rec), -1)?;
            let cigar_length = rec.cigar.read_length() as usize;
            if cigar_length != 0
                && cigar_length != rec.read_bases.len()
                && add(
                    &mut ret,
                    "MISMATCH_CIGAR_SEQ_LENGTH",
                    format!(
                        "CIGAR covers {cigar_length} bases but the sequence is {} read bases ",
                        rec.read_bases.len()
                    ),
                )
            {
                return Ok(ret);
            }
        }
        Ok(ret)
    }

    /// `SAMRecord.isValidReferenceIndexAndPosition`.
    fn reference_and_position(
        &self,
        rec: &BamRecord,
        index: i32,
        start: i32,
        is_mate: bool,
        first_only: bool,
    ) -> Vec<ValidationError> {
        let mut ret = Vec::new();
        let prefix = if is_mate { "Mate " } else { "" };
        let mut add = |message: String| {
            ret.push(ValidationError {
                // INVALID_REFERENCE_INDEX cannot happen: the index resolved when the record was built.
                kind: "INVALID_ALIGNMENT_START",
                read_name: rec.read_name.clone(),
                message: format!("{prefix}{message}"),
                record_number: 0,
            });
            first_only
        };
        if index == -1 {
            if start != 0 {
                add("Alignment start should be 0 because reference name = *.".into());
            }
            return ret;
        }
        if start == 0 && add("Alignment start should != 0 because reference name != *.".into()) {
            return ret;
        }
        if let Some((name, length)) = self.sequence(index) {
            if start > *length {
                add(format!(
                    "Alignment start ({start}) must be <= reference sequence length ({length}) on reference {name}"
                ));
            }
        }
        ret
    }

    /// `SAMRecord.validateCigar` as `BAMRecord.getCigar` runs it: mapped reads only.
    fn validate_cigar(&self, rec: &BamRecord) -> Vec<ValidationError> {
        if rec.flags & 0x4 != 0 {
            return Vec::new();
        }
        let error = |kind: &'static str, message: &str| ValidationError {
            kind,
            read_name: rec.read_name.clone(),
            message: message.to_string(),
            record_number: 0,
        };
        let mut ret: Vec<ValidationError> = cigar_errors(&rec.cigar)
            .into_iter()
            .map(|(kind, message)| error(kind, &message))
            .collect();
        if let Some((_, length)) = self.sequence(rec.reference_index) {
            // `SAMUtils.getAlignmentBlocks`: the first M, = or X block past the reference's end.
            let mut ref_base = i64::from(rec.alignment_start);
            for e in &rec.cigar.elements {
                let len = i64::from(e.length);
                match e.op {
                    Op::M | Op::Eq | Op::X => {
                        if ref_base + len - 1 > i64::from(*length) {
                            ret.push(error(
                                "CIGAR_MAPS_OFF_REFERENCE",
                                "Read CIGAR M operator maps off end of reference",
                            ));
                            break;
                        }
                        ref_base += len;
                    }
                    Op::N | Op::D => ref_base += len,
                    _ => {}
                }
            }
        }
        ret
    }
}

/// `Cigar.isValid`: each error's type and message, in htsjdk's order.
pub fn cigar_errors(cigar: &Cigar) -> Vec<(&'static str, String)> {
    let elements = &cigar.elements;
    if elements.is_empty() {
        return Vec::new();
    }
    let real = |op: Op| matches!(op, Op::M | Op::Eq | Op::X | Op::I | Op::D | Op::N);
    let indel = |op: Op| matches!(op, Op::I | Op::D);
    let last = elements.len() - 1;
    let mut ret = Vec::new();
    let invalid =
        |ret: &mut Vec<_>, message: &str| ret.push(("INVALID_CIGAR", message.to_string()));
    let mut seen_real = false;
    for (i, element) in elements.iter().enumerate() {
        if element.length == 0 {
            invalid(&mut ret, "CIGAR element with zero length");
        }
        let op = element.op;
        if op == Op::H {
            if i != 0 && i != last {
                invalid(
                    &mut ret,
                    "Hard clipping operator not at start or end of CIGAR",
                );
            }
        } else if op == Op::S {
            // At either end a soft clip is fine.
            if i == 0 || i == last {
                continue;
            }
            if i == 1 {
                if !(elements.len() == 3 && elements[2].op == Op::H) && elements[0].op != Op::H {
                    invalid(
                        &mut ret,
                        "Soft clipping CIGAR operator can only be inside of hard clipping operator",
                    );
                }
            } else if i == last - 1 {
                if elements[last].op != Op::H {
                    invalid(
                        &mut ret,
                        "Soft clipping CIGAR operator can only be inside of hard clipping operator",
                    );
                }
            } else {
                invalid(
                    &mut ret,
                    "Soft clipping CIGAR operator can at start or end of read, or be inside of hard clipping operator",
                );
            }
        } else if real(op) {
            seen_real = true;
            if indel(op) {
                for next in &elements[i + 1..] {
                    if (real(next.op) && !indel(next.op)) || next.op == Op::P {
                        break;
                    }
                    if next.op == op {
                        let name = if op == Op::I { "I" } else { "D" };
                        ret.push((
                            "ADJACENT_INDEL_IN_CIGAR",
                            format!("No M or N operator between pair of {name} operators in CIGAR"),
                        ));
                    }
                }
            }
        } else if op == Op::P {
            // A leading pad is allowed, for a read starting in the middle of a pad.
            if i == last && i != 0 {
                invalid(&mut ret, "Padding operator not valid at end of CIGAR");
            } else if i != 0 && (!real(elements[i - 1].op) || !real(elements[i + 1].op)) {
                invalid(
                    &mut ret,
                    "Padding operator not between real operators in CIGAR",
                );
            }
        }
    }
    if !seen_real {
        invalid(&mut ret, "No real operator (M|I|D|N) in CIGAR");
    }
    ret
}

/// A [`BamReader`] whose records are validated as htsjdk's reader validates them: a `STRICT`
/// failure ends the process as the uncaught exception ends Picard.
pub struct Validated<'a> {
    reader: BamReader<'a>,
    validator: Validator,
}

impl<'a> std::ops::Deref for Validated<'a> {
    type Target = BamReader<'a>;
    fn deref(&self) -> &BamReader<'a> {
        &self.reader
    }
}

impl<'a> std::ops::DerefMut for Validated<'a> {
    fn deref_mut(&mut self) -> &mut BamReader<'a> {
        &mut self.reader
    }
}

impl Iterator for Validated<'_> {
    type Item = Result<BamRecord, ReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        let item = self.reader.next()?;
        if let Ok(rec) = &item {
            if let Err(message) = self.validator.check(rec) {
                crate::metrics_cli::thrown(&message);
            }
        }
        Some(item)
    }
}

/// `BamReader::new`, its records validated at this process's stringency.
pub fn reader(data: &[u8]) -> Result<Validated<'_>, ReadError> {
    let reader = BamReader::new(data)?;
    let validator = Validator::new(&reader.header.text, &reader.header.binary_sequences);
    Ok(Validated { reader, validator })
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::cigar::CigarElement;
    use htsjdk_bam::header::SequenceRecord;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn the_last_stringency_given_in_either_form_wins() {
        assert_eq!(from_args(&args(&["I=x.bam"])), Stringency::Strict);
        assert_eq!(
            from_args(&args(&["VALIDATION_STRINGENCY=SILENT"])),
            Stringency::Silent
        );
        assert_eq!(
            from_args(&args(&[
                "--VALIDATION_STRINGENCY",
                "SILENT",
                "-VALIDATION_STRINGENCY",
                "LENIENT"
            ])),
            Stringency::Lenient
        );
    }

    #[test]
    fn a_mapped_mate_without_a_reference_is_refused_with_its_record_number() {
        let mut header = SamHeader::default();
        header.sequences.push(SequenceRecord::new("chr1", 1000));
        let mut validator = Validator::new(&header, &[]);
        validator.stringency = Stringency::Strict;
        let mut rec = BamRecord {
            read_name: "read3".into(),
            flags: 0x41,
            reference_index: 0,
            alignment_start: 10,
            read_bases: b"ACGT".to_vec(),
            ..BamRecord::default()
        };
        rec.cigar = Cigar::new(vec![CigarElement {
            length: 4,
            op: Op::M,
        }]);
        assert_eq!(validator.check(&rec), Err(
            "htsjdk.samtools.SAMFormatException: SAM validation error: ERROR::INVALID_FLAG_MATE_UNMAPPED:Record 1, Read name read3, Mapped mate should have mate reference name".into()
        ));
    }

    #[test]
    fn cigar_errors_follow_htsjdk_order() {
        let cigar = |ops: &[(u32, Op)]| {
            Cigar::new(
                ops.iter()
                    .map(|&(length, op)| CigarElement { length, op })
                    .collect(),
            )
        };
        assert_eq!(
            cigar_errors(&cigar(&[(2, Op::M), (1, Op::I), (1, Op::I)])),
            vec![(
                "ADJACENT_INDEL_IN_CIGAR",
                "No M or N operator between pair of I operators in CIGAR".to_string()
            )]
        );
        assert!(cigar_errors(&cigar(&[
            (1, Op::H),
            (1, Op::S),
            (2, Op::M),
            (1, Op::S),
            (1, Op::H)
        ]))
        .is_empty());
        assert_eq!(
            cigar_errors(&cigar(&[(4, Op::S)])),
            vec![(
                "INVALID_CIGAR",
                "No real operator (M|I|D|N) in CIGAR".to_string()
            )]
        );
    }
}
