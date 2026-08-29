//! The header `MergeBamAlignment` builds, which comes from three files rather than one.
//!
//! The sequences are the reference dictionary's, the read groups are the UNMAPPED bam's, and the
//! programs are the aligned one's, but only under a condition: the merger adopts a program record
//! when the aligned header holds EXACTLY one. None and it has nothing to adopt; two and it has no
//! way to choose, so the output carries no `@PG` at all.
//!
//! A program record named on the command line is set before the aligned files are opened, and the
//! merger refuses to set one twice, so the caller's record is not joined by the aligner's: it
//! replaces it.
//!
//! Ported from `picard.sam.AbstractAlignmentMerger`, `picard.sam.SamAlignmentMerger` and
//! `picard.sam.MergeBamAlignment` in Picard 3.4.0.

/// One `@PG` line, by its id and the rest of its fields in the order they were written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProgramRecord {
    pub id: String,
    pub fields: Vec<(String, String)>,
}

impl ProgramRecord {
    /// The line as a header writes it, the id first.
    pub fn line(&self) -> String {
        let mut text = format!("@PG\tID:{}", self.id);
        for (tag, value) in &self.fields {
            text.push('\t');
            text.push_str(&format!("{tag}:{value}"));
        }
        text
    }

    /// The record `--PROGRAM_RECORD_ID` and its fellows build.
    ///
    /// The four arguments come together or not at all, and the order the fields are written in is
    /// the order they are set: the version and the command line before the name.
    pub fn from_command_line(
        id: &str,
        name: Option<&str>,
        version: &str,
        command_line: &str,
    ) -> ProgramRecord {
        let mut fields = vec![
            ("VN".to_string(), version.to_string()),
            ("CL".to_string(), command_line.to_string()),
        ];
        if let Some(name) = name {
            fields.push(("PN".to_string(), name.to_string()));
        }
        ProgramRecord {
            id: id.to_string(),
            fields,
        }
    }
}

/// The order the output is written in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Coordinate,
    Queryname,
    Unsorted,
}

impl SortOrder {
    pub fn name(&self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
            SortOrder::Unsorted => "unsorted",
        }
    }
}

/// What a run was asked for, as far as the header is concerned.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Options {
    pub sort_order: SortOrder,
    /// The record the command line asks for, which is set before any file is opened.
    pub program_record: Option<ProgramRecord>,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            sort_order: SortOrder::Coordinate,
            program_record: None,
        }
    }
}

/// What a run is refused for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// The four arguments that make a program record were not all supplied.
    IncompleteProgramGroup,
    /// The unmapped bam already declares the id the command line asked for.
    ProgramRecordAlreadyInUse,
    /// The aligned header and the dictionary disagree about a sequence's length.
    SequenceLengthsDiffer {
        name: String,
        first: i64,
        second: i64,
    },
    /// They disagree about which sequences there are at all, and the message names both lists,
    /// the aligned header's first.
    DifferentSequences {
        aligned: Vec<String>,
        dictionary: Vec<String>,
    },
}

impl Refusal {
    /// The words the reference refuses in.
    pub fn message(&self) -> String {
        match self {
            Refusal::IncompleteProgramGroup => "PROGRAM_RECORD_ID, PROGRAM_GROUP_VERSION, and \
                 PROGRAM_GROUP_COMMAND_LINE must all be supplied or none should be included."
                .to_string(),
            Refusal::ProgramRecordAlreadyInUse => {
                "Program Record ID already in use in unmapped BAM file.".to_string()
            }
            Refusal::SequenceLengthsDiffer {
                name,
                first,
                second,
            } => format!(
                "Cannot merge the two dictionaries. Found sequence entry for which lengths \
                 differ: {name} has lengths {first} and {second}"
            ),
            Refusal::DifferentSequences {
                aligned,
                dictionary,
            } => format!(
                "Do not use this function to merge dictionaries with different sequences in them. \
                 Sequences must be in the same order as well. Found [{}] and [{}].",
                aligned.join(", "),
                dictionary.join(", ")
            ),
        }
    }
}

/// Which four arguments were given, and whether that is a whole program record.
///
/// The name is the one that may be left out; the other three come together.
pub fn command_line_program(
    id: Option<&str>,
    name: Option<&str>,
    version: Option<&str>,
    command_line: Option<&str>,
) -> Result<Option<ProgramRecord>, Refusal> {
    match (id, version, command_line) {
        (None, None, None) => Ok(None),
        (Some(id), Some(version), Some(command_line)) => Ok(Some(
            ProgramRecord::from_command_line(id, name, version, command_line),
        )),
        _ => Err(Refusal::IncompleteProgramGroup),
    }
}

/// The program record the merged header carries, if any.
///
/// A record from the command line wins outright: the merger refuses to set one twice, and it is
/// set before the aligned files are opened. Otherwise the aligned header's is adopted only when
/// there is exactly one of it.
pub fn adopted_program(
    aligned_programs: &[ProgramRecord],
    from_the_command_line: Option<&ProgramRecord>,
) -> Option<ProgramRecord> {
    if let Some(record) = from_the_command_line {
        return Some(record.clone());
    }
    if aligned_programs.len() == 1 {
        return Some(aligned_programs[0].clone());
    }
    None
}

/// Whether the aligned header's sequences agree with the dictionary's, and how they do not.
pub fn check_dictionary(
    dictionary: &[(String, i64)],
    aligned: &[(String, i64)],
) -> Result<(), Refusal> {
    // A sequence of another name is refused for the SET of sequences rather than for the name, and
    // the count is what the message names.
    if dictionary.len() != aligned.len()
        || dictionary
            .iter()
            .zip(aligned.iter())
            .any(|(left, right)| left.0 != right.0)
    {
        return Err(Refusal::DifferentSequences {
            aligned: aligned.iter().map(|entry| entry.0.clone()).collect(),
            dictionary: dictionary.iter().map(|entry| entry.0.clone()).collect(),
        });
    }
    for (left, right) in dictionary.iter().zip(aligned.iter()) {
        if left.1 != right.1 {
            return Err(Refusal::SequenceLengthsDiffer {
                name: left.0.clone(),
                first: right.1,
                second: left.1,
            });
        }
    }
    Ok(())
}

/// The merged header, line by line.
///
/// The sequence lines are the dictionary's own, copied whole: the `M5` and the `UR` come with
/// them, which is why the `UR` is the dictionary's canonical path and not the spelling the command
/// line used for the reference. The read groups are the unmapped bam's. The comments of neither
/// input are carried.
pub fn merged_header(
    dictionary_lines: &[String],
    unmapped_read_groups: &[String],
    aligned_programs: &[ProgramRecord],
    options: &Options,
) -> Vec<String> {
    let mut lines = vec![format!("@HD\tVN:1.6\tSO:{}", options.sort_order.name())];
    lines.extend(
        dictionary_lines
            .iter()
            .filter(|line| line.starts_with("@SQ"))
            .cloned(),
    );
    lines.extend(unmapped_read_groups.iter().cloned());
    if let Some(program) = adopted_program(aligned_programs, options.program_record.as_ref()) {
        lines.push(program.line());
    }
    lines
}
