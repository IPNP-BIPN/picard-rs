//! Conformance for `MergeBamAlignment` against Picard, through the port's own binary.
//!
//! Every run of the six dumps in `tools/mergebamalignment-conformance/` is replayed here with the
//! same inputs and the same arguments, and the rows the dump would print are rebuilt from what the
//! binary wrote: the record lines, the header, the exit status, the exception or the refusal. The
//! rebuilt text must equal the golden line for line, so a case the golden has and this file does
//! not (or the reverse) fails as well.
//!
//! What each golden measures:
//!
//!  * `merge_bam_alignment`: the merged record of a fragment and of a proper pair, its `NM`/`MD`/`UQ`
//!    and the coordinate order of the output;
//!  * `merge_bam_alignment_full`: the whole output file, header included (`UR` stripped, being the
//!    run's own temp path);
//!  * `merge_bam_alignment_header`: the header built from the dictionary, the unmapped bam's read
//!    groups and the aligned bams' programs, and the refusals;
//!  * `merge_bam_alignment_pair`: the mate fields, insert size, proper-pair flag, `MC`/`MQ` and the
//!    clipping of overlapping ends;
//!  * `merge_bam_alignment_clip`: the clip off the end of a contig, the adapter clip, and the
//!    contamination filter under its three unmapping strategies;
//!  * `merge_bam_alignment_hits`: which alignment the four strategies call primary.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

fn golden(name: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{name}.txt.gz"));
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("the golden"))
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

/// `ACGT` repeating, as every dump's reference is.
fn acgt(length: usize) -> String {
    "ACGT".chars().cycle().take(length).collect()
}

fn quals(length: usize) -> String {
    "I".repeat(length)
}

/// How a dump ran the tool and what it printed of the result.
#[derive(Clone, Copy, PartialEq)]
enum Shape {
    /// `unmapped`, `aligned`, `rc`, then the record lines.
    Records,
    /// `dict`, `unmapped`, `aligned`, `rc`, then the whole output.
    Full,
    /// `rc` and the record lines, or the exception.
    Outcome,
    /// `rc` and the header, or the refusal, or the exception.
    Header,
}

/// The rows of one dump, rebuilt from the binary's runs.
struct Dump {
    shape: Shape,
    contigs: Vec<(&'static str, usize)>,
    base: PathBuf,
    count: usize,
    text: String,
}

impl Dump {
    fn new(name: &str, shape: Shape, contigs: Vec<(&'static str, usize)>) -> Dump {
        let base = std::env::temp_dir().join(format!("mba-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        Dump {
            shape,
            contigs,
            base,
            count: 0,
            text: String::new(),
        }
    }

    fn emit(&mut self, kind: &str, case: &str, payload: &str) {
        let payload = payload
            .replace('\\', "\\\\")
            .replace('\t', "\\t")
            .replace('\n', "\\n");
        self.text.push_str(&format!("{kind}\t{case}\t{payload}\n"));
    }

    fn run(&mut self, case: &str, unmapped: &str, aligned: &str, extra: &[&str]) {
        self.run_from(case, unmapped, &[aligned], false, extra);
    }

    fn run_from(
        &mut self,
        case: &str,
        unmapped: &str,
        aligned: &[&str],
        long_way_round: bool,
        extra: &[&str],
    ) {
        self.count += 1;
        let dir = self.base.join(self.count.to_string());
        std::fs::create_dir_all(&dir).unwrap();
        let dir_name = dir.display().to_string();
        let reference = dir.join("ref.fasta");
        let fasta: String = self
            .contigs
            .iter()
            .map(|(name, length)| format!(">{name}\n{}\n", acgt(*length)))
            .collect();
        std::fs::write(&reference, fasta).unwrap();
        let dict = dir.join("ref.dict");
        let made = Command::new(env!("CARGO_BIN_EXE_CreateSequenceDictionary"))
            .arg(format!("R={}", reference.display()))
            .arg(format!("O={}", dict.display()))
            .output()
            .unwrap();
        assert!(made.status.success(), "{case}: no dictionary");

        let unmapped_file = dir.join("u.sam");
        std::fs::write(&unmapped_file, unmapped).unwrap();
        let mut argv = vec![format!("UNMAPPED_BAM={}", unmapped_file.display())];
        for (index, text) in aligned.iter().enumerate() {
            // The header dump numbers its aligned files; the others call the one file `a.sam`.
            let file = if self.shape == Shape::Header {
                dir.join(format!("a{index}.sam"))
            } else {
                dir.join("a.sam")
            };
            std::fs::write(&file, text).unwrap();
            argv.push(format!("ALIGNED_BAM={}", file.display()));
        }
        let reference_arg = if long_way_round {
            std::fs::create_dir_all(dir.join("sub")).unwrap();
            format!("{dir_name}/sub/../ref.fasta")
        } else {
            reference.display().to_string()
        };
        argv.push(format!("REFERENCE_SEQUENCE={reference_arg}"));
        let output = dir.join("o.sam");
        argv.push(format!("OUTPUT={}", output.display()));
        if matches!(self.shape, Shape::Outcome) {
            argv.push("SORT_ORDER=queryname".to_string());
        }
        argv.extend(extra.iter().map(|arg| arg.to_string()));

        let ran = Command::new(env!("CARGO_BIN_EXE_MergeBamAlignment"))
            .args(&argv)
            .output()
            .unwrap();
        let said = format!(
            "{}{}",
            String::from_utf8_lossy(&ran.stdout),
            String::from_utf8_lossy(&ran.stderr)
        );
        let code = ran.status.code().unwrap_or(-1);
        let written = std::fs::read_to_string(&output).ok();

        if matches!(self.shape, Shape::Outcome | Shape::Header) {
            // The dumps print the root cause as `class:message`, the binary as Java's uncaught
            // throwable, `class: message`.
            if let Some(thrown) = said
                .lines()
                .find_map(|line| line.strip_prefix("Exception in thread \"main\" "))
            {
                let thrown = thrown.replacen(": ", ":", 1).replace(&dir_name, "<dir>");
                self.emit("error", case, &thrown);
                return;
            }
        }
        match self.shape {
            Shape::Records => {
                self.emit("unmapped", case, unmapped);
                self.emit("aligned", case, aligned[0]);
                self.emit("rc", case, &code.to_string());
                self.records(case, written.as_deref());
            }
            Shape::Full => {
                let dict_text = std::fs::read_to_string(&dict).unwrap();
                self.emit("dict", case, &dict_text);
                self.emit("unmapped", case, unmapped);
                self.emit("aligned", case, aligned[0]);
                self.emit("rc", case, &code.to_string());
                self.emit("full", case, written.as_deref().unwrap_or(""));
            }
            Shape::Outcome => {
                self.emit("rc", case, &code.to_string());
                self.records(case, written.as_deref());
            }
            Shape::Header => {
                self.emit("rc", case, &code.to_string());
                if code != 0 {
                    // The parser's refusals are `ERROR` lines; the tool's own is the last line.
                    let lines: Vec<String> = said
                        .lines()
                        .filter(|line| !line.trim().is_empty())
                        .map(|line| line.replace(&dir_name, "<dir>"))
                        .collect();
                    let mut reasons: Vec<String> = lines
                        .iter()
                        .filter(|line| line.starts_with("ERROR"))
                        .cloned()
                        .collect();
                    if reasons.is_empty() {
                        reasons.extend(lines.last().cloned());
                    }
                    self.emit("refusal", case, &reasons.join("\n"));
                    return;
                }
                if let Some(written) = written {
                    let header: Vec<String> = written
                        .split('\n')
                        .filter(|line| line.starts_with('@'))
                        .map(|line| line.replace(&dir_name, "<dir>"))
                        .collect();
                    self.emit("header", case, &header.join("\n"));
                }
            }
        }
    }

    fn records(&mut self, case: &str, written: Option<&str>) {
        for line in written
            .unwrap_or("")
            .split('\n')
            .filter(|line| !line.starts_with('@') && !line.is_empty())
        {
            self.emit("record", case, line);
        }
    }

    /// The rebuilt dump against its golden, line by line.
    fn check(self, golden_name: &str) {
        let _ = std::fs::remove_dir_all(&self.base);
        let strip = |text: &str| -> Vec<String> {
            text.lines()
                .map(|line| {
                    if self.shape != Shape::Full {
                        return line.to_string();
                    }
                    // The `UR` is the run's own temp path, in the dictionary and in the header.
                    let mut out = String::new();
                    let mut rest = line;
                    while let Some(at) = rest.find("\\tUR:") {
                        out.push_str(&rest[..at]);
                        rest = &rest[at + 5..];
                        let end = [rest.find("\\t"), rest.find("\\n")]
                            .into_iter()
                            .flatten()
                            .min()
                            .unwrap_or(rest.len());
                        rest = &rest[end..];
                    }
                    out.push_str(rest);
                    out
                })
                .collect()
        };
        let ours = strip(&self.text);
        let picard = strip(&golden(golden_name));
        let differing: Vec<String> = (0..ours.len().max(picard.len()))
            .filter(|&index| ours.get(index) != picard.get(index))
            .map(|index| {
                format!(
                    "line {}:\n  ours:   {:?}\n  picard: {:?}",
                    index + 1,
                    ours.get(index),
                    picard.get(index)
                )
            })
            .collect();
        assert!(
            differing.is_empty(),
            "{golden_name}: {} lines differ\n{}",
            differing.len(),
            differing.join("\n")
        );
    }
}

const UNMAPPED_HEADER: &str =
    "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n";
const ALIGNED_HEADER_40: &str = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";

fn aligned_header_200() -> String {
    "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:200\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n".to_string()
}

/// `MergeBamAlignmentDump`: the merged record of a fragment, a reverse hit, three reads and a pair.
#[test]
fn every_merged_record_is_picards() {
    let (u, a) = (UNMAPPED_HEADER, ALIGNED_HEADER_40);
    let mut dump = Dump::new("records", Shape::Records, vec![("chr1", 40)]);
    dump.run(
        "forward_match",
        &format!("{u}r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tab:Z:keepme\n"),
        &format!(
            "{a}r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\tX0:i:1\n"
        ),
        &[],
    );
    dump.run(
        "forward_mismatch",
        &format!("{u}r2\t4\t*\t0\t0\t*\t*\t0\t0\tAAGTACGT\t#IIIIIII\tRG:Z:rg1\n"),
        &format!("{a}r2\t0\tchr1\t1\t60\t8M\t*\t0\t0\tAAGTACGT\t#IIIIIII\tRG:Z:rg1\n"),
        &[],
    );
    dump.run(
        "reverse",
        &format!("{u}r3\t4\t*\t0\t0\t*\t*\t0\t0\tGTACGT\tABCDEF\tRG:Z:rg1\n"),
        &format!("{a}r3\t16\tchr1\t1\t60\t6M\t*\t0\t0\tACGTAC\tFEDCBA\tRG:Z:rg1\n"),
        &[],
    );
    dump.run(
        "multi_read",
        &format!(
            "{u}a\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             b\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             c\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        ),
        &format!(
            "{a}a\t0\tchr1\t20\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             b\t0\tchr1\t5\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             c\t0\tchr1\t30\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        ),
        &[],
    );
    dump.run("proper_pair", &unmapped_p1(), &aligned_p1(), &[]);
    dump.check("merge_bam_alignment");
}

fn unmapped_p1() -> String {
    format!(
        "{UNMAPPED_HEADER}p1\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
         p1\t141\t*\t0\t0\t*\t*\t0\t0\tTTTTGGGG\tIIIIIIII\tRG:Z:rg1\n"
    )
}

fn aligned_p1() -> String {
    format!(
        "{ALIGNED_HEADER_40}p1\t99\tchr1\t1\t60\t8M\t=\t20\t27\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
         p1\t147\tchr1\t20\t60\t8M\t=\t1\t-27\tCCCCAAAA\tIIIIIIII\tRG:Z:rg1\n"
    )
}

/// `MergeBamAlignmentFullDump`: the whole output file, coordinate sorted.
#[test]
fn the_whole_file_is_picards() {
    let (u, a) = (UNMAPPED_HEADER, ALIGNED_HEADER_40);
    let mut dump = Dump::new("full", Shape::Full, vec![("chr1", 40)]);
    dump.run(
        "unpaired",
        &format!(
            "{u}a\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tab:Z:keepme\n\
             b\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        ),
        &format!(
            "{a}a\t0\tchr1\t20\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             b\t0\tchr1\t5\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        ),
        &[],
    );
    dump.run("paired", &unmapped_p1(), &aligned_p1(), &[]);
    dump.run(
        "unmapped_passthrough",
        &format!(
            "{u}a\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             z\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGCCCC\tJJJJJJJJ\tRG:Z:rg1\n"
        ),
        &format!(
            "{a}a\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             z\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGCCCC\tJJJJJJJJ\tRG:Z:rg1\n"
        ),
        &[],
    );
    dump.check("merge_bam_alignment_full");
}

/// `MergeBamAlignmentHeaderDump`: the header built from three files, and its refusals.
#[test]
fn the_header_is_picards() {
    let (u, a) = (UNMAPPED_HEADER, ALIGNED_HEADER_40);
    let unmapped_read = "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";
    let aligned_read = "r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";
    let unmapped = format!("{u}{unmapped_read}");
    let aligned = format!("{a}{aligned_read}");
    let mut dump = Dump::new("header", Shape::Header, vec![("chr1", 40)]);

    dump.run("the-whole-header", &unmapped, &aligned, &[]);
    dump.run_from(
        "a-reference-named-the-long-way-round",
        &unmapped,
        &[&aligned],
        true,
        &[],
    );
    dump.run(
        "a-read-group-the-aligner-rewrote",
        &unmapped,
        &format!(
            "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n\
             @RG\tID:rg1\tSM:other\tLB:otherlib\tPL:ILLUMINA\n\
             @PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n{aligned_read}"
        ),
        &[],
    );
    dump.run(
        "no-read-group-at-all",
        "@HD\tVN:1.6\tSO:queryname\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n\
         @PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n\
         r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        &[],
    );
    dump.run(
        "a-chain-of-programs",
        &unmapped,
        &format!(
            "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n\
             @RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n\
             @PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n\
             @PG\tID:samtools\tPN:samtools\tVN:1.9\tPP:bwa\tCL:samtools view\n{aligned_read}"
        ),
        &[],
    );
    dump.run_from(
        "two-aligned-files-with-the-same-program-id",
        &unmapped,
        &[
            &aligned,
            "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n\
             @RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n\
             @PG\tID:bwa\tPN:bwa\tVN:2.0\tCL:bwa mem -a\n\
             r2\t0\tchr1\t9\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n",
        ],
        false,
        &[],
    );
    dump.run(
        "a-program-from-the-command-line",
        &unmapped,
        &aligned,
        &[
            "PROGRAM_RECORD_ID=mine",
            "PROGRAM_GROUP_NAME=miner",
            "PROGRAM_GROUP_VERSION=3.0",
            "PROGRAM_GROUP_COMMAND_LINE=mine --do-it",
        ],
    );
    dump.run(
        "half-a-program-from-the-command-line",
        &unmapped,
        &aligned,
        &["PROGRAM_RECORD_ID=mine"],
    );
    dump.run(
        "comments-in-the-inputs",
        &format!(
            "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n\
             @CO\tfrom the unmapped\n{unmapped_read}"
        ),
        &format!(
            "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n\
             @RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n\
             @CO\tfrom the aligned\n{aligned_read}"
        ),
        &[],
    );
    dump.run(
        "a-queryname-output",
        &unmapped,
        &aligned,
        &["SORT_ORDER=queryname"],
    );
    dump.run(
        "an-unsorted-output",
        &unmapped,
        &aligned,
        &["SORT_ORDER=unsorted"],
    );
    dump.run(
        "an-aligned-header-that-disagrees",
        &unmapped,
        &format!(
            "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:41\n\
             @RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n\
             {aligned_read}"
        ),
        &[],
    );
    dump.run(
        "an-aligned-header-with-another-contig",
        &unmapped,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr2\tLN:40\n\
         @RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n\
         r1\t0\tchr2\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n",
        &[],
    );
    dump.check("merge_bam_alignment_header");
}

/// `MergeBamAlignmentPairDump`: what the merger does to a pair whose mate fields are wrong.
#[test]
fn a_pair_is_merged_as_picard_merges_it() {
    let u = UNMAPPED_HEADER;
    let a = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n@SQ\tSN:chr2\tLN:40\n\
             @RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";
    let unmapped = format!(
        "{u}p\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
         p\t141\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
    );
    let first = |start: i32, cigar: &str| {
        format!("p\t65\tchr1\t{start}\t60\t{cigar}\tchr2\t39\t999\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n")
    };
    let second = |start: i32, cigar: &str| {
        format!(
            "p\t129\tchr1\t{start}\t60\t{cigar}\tchr2\t39\t-999\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        )
    };
    let apart = format!("{a}{}{}", first(1, "8M"), second(25, "8M"));
    let mut dump = Dump::new("pair", Shape::Outcome, vec![("chr1", 40), ("chr2", 40)]);

    dump.run("mate-fields-from-the-other-end", &unmapped, &apart, &[]);
    dump.run(
        "the-second-end-first",
        &unmapped,
        &format!("{a}{}{}", first(25, "8M"), second(1, "8M")),
        &[],
    );
    let across = format!(
        "{a}p\t67\tchr1\t1\t60\t8M\tchr2\t1\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
         p\t131\tchr2\t1\t60\t8M\tchr1\t1\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
    );
    dump.run("across-two-contigs", &unmapped, &across, &[]);
    dump.run(
        "across-two-contigs-with-the-aligners-flags",
        &unmapped,
        &across,
        &["ALIGNER_PROPER_PAIR_FLAGS=true"],
    );
    dump.run(
        "without-the-mate-cigar",
        &unmapped,
        &apart,
        &["ADD_MATE_CIGAR=false"],
    );
    let overlapping = format!(
        "{a}p\t65\tchr1\t5\t60\t8M\tchr1\t3\t-10\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
         p\t145\tchr1\t3\t60\t8M\tchr1\t5\t10\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
    );
    dump.run("overlapping-ends", &unmapped, &overlapping, &[]);
    dump.run(
        "overlapping-ends-unclipped",
        &unmapped,
        &overlapping,
        &["CLIP_OVERLAPPING_READS=false"],
    );
    dump.run(
        "overlapping-ends-hard-clipped",
        &unmapped,
        &overlapping,
        &["HARD_CLIP_OVERLAPPING_READS=true"],
    );
    dump.run(
        "one-end-unaligned",
        &unmapped,
        &format!(
            "{a}p\t73\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n\
             p\t141\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        ),
        &[],
    );
    dump.run(
        "a-soft-clipped-end",
        &unmapped,
        &format!("{a}{}{}", first(1, "2S6M"), second(25, "8M")),
        &[],
    );
    dump.run(
        "the-deprecated-paired-run",
        &unmapped,
        &apart,
        &["PAIRED_RUN=false"],
    );
    dump.check("merge_bam_alignment_pair");
}

/// `MergeBamAlignmentClipDump`: the clips, and the contamination filter.
#[test]
fn the_clips_are_picards() {
    let u = UNMAPPED_HEADER;
    let a = aligned_header_200();
    let bases = &acgt(200)[40..80];
    let q = quals(40);
    let plain = format!("{u}r\t4\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{q}\tRG:Z:rg1\n");
    let aligned_at = |start: i32, cigar: &str| {
        format!("{a}r\t0\tchr1\t{start}\t60\t{cigar}\t*\t0\t0\t{bases}\t{q}\tRG:Z:rg1\n")
    };
    let mut dump = Dump::new("clip", Shape::Outcome, vec![("chr1", 200)]);

    dump.run(
        "off-the-end-of-the-contig",
        &plain,
        &aligned_at(171, "40M"),
        &["VALIDATION_STRINGENCY=SILENT"],
    );
    dump.run(
        "off-the-end-with-a-soft-clip",
        &plain,
        &aligned_at(171, "35M5S"),
        &["VALIDATION_STRINGENCY=SILENT"],
    );
    dump.run("inside-the-contig", &plain, &aligned_at(41, "40M"), &[]);
    let marked = format!("{u}r\t4\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{q}\tRG:Z:rg1\tXT:i:31\n");
    dump.run(
        "an-adapter-marked-in-the-unmapped-bam",
        &marked,
        &aligned_at(41, "40M"),
        &[],
    );
    dump.run(
        "an-adapter-left-alone",
        &marked,
        &aligned_at(41, "40M"),
        &["CLIP_ADAPTERS=false"],
    );
    let short = aligned_at(41, "15S10M15S");
    dump.run("a-short-alignment-kept", &plain, &short, &[]);
    dump.run(
        "clipped-at-one-end-only",
        &plain,
        &aligned_at(41, "10M30S"),
        &["UNMAP_CONTAMINANT_READS=true"],
    );
    dump.run(
        "a-short-alignment-unmapped",
        &plain,
        &short,
        &["UNMAP_CONTAMINANT_READS=true"],
    );
    dump.run(
        "a-short-alignment-under-a-lower-threshold",
        &plain,
        &short,
        &["UNMAP_CONTAMINANT_READS=true", "MIN_UNCLIPPED_BASES=8"],
    );
    for strategy in ["MOVE_TO_TAG", "COPY_TO_TAG", "DO_NOT_CHANGE_INVALID"] {
        dump.run(
            &format!("a-contaminant-{}", strategy.to_lowercase()),
            &plain,
            &short,
            &[
                "UNMAP_CONTAMINANT_READS=true",
                &format!("UNMAPPED_READ_STRATEGY={strategy}"),
            ],
        );
    }
    dump.check("merge_bam_alignment_clip");
}

/// `MergeBamAlignmentHitsDump`: which alignment each strategy calls primary.
#[test]
fn the_primary_is_picards() {
    let u = UNMAPPED_HEADER;
    let a = aligned_header_200();
    let bases = &acgt(200)[40..60];
    let q = quals(20);
    let unmapped = format!("{u}r\t4\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{q}\tRG:Z:rg1\n");
    let hit = |flag: u16, start: i32, mapq: i32, cigar: &str| {
        format!("r\t{flag}\tchr1\t{start}\t{mapq}\t{cigar}\t*\t0\t0\t{bases}\t{q}\tRG:Z:rg1\n")
    };
    let three = format!(
        "{a}{}{}{}",
        hit(256, 121, 60, "5S15M"),
        hit(256, 61, 30, "20M"),
        hit(256, 41, 10, "10S10M")
    );
    let mut dump = Dump::new("hits", Shape::Outcome, vec![("chr1", 200)]);

    dump.run("three-hits-best-mapq", &unmapped, &three, &[]);
    for (case, strategy) in [
        ("three-hits-earliest-fragment", "EarliestFragment"),
        ("three-hits-best-end-mapq", "BestEndMapq"),
        ("three-hits-most-distant", "MostDistant"),
    ] {
        dump.run(
            case,
            &unmapped,
            &three,
            &[&format!("PRIMARY_ALIGNMENT_STRATEGY={strategy}")],
        );
    }
    dump.run(
        "three-hits-without-the-secondaries",
        &unmapped,
        &three,
        &["INCLUDE_SECONDARY_ALIGNMENTS=false"],
    );
    let named = format!(
        "{a}{}{}{}",
        hit(256, 121, 60, "5S15M"),
        hit(0, 61, 30, "20M"),
        hit(256, 41, 10, "10S10M")
    );
    dump.run("the-aligners-own-primary", &unmapped, &named, &[]);
    dump.run(
        "the-aligners-own-primary-under-another-strategy",
        &unmapped,
        &named,
        &["PRIMARY_ALIGNMENT_STRATEGY=EarliestFragment"],
    );
    let unmapped_pair = format!(
        "{u}p\t77\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{q}\tRG:Z:rg1\n\
         p\t141\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{q}\tRG:Z:rg1\n"
    );
    let paired = format!(
        "{a}p\t321\tchr1\t41\t60\t20M\tchr1\t81\t60\t{bases}\t{q}\tRG:Z:rg1\tHI:i:1\n\
         p\t385\tchr1\t81\t60\t20M\tchr1\t41\t-60\t{bases}\t{q}\tRG:Z:rg1\tHI:i:1\n\
         p\t321\tchr1\t41\t30\t20M\tchr1\t161\t140\t{bases}\t{q}\tRG:Z:rg1\tHI:i:2\n\
         p\t385\tchr1\t161\t30\t20M\tchr1\t41\t-140\t{bases}\t{q}\tRG:Z:rg1\tHI:i:2\n"
    );
    dump.run("a-pair-with-two-hits", &unmapped_pair, &paired, &[]);
    dump.run(
        "a-pair-with-two-hits-most-distant",
        &unmapped_pair,
        &paired,
        &["PRIMARY_ALIGNMENT_STRATEGY=MostDistant"],
    );
    dump.run(
        "a-pair-under-earliest-fragment",
        &unmapped_pair,
        &paired,
        &["PRIMARY_ALIGNMENT_STRATEGY=EarliestFragment"],
    );
    dump.check("merge_bam_alignment_hits");
}
