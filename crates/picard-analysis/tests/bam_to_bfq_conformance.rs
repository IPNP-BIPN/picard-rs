//! Conformance for `BamToBfq` against Picard 3.4.0.
//!
//! Golden from `tools/bamtobfq-conformance/BamToBfqDump.java`, nineteen runs whose outputs are in
//! the golden twice over: as the bytes on disk, which are GZIP, and as the payload, which is what
//! this port reproduces.
//!
//! # What this suite is for
//!
//!  * **the encoding: one byte, two bits of base and six of quality**;
//!  * **an uncalled base, whose quality is not the read's**;
//!  * **the file names, and how many of them a run writes**;
//!  * **the two length arguments, which do different things to the record**;
//!  * **and the branch a BAM cannot reach.**

use std::io::Read;

use picard_analysis::bam_to_bfq::{
    base_code, encode, encode_base_and_quality, is_unknown_base, no_call_quality, output_file_name,
    output_file_prefix, FILE_BYTES_ARE_REPRODUCIBLE, MUTUALLY_EXCLUSIVE_MESSAGE,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bam_to_bfq.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| {
            line[prefix.len()..]
                .replace("\\t", "\t")
                .replace("\\n", "\n")
        })
}

fn files(text: &str, case: &str) -> Vec<String> {
    field(text, "files", case)
        .unwrap_or_else(|| panic!("files/{case}"))
        .split(' ')
        .map(str::to_string)
        .collect()
}

/// One file's payload, as bytes.
fn payload(text: &str, case: &str, file: &str) -> Vec<u8> {
    let hex = field(text, "plain", &format!("{case}/{file}"))
        .unwrap_or_else(|| panic!("plain/{case}/{file}"));
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).expect("a byte"))
        .collect()
}

/// One record as the format writes it: the name, the length, and the encoded bases.
fn record(name: &str, encoded: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    // The name is written with its length, NUL included, little-endian like everything else.
    out.extend((name.len() as u32 + 1).to_le_bytes());
    out.extend(name.as_bytes());
    out.push(0);
    out.extend((encoded.len() as u32).to_le_bytes());
    out.extend(encoded);
    out
}

/// One byte carries the base and the quality, and the port writes the golden's bytes.
#[test]
fn one_byte_carries_the_base_and_the_quality() {
    let text = corpus();
    assert_eq!(base_code(b'A'), Some(0));
    assert_eq!(base_code(b'C'), Some(1));
    assert_eq!(base_code(b'G'), Some(2));
    assert_eq!(base_code(b'T'), Some(3));
    // Lower case is the same code, which is what the `every-base` case shows.
    assert_eq!(base_code(b'a'), base_code(b'A'));
    assert_eq!(encode_base_and_quality(1, 40), 0x68);
    // `ACGT` at quality 40 is the first file of the plain pair, name and length included.
    let encoded = encode(b"ACGT", &[40; 4], 4).expect("the encoding");
    assert_eq!(encoded, vec![0x28, 0x68, 0xa8, 0xe8]);
    assert_eq!(
        payload(&text, "a-pair", "run.0.1.bfq"),
        record("read1/1", &encoded)
    );
    // And the second read of the pair is the reverse of those four bases, named `/2`.
    let second = encode(b"TGCA", &[40; 4], 4).expect("the encoding");
    assert_eq!(
        payload(&text, "a-pair", "run.0.2.bfq"),
        record("read1/2", &second)
    );
    // Every base in both cases is eight bytes, the lower-case half repeating the upper.
    let both = encode(b"ACGTacgt", &[40; 8], 8).expect("the encoding");
    assert_eq!(&both[..4], &both[4..]);
    assert_eq!(
        payload(&text, "every-base", "run.0.1.bfq"),
        record("read1/1", &both)
    );
}

/// An uncalled base is an `A` whose quality is not the read's.
#[test]
fn an_uncalled_base_carries_its_own_quality() {
    let text = corpus();
    assert_eq!(base_code(b'N'), base_code(b'A'));
    // Inside the seed region, the first two are written at one and the rest at zero.
    assert_eq!(no_call_quality(0, 0), (1, true));
    assert_eq!(no_call_quality(0, 2), (0, false));
    // Outside it, every one is written at one and none of them counts.
    assert_eq!(no_call_quality(30, 0), (1, false));
    let inside = encode(b"ANGT", &[40; 4], 4).expect("the encoding");
    assert_eq!(inside[1], encode_base_and_quality(0, 1));
    assert_eq!(
        payload(&text, "an-n-inside-the-seed", "run.0.1.bfq"),
        record("read1/1", &inside)
    );
    let mut bases = b"ACGT".repeat(10);
    bases[39] = b'N';
    let outside = encode(&bases, &[40; 40], 40).expect("the encoding");
    assert_eq!(outside[39], encode_base_and_quality(0, 1));
    assert_eq!(
        payload(&text, "an-n-outside-the-seed", "run.0.1.bfq"),
        record("read1/1", &outside)
    );
}

/// The file names, and how many of them a run writes.
#[test]
fn the_files_are_named_from_the_prefix() {
    let text = corpus();
    assert_eq!(output_file_name("run.", 0, 1), "run.0.1.bfq");
    assert_eq!(files(&text, "a-pair"), ["run.0.1.bfq", "run.0.2.bfq"]);
    // A single-end run writes one file per chunk rather than two.
    assert_eq!(files(&text, "single-end"), ["run.0.1.bfq"]);
    // A chunk size splits the output, and the run opens one more chunk than the reads need.
    assert_eq!(files(&text, "two-pairs").len(), 2);
    assert_eq!(files(&text, "two-pairs-chunked").len(), 6);
    // A flowcell and a lane build the prefix when no explicit one is given.
    assert_eq!(
        output_file_prefix(None, Some("30PYMAAXX"), Some(3)),
        "30PYMAAXX.3"
    );
    assert_eq!(
        files(&text, "a-flowcell-and-a-lane"),
        ["30PYMAAXX.3.0.1.bfq", "30PYMAAXX.3.0.2.bfq"]
    );
    assert_eq!(
        output_file_prefix(Some("run"), Some("30PYMAAXX"), Some(3)),
        "run"
    );
    // And the two ways of naming it may not be given together.
    let refusal = field(&text, "error", "both-ways-of-naming").expect("the refusal");
    assert!(refusal.contains(MUTUALLY_EXCLUSIVE_MESSAGE), "{refusal}");
}

/// The two length arguments do different things to the record.
#[test]
fn trimming_and_clipping_are_not_the_same() {
    let text = corpus();
    // `BASES_TO_WRITE` shortens the record: four bytes where the read had eight.
    let trimmed = payload(&text, "bases-to-write", "run.0.1.bfq");
    assert_eq!(
        trimmed,
        record("read1/1", &encode(b"ACGT", &[40; 4], 4).expect("it"))
    );
    // `CLIP_ADAPTERS` rewrites the tail as `A` at quality one and leaves the length alone, which
    // is why the clipped file is as long as the unclipped one.
    let clipped = payload(&text, "clip-adapters", "run.0.2.bfq");
    let untouched = payload(&text, "every-base", "run.0.2.bfq");
    assert_eq!(clipped.len(), untouched.len());
    // Nothing is clipped on this fixture, so the bytes are the read's own; what the case pins is
    // that the length did not move.
    let full = encode(b"TGCATGCA", &[40; 8], 8).expect("the encoding");
    assert_eq!(clipped, record("read1/2", &full));
    // And the rewrite itself is the port's, exercised where a clip does land.
    let rewritten = encode(b"ACGTACGT", &[40; 8], 4).expect("the encoding");
    assert_eq!(&rewritten[4..], &[encode_base_and_quality(0, 1); 4]);
}

/// The branch a BAM cannot reach, and the bytes this port does not claim.
#[test]
fn the_unknown_base_branch_cannot_be_reached_through_a_bam() {
    let text = corpus();
    assert!(is_unknown_base(b'X'));
    assert!(encode(b"ACXT", &[40; 4], 4).is_none());
    // The fixture could not even be built: htsjdk refuses the record at write time, so the tool
    // never sees a base it has no code for.
    let refusal = field(&text, "error", "an-unknown-base").expect("the refusal");
    assert!(
        refusal.starts_with("fixture java.lang.IllegalStateException"),
        "{refusal}"
    );
    assert!(
        refusal.contains("Bad base passed to charToCompressedBaseHigh"),
        "{refusal}"
    );
    // And the bytes on disk are GZIP, which this port does not reproduce: the golden holds both
    // forms, and the constant says which of the two is a claim.
    let raw = field(&text, "bfq", "a-pair/run.0.1.bfq").expect("the file");
    assert!(raw.starts_with("1f8b08"), "{raw}");
    assert_eq!(
        FILE_BYTES_ARE_REPRODUCIBLE,
        raw == field(&text, "plain", "a-pair/run.0.1.bfq").expect("the payload")
    );
}
