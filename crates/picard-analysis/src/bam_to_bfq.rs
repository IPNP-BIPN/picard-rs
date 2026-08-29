//! `BamToBfq`: Maq's binary fastq, one byte per base.
//!
//! A record is a name, a length and one byte per base, the base in the top two bits and the
//! quality in the bottom six. Everything is fixed-width and little-endian, so the payload is
//! byte-comparable; the file on disk is GZIP, which `IOUtil.openFileForWriting` applies to a
//! `.bfq` the way it applies it to a `.gz`, and reproducing THOSE bytes means reproducing a
//! DEFLATE stream, which this port does not do.
//!
//! Ported from `picard.fastq.BamToBfq` and `picard.fastq.SamToBfqWriter`.

/// `SEED_REGION_LENGTH`, inside which an uncalled base is treated differently.
pub const SEED_REGION_LENGTH: usize = 28;
/// `MAX_SEED_REGION_NOCALL_FIXES`, after which an uncalled base in the seed region loses even the
/// quality of one.
pub const MAX_SEED_REGION_NOCALL_FIXES: usize = 2;

/// The two-bit code of a base, which is what the top of the byte carries.
///
/// `A` is zero and so is an uncalled base, which is why an `N` is written as an `A` and told apart
/// only by its quality.
pub fn base_code(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        b'N' | b'.' => Some(0),
        _ => None,
    }
}

/// Whether a base is one the encoder has no code for.
///
/// The reference throws `Unknown base when writing bfq file` here, and that branch cannot be
/// reached through a BAM at all: htsjdk refuses such a record at write time, so the file the tool
/// would have read cannot be built.
pub fn is_unknown_base(base: u8) -> bool {
    base_code(base).is_none()
}

/// `encodeBaseAndQuality`, which is the whole of the format's payload.
pub fn encode_base_and_quality(base: u8, quality: u8) -> u8 {
    (base << 6) | quality
}

/// The quality an uncalled base is written with, which is not the read's own.
///
/// Inside the seed region the first two are written at one and the rest at zero; outside it every
/// one is written at one. The counter is per READ and not per file.
pub fn no_call_quality(index: usize, fixes_so_far: usize) -> (u8, bool) {
    if index < SEED_REGION_LENGTH {
        if fixes_so_far < MAX_SEED_REGION_NOCALL_FIXES {
            (1, true)
        } else {
            (0, false)
        }
    } else {
        (1, false)
    }
}

/// One read's bases and qualities, encoded.
///
/// `retained_length` is where a clipped adapter starts: the tail is rewritten as `A` at quality
/// one rather than dropped, so the record's length does not move. `--BASES_TO_WRITE` is a
/// different thing entirely and shortens the arrays before this is called.
pub fn encode(bases: &[u8], qualities: &[u8], retained_length: usize) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(bases.len());
    let mut fixes = 0;
    for (index, base) in bases.iter().enumerate() {
        let code = base_code(*base)?;
        let mut quality = qualities[index];
        if matches!(base.to_ascii_uppercase(), b'N' | b'.') {
            let (written, counted) = no_call_quality(index, fixes);
            quality = written;
            if counted {
                fixes += 1;
            }
        }
        out.push(encode_base_and_quality(code, quality));
    }
    for byte in out.iter_mut().skip(retained_length) {
        *byte = encode_base_and_quality(0, 1);
    }
    Some(out)
}

/// The name of one output file: `<prefix><index>.<read>.bfq`.
///
/// The index counts chunks from ZERO, and the read is 1 or 2, so a single-end run writes half as
/// many files as a paired one.
pub fn output_file_name(prefix: &str, index: usize, read: usize) -> String {
    format!("{prefix}{index}.{read}.bfq")
}

/// `customCommandLineValidation`: the prefix a flowcell and a lane build when no prefix is given.
pub fn output_file_prefix(
    explicit: Option<&str>,
    flowcell_barcode: Option<&str>,
    lane: Option<i32>,
) -> String {
    match explicit {
        Some(prefix) => prefix.to_string(),
        None => format!(
            "{}.{}",
            flowcell_barcode.unwrap_or("null"),
            lane.map(|value| value.to_string())
                .unwrap_or_else(|| "null".to_string())
        ),
    }
}

/// The message the parser writes when both ways of naming the output are given.
pub const MUTUALLY_EXCLUSIVE_MESSAGE: &str =
    "ERROR: Option 'FLOWCELL_BARCODE' cannot be used in conjunction with option(s) OUTPUT_FILE_PREFIX";

/// Whether the bytes on disk are a claim this port makes, which they are not.
///
/// The file is GZIP, so its bytes are a DEFLATE stream. The payload beside it in the golden is
/// what this port reproduces, and the compressed form waits on the deflater Milestone X.4 is for.
pub const FILE_BYTES_ARE_REPRODUCIBLE: bool = false;
