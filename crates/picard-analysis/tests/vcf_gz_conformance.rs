//! Conformance for every VCF writer asked for a `.vcf.gz`, against Picard 3.5.0.
//!
//! Golden from `tools/vcfgz-conformance/`: twenty-nine runs over eleven tools, each output's MD5
//! and size, and every `.tbi` inflated to hex. This test repeats `port_side.py`: it writes the
//! inputs by the dump's formula, runs each binary under the dump's arguments, and compares the
//! rows to the golden one by one.
//!
//! # What this suite is for
//!
//!  * **the BGZF bytes at every installed level and deflater**, GKL and libdeflate;
//!  * **the last tabix chunk ending at the flushed block**, offset zero, as `close()` leaves it;
//!  * **an index without a dictionary**, which the tabix creator does not need;
//!  * **the two refusals only the tabix creator makes**: a record out of order and a contig the
//!    dictionary lacks;
//!  * **and GatherVcfs copying blocks** when its inputs are already compressed.

use std::io::Read;
use std::path::Path;
use std::process::Command;

use md5::{Digest, Md5};

/// `VcfGzDump.DIR`. LiftoverVcf writes R's absolute path into its `##reference` line, so this
/// directory is part of the golden's bytes.
const DIR: &str = "/tmp/picard-rs-vcfgz";
const CONTIGS: [&str; 3] = ["chr1", "chr2", "chr3"];
const LENGTH: usize = 60000;
const PER_CONTIG: usize = 4000;
const STEP: usize = 13;
const BASES: &[u8; 4] = b"ACGT";

fn golden() -> Vec<String> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcfgz.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

fn emit(rows: &mut Vec<String>, kind: &str, name: &str, payload: &str) {
    let payload = payload
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n");
    rows.push(format!("{kind}\t{name}\t{payload}"));
}

fn md5(data: &[u8]) -> String {
    Md5::digest(data)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn base(contig: usize, position: usize) -> char {
    BASES[(position * 31 + position / 7 + contig * 5) % 4] as char
}

fn next(b: char) -> char {
    let i = BASES.iter().position(|&x| x as char == b).expect("a base");
    BASES[(i + 1) % 4] as char
}

fn record(i: usize, contig: usize, position: usize, samples: usize) -> String {
    let indel = i.is_multiple_of(5);
    let first = if i.is_multiple_of(17) {
        next(base(contig, position))
    } else {
        base(contig, position)
    };
    let (reference, alternate) = if indel {
        (
            format!("{first}{}", base(contig, position + 1)),
            first.to_string(),
        )
    } else {
        (first.to_string(), next(first).to_string())
    };
    let id = if i.is_multiple_of(3) {
        format!("rs{i}")
    } else {
        ".".to_string()
    };
    let filter = if i.is_multiple_of(11) {
        "LowQual"
    } else {
        "PASS"
    };
    let mut line = format!(
        "{}\t{position}\t{id}\t{reference}\t{alternate}\t{}.{}\t{filter}\t\
         DP={};AC=1;AF=0.167;AN=6\tGT:AD:DP:GQ:PL",
        CONTIGS[contig],
        i % 997,
        i % 10,
        i * 7 % 200
    );
    for s in 0..samples {
        let g = (i + s) % 4;
        line += "\t";
        line += ["0/0", "0/1", "1/1", "./."][g];
        if g != 3 {
            let (ref_depth, alt_depth) = (10 + i % 30, i % 20);
            line += &format!(
                ":{ref_depth},{alt_depth}:{}:{}:0,{},{}",
                ref_depth + alt_depth,
                i * 13 % 99,
                i % 50,
                i % 90 + 10
            );
        }
    }
    line + "\n"
}

fn header(contig_lines: &[&str], samples: usize) -> String {
    let mut h = String::from("##fileformat=VCFv4.2\n");
    for c in contig_lines {
        h += &format!("##contig=<ID={c},length={LENGTH}>\n");
    }
    h += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
          ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"AC\">\n\
          ##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
          ##INFO=<ID=AN,Number=1,Type=Integer,Description=\"AN\">\n\
          ##FILTER=<ID=LowQual,Description=\"Low quality\">\n\
          ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
          ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
          ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
          ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n\
          ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for name in ["child", "father", "mother"].iter().take(samples) {
        h += "\t";
        h += name;
    }
    h + "\n"
}

/// The records of `contigs` in order; `parity` keeps every other one, `swap` exchanges a record
/// with the next so the file is out of order.
fn body(contigs: &[usize], samples: usize, parity: Option<usize>, swap: Option<usize>) -> String {
    let mut lines = Vec::new();
    let mut i = 0;
    for &c in contigs {
        for k in 0..PER_CONTIG {
            if parity.is_none_or(|p| i % 2 == p) {
                lines.push(record(i, c, 11 + k * STEP, samples));
            }
            i += 1;
        }
    }
    if let Some(s) = swap {
        lines.swap(s, s + 1);
    }
    lines.concat()
}

fn inputs(d: &Path, rows: &mut Vec<String>) {
    let mut put = |name: &str, text: String| {
        std::fs::write(d.join(name), &text).expect("an input");
        let payload = format!("{}\t{}", md5(text.as_bytes()), text.len());
        emit(rows, "input", name, &payload);
    };
    let order = [0, 1, 2];
    put("in.vcf", header(&CONTIGS, 3) + &body(&order, 3, None, None));
    put(
        "unsorted.vcf",
        header(&CONTIGS, 3) + &body(&[2, 0, 1], 3, None, None),
    );
    put(
        "swapped.vcf",
        header(&CONTIGS, 3) + &body(&order, 3, None, Some(100)),
    );
    put("nodict.vcf", header(&[], 3) + &body(&order, 3, None, None));
    put(
        "missing.vcf",
        header(&CONTIGS[..2], 3) + &body(&order, 3, None, None),
    );
    put(
        "single.vcf",
        header(&CONTIGS, 1) + &body(&order, 1, None, None),
    );
    put(
        "single-nodict.vcf",
        header(&[], 1) + &body(&order, 1, None, None),
    );
    put(
        "m1.vcf",
        header(&CONTIGS, 3) + &body(&order, 3, Some(0), None),
    );
    put(
        "m2.vcf",
        header(&CONTIGS, 3) + &body(&order, 3, Some(1), None),
    );
    for c in 0..3 {
        let name = format!("p{}.vcf", c + 1);
        put(&name, header(&CONTIGS, 3) + &body(&[c], 3, None, None));
    }
    let mut fasta = String::new();
    let mut dictionary = String::from("@HD\tVN:1.6\n");
    let mut other = String::from("@HD\tVN:1.6\n");
    let mut chain = String::new();
    for (c, name) in CONTIGS.iter().enumerate() {
        fasta += &format!(">{name}\n");
        let sequence: String = (1..=LENGTH).map(|p| base(c, p)).collect();
        for line in sequence.as_bytes().chunks(60) {
            fasta += std::str::from_utf8(line).expect("ASCII");
            fasta += "\n";
        }
        dictionary += &format!("@SQ\tSN:{name}\tLN:{LENGTH}\n");
        other += &format!("@SQ\tSN:{name}\tLN:{}\n", 2 * LENGTH);
        chain += &format!(
            "chain 1000 {name} {LENGTH} + 0 {LENGTH} {name} {LENGTH} + 0 {LENGTH} {}\n{LENGTH}\n\n",
            c + 1
        );
    }
    put("ref.fasta", fasta);
    put("ref.dict", dictionary);
    put("other.dict", other);
    put("id.chain", chain);
}

/// The binary Cargo built for `tool`; every tool here is a binary of this crate.
fn binary(tool: &str) -> &'static str {
    match tool {
        "SortVcf" => env!("CARGO_BIN_EXE_SortVcf"),
        "MakeSitesOnlyVcf" => env!("CARGO_BIN_EXE_MakeSitesOnlyVcf"),
        "RenameSampleInVcf" => env!("CARGO_BIN_EXE_RenameSampleInVcf"),
        "UpdateVcfSequenceDictionary" => env!("CARGO_BIN_EXE_UpdateVcfSequenceDictionary"),
        "FilterVcf" => env!("CARGO_BIN_EXE_FilterVcf"),
        "FixVcfHeader" => env!("CARGO_BIN_EXE_FixVcfHeader"),
        "MergeVcfs" => env!("CARGO_BIN_EXE_MergeVcfs"),
        "SplitVcfs" => env!("CARGO_BIN_EXE_SplitVcfs"),
        "LiftoverVcf" => env!("CARGO_BIN_EXE_LiftoverVcf"),
        "VcfFormatConverter" => env!("CARGO_BIN_EXE_VcfFormatConverter"),
        "GatherVcfs" => env!("CARGO_BIN_EXE_GatherVcfs"),
        _ => panic!("no binary for {tool}"),
    }
}

/// One run in its own directory: `@` in an argument is the inputs' directory, `%` this case's.
/// The exit row is the uncaught exception's `Class: message` when there is one, as the dump's
/// catch reports it, and the exit status otherwise.
fn run(d: &Path, rows: &mut Vec<String>, name: &str, tool: &str, args: &[&str]) {
    let out = d.join(format!("case-{name}"));
    std::fs::create_dir(&out).expect("the case directory");
    let argv: Vec<String> = args
        .iter()
        .map(|a| {
            a.replace('@', d.to_str().expect("UTF-8"))
                .replace('%', out.to_str().expect("UTF-8"))
        })
        .collect();
    let done = Command::new(binary(tool))
        .args(&argv)
        .current_dir(&out)
        .output()
        .expect("the binary runs");
    let stderr = String::from_utf8_lossy(&done.stderr);
    let exception = stderr
        .lines()
        .find_map(|line| line.strip_prefix("Exception in thread \"main\" "));
    let status = match exception {
        Some(exception) => exception.to_string(),
        None => done.status.code().expect("an exit code").to_string(),
    };
    emit(rows, "exit", name, &status);
    let mut files: Vec<_> = std::fs::read_dir(&out)
        .expect("the case directory")
        .map(|entry| entry.expect("an entry").path())
        .collect();
    files.sort();
    for f in files {
        let file_name = f.file_name().expect("a name").to_str().expect("UTF-8");
        let key = format!("{name}/{file_name}");
        if file_name.ends_with(".idx") {
            emit(rows, "file", &key, "present");
            continue;
        }
        let data = std::fs::read(&f).expect("an output");
        emit(
            rows,
            "file",
            &key,
            &format!("{}\t{}", md5(&data), data.len()),
        );
        if file_name.ends_with(".tbi") {
            let mut inflated = Vec::new();
            flate2::read::MultiGzDecoder::new(&data[..])
                .read_to_end(&mut inflated)
                .expect("the index inflates");
            let hex: String = inflated.iter().map(|b| format!("{b:02x}")).collect();
            emit(rows, "tbi", &key, &hex);
        }
    }
}

fn rows() -> Vec<String> {
    let d = Path::new(DIR);
    // A leftover from an earlier run, or from `port_side.py`, which uses the same directory.
    let _ = std::fs::remove_dir_all(d);
    std::fs::create_dir_all(d).expect("the dump's directory");
    let mut rows = Vec::new();
    inputs(d, &mut rows);
    let r = &mut rows;
    for (name, extra) in [
        ("sort", &[][..]),
        ("sort-no-index", &["CREATE_INDEX=false"][..]),
        ("sort-level-3", &["COMPRESSION_LEVEL=3"][..]),
        ("sort-level-9", &["COMPRESSION_LEVEL=9"][..]),
        ("sort-jdk-deflater", &["USE_JDK_DEFLATER=true"][..]),
        (
            "sort-jdk-deflater-level-1",
            &["USE_JDK_DEFLATER=true", "COMPRESSION_LEVEL=1"][..],
        ),
    ] {
        let args: Vec<&str> = ["I=@/unsorted.vcf", "O=%/out.vcf.gz"]
            .into_iter()
            .chain(extra.iter().copied())
            .collect();
        run(d, r, name, "SortVcf", &args);
    }
    run(
        d,
        r,
        "sort-bgz",
        "SortVcf",
        &["I=@/unsorted.vcf", "O=%/out.vcf.bgz"],
    );
    run(
        d,
        r,
        "sort-gz-input",
        "SortVcf",
        &["I=@/case-sort/out.vcf.gz", "O=%/out.vcf.gz"],
    );

    run(
        d,
        r,
        "sites-only",
        "MakeSitesOnlyVcf",
        &["I=@/in.vcf", "O=%/out.vcf.gz"],
    );
    run(
        d,
        r,
        "sites-only-out-of-order",
        "MakeSitesOnlyVcf",
        &["I=@/swapped.vcf", "O=%/out.vcf.gz"],
    );
    run(
        d,
        r,
        "sites-only-missing-contig",
        "MakeSitesOnlyVcf",
        &["I=@/missing.vcf", "O=%/out.vcf.gz"],
    );

    let rename = [
        "I=@/single.vcf",
        "O=%/out.vcf.gz",
        "NEW_SAMPLE_NAME=renamed",
    ];
    run(d, r, "rename", "RenameSampleInVcf", &rename);
    let rename = [
        "I=@/single-nodict.vcf",
        "O=%/out.vcf.gz",
        "NEW_SAMPLE_NAME=renamed",
        "CREATE_INDEX=true",
    ];
    run(d, r, "rename-no-dictionary", "RenameSampleInVcf", &rename);

    let update = [
        "I=@/in.vcf",
        "O=%/out.vcf.gz",
        "SEQUENCE_DICTIONARY=@/other.dict",
        "CREATE_INDEX=true",
    ];
    run(
        d,
        r,
        "update-dictionary",
        "UpdateVcfSequenceDictionary",
        &update,
    );

    run(
        d,
        r,
        "filter",
        "FilterVcf",
        &["I=@/in.vcf", "O=%/out.vcf.gz", "MIN_DP=30"],
    );
    let filter = ["I=@/nodict.vcf", "O=%/out.vcf.gz", "MIN_DP=30"];
    run(d, r, "filter-no-dictionary", "FilterVcf", &filter);

    run(
        d,
        r,
        "fix-header",
        "FixVcfHeader",
        &["I=@/in.vcf", "O=%/out.vcf.gz"],
    );

    run(
        d,
        r,
        "merge",
        "MergeVcfs",
        &["I=@/m1.vcf", "I=@/m2.vcf", "O=%/out.vcf.gz"],
    );

    let split = [
        "I=@/in.vcf",
        "SNP_OUTPUT=%/snp.vcf.gz",
        "INDEL_OUTPUT=%/indel.vcf.gz",
        "STRICT=false",
    ];
    run(d, r, "split", "SplitVcfs", &split);

    let liftover = [
        "I=@/in.vcf",
        "O=%/out.vcf.gz",
        "REJECT=%/reject.vcf.gz",
        "CHAIN=@/id.chain",
        "R=@/ref.fasta",
    ];
    run(d, r, "liftover", "LiftoverVcf", &liftover);

    let convert = ["I=@/in.vcf", "O=%/out.vcf.gz", "REQUIRE_INDEX=false"];
    run(d, r, "convert", "VcfFormatConverter", &convert);
    let convert = ["I=@/case-sort/out.vcf.gz", "O=%/out.vcf.gz"];
    run(d, r, "convert-gz-to-gz", "VcfFormatConverter", &convert);
    let convert = ["I=@/case-sort/out.vcf.gz", "O=%/out.vcf"];
    run(d, r, "convert-gz-to-vcf", "VcfFormatConverter", &convert);

    for c in 1..=3 {
        let (input, output) = (format!("I=@/p{c}.vcf"), format!("O=%/p{c}.vcf.gz"));
        run(d, r, &format!("part{c}"), "SortVcf", &[&input, &output]);
    }
    let gather = [
        "I=@/case-part1/p1.vcf.gz",
        "I=@/case-part2/p2.vcf.gz",
        "I=@/case-part3/p3.vcf.gz",
        "O=%/out.vcf.gz",
    ];
    run(d, r, "gather-block-copy", "GatherVcfs", &gather);
    let gather = ["I=@/p1.vcf", "I=@/p2.vcf", "I=@/p3.vcf", "O=%/out.vcf.gz"];
    run(d, r, "gather-conventional", "GatherVcfs", &gather);
    rows
}

/// Every row the reference wrote, in its order, and nothing else. The levels GKL routes through
/// ISA-L (1 and 2) are not among the cases, so the rows do not depend on the build's assembler.
#[test]
fn every_vcf_gz_and_tbi_is_the_reference_s() {
    let expected = golden();
    let actual = rows();
    for (i, (want, got)) in expected.iter().zip(&actual).enumerate() {
        assert_eq!(got, want, "row {}", i + 1);
    }
    assert_eq!(actual.len(), expected.len(), "the number of rows");
}
