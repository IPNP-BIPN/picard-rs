//! What the Illumina dumps in `tools/illumina-conformance/` share: the run directory
//! `MakeIlluminaFixtures.java` writes byte by byte, the binary's run, and the rows the dumps print.
//!
//! Each test file uses part of it.
#![allow(dead_code)]

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The bases a cycle's file carries, one character per cluster.
const CYCLES: [&str; 4] = ["ACGT", "ACGT", "AACC", "GGTT"];
const CLUSTERS: u32 = 4;

fn put(file: &Path, bytes: &[u8]) {
    std::fs::create_dir_all(file.parent().unwrap()).unwrap();
    std::fs::write(file, bytes).unwrap();
}

fn locs(file: &Path) {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&1u32.to_le_bytes());
    bytes.extend_from_slice(&1.0f32.to_le_bytes());
    bytes.extend_from_slice(&CLUSTERS.to_le_bytes());
    for cluster in 0..CLUSTERS {
        bytes.extend_from_slice(&(100.0f32 * (cluster + 1) as f32).to_le_bytes());
        bytes.extend_from_slice(&(200.0f32 * (cluster + 1) as f32).to_le_bytes());
    }
    put(file, &bytes);
}

/// `MakeIlluminaFixtures.write(run, cycles)`: one lane, tile 1101, four clusters of which three
/// pass the filter, every base at quality thirty. Returns the basecalls directory.
pub fn write_run(run: &Path, cycles: usize) -> PathBuf {
    let intensities = run.join("Data/Intensities");
    let root = intensities.join("BaseCalls");
    // Version two, ten-byte records: lane 1, tile 1101, code 100 (the cluster count), 4.0.
    let mut tile_metrics = vec![2u8, 10];
    for value in [1u16, 1101, 100] {
        tile_metrics.extend_from_slice(&value.to_le_bytes());
    }
    tile_metrics.extend_from_slice(&(CLUSTERS as f32).to_le_bytes());
    put(&run.join("InterOp/TileMetricsOut.bin"), &tile_metrics);
    locs(&intensities.join("s.locs"));
    let lane = root.join("L001");
    for cycle in 1..=cycles {
        let mut bcl = (CLUSTERS).to_le_bytes().to_vec();
        bcl.extend(CYCLES[(cycle - 1) % CYCLES.len()].bytes().map(|base| {
            let code = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
            code | (30 << 2)
        }));
        put(&lane.join(format!("C{cycle}.1/s_1_1101.bcl")), &bcl);
    }
    let mut filter = Vec::new();
    for value in [0u32, 3, CLUSTERS] {
        filter.extend_from_slice(&value.to_le_bytes());
    }
    filter.extend_from_slice(&[1, 1, 1, 0]);
    put(&lane.join("s_1_1101.filter"), &filter);
    locs(&lane.join("s_1_1101.locs"));
    std::fs::write(root.join("config.xml"), "<BaseCallAnalysis/>\n").unwrap();
    root
}

/// The dumps' multiplexed pipeline: `ExtractIlluminaBarcodes` over `2T2B` with `AG` and `CT`,
/// writing its barcode files into the basecalls directory. Returns the barcode file.
pub fn extract_barcodes(dir: &Path, basecalls: &Path) -> PathBuf {
    let barcodes = dir.join("barcodes.tsv");
    std::fs::write(
        &barcodes,
        "barcode_sequence_1\tbarcode_name\tlibrary_name\nAG\tfirst\tlibraryA\nCT\tsecond\tlibraryB\n",
    )
    .unwrap();
    let ran = Command::new(env!("CARGO_BIN_EXE_ExtractIlluminaBarcodes"))
        .arg(format!("BASECALLS_DIR={}", basecalls.display()))
        .arg("LANE=1")
        .arg("READ_STRUCTURE=2T2B")
        .arg(format!("BARCODE_FILE={}", barcodes.display()))
        .arg(format!(
            "METRICS_FILE={}",
            dir.join("barcode-metrics.txt").display()
        ))
        .arg(format!("OUTPUT_DIR={}", basecalls.display()))
        .output()
        .unwrap();
    assert!(
        ran.status.success(),
        "ExtractIlluminaBarcodes: {}",
        String::from_utf8_lossy(&ran.stderr)
    );
    barcodes
}

/// A binary's exit status and everything it said, stdout then stderr.
pub struct Ran {
    pub code: i32,
    pub stdout: String,
    pub said: String,
}

pub fn run(binary: &str, argv: &[String]) -> Ran {
    let ran = Command::new(binary).args(argv).output().unwrap();
    Ran {
        code: ran.status.code().unwrap_or(-1),
        stdout: String::from_utf8_lossy(&ran.stdout).into_owned(),
        said: format!(
            "{}{}",
            String::from_utf8_lossy(&ran.stdout),
            String::from_utf8_lossy(&ran.stderr)
        ),
    }
}

/// `table(file)`: the metrics file without its comments and blank lines.
pub fn table(file: &Path) -> String {
    std::fs::read_to_string(file)
        .unwrap()
        .split('\n')
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect::<Vec<_>>()
        .join("\n")
}

/// The names of the files a run wrote, sorted.
pub fn written(out: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(out)
        .unwrap()
        .map(|entry| entry.unwrap())
        .filter(|entry| entry.path().is_file())
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// The rows of one dump, rebuilt from the binary's runs.
pub struct Dump {
    base: PathBuf,
    count: usize,
    text: String,
}

impl Dump {
    pub fn new(name: &str) -> Dump {
        let base = std::env::temp_dir().join(format!("illumina-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        Dump {
            base,
            count: 0,
            text: String::new(),
        }
    }

    /// A fresh directory for one run, as the dump's `createTempDirectory`.
    pub fn dir(&mut self) -> PathBuf {
        self.count += 1;
        let dir = self.base.join(self.count.to_string());
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    pub fn emit(&mut self, kind: &str, case: &str, payload: &str) {
        let payload = payload
            .replace('\\', "\\\\")
            .replace('\t', "\\t")
            .replace('\n', "\\n");
        self.text.push_str(&format!("{kind}\t{case}\t{payload}\n"));
    }

    /// The dump's `error` row, when the run threw or exited non-zero; true when it did.
    ///
    /// A thrown exception is `class:message` with the run's directory as `<dir>` (the binary
    /// prints Java's uncaught form, `class: message`). A non-zero exit is `exit <code>`, followed
    /// by the `ERROR:` lines when the dump `with_reasons` collects them.
    pub fn failed(&mut self, case: &str, dir: &Path, ran: &Ran, with_reasons: bool) -> bool {
        if let Some(thrown) = ran
            .said
            .lines()
            .find_map(|line| line.strip_prefix("Exception in thread \"main\" "))
        {
            let thrown = thrown
                .replacen(": ", ":", 1)
                .replace(&dir.display().to_string(), "<dir>");
            self.emit("error", case, &thrown);
            return true;
        }
        if ran.code == 0 {
            return false;
        }
        let mut row = format!("exit {}", ran.code);
        if with_reasons {
            let reasons: Vec<&str> = ran
                .said
                .split('\n')
                .filter(|line| line.starts_with("ERROR:"))
                .collect();
            if !reasons.is_empty() {
                row = format!("{row} {}", reasons.join(" "));
            }
        }
        self.emit("error", case, &row);
        true
    }

    /// The rebuilt dump against its golden, line for line.
    pub fn check(self, golden: &str) {
        let _ = std::fs::remove_dir_all(&self.base);
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("tests/data/{golden}.txt.gz"));
        let mut expected = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("the golden"))
            .read_to_string(&mut expected)
            .expect("the golden decompresses");
        let expected: Vec<&str> = expected.lines().collect();
        let rebuilt: Vec<&str> = self.text.lines().collect();
        for (index, (want, got)) in expected.iter().zip(&rebuilt).enumerate() {
            assert_eq!(got, want, "{golden}: line {}", index + 1);
        }
        assert_eq!(
            rebuilt.len(),
            expected.len(),
            "{golden}: line count\n{}",
            self.text
        );
    }
}
