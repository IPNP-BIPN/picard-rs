//! The Illumina run directory as Picard's basecalling tools see it: the read structure, the
//! output mapping, the tile metrics, and which files of which format a lane has.
//!
//! [`crate::illumina_files`] holds a reading of the record layouts; this is the layer above it that
//! the binaries need, ported line by line because almost everything a tool reports about a
//! directory is decided here rather than by the records:
//!
//!  * **file names are regular expressions with unescaped dots.** `ParameterizedFileUtil.escapePeriods`
//!    is `replaceAll("\\.", "\\.")`, which replaces a dot with a dot, so `.bcl` in a pattern is
//!    "any character, then bcl" and `s.locs` matches `sXlocs`. The matchers below keep that;
//!  * **a format is chosen before it is checked**: `IlluminaDataProviderFactory.determineFormats`
//!    asks each data type's preferred formats, in order, whether files are available, and only
//!    the first available one is verified or read;
//!  * **an empty file is not a file** for every per-tile format except the barcode files, which
//!    are built with `skipEmptyFiles = false`.
//!
//! What is not ported: the multi-tile formats need a `s_<lane>.bci` tile index, and the CBCL
//! format its own reader. A directory that has either is refused by name ([`Unsupported`])
//! rather than read approximately.
//!
//! Ported from `picard.illumina.parser.ReadStructure`, `OutputMapping`, `IlluminaFileUtil`,
//! `ParameterizedFileUtil`, `PerTileFileUtil`, `PerTileOrPerRunFileUtil`,
//! `PerTilePerCycleFileUtil`, `IlluminaDataProviderFactory` and
//! `readers.TileMetricsOutReader` (Picard 3.4.0), and `htsjdk.samtools.util.IOUtil` (htsjdk 4.2.0).

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

/// An uncaught Java throwable: its class and its message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Thrown {
    pub class: &'static str,
    pub message: String,
}

impl Thrown {
    pub fn new(class: &'static str, message: impl Into<String>) -> Self {
        Thrown {
            class,
            message: message.into(),
        }
    }
    pub fn picard(message: impl Into<String>) -> Self {
        Self::new("picard.PicardException", message)
    }
    pub fn sam(message: impl Into<String>) -> Self {
        Self::new("htsjdk.samtools.SAMException", message)
    }
    pub fn illegal_argument(message: impl Into<String>) -> Self {
        Self::new("java.lang.IllegalArgumentException", message)
    }
    /// What the JVM prints for it.
    pub fn render(&self) -> String {
        format!("{}: {}", self.class, self.message)
    }
}

/// A port limit, said out loud: the directory holds a format this port does not read.
pub fn unsupported(what: &str) -> Thrown {
    Thrown::new(
        "picard_rs.Unsupported",
        format!("{what} is not supported by this port"),
    )
}

/// `ReadType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ReadType {
    T,
    B,
    M,
    S,
}

impl ReadType {
    pub fn letter(self) -> char {
        match self {
            ReadType::T => 'T',
            ReadType::B => 'B',
            ReadType::M => 'M',
            ReadType::S => 'S',
        }
    }
}

/// `ReadDescriptor`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadDescriptor {
    pub length: i32,
    pub kind: ReadType,
}

/// `ReadStructure`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadStructure {
    pub descriptors: Vec<ReadDescriptor>,
    pub total_cycles: i32,
    pub read_lengths: Vec<i32>,
    /// Each descriptor's zero-based inclusive cycle range.
    ranges: Vec<(i32, i32)>,
}

const READ_STRUCTURE_MSG: &str = "Read structure must be formatted as follows: <number of bases><type><number of bases><type>...<number of bases> where number of bases is a positive (NON-ZERO) integer and type is one of the following characters T,B,M,S (e.g. 76T8B68T would denote a paired-end run with a 76 base first end an 8 base barcode followed by a 68 base second end).";

impl ReadStructure {
    /// `new ReadStructure(List<ReadDescriptor>)`.
    pub fn new(descriptors: Vec<ReadDescriptor>) -> Result<Self, Thrown> {
        if descriptors.is_empty() {
            return Err(Thrown::illegal_argument(
                "ReadStructure does not support 0 length clusters!",
            ));
        }
        let mut ranges = Vec::new();
        let mut cycle = 0;
        let mut total = 0;
        for descriptor in &descriptors {
            if descriptor.length <= 0 {
                return Err(Thrown::illegal_argument(format!(
                    "ReadStructure only supports ReadDescriptor lengths > 0, found({})",
                    descriptor.length
                )));
            }
            let end = cycle + descriptor.length - 1;
            ranges.push((cycle, end));
            cycle = end + 1;
            total += descriptor.length;
        }
        Ok(ReadStructure {
            read_lengths: descriptors.iter().map(|d| d.length).collect(),
            descriptors,
            total_cycles: total,
            ranges,
        })
    }

    /// `new ReadStructure(String)`: `^((\d+)([TBMS]{1}))+$`, each count through
    /// `Integer.parseInt`.
    pub fn parse(text: &str) -> Result<Self, Thrown> {
        let refuse = || {
            Thrown::illegal_argument(format!(
                "{text} cannot be parsed as a ReadStructure! {READ_STRUCTURE_MSG}"
            ))
        };
        let bytes = text.as_bytes();
        let mut descriptors = Vec::new();
        let mut at = 0;
        if bytes.is_empty() {
            return Err(refuse());
        }
        while at < bytes.len() {
            let start = at;
            while at < bytes.len() && bytes[at].is_ascii_digit() {
                at += 1;
            }
            if at == start || at >= bytes.len() {
                return Err(refuse());
            }
            let kind = match bytes[at] {
                b'T' => ReadType::T,
                b'B' => ReadType::B,
                b'M' => ReadType::M,
                b'S' => ReadType::S,
                _ => return Err(refuse()),
            };
            let digits = &text[start..at];
            let length: i32 = digits.parse().map_err(|_| {
                Thrown::new(
                    "java.lang.NumberFormatException",
                    format!("For input string: \"{digits}\""),
                )
            })?;
            descriptors.push(ReadDescriptor { length, kind });
            at += 1;
        }
        Self::new(descriptors)
    }

    /// The indices of the descriptors of the given types.
    pub fn indices(&self, kinds: &[ReadType]) -> Vec<usize> {
        (0..self.descriptors.len())
            .filter(|i| kinds.contains(&self.descriptors[*i].kind))
            .collect()
    }

    /// `Substructure.getCycles()`: one-based cycles of the given descriptors, in order.
    pub fn cycles_of(&self, indices: &[usize]) -> Vec<i32> {
        let mut out = Vec::new();
        for index in indices {
            let (start, end) = self.ranges[*index];
            out.extend((start..=end).map(|c| c + 1));
        }
        out
    }

    /// `ReadStructure.toString()`.
    pub fn render(&self) -> String {
        self.descriptors
            .iter()
            .map(|d| format!("{}{}", d.length, d.kind.letter()))
            .collect()
    }
}

/// `OutputMapping`: the non-skip descriptors, which are what a cluster is read into.
#[derive(Debug, Clone)]
pub struct OutputMapping {
    pub output_cycles: Vec<i32>,
    pub output_read_lengths: Vec<i32>,
    pub output_descriptors: Vec<ReadDescriptor>,
}

impl OutputMapping {
    pub fn new(structure: &ReadStructure) -> Result<Self, Thrown> {
        let non_skips = structure.indices(&[ReadType::T, ReadType::B, ReadType::M]);
        let descriptors: Vec<ReadDescriptor> = non_skips
            .iter()
            .map(|i| structure.descriptors[*i])
            .collect();
        // `outputSubstructure.toReadStructure()` refuses an empty list.
        let output = ReadStructure::new(descriptors.clone())?;
        Ok(OutputMapping {
            output_cycles: structure.cycles_of(&non_skips),
            output_read_lengths: output.read_lengths,
            output_descriptors: descriptors,
        })
    }

    pub fn total_output_cycles(&self) -> i32 {
        self.output_cycles.len() as i32
    }
}

/// `IlluminaFileUtil.longLaneStr`.
pub fn long_lane_str(lane: i32) -> String {
    let digits = lane.to_string();
    let zeros = 3usize.saturating_sub(digits.len());
    format!("L{}{}", "0".repeat(zeros), digits)
}

/// `File.getAbsolutePath()`.
pub fn absolute(path: &Path) -> String {
    if path.is_absolute() {
        return path.display().to_string();
    }
    let cwd = std::env::current_dir().unwrap_or_default();
    cwd.join(path).display().to_string()
}

/// `Path.toUri().toString()` for an absolute path of plain characters.
pub fn file_uri(path: &Path) -> String {
    format!("file://{}", absolute(path))
}

/// `IOUtil.assertDirectoryIsReadable`.
pub fn assert_directory_is_readable(dir: &Path) -> Result<(), Thrown> {
    if !dir.exists() {
        return Err(Thrown::sam(format!(
            "Directory does not exist: {}",
            absolute(dir)
        )));
    }
    if !dir.is_dir() {
        return Err(Thrown::sam(format!(
            "Cannot read from directory because it is not a directory: {}",
            absolute(dir)
        )));
    }
    Ok(())
}

/// `IOUtil.assertFileIsReadable`.
pub fn assert_file_is_readable(file: &Path) -> Result<(), Thrown> {
    if !file.exists() {
        return Err(Thrown::sam(format!(
            "Cannot read non-existent file: {}",
            file_uri(file)
        )));
    }
    if file.is_dir() {
        return Err(Thrown::sam(format!(
            "Cannot read file because it is a directory: {}",
            file_uri(file)
        )));
    }
    Ok(())
}

/// `File.length()`, which is zero for a file that is not there.
pub fn file_length(path: &Path) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

/// `File.listFiles(filter)`: `None` where the directory is not there, the matching entries
/// otherwise, in the directory's own order.
pub fn list_matching(dir: &Path, matches: &dyn Fn(&str) -> bool) -> Option<Vec<PathBuf>> {
    let entries = std::fs::read_dir(dir).ok()?;
    let mut out = Vec::new();
    for entry in entries.flatten() {
        let name = entry.file_name().to_string_lossy().to_string();
        if matches(&name) {
            out.push(entry.path());
        }
    }
    Some(out)
}

/// A pattern literal as `Pattern` reads it after `escapePeriods`: a dot is any character.
fn literal_matches(text: &[char], pattern: &[char]) -> bool {
    text.len() == pattern.len()
        && text
            .iter()
            .zip(pattern)
            .all(|(t, p)| *p == '.' && !matches!(t, '\n' | '\r') || t == p)
}

/// `^<prefix>(\d{1,5})<suffix>$` with an optional `(\.gz|\.bz2)?` after it, where `prefix` and
/// `suffix` are literals whose dots match anything. Returns the captured number, taking the
/// longest run of digits that lets the rest match, which is the order a greedy quantifier tries.
fn match_numbered(name: &str, prefix: &str, suffix: &str, compression: bool) -> Option<i32> {
    let chars: Vec<char> = name.chars().collect();
    let prefix: Vec<char> = prefix.chars().collect();
    let suffix: Vec<char> = suffix.chars().collect();
    if chars.len() < prefix.len() || !literal_matches(&chars[..prefix.len()], &prefix) {
        return None;
    }
    let rest = &chars[prefix.len()..];
    let max_digits = rest
        .iter()
        .take(5)
        .take_while(|c| c.is_ascii_digit())
        .count();
    for digits in (1..=max_digits).rev() {
        let tail = &rest[digits..];
        let tails: Vec<&[char]> = if compression {
            let mut options = vec![tail];
            for ext in [".gz", ".bz2"] {
                let ext: Vec<char> = ext.chars().collect();
                if tail.len() >= ext.len()
                    && tail[tail.len() - ext.len()..] == ext[..]
                    && tail.len() - ext.len() >= suffix.len()
                {
                    options.push(&tail[..tail.len() - ext.len()]);
                }
            }
            options
        } else {
            vec![tail]
        };
        if tails.iter().any(|t| literal_matches(t, &suffix)) {
            let number: String = rest[..digits].iter().collect();
            return number.parse().ok();
        }
    }
    None
}

/// `IlluminaFileUtil.CYCLE_SUBDIRECTORY_PATTERN`, `^C(\d+)\.1$`, with a real dot.
pub fn cycle_of_dir(name: &str) -> Option<i32> {
    let rest = name.strip_prefix('C')?.strip_suffix(".1")?;
    if rest.is_empty() || !rest.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    rest.parse().ok()
}

/// `IlluminaFileUtil.hasCbcls`.
pub fn has_cbcls(basecall_dir: &Path, lane: i32) -> bool {
    let lane_dir = basecall_dir.join(long_lane_str(lane));
    let Some(cycle_dirs) = list_matching(&lane_dir, &|name| cycle_of_dir(name).is_some()) else {
        return false;
    };
    let prefix = format!("{}_", long_lane_str(lane));
    cycle_dirs.iter().any(|dir| {
        list_matching(dir, &|name| {
            match_numbered(name, &prefix, ".cbcl", false).is_some()
        })
        .is_some_and(|files| !files.is_empty())
    })
}

/// One `TileMetricsOut.bin` record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TileMetricRecord {
    pub lane: i32,
    pub tile: i32,
    pub code: i32,
    pub value: f32,
    /// Version three's record type byte (`'t'` for a cluster record).
    pub kind: u8,
}

/// `TileMetricsOutReader`: the version, version three's density, and the whole records.
#[derive(Debug, Clone)]
pub struct TileMetricsOut {
    pub version: i32,
    pub density: f32,
    pub records: Vec<TileMetricRecord>,
}

/// `MMapBackedIteratorFactory.checkFactoryVars`.
fn check_factory_vars(header_size: usize, file: &Path) -> Result<Vec<u8>, Thrown> {
    assert_file_is_readable(file)?;
    let bytes = std::fs::read(file).map_err(|e| Thrown::picard(e.to_string()))?;
    if header_size > bytes.len() {
        return Err(Thrown::picard(format!(
            "Header size({header_size}) is greater than file size({}) for file {}",
            bytes.len(),
            absolute(file)
        )));
    }
    Ok(bytes)
}

pub fn read_tile_metrics(file: &Path) -> Result<TileMetricsOut, Thrown> {
    let bytes = check_factory_vars(1, file)?;
    let version = bytes[0] as i32;
    let (header_size, record_size) = match version {
        2 => (2usize, 10usize),
        3 => (6, 15),
        _ => {
            return Err(Thrown::new(
                "java.lang.NullPointerException",
                "Cannot read field \"headerSize\" because \"this.version\" is null",
            ))
        }
    };
    let bytes = check_factory_vars(header_size, file)?;
    let record_size_in_header = bytes[1] as usize;
    if record_size != record_size_in_header {
        return Err(Thrown::picard(format!(
            "TileMetricsOutReader expects the record size to be {record_size}.  Actual Record Size in Header({record_size_in_header})"
        )));
    }
    let density = if version == 3 {
        // `UnsignedTypeUtil.uIntToFloat`: the unsigned VALUE as a float, not its bits.
        u32::from_le_bytes(bytes[2..6].try_into().unwrap()) as f32
    } else {
        0.0
    };
    let u16_at = |at: usize| u16::from_le_bytes([bytes[at], bytes[at + 1]]) as i32;
    let f32_at = |at: usize| f32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
    let mut records = Vec::new();
    let mut at = header_size;
    while bytes.len() - at >= record_size {
        records.push(if version == 3 {
            TileMetricRecord {
                lane: u16_at(at),
                tile: i32::from_le_bytes(bytes[at + 2..at + 6].try_into().unwrap()),
                code: 0,
                kind: bytes[at + 6],
                value: f32_at(at + 7),
            }
        } else {
            TileMetricRecord {
                lane: u16_at(at),
                tile: u16_at(at + 2),
                code: u16_at(at + 4),
                value: f32_at(at + 6),
                kind: 0,
            }
        });
        at += record_size;
    }
    Ok(TileMetricsOut {
        version,
        density,
        records,
    })
}

/// `IlluminaDataType`, in declaration order (which is what a `TreeSet` of them iterates in).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DataType {
    Position,
    BaseCalls,
    QualityScores,
    PF,
    Barcodes,
}

impl DataType {
    pub fn parse(name: &str) -> Option<Self> {
        Some(match name {
            "Position" => DataType::Position,
            "BaseCalls" => DataType::BaseCalls,
            "QualityScores" => DataType::QualityScores,
            "PF" => DataType::PF,
            "Barcodes" => DataType::Barcodes,
            _ => return None,
        })
    }
    pub fn name(self) -> &'static str {
        match self {
            DataType::Position => "Position",
            DataType::BaseCalls => "BaseCalls",
            DataType::QualityScores => "QualityScores",
            DataType::PF => "PF",
            DataType::Barcodes => "Barcodes",
        }
    }
}

/// `IlluminaFileUtil.SupportedIlluminaFormat`, in declaration order (an `EnumMap`'s order).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Format {
    Bcl,
    Cbcl,
    Locs,
    Clocs,
    Pos,
    Filter,
    Barcode,
    MultiTileFilter,
    MultiTileLocs,
    MultiTileBcl,
}

/// `IlluminaDataProviderFactory.DATA_TYPE_TO_PREFERRED_FORMATS`.
pub fn preferred_formats(data_type: DataType) -> &'static [Format] {
    match data_type {
        DataType::BaseCalls | DataType::QualityScores => {
            &[Format::MultiTileBcl, Format::Bcl, Format::Cbcl]
        }
        DataType::PF => &[Format::MultiTileFilter, Format::Filter],
        DataType::Position => &[
            Format::MultiTileLocs,
            Format::Locs,
            Format::Clocs,
            Format::Pos,
        ],
        DataType::Barcodes => &[Format::Barcode],
    }
}

/// A `PerTileFileUtil` (and the per-run `s.locs` variant).
#[derive(Debug, Clone)]
pub struct PerTileUtil {
    pub extension: String,
    pub base: PathBuf,
    /// tile -> file, a `TreeMap`.
    pub files: BTreeMap<i32, PathBuf>,
    pub tiles: Vec<i32>,
    /// `PerTileOrPerRunFileUtil`'s run file.
    pub run_file: Option<PathBuf>,
    per_run: bool,
}

impl PerTileUtil {
    fn new(
        extension: &str,
        base: &Path,
        lane: i32,
        skip_empty: bool,
        per_run: bool,
    ) -> Result<Self, Thrown> {
        let compression = extension.ends_with(".txt");
        let prefix = format!("s_{lane}_");
        let mut files = BTreeMap::new();
        if base.exists() {
            assert_directory_is_readable(base)?;
            let matched = list_matching(base, &|name| {
                match_numbered(name, &prefix, extension, compression).is_some()
            })
            .unwrap_or_default();
            for file in matched {
                if !skip_empty || file_length(&file) > 0 {
                    let name = file.file_name().unwrap().to_string_lossy().to_string();
                    let tile = match_numbered(&name, &prefix, extension, compression).unwrap();
                    files.insert(tile, file);
                }
            }
        }
        let tiles: Vec<i32> = files.keys().copied().collect();
        let mut util = PerTileUtil {
            extension: extension.to_string(),
            base: base.to_path_buf(),
            files,
            tiles,
            run_file: None,
            per_run,
        };
        if per_run {
            // `getRunFile(base.getParentFile(), ^s<extension>$)`: exactly one match, or none.
            let parent = base.parent().map(Path::to_path_buf).unwrap_or_default();
            if parent.exists() {
                assert_directory_is_readable(&parent)?;
                let pattern: Vec<char> = format!("s{extension}").chars().collect();
                let matched = list_matching(&parent, &|name| {
                    literal_matches(&name.chars().collect::<Vec<_>>(), &pattern)
                })
                .unwrap_or_default();
                if matched.len() == 1 {
                    util.run_file = Some(matched[0].clone());
                }
            }
        }
        Ok(util)
    }

    pub fn files_available(&self) -> bool {
        !self.files.is_empty() || self.run_file.is_some()
    }

    /// `setTilesForPerRunFile`: the run file stands for every tile.
    pub fn set_tiles_for_per_run_file(&mut self, tiles: &[i32]) {
        if let Some(run) = &self.run_file {
            for tile in tiles {
                self.files.insert(*tile, run.clone());
            }
            self.tiles = tiles.to_vec();
        }
    }

    pub fn check_tile_count(&self) -> bool {
        self.run_file.is_none()
    }

    /// `PerTileFileUtil.verify` / `PerTileOrPerRunFileUtil.verify`.
    pub fn verify(&self, expected_tiles: &[i32]) -> Vec<String> {
        let base = match (&self.run_file, self.per_run) {
            (Some(run), true) => run.parent().map(Path::to_path_buf).unwrap_or_default(),
            _ => self.base.clone(),
        };
        let mut failures = Vec::new();
        if !base.exists() {
            failures.push(format!(
                "Base directory({}) does not exist!",
                absolute(&base)
            ));
        } else if !expected_tiles.iter().all(|t| self.tiles.contains(t)) {
            let missing: Vec<String> = expected_tiles
                .iter()
                .filter(|t| !self.tiles.contains(t))
                .map(|t| t.to_string())
                .collect();
            failures.push(format!(
                "Missing tile [{}] for file type {}.",
                missing.join(", "),
                self.extension
            ));
        }
        failures
    }
}

/// A `PerTilePerCycleFileUtil`.
#[derive(Debug, Clone)]
pub struct PerTilePerCycleUtil {
    pub extension: String,
    pub base: PathBuf,
    /// cycle -> tile -> file.
    pub cycles: BTreeMap<i32, BTreeMap<i32, PathBuf>>,
    pub detected_cycles: BTreeSet<i32>,
    /// `null` until a cycle directory is found, as in the reference.
    pub tiles: Option<Vec<i32>>,
}

impl PerTilePerCycleUtil {
    fn new(extension: &str, base: &Path, lane: i32) -> Result<Self, Thrown> {
        let mut util = PerTilePerCycleUtil {
            extension: extension.to_string(),
            base: base.to_path_buf(),
            cycles: BTreeMap::new(),
            detected_cycles: BTreeSet::new(),
            tiles: None,
        };
        let Some(cycle_dirs) = list_matching(base, &|name| cycle_of_dir(name).is_some()) else {
            return Ok(util);
        };
        if cycle_dirs.is_empty() {
            return Ok(util);
        }
        let prefix = format!("s_{lane}_");
        let mut unique: BTreeSet<i32> = BTreeSet::new();
        for dir in &cycle_dirs {
            let cycle = cycle_of_dir(&dir.file_name().unwrap().to_string_lossy()).unwrap();
            util.detected_cycles.insert(cycle);
            let mut files = BTreeMap::new();
            if dir.exists() {
                assert_directory_is_readable(dir)?;
                for file in list_matching(dir, &|name| {
                    match_numbered(name, &prefix, extension, false).is_some()
                })
                .unwrap_or_default()
                {
                    if file_length(&file) > 0 {
                        let name = file.file_name().unwrap().to_string_lossy().to_string();
                        files.insert(
                            match_numbered(&name, &prefix, extension, false).unwrap(),
                            file,
                        );
                    }
                }
            }
            unique.extend(files.keys().copied());
            util.cycles.insert(cycle, files);
        }
        // A `HashSet<Integer>` of tile numbers; only its membership is ever read.
        util.tiles = Some(unique.into_iter().collect());
        Ok(util)
    }

    pub fn files_available(&self) -> bool {
        self.cycles.values().any(|files| !files.is_empty())
    }

    /// `getFiles(tiles, cycles)`: the cycles that exist, each with the tiles asked for.
    pub fn files(&self, tiles: &[i32], cycles: &[i32]) -> BTreeMap<i32, BTreeMap<i32, PathBuf>> {
        let wanted: BTreeSet<i32> = cycles
            .iter()
            .copied()
            .filter(|c| self.detected_cycles.contains(c))
            .collect();
        let mut out = BTreeMap::new();
        for cycle in wanted {
            if let Some(files) = self.cycles.get(&cycle) {
                let kept: BTreeMap<i32, PathBuf> = tiles
                    .iter()
                    .filter_map(|t| files.get(t).map(|f| (*t, f.clone())))
                    .collect();
                out.insert(cycle, kept);
            }
        }
        out
    }

    /// `PerTilePerCycleFileUtil.verify`.
    pub fn verify(&self, expected_tiles: &[i32], expected_cycles: &[i32]) -> Vec<String> {
        let mut failures = Vec::new();
        if !self.base.exists() {
            failures.push(format!(
                "Base directory({}) does not exist!",
                absolute(&self.base)
            ));
            return failures;
        }
        let files = self.files(expected_tiles, expected_cycles);
        let mut first_length: BTreeMap<i32, u64> = BTreeMap::new();
        for cycle in expected_cycles {
            match files.get(cycle) {
                Some(tiles) => {
                    for tile in expected_tiles {
                        match tiles.get(tile) {
                            Some(file) => {
                                let length = file_length(file);
                                match first_length.get(tile) {
                                    None => {
                                        first_length.insert(*tile, length);
                                    }
                                    Some(first) => {
                                        if self.extension != ".bcl.gz" && *first != length {
                                            failures.push(format!(
                                                "File type {} has cycles files of different length.  Current cycle ({cycle}) Length of first non-empty file ({first}) length of current cycle ({length}) File({})",
                                                self.extension,
                                                absolute(file)
                                            ));
                                        }
                                    }
                                }
                            }
                            None => failures.push(format!(
                                "File type {} is missing a file for cycle {cycle} and tile {tile}",
                                self.extension
                            )),
                        }
                    }
                }
                None => failures.push(format!(
                    "Missing file for cycle {cycle} in directory {} for file type {}",
                    absolute(&self.base),
                    self.extension
                )),
            }
        }
        failures
    }
}

/// One format's file utility, as `IlluminaFileUtil.getUtil` builds it.
#[derive(Debug, Clone)]
pub enum Util {
    PerTile(PerTileUtil),
    PerTilePerCycle(PerTilePerCycleUtil),
    /// A multi-tile format with no `.bci`, which reports no files.
    Absent,
    /// CBCL, whose availability is `hasCbcls`.
    Cbcl(bool),
}

impl Util {
    pub fn files_available(&self) -> bool {
        match self {
            Util::PerTile(u) => u.files_available(),
            Util::PerTilePerCycle(u) => u.files_available(),
            Util::Absent => false,
            Util::Cbcl(available) => *available,
        }
    }
}

/// `IlluminaFileUtil`: the directories of one lane and the utilities it has built so far.
pub struct IlluminaFileUtil {
    pub lane: i32,
    pub basecall_dir: PathBuf,
    pub barcode_dir: Option<PathBuf>,
    pub basecall_lane_dir: PathBuf,
    pub intensity_dir: PathBuf,
    pub intensity_lane_dir: PathBuf,
    pub tile_metrics_out: PathBuf,
    utils: BTreeMap<Format, Util>,
}

fn parent(path: &Path) -> PathBuf {
    path.parent().map(Path::to_path_buf).unwrap_or_default()
}

impl IlluminaFileUtil {
    pub fn new(basecall_dir: &Path, barcode_dir: Option<&Path>, lane: i32) -> Self {
        let intensity_dir = parent(basecall_dir);
        let data_dir = parent(&intensity_dir);
        IlluminaFileUtil {
            lane,
            basecall_dir: basecall_dir.to_path_buf(),
            barcode_dir: barcode_dir.map(Path::to_path_buf),
            basecall_lane_dir: basecall_dir.join(long_lane_str(lane)),
            intensity_lane_dir: intensity_dir.join(long_lane_str(lane)),
            tile_metrics_out: parent(&data_dir).join("InterOp").join("TileMetricsOut.bin"),
            intensity_dir,
            utils: BTreeMap::new(),
        }
    }

    /// `getUtil`, building and caching the utility on first use.
    pub fn util(&mut self, format: Format) -> Result<&mut Util, Thrown> {
        if !self.utils.contains_key(&format) {
            let lane = self.lane;
            let util = match format {
                Format::Bcl => {
                    let plain = PerTilePerCycleUtil::new(".bcl", &self.basecall_lane_dir, lane)?;
                    let gz = PerTilePerCycleUtil::new(".bcl.gz", &self.basecall_lane_dir, lane)?;
                    match (plain.files_available(), gz.files_available()) {
                        (true, false) | (false, false) => Util::PerTilePerCycle(plain),
                        (false, true) => Util::PerTilePerCycle(gz),
                        (true, true) => {
                            return Err(Thrown::picard(format!(
                                "Not all BCL files in {} have the same extension!",
                                absolute(&self.basecall_lane_dir)
                            )))
                        }
                    }
                }
                Format::Cbcl => Util::Cbcl(has_cbcls(&self.basecall_dir, lane)),
                Format::Locs => Util::PerTile(PerTileUtil::new(
                    ".locs",
                    &self.intensity_lane_dir,
                    lane,
                    true,
                    true,
                )?),
                Format::Clocs => Util::PerTile(PerTileUtil::new(
                    ".clocs",
                    &self.intensity_lane_dir,
                    lane,
                    true,
                    false,
                )?),
                Format::Pos => Util::PerTile(PerTileUtil::new(
                    "_pos.txt",
                    &self.intensity_dir,
                    lane,
                    true,
                    false,
                )?),
                Format::Filter => Util::PerTile(PerTileUtil::new(
                    ".filter",
                    &self.basecall_lane_dir,
                    lane,
                    true,
                    false,
                )?),
                Format::Barcode => {
                    let dir = self
                        .barcode_dir
                        .clone()
                        .unwrap_or_else(|| self.basecall_dir.clone());
                    Util::PerTile(PerTileUtil::new("_barcode.txt", &dir, lane, false, false)?)
                }
                Format::MultiTileFilter | Format::MultiTileLocs | Format::MultiTileBcl => {
                    // All three need the lane's tile index; without it there are no files.
                    let bci = self.basecall_lane_dir.join(format!("s_{lane}.bci"));
                    if bci.exists() {
                        return Err(unsupported("A multi-tile run directory (s_<lane>.bci)"));
                    }
                    Util::Absent
                }
            };
            self.utils.insert(format, util);
        }
        Ok(self.utils.get_mut(&format).unwrap())
    }

    /// `getExpectedTiles`: the lane's tiles in the tile metrics, sorted.
    pub fn expected_tiles(&self) -> Result<Vec<i32>, Thrown> {
        assert_file_is_readable(&self.tile_metrics_out)?;
        let metrics = read_tile_metrics(&self.tile_metrics_out)?;
        let tiles: BTreeSet<i32> = metrics
            .records
            .iter()
            .filter(|r| r.lane == self.lane)
            .map(|r| r.tile)
            .collect();
        Ok(tiles.into_iter().collect())
    }
}

/// `IlluminaDataProviderFactory.determineFormats`: each requested type's first available
/// format, grouped by format in `EnumMap` order.
pub fn determine_formats(
    requested: &BTreeSet<DataType>,
    util: &mut IlluminaFileUtil,
) -> Result<BTreeMap<Format, BTreeSet<DataType>>, Thrown> {
    let mut out: BTreeMap<Format, BTreeSet<DataType>> = BTreeMap::new();
    for data_type in requested {
        for format in preferred_formats(*data_type) {
            if util.util(*format)?.files_available() {
                out.entry(*format).or_default().insert(*data_type);
                break;
            }
        }
    }
    Ok(out)
}

/// `findUnmatchedTypes`.
pub fn unmatched_types(
    requested: &BTreeSet<DataType>,
    formats: &BTreeMap<Format, BTreeSet<DataType>>,
) -> Vec<DataType> {
    requested
        .iter()
        .copied()
        .filter(|t| !formats.values().any(|types| types.contains(t)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dots_in_extensions_match_any_character() {
        assert_eq!(
            match_numbered("s_1_1101.bcl", "s_1_", ".bcl", false),
            Some(1101)
        );
        assert_eq!(
            match_numbered("s_1_1101xbcl", "s_1_", ".bcl", false),
            Some(1101)
        );
        assert_eq!(
            match_numbered("s_1_11012bcl", "s_1_", ".bcl", false),
            Some(1101)
        );
        assert_eq!(match_numbered("s_2_1101.bcl", "s_1_", ".bcl", false), None);
        assert_eq!(
            match_numbered("s_1_1101_barcode.txt.gz", "s_1_", "_barcode.txt", true),
            Some(1101)
        );
    }

    #[test]
    fn read_structures_parse_and_refuse() {
        let rs = ReadStructure::parse("4T4B4T").unwrap();
        assert_eq!(rs.total_cycles, 12);
        assert_eq!(rs.cycles_of(&[1]), vec![5, 6, 7, 8]);
        assert!(ReadStructure::parse("4t").is_err());
        assert!(ReadStructure::parse("0T").is_err());
        assert_eq!(long_lane_str(1), "L001");
    }
}
