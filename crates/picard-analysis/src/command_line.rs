//! The lines a Picard tool writes about its own run: the command line, `Started on:`, and `@PG`.
//!
//! `CommandLineProgram.instanceMain` builds two default headers right after parsing,
//! `new StringHeader(commandLine)` and `new StringHeader("Started on: " + startDate)`, and
//! `getMetricsFile()` puts both at the head of every metrics file the tool writes. The tools that
//! stamp their provenance into a SAM header (`getPGRecord`, `getChainedPgIds`) put the same
//! command line in the record's `CL`.
//!
//! # Which command line
//!
//! picard-rs takes `NAME=value`, which is what `CommandLineSyntaxTranslater.isLegacyPicardStyle`
//! recognises, so the parser Picard uses for the same command line is Barclay's
//! `LegacyCommandLineArgumentParser`, and the line is the one its `checkNumArguments` builds:
//!
//! * the tool's simple name;
//! * every positional argument, ` <value>`;
//! * every option that was SET, ` NAME=<the field's value>`, in the parser's declaration order
//!   rather than the order typed, and by its long name whatever name was typed;
//! * a separator of three spaces;
//! * every option left unset whose default string is not `null`, ` NAME=<default>`, in the same
//!   order. An empty collection's default string is `null`, so an empty list is left out.
//!
//! The field's value is what the type's `String` constructor made of what was typed, printed with
//! `toString()`, so `0.250` comes back `0.25`, `+0042` comes back `42`, `TRUE` comes back `true`
//! and `dir//x/` comes back `dir/x`. A collection prints as a Java list, `[a, b]`, and a value of
//! `null` empties it to `[]`; and the legacy parser ADDS to the collection the field already holds,
//! so a typed value lands after the default elements rather than replacing them.
//!
//! What the parser declares, in which order and with which default, is data read off the parser
//! itself by `tools/headerlines-conformance/HeaderLinesDump.java`: see `tool_options`.

use crate::tool_options::{Kind, Order, ToolOptions, TOOLS};

/// `CommandLineProgram.getVersion()` under the legacy parser: the jar's implementation version.
pub const PICARD_VERSION: &str = "3.5.0";

/// What one option holds once the command line is parsed.
#[derive(Clone, Debug)]
enum Value {
    /// A scalar the user named, or `null`.
    Scalar(Option<String>),
    /// A collection's elements, defaults first.
    List(Vec<String>),
}

/// The declarations of one tool, by its simple name.
pub fn tool_options(tool: &str) -> Option<&'static ToolOptions> {
    TOOLS.iter().find(|t| t.name == tool)
}

/// `java.io.UnixFileSystem.normalize`: repeated slashes collapse and a trailing one goes, except
/// for the root itself.
fn normalize_file(path: &str) -> String {
    let mut out = String::with_capacity(path.len());
    let mut previous_slash = false;
    for c in path.chars() {
        if c == '/' {
            if previous_slash {
                continue;
            }
            previous_slash = true;
        } else {
            previous_slash = false;
        }
        out.push(c);
    }
    if out.len() > 1 && out.ends_with('/') {
        out.pop();
    }
    out
}

/// `Integer.parseInt`'s grammar: an optional sign and decimal digits.
fn java_integer(text: &str) -> Option<i64> {
    let digits = text.strip_prefix(['+', '-']).unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    text.strip_prefix('+').unwrap_or(text).parse().ok()
}

/// `Double.valueOf`, close enough for the values a command line carries: Rust's grammar plus
/// Java's trailing type suffix and its spelling of the two specials.
fn java_double(text: &str) -> Option<f64> {
    let trimmed = text.trim();
    let body = trimmed
        .strip_suffix(['d', 'D', 'f', 'F'])
        .filter(|b| !b.is_empty() && !b.ends_with(|c: char| c.is_ascii_alphabetic()))
        .unwrap_or(trimmed);
    match body.strip_prefix(['+', '-']).unwrap_or(body) {
        "Infinity" | "NaN" => {
            let value = if body.ends_with("NaN") {
                f64::NAN
            } else {
                f64::INFINITY
            };
            Some(if body.starts_with('-') { -value } else { value })
        }
        rest if rest
            .chars()
            .any(|c| c.is_ascii_alphabetic() && c != 'e' && c != 'E') =>
        {
            None
        }
        _ => body.parse().ok(),
    }
}

/// What the option's type makes of a typed value, printed as `toString()` prints it. A value the
/// type refuses is kept as typed: Picard stops on it, so the line is never written.
fn render(kind: Kind, text: &str) -> String {
    match kind {
        Kind::Text => text.to_string(),
        Kind::File => normalize_file(text),
        Kind::Integer => java_integer(text).map_or_else(|| text.to_string(), |v| v.to_string()),
        Kind::Double => {
            java_double(text).map_or_else(|| text.to_string(), jmath::to_string::double_to_string)
        }
        Kind::Float => java_double(text).map_or_else(
            || text.to_string(),
            |v| jmath::to_string::float_to_string(v as f32),
        ),
        // `new Boolean(s)` is `parseBoolean`: true for `true` in any case, false for anything else.
        Kind::Boolean => text.eq_ignore_ascii_case("true").to_string(),
    }
}

/// A collection's elements in the order its class iterates them.
fn ordered(order: Order, constants: &[&str], list: &[String]) -> Vec<String> {
    match order {
        Order::Enum => {
            let mut sorted = list.to_vec();
            sorted.sort_by_key(|c| constants.iter().position(|k| k == c).unwrap_or(usize::MAX));
            sorted
        }
        Order::List | Order::Linked => list.to_vec(),
        Order::Tree => {
            let mut sorted = list.to_vec();
            sorted.sort();
            sorted
        }
        Order::Hash => {
            let keyed = list
                .iter()
                .map(|s| (crate::java_hash_map::string_hash_code(s), s.clone()))
                .collect();
            crate::barcode_extractor::hash_set_order(keyed, 16)
        }
    }
}

/// A non-empty collection's default string back into its elements: `AbstractCollection.toString`
/// is `[` + the elements joined with `, ` + `]`.
fn default_elements(default: &str) -> Vec<String> {
    match default.strip_prefix('[').and_then(|d| d.strip_suffix(']')) {
        Some("") => Vec::new(),
        Some(inner) => inner.split(", ").map(str::to_string).collect(),
        None => Vec::new(),
    }
}

/// `LegacyCommandLineArgumentParser.getCommandLine()` for `args` (the arguments after the program
/// name), as `parseArguments` and `checkNumArguments` build it.
///
/// A name the tool does not declare is left out: Picard refuses the whole command line there, so
/// whatever this returns for it is never compared with anything.
/// `parseArguments`: what each option holds, and the positional values, after `args`.
fn parse(options: &ToolOptions, args: &[String]) -> (Vec<String>, Vec<Option<Value>>) {
    let mut values: Vec<Option<Value>> = vec![None; options.options.len()];
    let mut positional = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        i += 1;
        // `arg.split("=", 2)`: a NAME= with nothing after it takes the next word as its value.
        let Some((key, value)) = arg.split_once('=') else {
            if let Some(kind) = options.positional {
                positional.push(render(kind, arg));
            }
            continue;
        };
        let value = if value.is_empty() && i < args.len() {
            i += 1;
            args[i - 1].as_str()
        } else {
            value
        };
        let key = key.to_ascii_uppercase();
        let Some(index) = options
            .options
            .iter()
            .position(|o| o.name == key || (!o.short.is_empty() && o.short == key))
        else {
            continue;
        };
        let option = &options.options[index];
        let slot = &mut values[index];
        if option.collection {
            let mut list = match slot.take() {
                Some(Value::List(list)) => list,
                _ => default_elements(option.default),
            };
            if value == "null" {
                list.clear();
            } else {
                let element = render(option.kind, value);
                if option.order == Order::List || !list.contains(&element) {
                    list.push(element);
                }
            }
            *slot = Some(Value::List(list));
        } else if value == "null" {
            *slot = Some(Value::Scalar(None));
        } else {
            *slot = Some(Value::Scalar(Some(render(option.kind, value))));
        }
    }

    (positional, values)
}

pub fn legacy_command_line(options: &ToolOptions, args: &[String]) -> String {
    let (positional, values) = parse(options, args);
    let mut line = String::from(options.name);
    for value in &positional {
        line.push(' ');
        line.push_str(value);
    }
    for (option, value) in options.options.iter().zip(&values) {
        let Some(value) = value else { continue };
        line.push(' ');
        line.push_str(option.name);
        line.push('=');
        match value {
            Value::Scalar(Some(v)) => line.push_str(v),
            Value::Scalar(None) => line.push_str("null"),
            Value::List(list) => {
                line.push('[');
                line.push_str(&ordered(option.order, option.constants, list).join(", "));
                line.push(']');
            }
        }
    }
    line.push_str("   ");
    for (option, value) in options.options.iter().zip(&values) {
        if value.is_none() && option.default != "null" {
            line.push(' ');
            line.push_str(option.name);
            line.push('=');
            match heap_dependent_default(options.name, option.name) {
                Some(default) => line.push_str(&default),
                None => line.push_str(option.default),
            }
        }
    }
    line
}

/// `PairedReadSequence.getSizeInBytes()`: EstimateLibraryComplexity's guess at one pair's footprint,
/// `16 + 4 + 2*4 + 1 + 2*(24 + 8 + NUMBER_BASES_IN_READ) + 2 + 2*(24 + 8) + 8 + 4` with 150 bases.
const PAIRED_READ_SEQUENCE_BYTES: i64 = 471;

/// The `MAX_RECORDS_IN_RAM` default EstimateLibraryComplexity's constructor computes from the heap:
/// `Math.min((int) (maxMemory / sizeInBytes) / 2, Integer.MAX_VALUE - 32)`. The constructor runs
/// before any argument is parsed, so the barcode variant's size never applies.
pub fn library_complexity_records_in_ram(max_memory: i64) -> i64 {
    ((max_memory / PAIRED_READ_SEQUENCE_BYTES) as i32 as i64 / 2).min(i64::from(i32::MAX - 32))
}

/// A default no table can hold, because the reference computes it from its own JVM.
///
/// Only one tool here has one: EstimateLibraryComplexity's `MAX_RECORDS_IN_RAM`, which follows
/// `Runtime.maxMemory()`. Without a `-Xmx` the JVM's heap is a quarter of physical memory
/// (`MaxRAMPercentage=25`), which is what this takes `maxMemory` to be; a Picard run with another
/// heap prints another number.
fn heap_dependent_default(tool: &str, option: &str) -> Option<String> {
    if tool != "EstimateLibraryComplexity" || option != "MAX_RECORDS_IN_RAM" {
        return None;
    }
    // SAFETY: `sysconf` only reads system configuration.
    let (pages, page_size) = unsafe {
        (
            libc::sysconf(libc::_SC_PHYS_PAGES),
            libc::sysconf(libc::_SC_PAGESIZE),
        )
    };
    if pages <= 0 || page_size <= 0 {
        return None;
    }
    // `c_long` is 64 bits on every target the release builds.
    let max_memory: i64 = pages * page_size / 4;
    Some(library_complexity_records_in_ram(max_memory).to_string())
}

/// One scalar option of `tool` as `args` set it: `None` when unset, `Some(None)` when set to
/// `null`, otherwise the value as its type prints it.
pub fn argument(tool: &str, args: &[String], name: &str) -> Option<Option<String>> {
    let options = tool_options(tool)?;
    let index = options.options.iter().position(|o| o.name == name)?;
    match parse(options, args).1.swap_remove(index)? {
        Value::Scalar(value) => Some(value),
        Value::List(list) => Some(Some(format!("[{}]", list.join(", ")))),
    }
}

/// The command line of this process, as Picard would record it for `tool`.
///
/// The arguments are this process's own. A tool picard-rs has no declarations for is a bug in the
/// generated table, and it is reported as one rather than papered over.
pub fn command_line(tool: &str) -> String {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let options = tool_options(tool)
        .unwrap_or_else(|| panic!("no option declarations for {tool} in tool_options"));
    legacy_command_line(options, &args)
}

/// `java.util.Date.toString()`: `EEE MMM dd HH:mm:ss zzz yyyy`, in English whatever the locale, in
/// the default time zone, the zone written as its short name.
pub fn java_date(seconds_since_epoch: i64) -> String {
    const DAYS: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
    const MONTHS: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let local = local_time(seconds_since_epoch);
    format!(
        "{} {} {:02} {:02}:{:02}:{:02} {} {}",
        DAYS[local.weekday as usize],
        MONTHS[local.month as usize],
        local.day,
        local.hour,
        local.minute,
        local.second,
        local.zone,
        local.year
    )
}

struct LocalTime {
    year: i64,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: u32,
    weekday: u32,
    zone: String,
}

/// The C library's `localtime_r`, which reads `TZ` and `/etc/localtime` the way the JVM picks its
/// default zone, and names the zone by its abbreviation (`UTC`, `CEST`, `EST`).
fn local_time(seconds: i64) -> LocalTime {
    // SAFETY: `localtime_r` writes only into the `tm` it is given, and `tm_zone` points at a
    // static string the C library owns, read before anything else can change the zone.
    unsafe {
        let time: libc::time_t = seconds as libc::time_t;
        let mut tm: libc::tm = std::mem::zeroed();
        if libc::localtime_r(&time, &mut tm).is_null() {
            return utc(seconds);
        }
        let zone = if tm.tm_zone.is_null() {
            "UTC".to_string()
        } else {
            std::ffi::CStr::from_ptr(tm.tm_zone)
                .to_string_lossy()
                .into_owned()
        };
        LocalTime {
            year: i64::from(tm.tm_year) + 1900,
            month: tm.tm_mon as u32,
            day: tm.tm_mday as u32,
            hour: tm.tm_hour as u32,
            minute: tm.tm_min as u32,
            second: tm.tm_sec as u32,
            weekday: tm.tm_wday as u32,
            zone,
        }
    }
}

/// The same fields in UTC, for when the C library cannot say.
fn utc(seconds: i64) -> LocalTime {
    let days = seconds.div_euclid(86_400);
    let rem = seconds.rem_euclid(86_400);
    // Howard Hinnant's civil_from_days.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let month = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let year = yoe + era * 400 + i64::from(month <= 2);
    LocalTime {
        year,
        month: month - 1,
        day,
        hour: (rem / 3600) as u32,
        minute: (rem % 3600 / 60) as u32,
        second: (rem % 60) as u32,
        weekday: (days + 4).rem_euclid(7) as u32,
        zone: "UTC".to_string(),
    }
}

/// The start of the run: `instanceMain`'s `startDate`, taken once, so every metrics file one run
/// writes says the same.
pub fn started_on() -> &'static str {
    static STARTED: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    STARTED.get_or_init(|| {
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_secs() as i64);
        format!("Started on: {}", java_date(now))
    })
}

/// `getMetricsFile()`'s two headers for `tool`: the command line and the start of the run.
pub fn default_headers(tool: &str) -> [String; 2] {
    [command_line(tool), started_on().to_string()]
}

/// `CommandLineProgram.getMetricsFile()`: a metrics file opened with the default headers.
pub fn metrics_file(tool: &str) -> htsjdk_metrics::file::MetricsFile {
    let mut file = htsjdk_metrics::file::MetricsFile::new();
    for header in default_headers(tool) {
        file.add_header(&header);
    }
    file
}

/// `SAMFileHeader.PgIdGenerator`, one generator over a header's `@PG` ids.
///
/// The counter starts at the number of distinct ids already taken and moves on every id handed
/// out, a fresh one included, which is what spaces two collisions apart.
pub struct PgIdGenerator {
    taken: std::collections::HashSet<String>,
    counter: u32,
}

impl PgIdGenerator {
    pub fn new(header: &htsjdk_bam::header::SamHeader) -> Self {
        let taken: std::collections::HashSet<String> =
            header.programs.iter().map(|pg| pg.id.clone()).collect();
        let counter = taken.len() as u32;
        PgIdGenerator { taken, counter }
    }

    /// `getNonCollidingId`.
    pub fn non_colliding_id(&mut self, id: &str) -> String {
        if self.taken.insert(id.to_string()) {
            self.counter += 1;
            return id.to_string();
        }
        loop {
            let candidate = format!("{id}.{}", crate::merge_sam_files::base36(self.counter));
            self.counter += 1;
            if self.taken.insert(candidate.clone()) {
                return candidate;
            }
        }
    }
}

/// `CommandLineProgram.getPGRecord`: the tool's own `@PG`, id made unique against the header,
/// with `PN`, `CL` and `VN` in that order. The caller adds it to the header.
pub fn pg_record(
    header: &htsjdk_bam::header::SamHeader,
    tool: &str,
) -> htsjdk_bam::header::ProgramRecord {
    let id = PgIdGenerator::new(header).non_colliding_id(tool);
    let mut pg = htsjdk_bam::header::ProgramRecord::new(&id);
    pg.attributes.set("PN", tool);
    pg.attributes.set("CL", &command_line(tool));
    pg.attributes.set("VN", PICARD_VERSION);
    pg
}

/// `java.util.HashSet` iteration order over program ids, `null` included (it hashes to bucket 0).
pub fn hash_set_order(ids: Vec<Option<String>>) -> Vec<Option<String>> {
    let keyed = ids
        .into_iter()
        .map(|id| {
            let hash = id
                .as_deref()
                .map_or(0, crate::java_hash_map::string_hash_code);
            (hash, id)
        })
        .collect();
    crate::barcode_extractor::hash_set_order(keyed, 16)
}

/// The `PROGRAM_*` arguments of `AbstractMarkDuplicatesCommandLineProgram`.
pub struct ProgramGroup {
    /// `PROGRAM_RECORD_ID`; `None` when set to `null`, which turns the whole chain off.
    pub record_id: Option<String>,
    /// `PROGRAM_GROUP_VERSION`, `getVersion()` when unset.
    pub version: Option<String>,
    /// `PROGRAM_GROUP_COMMAND_LINE`, `getCommandLine()` when unset.
    pub command_line: Option<String>,
    /// `PROGRAM_GROUP_NAME`, the tool's simple name when unset; `None` when set to `null`.
    pub name: Option<String>,
}

impl ProgramGroup {
    /// The four arguments as this process's command line sets them for `tool`.
    pub fn of(tool: &str) -> Self {
        let args: Vec<String> = std::env::args().skip(1).collect();
        let get = |name: &str, default: Option<String>| match argument(tool, &args, name) {
            None => default,
            Some(value) => value,
        };
        ProgramGroup {
            record_id: get("PROGRAM_RECORD_ID", Some("MarkDuplicates".to_string())),
            version: get("PROGRAM_GROUP_VERSION", None),
            command_line: get("PROGRAM_GROUP_COMMAND_LINE", None),
            name: get("PROGRAM_GROUP_NAME", Some(tool.to_string())),
        }
    }
}

/// `AbstractMarkDuplicatesCommandLineProgram.getChainedPgIds`: one new `@PG` per program id the
/// run saw, each chained to it by `PP`, appended to `header` in `seen`'s order (the order of the
/// `HashSet` Picard keeps them in, see `hash_set_order`). Returns the map from a read's old `PG`
/// (or none) to its new one, or `None` when `PROGRAM_RECORD_ID` is `null`.
pub fn chained_pg_ids(
    header: &mut htsjdk_bam::header::SamHeader,
    seen: &[Option<String>],
    group: &ProgramGroup,
    tool: &str,
) -> Option<Vec<(Option<String>, String)>> {
    let record_id = group.record_id.as_ref()?;
    let mut generator = PgIdGenerator::new(header);
    let version = group
        .version
        .clone()
        .unwrap_or_else(|| PICARD_VERSION.to_string());
    let command_line = group
        .command_line
        .clone()
        .unwrap_or_else(|| command_line(tool));
    let mut chained = Vec::new();
    for existing in seen {
        let id = generator.non_colliding_id(record_id);
        let mut pg = htsjdk_bam::header::ProgramRecord::new(&id);
        pg.attributes.set("VN", &version);
        pg.attributes.set("CL", &command_line);
        if let Some(name) = &group.name {
            pg.attributes.set("PN", name);
        }
        if let Some(previous) = existing {
            pg.attributes.set("PP", previous);
        }
        header.programs.push(pg);
        chained.push((existing.clone(), id));
    }
    Some(chained)
}

/// The new id a read whose `PG` is `old` is given: `chainedPgIds.get(old)`, which is `null` for an
/// id the chain does not hold.
pub fn chained_id<'a>(
    chained: &'a [(Option<String>, String)],
    old: Option<&str>,
) -> Option<&'a str> {
    chained
        .iter()
        .find(|(from, _)| from.as_deref() == old)
        .map(|(_, to)| to.as_str())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_file_is_normalized_as_java_io_file_does() {
        assert_eq!(normalize_file("dir//x/"), "dir/x");
        assert_eq!(normalize_file("/abs//x/"), "/abs/x");
        assert_eq!(normalize_file("/"), "/");
        assert_eq!(normalize_file("./b"), "./b");
    }

    #[test]
    fn numbers_print_as_their_java_type_prints_them() {
        assert_eq!(render(Kind::Integer, "+0042"), "42");
        assert_eq!(render(Kind::Integer, "-3"), "-3");
        assert_eq!(render(Kind::Double, "0.250"), "0.25");
        assert_eq!(render(Kind::Double, "1E7"), "1.0E7");
        assert_eq!(render(Kind::Double, "1e-3"), "0.001");
        assert_eq!(render(Kind::Boolean, "TRUE"), "true");
        assert_eq!(render(Kind::Boolean, "T"), "false");
    }

    #[test]
    fn the_heap_default_is_what_a_two_gigabyte_heap_gives_picard() {
        // The header-lines dump runs with -Xmx2g, and EstimateLibraryComplexity's default there
        // is the number this has to reproduce for that heap.
        assert_eq!(
            library_complexity_records_in_ram(2 * 1024 * 1024 * 1024),
            2_279_706
        );
    }

    #[test]
    fn the_epoch_is_a_thursday_in_utc() {
        let t = utc(0);
        assert_eq!((t.year, t.month, t.day, t.weekday), (1970, 0, 1, 4));
        let t = utc(1_792_168_424);
        assert_eq!((t.year, t.month, t.day, t.hour), (2026, 9, 16, 16));
    }

    #[test]
    fn a_header_collision_takes_the_counter_after_the_taken_ids() {
        let header =
            htsjdk_bam::sam_file::read_sam("@HD\tVN:1.6\n@PG\tID:bwa\n@PG\tID:MarkDuplicates\n")
                .unwrap()
                .0;
        let mut generator = PgIdGenerator::new(&header);
        assert_eq!(
            generator.non_colliding_id("MarkDuplicates"),
            "MarkDuplicates.2"
        );
        assert_eq!(
            generator.non_colliding_id("MarkDuplicates"),
            "MarkDuplicates.3"
        );
        assert_eq!(generator.non_colliding_id("Fresh"), "Fresh");
        assert_eq!(generator.non_colliding_id("Fresh"), "Fresh.5");
    }

    #[test]
    fn null_comes_first_in_a_hash_set() {
        let order = hash_set_order(vec![
            Some("bwa".to_string()),
            Some("samtools".to_string()),
            None,
        ]);
        assert_eq!(
            order,
            vec![None, Some("samtools".to_string()), Some("bwa".to_string())]
        );
    }
}
