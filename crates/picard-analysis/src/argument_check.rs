//! What Picard does with a command line before the tool runs: refuse an argument the tool does not
//! declare, and print the usage for `-h`, `--help` and `--version`.
//!
//! Ported from Picard 3.5.0 (`CommandLineProgram.instanceMain`, `getCommandLineParserForArgs`,
//! `CommandLineSyntaxTranslater`) and Barclay 5.0.0 (`LegacyCommandLineArgumentParser`,
//! `CommandLineArgumentParser`), as far as the refusals and the help go.
//!
//! Picard picks a parser per command line. A word with `=` that does not start with `-` makes it
//! legacy syntax, which `LegacyCommandLineArgumentParser` reads after htsjdk's `Log` prints a
//! banner saying the syntax is changing; any `-` word makes it Barclay syntax, with a warning when
//! both kinds are present. The legacy parser upper-cases a name before it looks it up, so
//! `input=x.bam` is `INPUT`, and refuses the first name it does not hold, or the first bare word
//! of a tool without positional arguments, with the tool's whole usage after the message.
//!
//! Nothing here is written from the documentation. The names each parser accepts are dumped from
//! the parsers' own tables (`tools/commandline-conformance/CommandLineDeclarationsDump.java`), and
//! every text a refusal or a help request prints is the reference's own output for a probe command
//! line (`CommandLineDump.java`), with only the user's words put in. Both goldens are embedded, so
//! a binary prints the bytes Picard 3.5.0 printed, and the second one is decompressed only when a
//! command line is refused.
//!
//! The port reads `NAME=value` only. A Barclay-syntax argument the tool declares, `--INPUT x`, is
//! left to the binary, which does not read it.

use std::collections::HashSet;
use std::io::Read;
use std::sync::{Mutex, OnceLock};

static DECLARATIONS: &[u8] = include_bytes!("../tests/data/command_line_declarations.txt.gz");
static PROBES: &[u8] = include_bytes!("../tests/data/command_line.txt.gz");

/// Names the legacy parser refuses that the covering arrays pass in `NAME=value` form all the same.
/// The flow arguments of `CollectQualityYieldMetricsFlow` are in a collection only Barclay's parser
/// sees, and the runner strips the dashes off every argument it hands a port, so refusing them
/// would refuse rows the reference ran.
const ARRAY_ONLY: &[(&str, &str)] = &[
    (
        "CollectQualityYieldMetricsFlow",
        "FLOW-FILL-EMPTY-BINS-VALUE",
    ),
    ("CollectQualityYieldMetricsFlow", "FLOW-IGNORE-T0-TAG"),
];

const SYNTAX_TRANSITION_URL: &str =
    "https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)";

/// What a refused or answered command line prints, and the status the process exits with.
#[derive(Debug, PartialEq, Eq)]
pub struct Answer {
    pub exit: i32,
    pub stdout: String,
    pub stderr: String,
}

/// The names one tool's two parsers accept.
pub struct Declared {
    /// The legacy parser's `optionMap` keys, upper-cased as it holds them.
    pub legacy: HashSet<String>,
    /// Barclay's full and short names.
    pub barclay: HashSet<String>,
    /// Whether the tool takes positional words.
    pub positional: bool,
}

/// Refuses the process's command line the way Picard would, or prints the help it asked for, and
/// exits; returns when the tool should run. Every binary calls this first.
pub fn check(tool: &str) {
    crate::args::set_tool(tool);
    let args: Vec<String> = crate::args::process().into_iter().skip(1).collect();
    if let Some(answer) = answer(tool, &args) {
        print!("{}", answer.stdout);
        eprint!("{}", answer.stderr);
        std::process::exit(answer.exit);
    }
}

/// An embedded golden, inflated once per process.
fn gunzip(bytes: &'static [u8]) -> &'static str {
    static INFLATED: OnceLock<Mutex<Vec<(usize, &'static str)>>> = OnceLock::new();
    let mut cache = INFLATED.get_or_init(Default::default).lock().unwrap();
    let key = bytes.as_ptr() as usize;
    if let Some((_, text)) = cache.iter().find(|(k, _)| *k == key) {
        return text;
    }
    let mut text = String::new();
    flate2::read::GzDecoder::new(bytes)
        .read_to_string(&mut text)
        .expect("an embedded golden is gzip");
    let text: &'static str = Box::leak(text.into_boxed_str());
    cache.push((key, text));
    text
}

/// What `CommandLineDeclarationsDump` recorded for `tool`, or `None` for a name Picard does not have.
pub fn declared(tool: &str) -> Option<Declared> {
    let mut declared = Declared {
        legacy: HashSet::new(),
        barclay: HashSet::new(),
        positional: false,
    };
    let mut found = false;
    for line in gunzip(DECLARATIONS).lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 || fields[1] != tool {
            continue;
        }
        found = true;
        match fields[0] {
            "legacy" => drop(declared.legacy.insert(fields[2].to_string())),
            "barclay" => drop(declared.barclay.insert(fields[2].to_string())),
            "positional" => declared.positional = true,
            _ => {}
        }
    }
    for (owner, name) in ARRAY_ONLY {
        if *owner == tool {
            declared.legacy.insert(name.to_string());
        }
    }
    found.then_some(declared)
}

/// The full name the legacy parser resolves `name` to for `tool`: a short name or the full name
/// itself, upper-cased as the parser upper-cases it before the lookup.
pub fn full_name(tool: &str, name: &str) -> Option<&'static str> {
    gunzip(DECLARATIONS).lines().find_map(|line| {
        let fields: Vec<&str> = line.splitn(4, '\t').collect();
        (fields.len() == 4 && fields[0] == "legacy" && fields[1] == tool && fields[2] == name)
            .then(|| fields[3])
    })
}

/// One probe `CommandLineDump` ran, as the reference answered it.
pub fn probe(tool: &str, name: &str) -> Option<Answer> {
    gunzip(PROBES).lines().find_map(|line| {
        let fields: Vec<&str> = line.split('\t').collect();
        (fields.len() == 6 && fields[0] == "run" && fields[1] == tool && fields[2] == name).then(
            || Answer {
                exit: fields[3].parse().expect("an exit status"),
                stdout: unescape(fields[4]),
                stderr: unescape(fields[5]),
            },
        )
    })
}

/// The dump's escaping undone: `\\`, `\n`, `\t` and `\r`.
pub fn unescape(field: &str) -> String {
    let mut text = String::with_capacity(field.len());
    let mut chars = field.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            text.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => text.push('\n'),
            Some('t') => text.push('\t'),
            Some('r') => text.push('\r'),
            Some(other) => text.push(other),
            None => text.push('\\'),
        }
    }
    text
}

/// `CommandLineSyntaxTranslater.isLegacyPicardStyle`, without its warning: some word is
/// `NAME=value` and none starts with `-`.
fn legacy_words(args: &[String]) -> bool {
    args.iter().any(|a| !a.starts_with('-') && a.contains('='))
}

/// The banner `instanceMain` logs before a legacy command line is parsed, the command line put in
/// the new syntax by `convertPicardStyleToPosixStyle`.
pub fn banner(tool: &str, args: &[String]) -> String {
    let syntax: Vec<String> = args
        .iter()
        .map(|arg| match arg.split_once('=') {
            None => arg.clone(),
            Some(("EXTRA_ARGUMENT", value)) => {
                format!(
                    "-EXTRA_ARGUMENT {}",
                    value.replace("::", "::--").replace('=', " ")
                )
            }
            Some((name, value)) => format!("-{name} {value}"),
        })
        .collect();
    let message = [
        "",
        "",
        "********** NOTE: Picard's command line syntax is changing.",
        "**********",
        "********** For more information, please see:",
        "********** ",
        SYNTAX_TRANSITION_URL,
        "**********",
        "********** The command line looks like this in the new syntax:",
        "**********",
        &format!("**********    {tool} {}", syntax.join(" ")),
        "**********",
        "",
        "",
    ]
    .join("\n");
    format!("INFO\t1970-01-01 00:00:00\t{tool}\t{message}\n")
}

/// A recorded legacy probe's answer for `args`: the banner rebuilt for them, the message replaced,
/// and everything else as the reference printed it. Between the banner and the message there can be
/// log4j's warnings from building the parser (`Hidden annotation is not honored for
/// --DROP_MISSING_CONTIGS`), which are kept.
fn reworded(tool: &str, name: &str, message: &str, args: &[String], replacement: &str) -> Answer {
    let recorded = probe(tool, name).expect("every tool's probes are in the golden");
    let (before, after) = recorded
        .stderr
        .split_once(message)
        .expect("the probe's message is in what the reference printed");
    let words: Vec<String> = probe_args(name).iter().map(|w| w.to_string()).collect();
    let between = before
        .strip_prefix(banner(tool, &words).as_str())
        .expect("the probe's banner is the one this module builds");
    Answer {
        exit: recorded.exit,
        stdout: recorded.stdout,
        stderr: format!("{}{between}{replacement}{after}", banner(tool, args)),
    }
}

/// The words `CommandLineDump` put after the tool's name for each probe.
pub fn probe_args(name: &str) -> &'static [&'static str] {
    match name {
        "legacy_unknown" => &["PICARD_RS_PROBE=1"],
        "legacy_lowercase" => &["picard_rs_probe=1"],
        "legacy_word" => &["VERBOSITY=INFO", "picard_rs_probe"],
        "word" => &["picard_rs_probe"],
        "help_short" => &["-h"],
        "help_long" => &["--help"],
        "version" => &["--version"],
        "barclay_unknown" => &["--PICARD_RS_PROBE", "1"],
        "barclay_embedded" => &["--PICARD_RS_PROBE=1"],
        "mixed_help" => &["PICARD_RS_PROBE=1", "-h"],
        _ => &[],
    }
}

/// What Picard 3.5.0 does with `args` for `tool` before the tool runs, or `None` when it goes on to
/// run the tool.
pub fn answer(tool: &str, args: &[String]) -> Option<Answer> {
    let declared = declared(tool)?;
    let dashed = args.iter().any(|a| a.starts_with('-'));
    if legacy_words(args) && !dashed {
        return legacy(tool, &declared, args);
    }
    if !dashed {
        // No `-` word and no `NAME=value`: Barclay's parser, and only bare words are left for it to
        // refuse.
        let word = args.first().filter(|_| !declared.positional)?;
        let recorded = probe(tool, "word")?;
        let message = "picard_rs_probe";
        let (before, after) = recorded.stderr.rsplit_once(message)?;
        return Some(Answer {
            exit: recorded.exit,
            stdout: recorded.stdout,
            stderr: format!("{before}{word}{after}"),
        });
    }
    barclay(tool, &declared, args)
}

fn legacy(tool: &str, declared: &Declared, args: &[String]) -> Option<Answer> {
    for arg in args {
        match arg.split_once('=') {
            Some((name, _)) => {
                let key = name.to_uppercase();
                // An options file is read by the parser itself, and the binary is left to it.
                if key == "OPTIONS_FILE" || declared.legacy.contains(&key) {
                    continue;
                }
                let message = "ERROR: Unrecognized option: PICARD_RS_PROBE";
                return Some(reworded(
                    tool,
                    "legacy_unknown",
                    message,
                    args,
                    &format!("ERROR: Unrecognized option: {key}"),
                ));
            }
            None if declared.positional => {}
            None => {
                let message = "ERROR: Invalid argument 'picard_rs_probe'.";
                return Some(reworded(
                    tool,
                    "legacy_word",
                    message,
                    args,
                    &format!("ERROR: Invalid argument '{arg}'."),
                ));
            }
        }
    }
    None
}

fn barclay(tool: &str, declared: &Declared, args: &[String]) -> Option<Answer> {
    // `isLegacyPicardStyle` runs twice on a mixed command line, in `instanceMain` and in
    // `getCommandLineParserForArgs`, and warns each time.
    let warning = if legacy_words(args) {
        let mixed = probe(tool, "mixed_help")?;
        let help = probe(tool, "help_short")?;
        mixed.stderr.strip_suffix(help.stderr.as_str())?.to_string()
    } else {
        String::new()
    };
    let mut help = None;
    let mut version = false;
    for arg in args.iter().filter(|a| a.starts_with('-')) {
        let word = arg.trim_start_matches('-');
        let name = word.split_once('=').map_or(word, |(name, _)| name);
        if !declared.barclay.contains(name) {
            let (probe_name, message, replacement) = if word.contains('=') {
                (
                    "barclay_embedded",
                    "Can't parse option name containing an embedded '=' (PICARD_RS_PROBE=1)",
                    format!("Can't parse option name containing an embedded '=' ({word})"),
                )
            } else {
                (
                    "barclay_unknown",
                    "PICARD_RS_PROBE is not a recognized option",
                    format!("{word} is not a recognized option"),
                )
            };
            let recorded = probe(tool, probe_name)?;
            let (before, after) = recorded.stderr.rsplit_once(message)?;
            return Some(Answer {
                exit: recorded.exit,
                stdout: recorded.stdout,
                stderr: format!("{warning}{before}{replacement}{after}"),
            });
        }
        match name {
            "h" if help.is_none() => help = Some("help_short"),
            "help" if help.is_none() => help = Some("help_long"),
            "version" => version = true,
            _ => {}
        }
    }
    let shown = help.or(version.then_some("version"))?;
    let recorded = probe(tool, shown)?;
    Some(Answer {
        exit: recorded.exit,
        stdout: recorded.stdout,
        stderr: format!("{warning}{}", recorded.stderr),
    })
}
