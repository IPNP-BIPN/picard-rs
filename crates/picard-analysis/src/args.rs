//! The `NAME=value` words a tool reads, under the rules Picard's legacy parser applies to them.
//!
//! A command line written `NAME=value` goes through Barclay's `LegacyCommandLineArgumentParser`,
//! which builds each value with the option type's `String` constructor. These helpers are the
//! ones every binary shares, so one rule holds for every tool.

static TOOL: std::sync::OnceLock<String> = std::sync::OnceLock::new();

/// The tool whose names [`process`] resolves, set by `argument_check::check`.
pub fn set_tool(tool: &str) {
    let _ = TOOL.set(tool.to_string());
}

/// The process's arguments as a binary reads them: each `NAME=value` word under the full name
/// Picard's legacy parser resolves NAME to, in any case and short or full, so a tool that asks for
/// `PROBABILITY=` reads `P=0.5`, and one that asks for `INPUT=` reads `input=x`. A word the tool
/// does not declare, or a Barclay `-` word, is left as it is. A `NAME=` with nothing after the `=`
/// takes the next word as its value, as the legacy parser does. The first word is the program.
pub fn process() -> Vec<String> {
    let Some(tool) = TOOL.get() else {
        return std::env::args().collect();
    };
    joined(std::env::args().collect())
        .into_iter()
        .enumerate()
        .map(|(i, word)| if i == 0 { word } else { canonical(tool, &word) })
        .collect()
}

/// Each `NAME=` word, after the program, joined with the word after it: `I= x.bam` is `I=x.bam`.
fn joined(raw: Vec<String>) -> Vec<String> {
    let mut words = Vec::with_capacity(raw.len());
    let mut raw = raw.into_iter();
    words.extend(raw.next());
    while let Some(word) = raw.next() {
        let open =
            !word.starts_with('-') && word.split_once('=').is_some_and(|(_, v)| v.is_empty());
        match if open { raw.next() } else { None } {
            Some(value) => words.push(word + &value),
            None => words.push(word),
        }
    }
    words
}

/// One word under its full name, when `tool` declares the name it starts with.
fn canonical(tool: &str, word: &str) -> String {
    if word.starts_with('-') {
        return word.to_string();
    }
    match word.split_once('=') {
        Some((name, value)) => match crate::argument_check::full_name(tool, &name.to_uppercase()) {
            Some(full) => format!("{full}={value}"),
            None => word.to_string(),
        },
        None => word.to_string(),
    }
}

/// The value of the first word that starts with `key`, which carries its `=` (`"INPUT="`).
/// `null` is no value, as Barclay sets the field to `null`.
pub fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key))
        .filter(|v| *v != "null")
        .map(str::to_string)
}

/// A `Boolean` option's value as `new Boolean(value)` builds it, which is what the legacy parser's
/// `constructFromString` calls: `true` in any case is true, and every other word is false, without
/// an error. `CREATE_INDEX=TRUE` is an index, as it is in Picard.
pub fn is_true(value: &str) -> bool {
    value.eq_ignore_ascii_case("true")
}

/// Every value of a list option without a default, in command-line order, under its long and
/// short names (`"INPUT="` and `"I="`). `null` empties the list, as Barclay clears it.
pub fn all(args: &[String], long: &str, short: &str) -> Vec<String> {
    let mut values = Vec::new();
    for value in args
        .iter()
        .filter_map(|a| a.strip_prefix(long).or_else(|| a.strip_prefix(short)))
    {
        if value == "null" {
            values.clear();
        } else {
            values.push(value.to_string());
        }
    }
    values
}

/// The files a tool reads from `INPUT` (or `I`), a list in Picard, at least one of them.
pub fn inputs(args: &[String]) -> Result<Vec<String>, String> {
    let inputs = all(args, "INPUT=", "I=");
    if inputs.is_empty() {
        return Err("INPUT= is required".to_string());
    }
    Ok(inputs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_short_or_lower_case_name_reads_as_the_full_one() {
        assert_eq!(canonical("DownsampleSam", "P=0.5"), "PROBABILITY=0.5");
        assert_eq!(
            canonical("DownsampleSam", "probability=0.5"),
            "PROBABILITY=0.5"
        );
        assert_eq!(canonical("DownsampleSam", "I=a=b.bam"), "INPUT=a=b.bam");
        assert_eq!(
            canonical("DownsampleSam", "NOT_AN_OPTION=1"),
            "NOT_AN_OPTION=1"
        );
        assert_eq!(canonical("DownsampleSam", "--P"), "--P");
    }

    fn words(line: &str) -> Vec<String> {
        line.split(' ').map(str::to_string).collect()
    }

    #[test]
    fn a_boolean_is_true_in_any_case_and_false_otherwise() {
        assert!(is_true("true"));
        assert!(is_true("TRUE"));
        assert!(is_true("True"));
        assert!(!is_true("false"));
        assert!(!is_true("T"));
        assert!(!is_true("yes"));
        assert!(!is_true(""));
    }

    #[test]
    fn the_first_word_wins() {
        let args = words("tool CREATE_INDEX=true CREATE_INDEX=false");
        assert_eq!(arg(&args, "CREATE_INDEX=").as_deref(), Some("true"));
        assert_eq!(arg(&args, "OUTPUT="), None);
    }

    #[test]
    fn null_is_no_value() {
        let args = words("tool R=null I=a.bam I=null I=b.bam");
        assert_eq!(arg(&args, "R="), None);
        assert_eq!(all(&args, "INPUT=", "I="), vec!["b.bam".to_string()]);
    }

    #[test]
    fn an_empty_value_takes_the_next_word() {
        assert_eq!(
            joined(words("tool I= in.bam O=x.bam -X= y A=")),
            words("tool I=in.bam O=x.bam -X= y A=")
        );
    }

    #[test]
    fn every_input_is_read_in_order() {
        assert_eq!(
            inputs(&words("tool I=a.bam")),
            Ok(vec!["a.bam".to_string()])
        );
        assert_eq!(
            inputs(&words("tool I=a.bam INPUT=b.bam")),
            Ok(vec!["a.bam".to_string(), "b.bam".to_string()])
        );
        assert!(inputs(&words("tool O=x.bam")).is_err());
    }
}
