//! The `NAME=value` words a tool reads, under the rules Picard's legacy parser applies to them.
//!
//! A command line written `NAME=value` goes through Barclay's `LegacyCommandLineArgumentParser`,
//! which builds each value with the option type's `String` constructor. These helpers are the
//! ones every binary shares, so one rule holds for every tool.

/// The value of the first word that starts with `key`, which carries its `=` (`"INPUT="`).
pub fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

/// A `Boolean` option's value as `new Boolean(value)` builds it, which is what the legacy parser's
/// `constructFromString` calls: `true` in any case is true, and every other word is false, without
/// an error. `CREATE_INDEX=TRUE` is an index, as it is in Picard.
pub fn is_true(value: &str) -> bool {
    value.eq_ignore_ascii_case("true")
}

/// Every value of a list option, in command-line order, under its long and short names (`"INPUT="`
/// and `"I="`).
pub fn all(args: &[String], long: &str, short: &str) -> Vec<String> {
    args.iter()
        .filter_map(|a| a.strip_prefix(long).or_else(|| a.strip_prefix(short)))
        .map(str::to_string)
        .collect()
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
