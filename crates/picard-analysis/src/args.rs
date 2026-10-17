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

/// The one file a tool reads from `INPUT` (or `I`), which Picard declares as a list. Reading
/// several inputs as one is not ported, so a second one is refused rather than silently dropped.
pub fn single_input(args: &[String]) -> Result<String, String> {
    let inputs: Vec<String> = args
        .iter()
        .filter_map(|a| a.strip_prefix("INPUT=").or_else(|| a.strip_prefix("I=")))
        .map(str::to_string)
        .collect();
    match inputs.as_slice() {
        [] => Err("INPUT= is required".to_string()),
        [one] => Ok(one.clone()),
        several => Err(format!(
            "more than one INPUT is not ported: {} were given",
            several.len()
        )),
    }
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
    fn a_second_input_is_refused() {
        assert_eq!(single_input(&words("tool I=a.bam")).as_deref(), Ok("a.bam"));
        assert_eq!(
            single_input(&words("tool INPUT=a.bam")).as_deref(),
            Ok("a.bam")
        );
        assert!(single_input(&words("tool O=x.bam")).is_err());
        let refused = single_input(&words("tool I=a.bam INPUT=b.bam")).unwrap_err();
        assert!(refused.contains("2 were given"), "{refused}");
    }
}
