//! `Float.parseFloat`, `Double.parseDouble` and `Integer.parseInt`, with the messages their
//! `NumberFormatException`s carry.
//!
//! The double is htsjdk-rs's `parse_java_double`, which already takes Java's spellings (`1.5f`,
//! surrounding control characters, `Infinity`, hexadecimal) and refuses Rust's (`inf`, `nan`). The
//! float goes through the same gate and is then read AS a float: Java parses straight to the
//! nearest `float`, and rounding through a double first can land on the other neighbour.

use htsjdk_vcf::genotype_likelihoods::parse_java_double;

/// `Double.parseDouble`.
pub fn parse_double(text: &str) -> Result<f64, String> {
    parse_java_double(text).ok_or_else(|| number_format_message(text))
}

/// `Float.parseFloat`.
pub fn parse_float(text: &str) -> Result<f32, String> {
    let value = parse_double(text)?;
    if !value.is_finite() {
        return Ok(value as f32);
    }
    let trimmed = text.trim_matches(|c: char| c <= ' ');
    let body = trimmed.trim_end_matches(['f', 'F', 'd', 'D']);
    Ok(body.parse::<f32>().unwrap_or(value as f32))
}

/// `Integer.parseInt(text)` in radix ten: an optional sign and ASCII digits, nothing else.
pub fn parse_int(text: &str) -> Result<i32, String> {
    let digits = text
        .strip_prefix(['-', '+'])
        .filter(|rest| !rest.is_empty())
        .unwrap_or(text);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!("For input string: \"{text}\""));
    }
    text.parse::<i32>()
        .map_err(|_| format!("For input string: \"{text}\""))
}

/// `FloatingDecimal.readJavaFormatString`'s message: an empty string is named as such, anything
/// else is quoted.
pub fn number_format_message(text: &str) -> String {
    if text.trim_matches(|c: char| c <= ' ').is_empty() {
        "empty String".to_string()
    } else {
        format!("For input string: \"{text}\"")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_spellings() {
        assert_eq!(parse_double("1.5d"), Ok(1.5));
        assert!(parse_double("inf").is_err());
        assert_eq!(parse_float(" 0.1f"), Ok(0.1f32));
        assert_eq!(parse_int("+7"), Ok(7));
        assert_eq!(
            parse_int("1.0"),
            Err("For input string: \"1.0\"".to_string())
        );
        assert_eq!(parse_double(""), Err("empty String".to_string()));
    }
}
