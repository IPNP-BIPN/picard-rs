//! `SimpleDateFormat.parse` for the two patterns an array VCF's header dates are read with, and
//! `Iso8601Date.toString` for writing them back.
//!
//! `CollectArraysVariantCallingMetrics` reads `autocallDate` as `MM/dd/yyyy HH:mm` in the JVM's
//! zone (UTC in the oracle image) and `imagingDate` as `MM/dd/yyyy hh:mm:ss a` in
//! `America/New_York`, and writes both through `yyyy-MM-dd'T'HH:mm:ssZ` in the JVM's zone. What
//! the parse does is the default `SimpleDateFormat`'s, which is LENIENT:
//!
//! * a numeric field takes as many digits as are there (`4/2/2019 0:07` is a date), after any
//!   spaces or tabs, and an out-of-range value rolls over into the next field rather than failing;
//! * a literal in the pattern must be the next character exactly;
//! * `hh` is the 1-12 hour, and 12 means 0 before the AM/PM marker is added; the marker matches
//!   `AM` or `PM` in any case;
//! * text after the last field is ignored: `DateFormat.parse(String)` fails only when nothing at
//!   all was consumed.
//!
//! The New York offset follows the US rules since 1987: daylight time from the first Sunday of
//! April (the second of March from 2007) at 02:00 to the last Sunday of October (the first of
//! November from 2007) at 02:00. A wall time in the spring gap or the autumn overlap is resolved
//! as standard time, which is not proven against the reference; the corpus has none.

/// The fields `SimpleDateFormat` fills, before the calendar resolves them.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
struct Fields {
    year: i64,
    month: i64,
    day: i64,
    hour: i64,
    minute: i64,
    second: i64,
}

enum Token {
    Number(char),
    Literal(char),
    AmPm,
}

fn compile(pattern: &str) -> Vec<Token> {
    let chars: Vec<char> = pattern.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_ascii_alphabetic() {
            while i < chars.len() && chars[i] == c {
                i += 1;
            }
            tokens.push(if c == 'a' {
                Token::AmPm
            } else {
                Token::Number(c)
            });
        } else {
            tokens.push(Token::Literal(c));
            i += 1;
        }
    }
    tokens
}

/// The pattern's fields, or `None` for a `ParseException`.
fn parse_fields(pattern: &str, text: &str) -> Option<Fields> {
    let chars: Vec<char> = text.chars().collect();
    let mut at = 0usize;
    let mut fields = Fields {
        year: 1970,
        month: 1,
        day: 1,
        ..Fields::default()
    };
    let mut hour12: Option<i64> = None;
    let mut pm = false;
    for token in compile(pattern) {
        match token {
            Token::Literal(c) => {
                if chars.get(at) != Some(&c) {
                    return None;
                }
                at += 1;
            }
            Token::Number(letter) => {
                while matches!(chars.get(at), Some(' ') | Some('\t')) {
                    at += 1;
                }
                let negative = chars.get(at) == Some(&'-');
                let begin = if negative { at + 1 } else { at };
                let mut end = begin;
                while chars.get(end).is_some_and(char::is_ascii_digit) {
                    end += 1;
                }
                if end == begin {
                    return None;
                }
                let digits: String = chars[begin..end].iter().collect();
                let mut value: i64 = digits.parse().ok()?;
                if negative {
                    value = -value;
                }
                at = end;
                match letter {
                    'y' => fields.year = value,
                    'M' => fields.month = value,
                    'd' => fields.day = value,
                    'H' => fields.hour = value,
                    'h' => hour12 = Some(if value == 12 { 0 } else { value }),
                    'm' => fields.minute = value,
                    's' => fields.second = value,
                    _ => return None,
                }
            }
            Token::AmPm => {
                let rest: String = chars[at.min(chars.len())..].iter().collect();
                let upper = rest.to_ascii_uppercase();
                if upper.starts_with("AM") {
                    pm = false;
                } else if upper.starts_with("PM") {
                    pm = true;
                } else {
                    return None;
                }
                at += 2;
            }
        }
    }
    if let Some(hour) = hour12 {
        fields.hour = hour + if pm { 12 } else { 0 };
    }
    Some(fields)
}

/// Days from 1970-01-01 to a proleptic Gregorian date (Howard Hinnant's `days_from_civil`).
fn days_from_civil(year: i64, month: i64, day: i64) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = (month + 9) % 12;
    let doy = (153 * mp + 2) / 5 + day - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// `civil_from_days`, the inverse.
fn civil_from_days(days: i64) -> (i64, i64, i64) {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    (yoe + era * 400 + i64::from(month <= 2), month, day)
}

/// The lenient calendar: the month rolls into the year, everything else is plain arithmetic on
/// seconds, so a day 32 or an hour 25 is simply later.
fn local_seconds(fields: &Fields) -> i64 {
    let months = fields.month - 1;
    let year = fields.year + months.div_euclid(12);
    let month = months.rem_euclid(12) + 1;
    let days = days_from_civil(year, month, 1) + fields.day - 1;
    days * 86_400 + fields.hour * 3_600 + fields.minute * 60 + fields.second
}

/// The day of the month of the nth (1-based) Sunday of a month, or the last one with `n == 0`.
fn sunday(year: i64, month: i64, n: i64) -> i64 {
    let first = days_from_civil(year, month, 1);
    // 1970-01-01 was a Thursday; Sunday is weekday 0.
    let weekday = (first + 4).rem_euclid(7);
    let first_sunday = 1 + (7 - weekday) % 7;
    if n > 0 {
        return first_sunday + 7 * (n - 1);
    }
    let next = if month == 12 {
        days_from_civil(year + 1, 1, 1)
    } else {
        days_from_civil(year, month + 1, 1)
    };
    let length = next - first;
    let mut day = first_sunday;
    while day + 7 <= length {
        day += 7;
    }
    day
}

/// The UTC offset of a New York wall time, in seconds.
fn new_york_offset(local: i64) -> i64 {
    let (year, _, _) = civil_from_days(local.div_euclid(86_400));
    let (start, end) = if year >= 2007 {
        ((3, sunday(year, 3, 2)), (11, sunday(year, 11, 1)))
    } else if year >= 1987 {
        ((4, sunday(year, 4, 1)), (10, sunday(year, 10, 0)))
    } else {
        return -5 * 3_600;
    };
    let begins = days_from_civil(year, start.0, start.1) * 86_400 + 2 * 3_600;
    // 02:00 daylight time, which is 01:00 standard: the overlap reads as standard time.
    let ends = days_from_civil(year, end.0, end.1) * 86_400 + 3_600;
    if local >= begins + 3_600 && local < ends {
        -4 * 3_600
    } else {
        -5 * 3_600
    }
}

/// The zone a pattern is read in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zone {
    Utc,
    NewYork,
}

/// `new SimpleDateFormat(pattern)` (with its zone set) `.parse(text)`, as seconds since the epoch.
pub fn parse(pattern: &str, zone: Zone, text: &str) -> Option<i64> {
    let fields = parse_fields(pattern, text)?;
    let local = local_seconds(&fields);
    Some(match zone {
        Zone::Utc => local,
        Zone::NewYork => local - new_york_offset(local),
    })
}

/// `Iso8601Date.toString` in UTC: `yyyy-MM-dd'T'HH:mm:ssZ`.
pub fn iso8601(seconds: i64) -> String {
    let (year, month, day) = civil_from_days(seconds.div_euclid(86_400));
    let rest = seconds.rem_euclid(86_400);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}+0000",
        rest / 3_600,
        rest % 3_600 / 60,
        rest % 60
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_oracle_s_dates() {
        let autocall = parse("MM/dd/yyyy HH:mm", Zone::Utc, "4/2/2019 0:07").unwrap();
        assert_eq!(iso8601(autocall), "2019-04-02T00:07:00+0000");
        let imaging = "MM/dd/yyyy hh:mm:ss a";
        let summer = parse(imaging, Zone::NewYork, "04/10/2019 3:12:31 PM").unwrap();
        assert_eq!(iso8601(summer), "2019-04-10T19:12:31+0000");
        let winter = parse(imaging, Zone::NewYork, "12/05/2019 12:05:00 AM").unwrap();
        assert_eq!(iso8601(winter), "2019-12-05T05:05:00+0000");
        assert!(parse(imaging, Zone::NewYork, "2019-12-05T00:05:00").is_none());
    }
}
