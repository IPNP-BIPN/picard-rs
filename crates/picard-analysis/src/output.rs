//! How every block-compressed file a tool writes is deflated: `COMPRESSION_LEVEL` and
//! `USE_JDK_DEFLATER`.
//!
//! Picard sets both once, in `CommandLineProgram.instanceMain`, before `doWork` runs
//! (Picard 3.5.0, `src/main/java/picard/cmdline/CommandLineProgram.java:221-222`):
//!
//! ```java
//! BlockCompressedOutputStream.setDefaultCompressionLevel(COMPRESSION_LEVEL);
//! IOUtil.setCompressionLevel(COMPRESSION_LEVEL);
//! ```
//!
//! and, further down the same method (`:240-242`):
//!
//! ```java
//! if (!USE_JDK_DEFLATER) {
//!     BlockCompressedOutputStream.setDefaultDeflaterFactory(new IntelDeflaterFactory());
//! }
//! ```
//!
//! The second line of the first quote is the level of htsjdk's plain gzip writer, which
//! [`crate::java_writer::gzip_writer`] takes as an argument; this file serves the BGZF writers.
//!
//! They are htsjdk statics, so every BAM, BGZF VCF, FASTQ and tabix index the tool opens afterwards
//! takes them, whichever class opens it. [`install`] is that call, made by each binary;
//! [`bam_writer`], [`bgzf_writer`] and [`current`] are what the writers read.
//!
//! `USE_JDK_DEFLATER=true` does not mean `java.util.zip`. It only stops Picard installing GKL, and
//! the factory left in place is htsjdk 5.0.0's own, which returns libdeflate whenever its native
//! library loads (`samjdk.use_libdeflate`, true by default, and jlibdeflate ships one for every
//! platform Picard runs on). `DeflaterFactory.java:56` and `Defaults.java:140`:
//!
//! ```java
//! if (Defaults.USE_LIBDEFLATE && isLibdeflateAvailable()) {
//! ```
//!
//! ```java
//! USE_LIBDEFLATE = getBooleanProperty("use_libdeflate", true);
//! ```
//!
//! A caller that never runs a binary, a test among them, keeps htsjdk-rs's defaults, the JDK
//! deflater at level five, which is what every committed golden was measured with.
//!
//! # Reading this file without reading Rust
//!
//! - `static SETTING: OnceLock<..>` is a global that can be set once; `get` returns `None`
//!   until it is. It stands for htsjdk's two static fields.
//! - `(u32, Deflater)` is a tuple, a two-field value with no class of its own.
//! - `Result<T, E>` is "a value or an error"; `match` on `Ok(v) if ...` takes the success case
//!   only when the condition after `if` (a guard) holds too.
//! - `W: Write` makes a function generic over any output stream, Java's `OutputStream`
//!   parameter.

use std::io::{self, Write};
use std::sync::OnceLock;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::writer::BamWriter;
use htsjdk_bgzf::{BgzfWriter, Deflater, DEFAULT_COMPRESSION_LEVEL};

/// `CommandLineProgram.COMPRESSION_LEVEL`'s default, which is htsjdk's
/// `Defaults.COMPRESSION_LEVEL`.
///
/// `public int COMPRESSION_LEVEL = Defaults.COMPRESSION_LEVEL;` (`CommandLineProgram.java:116`)
/// and `COMPRESSION_LEVEL = getIntProperty("compression_level", 5);` (`Defaults.java:125`).
const PICARD_COMPRESSION_LEVEL: u32 = 5;

/// The level and deflater [`install`] chose, unset in a process that never calls it.
static SETTING: OnceLock<(u32, Deflater)> = OnceLock::new();

/// Reads `COMPRESSION_LEVEL` and `USE_JDK_DEFLATER` (short name `use_jdk_deflater`) from the
/// process's command line and makes them the level and deflater of every writer opened after it.
///
/// Without `USE_JDK_DEFLATER=true` the deflater is Intel's GKL, which is Picard's default and
/// therefore what an unchanged Picard command line writes; with it, libdeflate.
pub fn install() {
    // `skip(1)` drops the program's own name. A second `set` would be refused; `let _ =`
    // discards that result, so the first call wins.
    let args: Vec<String> = crate::args::process().into_iter().skip(1).collect();
    let _ = SETTING.set(setting(&args));
    crate::bam_stream::budget_threads();
}

/// The setting a command line asks for, the last value of each argument winning as Barclay's does.
///
/// Both `NAME=value` and `--NAME value` are read, with `USE_JDK_DEFLATER`'s short name
/// (`shortName = "use_jdk_deflater"`, `CommandLineProgram.java:141`). Every other argument is
/// skipped. A level of `null` is the default; any other level that is not an integer from 0 to 9
/// ends the tool with the message below.
pub fn setting(args: &[String]) -> (u32, Deflater) {
    let mut level = PICARD_COMPRESSION_LEVEL;
    let mut jdk = false;
    let mut words = args.iter();
    // `while let Some(word) = words.next()` loops until the iterator is exhausted.
    while let Some(word) = words.next() {
        let bare = word.trim_start_matches('-');
        let (name, value) = match bare.split_once('=') {
            Some((name, value)) => (name, value.to_string()),
            // `--NAME value`, Barclay's own form; a flag given alone is `true`.
            // `words.clone().next()` looks at the next word without consuming it; it is consumed
            // below only when it is the value.
            None if word.starts_with("--") || word.starts_with('-') => {
                let value = words
                    .clone()
                    .next()
                    .filter(|next| !next.starts_with('-'))
                    .cloned();
                if value.is_some() {
                    words.next();
                }
                (bare, value.unwrap_or_else(|| "true".to_string()))
            }
            None => continue,
        };
        match name {
            "COMPRESSION_LEVEL" => {
                level = match value.parse::<u32>() {
                    Ok(v) if v <= 9 => v,
                    _ if value == "null" => PICARD_COMPRESSION_LEVEL,
                    _ => crate::metrics_cli::fail(&format!(
                        "Argument 'COMPRESSION_LEVEL' cannot be set to '{value}'"
                    )),
                }
            }
            "USE_JDK_DEFLATER" | "use_jdk_deflater" => {
                jdk = value.eq_ignore_ascii_case("true");
            }
            _ => {}
        }
    }
    // `USE_JDK_DEFLATER=true` leaves htsjdk's own factory, libdeflate; the default installs GKL.
    (
        level,
        if jdk {
            Deflater::Libdeflate
        } else {
            Deflater::Gkl
        },
    )
}

/// The level and deflater [`install`] set, or htsjdk-rs's defaults when no binary installed any.
///
/// `copied` turns the borrowed `&(u32, Deflater)` into a copy of the pair.
pub fn current() -> (u32, Deflater) {
    SETTING
        .get()
        .copied()
        .unwrap_or((DEFAULT_COMPRESSION_LEVEL, Deflater::Jdk))
}

/// `SAMFileWriterFactory.makeBAMWriter` under the installed setting.
pub fn bam_writer<W: Write>(inner: W, header: &SamHeader) -> io::Result<BamWriter<W>> {
    let (level, deflater) = current();
    BamWriter::with_compression(inner, header, level, deflater)
}

/// `new BlockCompressedOutputStream(...)` under the installed setting.
pub fn bgzf_writer<W: Write>(inner: W) -> BgzfWriter<W> {
    let (level, deflater) = current();
    BgzfWriter::with_deflater(inner, level, deflater)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(line: &str) -> Vec<String> {
        line.split_whitespace().map(str::to_string).collect()
    }

    #[test]
    fn picard_defaults_to_gkl_at_level_five() {
        assert_eq!(setting(&words("I=a.bam O=b.bam")), (5, Deflater::Gkl));
    }

    #[test]
    fn every_spelling_of_the_two_arguments_is_read() {
        assert_eq!(
            setting(&words("COMPRESSION_LEVEL=1 USE_JDK_DEFLATER=true")),
            (1, Deflater::Libdeflate)
        );
        assert_eq!(
            setting(&words("--COMPRESSION_LEVEL 9 --use_jdk_deflater true")),
            (9, Deflater::Libdeflate)
        );
        assert_eq!(
            setting(&words("--USE_JDK_DEFLATER --COMPRESSION_LEVEL=0")),
            (0, Deflater::Libdeflate)
        );
        assert_eq!(
            setting(&words("USE_JDK_DEFLATER=true USE_JDK_DEFLATER=false")),
            (5, Deflater::Gkl)
        );
    }

    #[test]
    fn a_library_caller_keeps_htsjdks_defaults() {
        assert_eq!(current(), (DEFAULT_COMPRESSION_LEVEL, Deflater::Jdk));
    }
}
