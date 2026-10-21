//! How every block-compressed file a tool writes is deflated: `COMPRESSION_LEVEL` and
//! `USE_JDK_DEFLATER`.
//!
//! Picard sets both once, in `CommandLineProgram.instanceMain`, before `doWork` runs:
//!
//! ```java
//! BlockCompressedOutputStream.setDefaultCompressionLevel(COMPRESSION_LEVEL);
//! if (!USE_JDK_DEFLATER) {
//!     BlockCompressedOutputStream.setDefaultDeflaterFactory(new IntelDeflaterFactory());
//! }
//! ```
//!
//! They are htsjdk statics, so every BAM, BGZF VCF, FASTQ and tabix index the tool opens afterwards
//! takes them, whichever class opens it. [`install`] is that call, made by each binary;
//! [`bam_writer`], [`bgzf_writer`] and [`current`] are what the writers read.
//!
//! `USE_JDK_DEFLATER=true` does not mean `java.util.zip`. It only stops Picard installing GKL, and
//! the factory left in place is htsjdk 5.0.0's own, which returns libdeflate whenever its native
//! library loads (`samjdk.use_libdeflate`, true by default, and jlibdeflate ships one for every
//! platform Picard runs on).
//!
//! A caller that never runs a binary, a test among them, keeps htsjdk-rs's defaults, the JDK
//! deflater at level five, which is what every committed golden was measured with.

use std::io::{self, Write};
use std::sync::OnceLock;

use htsjdk_bam::header::SamHeader;
use htsjdk_bam::writer::BamWriter;
use htsjdk_bgzf::{BgzfWriter, Deflater, DEFAULT_COMPRESSION_LEVEL};

/// `CommandLineProgram.COMPRESSION_LEVEL`'s default, which is htsjdk's
/// `Defaults.COMPRESSION_LEVEL`.
const PICARD_COMPRESSION_LEVEL: u32 = 5;

static SETTING: OnceLock<(u32, Deflater)> = OnceLock::new();

/// Reads `COMPRESSION_LEVEL` and `USE_JDK_DEFLATER` (short name `use_jdk_deflater`) from the
/// process's command line and makes them the level and deflater of every writer opened after it.
///
/// Without `USE_JDK_DEFLATER=true` the deflater is Intel's GKL, which is Picard's default and
/// therefore what an unchanged Picard command line writes; with it, libdeflate.
pub fn install() {
    let args: Vec<String> = crate::args::process().into_iter().skip(1).collect();
    let _ = SETTING.set(setting(&args));
    crate::bam_stream::budget_threads();
}

/// The setting a command line asks for, the last value of each argument winning as Barclay's does.
pub fn setting(args: &[String]) -> (u32, Deflater) {
    let mut level = PICARD_COMPRESSION_LEVEL;
    let mut jdk = false;
    let mut words = args.iter();
    while let Some(word) = words.next() {
        let bare = word.trim_start_matches('-');
        let (name, value) = match bare.split_once('=') {
            Some((name, value)) => (name, value.to_string()),
            // `--NAME value`, Barclay's own form; a flag given alone is `true`.
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
