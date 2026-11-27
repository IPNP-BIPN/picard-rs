//! CRAM outputs, for every tool whose writer comes from `SAMFileWriterFactory.makeWriter`.
//!
//! `makeWriter(header, presorted, output, REFERENCE_SEQUENCE)` returns a `CRAMFileWriter` when the
//! output's name ends in `.cram`, and a SAM or BAM writer otherwise. The tools here build their
//! BAM in memory and write it at the end, so a CRAM output is that BAM's header and records
//! written again as CRAM: the header the BAM writer stamped is the one the CRAM writer stamps,
//! and the records are the same records.
//!
//! The CRAM is htsjdk 5.0.0's, the version Picard 3.5.0 bundles: CRAM 3.1 under the NORMAL
//! profile, with the output's `file:` URI as the file definition's id.
//!
//! # A missing reference leaves an empty file
//!
//! `Files.newOutputStream` creates the output when the writer is made, and the header the writer
//! writes first sits in a buffer that is never flushed when the first mapped record finds no
//! reference. Picard exits with the file there and empty, and so does this.

use htsjdk_bam::reader::BamReader;

/// `CRAM_TYPE.hasValidFileExtension`: the name ends in `.cram`.
pub fn is_cram(path: &str) -> bool {
    path.ends_with(".cram")
}

/// The id `SAMFileWriterFactory` gives a CRAM it writes to a file: `outputFile.toUri()`.
pub fn file_id(path: &str) -> String {
    crate::fingerprint::file_uri(path)
}

/// A CRAM of `header` and `records`, against this process's reference, or the line an uncaught
/// exception prints when it fails.
pub fn encode(
    path: &str,
    header: &htsjdk_bam::header::SamHeader,
    records: &[htsjdk_bam::record::BamRecord],
) -> Result<Vec<u8>, String> {
    let mut source = writer_reference();
    htsjdk_cram::writer::write_cram(header, records, &file_id(path), source.as_mut())
        .map_err(|error| error.java_line())
}

/// The reference `makeWriter(..., REFERENCE_SEQUENCE)` builds a CRAM writer with: the FASTA, or,
/// with none, `CRAMLazyReferenceSource`, which refuses the first mapped record with its own
/// message rather than the one a reader without a reference gives.
fn writer_reference() -> Box<dyn htsjdk_cram::reader::ReferenceSource> {
    match crate::cram_input::argv_reference() {
        Some(path) => Box::new(htsjdk_cram::reader::FastaReference::open(
            std::path::Path::new(&path),
        )),
        None => Box::new(htsjdk_cram::reader::LazyReference),
    }
}

/// The finished BAM `bam` written to `path`: as it is, or, for a `.cram` path, as the CRAM of the
/// same header and records, indexed when `CREATE_INDEX` is set, which is the default every
/// `new SAMFileWriterFactory()` takes from the command line. A CRAM that cannot be written leaves
/// the file empty and exits as Picard does.
pub fn write(path: &str, bam: Vec<u8>) -> std::io::Result<()> {
    if !is_cram(path) {
        return std::fs::write(path, bam);
    }
    write_from_bam(path, &bam, create_index_argument())
}

/// `CREATE_INDEX` as this process was given it, in either spelling, the last one winning.
pub fn create_index_argument() -> bool {
    boolean_argument("CREATE_INDEX=")
}

fn boolean_argument(prefix: &str) -> bool {
    let mut value = false;
    for raw in crate::args::process().into_iter().skip(1) {
        if let Some(v) = raw.trim_start_matches('-').strip_prefix(prefix) {
            value = v.eq_ignore_ascii_case("true");
        }
    }
    value
}

/// The CRAM, and with `CREATE_MD5_FILE` its digest beside it as `Md5CalculatingOutputStream`
/// leaves it: thirty-two lower-case hex digits and no newline.
fn write_with_md5(path: &str, bytes: &[u8]) -> std::io::Result<()> {
    use md5::{Digest, Md5};
    std::fs::write(path, bytes)?;
    if boolean_argument("CREATE_MD5_FILE=") {
        let digest: String = Md5::digest(bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        std::fs::write(format!("{path}.md5"), digest)?;
    }
    Ok(())
}

/// `header` and `records` written to `path` as a CRAM.
pub fn write_records(
    path: &str,
    header: &htsjdk_bam::header::SamHeader,
    records: &[htsjdk_bam::record::BamRecord],
) -> std::io::Result<()> {
    match encode(path, header, records) {
        Ok(bytes) => write_with_md5(path, &bytes),
        Err(line) => {
            std::fs::write(path, [])?;
            crate::metrics_cli::thrown(&line)
        }
    }
}

/// [`write`] for a tool that also honours `CREATE_INDEX`, given the BAM it built: a CRAM output is
/// indexed by `CRAMBAIIndexer` into `<output>.bai`, which refuses a header that is not
/// coordinate-sorted where a BAM writer would quietly skip the index.
pub fn write_from_bam(path: &str, bam: &[u8], create_index: bool) -> std::io::Result<()> {
    let plain =
        crate::bgzf::decompress_all(bam).map_err(|e| std::io::Error::other(format!("{e:?}")))?;
    let reader = BamReader::new(&plain).map_err(|e| std::io::Error::other(format!("{e:?}")))?;
    let header = reader.header.text.clone();
    let records = reader
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| std::io::Error::other(format!("{e:?}")))?;
    write_records_indexed(path, &header, &records, create_index)
}

/// A CRAM written to its file as each container completes, as `CRAMFileWriter` writes it: memory
/// holds the container being filled, not the file.
pub struct CramFile {
    writer: htsjdk_cram::writer::CramWriter<'static>,
    file: std::io::BufWriter<std::fs::File>,
    path: String,
    coordinate: bool,
    md5: Option<md5::Md5>,
}

impl CramFile {
    /// Creates `path` for `header`, as the BAM writer would have stamped it.
    pub fn create(path: &str, header: &htsjdk_bam::header::SamHeader) -> std::io::Result<CramFile> {
        use md5::Digest;
        // The header the whole-file path read back from the BAM it wrote.
        let header = htsjdk_bam::reader::parse_header_text(&header.encode());
        // shortcut: the reference source lives until exit, as `CramWriter` borrows it; one per
        // CRAM output.
        let source = Box::leak(writer_reference());
        let coordinate = crate::metrics_cli::sort_order(&header) == "coordinate";
        let mut writer = htsjdk_cram::writer::CramWriter::new(&header, &file_id(path), source);
        // Indexed as it goes whenever it can be, the index kept only when `CREATE_INDEX` asks.
        if coordinate {
            writer = writer.with_index();
        }
        Ok(CramFile {
            writer,
            file: std::io::BufWriter::new(std::fs::File::create(path)?),
            path: path.to_string(),
            coordinate,
            md5: boolean_argument("CREATE_MD5_FILE=").then(md5::Md5::new),
        })
    }

    /// Writes `record`, and the container it completes; a record the CRAM cannot carry fails the
    /// run as the whole-file write did.
    pub fn write(&mut self, record: htsjdk_bam::record::BamRecord) -> std::io::Result<()> {
        if let Err(error) = self.writer.write(record) {
            self.fail(&error.java_line());
        }
        let bytes = self.writer.take_output();
        self.put(&bytes)
    }

    fn put(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        use md5::Digest;
        if let Some(md5) = &mut self.md5 {
            md5.update(bytes);
        }
        std::io::Write::write_all(&mut self.file, bytes)
    }

    /// The last container, the EOF container, and with `create_index` the `.bai`, which
    /// `CRAMBAIIndexer` refuses for a header that is not coordinate-sorted.
    pub fn finish(self, create_index: bool) -> std::io::Result<()> {
        use md5::Digest;
        let CramFile {
            writer,
            mut file,
            path,
            coordinate,
            md5,
        } = self;
        if create_index && !coordinate {
            fail(
                &mut file,
                &path,
                "htsjdk.samtools.SAMException: CRAM file mut be coordinate-sorted for indexing.",
            );
        }
        let (rest, index) = writer
            .finish()
            .unwrap_or_else(|error| fail(&mut file, &path, &error.java_line()));
        std::io::Write::write_all(&mut file, &rest)?;
        std::io::Write::flush(&mut file)?;
        if let Some(mut md5) = md5 {
            md5.update(&rest);
            let digest: String = md5.finalize().iter().map(|b| format!("{b:02x}")).collect();
            std::fs::write(format!("{path}.md5"), digest)?;
        }
        if create_index {
            std::fs::write(format!("{path}.bai"), index.unwrap_or_default())?;
        }
        Ok(())
    }

    /// The file as an unclosed writer leaves it: the containers completed so far.
    pub fn leave_unclosed(mut self) -> std::io::Result<()> {
        std::io::Write::flush(&mut self.file)
    }

    fn fail(&mut self, line: &str) -> ! {
        fail(&mut self.file, &self.path, line)
    }
}

/// What a CRAM that cannot be written leaves: the file there and empty, and with `CREATE_INDEX`
/// the index's too, then the uncaught exception.
fn fail(file: &mut std::io::BufWriter<std::fs::File>, path: &str, line: &str) -> ! {
    let _ = std::io::Write::flush(file);
    let _ = file.get_ref().set_len(0);
    if create_index_argument() {
        let _ = std::fs::write(format!("{path}.bai"), []);
    }
    // Left empty, not removed as a failing output is.
    crate::stream_output::forget(path);
    crate::metrics_cli::thrown(line)
}

/// [`write_records`], with `CREATE_INDEX`.
pub fn write_records_indexed(
    path: &str,
    header: &htsjdk_bam::header::SamHeader,
    records: &[htsjdk_bam::record::BamRecord],
    create_index: bool,
) -> std::io::Result<()> {
    if !create_index {
        return write_records(path, header, records);
    }
    let index_path = format!("{path}.bai");
    // The index stream is opened, then the CRAM's, then the indexer refuses.
    if crate::metrics_cli::sort_order(header) != "coordinate" {
        std::fs::write(&index_path, [])?;
        std::fs::write(path, [])?;
        crate::metrics_cli::thrown(
            "htsjdk.samtools.SAMException: CRAM file mut be coordinate-sorted for indexing.",
        );
    }
    let mut source = writer_reference();
    let written =
        htsjdk_cram::writer::write_cram_indexed(header, records, &file_id(path), source.as_mut())
            .map_err(|error| error.java_line());
    match written {
        Ok((bytes, index)) => {
            write_with_md5(path, &bytes)?;
            std::fs::write(&index_path, index)
        }
        Err(line) => {
            std::fs::write(&index_path, [])?;
            std::fs::write(path, [])?;
            crate::metrics_cli::thrown(&line)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_cram_written_as_it_goes_is_the_one_encoded_whole() {
        // Unplaced reads need no reference; 25,000 of them make three containers.
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n");
        for i in 0..25_000 {
            sam.push_str(&format!("r{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGTN\tIIIII\n"));
        }
        let (header, records) = htsjdk_bam::sam_file::read_sam(&sam).unwrap();
        let dir = std::env::temp_dir().join(format!("cram-file-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("x.cram").display().to_string();
        let whole = encode(
            &path,
            &htsjdk_bam::reader::parse_header_text(&header.encode()),
            &records,
        )
        .unwrap();

        let mut cram = CramFile::create(&path, &header).unwrap();
        let mut sizes = Vec::new();
        for record in records {
            cram.write(record).unwrap();
            sizes.push(std::fs::metadata(&path).unwrap().len());
        }
        cram.finish(false).unwrap();
        let streamed = std::fs::read(&path).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(streamed, whole);
        // The file grew before the end: containers reached it as they completed.
        assert!(sizes.last().unwrap() > &0, "{:?}", sizes.last());
    }
}
