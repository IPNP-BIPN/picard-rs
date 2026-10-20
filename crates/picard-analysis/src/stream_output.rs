//! An output written record by record, with the files the whole-file writers left.
//!
//! The tools used to build the whole SAM text or BAM stream in memory and write it at the end, so a
//! record that failed to encode, or a reader that threw halfway, left either no file or the empty
//! one the tool had opened first. This writer streams to the file instead and, on any of those
//! failures, leaves what the whole-file writer left: the file removed, or emptied.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::sync::Mutex;

use htsjdk_bam::build_index::build_bam_index;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::writer::BamWriter;

use crate::java_writer;

/// The outputs being written, each with whether a refusal leaves it empty rather than absent, for
/// [`crate::metrics_cli::thrown`] to clean up after a reader that throws mid-stream.
static OPEN: Mutex<Vec<(String, bool)>> = Mutex::new(Vec::new());

/// What a throw while an [`Output`] is open leaves on disk: called by
/// [`crate::metrics_cli::thrown`] before it exits.
pub fn on_thrown() {
    let open = OPEN.lock().map(|mut open| std::mem::take(&mut *open));
    // Last opened first, so that a directory goes after the files in it.
    for (path, empty) in open.unwrap_or_default().into_iter().rev() {
        leave(&path, empty);
    }
}

fn leave(path: &str, empty: bool) {
    // Best effort: the run is failing already, with its own message.
    let _ = if empty {
        File::create(path).map(drop)
    } else {
        // An [`OutputDir`] goes with the outputs already finished in it.
        std::fs::remove_file(path).or_else(|_| std::fs::remove_dir_all(path))
    };
}

/// An output directory the tool creates, removed again with everything in it if the run fails
/// before [`OutputDir::keep`], as the whole-file writers had not created it yet.
pub struct OutputDir(Option<String>);

impl OutputDir {
    pub fn create(path: &str) -> io::Result<OutputDir> {
        if std::path::Path::new(path).exists() {
            return Ok(OutputDir(None));
        }
        std::fs::create_dir_all(path)?;
        OPEN.lock()
            .expect("not poisoned")
            .push((path.to_string(), false));
        Ok(OutputDir(Some(path.to_string())))
    }

    pub fn keep(mut self) {
        if let Some(path) = self.0.take() {
            close(&path);
        }
    }
}

impl Drop for OutputDir {
    fn drop(&mut self) {
        if let Some(path) = self.0.take() {
            close(&path);
            let _ = std::fs::remove_dir_all(path);
        }
    }
}

/// Puts a file the tool writes itself on the open list: a throw before [`untrack`] removes it, as
/// a writer that held the file in memory had not created it yet.
pub fn track(path: &str) {
    OPEN.lock()
        .expect("not poisoned")
        .push((path.to_string(), false));
}

/// Takes a file [`track`] listed off the open list, once the tool is done with it.
pub fn untrack(path: &str) {
    close(path);
}

/// Takes `path` off the open list: whether it was on it, and how a failure leaves it.
fn close(path: &str) -> Option<bool> {
    let mut open = OPEN.lock().ok()?;
    let at = open.iter().position(|(p, _)| p == path)?;
    Some(open.remove(at).1)
}

enum Sink {
    /// The text writer, and the header length that [`java_writer::sam_on_disk`] keeps.
    Sam(BufWriter<File>, SamHeader, usize),
    Bam(BamWriter<BufWriter<File>>),
    /// shortcut: a CRAM is still encoded from the whole BAM at the end, as before; stream it once
    /// a CRAM output on a large input matters.
    Cram(BamWriter<Vec<u8>>),
}

/// `SAMFileWriterFactory.makeWriter`: SAM for a `.sam` name, CRAM for a CRAM one, BAM otherwise.
/// Dropped before [`Output::finish`], as by a `?` in the tool, it is abandoned.
pub struct Output {
    path: String,
    sink: Option<Sink>,
    checker: Option<Checker>,
}

/// `SAMSortOrderChecker`, which a writer opened presorted runs on every record when its header's
/// order is coordinate or queryname, as almost every Picard tool opens its writer: the previous
/// record's `fileOrderCompare` key, and the contig names its coordinate key is printed with.
struct Checker {
    queryname: bool,
    names: Vec<String>,
    previous: Option<Key>,
}

/// What `fileOrderCompare` reads: the reference index and start, or the read name.
#[derive(Clone, PartialEq)]
enum Key {
    Coordinate(i32, i32),
    Queryname(String),
}

impl Checker {
    fn for_header(header: &SamHeader) -> Option<Checker> {
        let queryname = match header.attributes.get("SO") {
            Some("coordinate") => false,
            Some("queryname") => true,
            _ => return None,
        };
        Some(Checker {
            queryname,
            names: header.sequences.iter().map(|s| s.name.clone()).collect(),
            previous: None,
        })
    }

    fn key(&self, record: &BamRecord) -> Key {
        if self.queryname {
            Key::Queryname(record.read_name.clone())
        } else {
            Key::Coordinate(record.reference_index, record.alignment_start)
        }
    }

    /// `SAMRecordCoordinateComparator` / `SAMRecordQueryNameComparator.fileOrderCompare`.
    fn after(a: &Key, b: &Key) -> bool {
        match (a, b) {
            (Key::Queryname(a), Key::Queryname(b)) => a > b,
            (Key::Coordinate(r1, s1), Key::Coordinate(r2, s2)) => match (*r1 == -1, *r2 == -1) {
                (true, true) => false,
                (true, false) => true,
                (false, true) => false,
                (false, false) => (r1, s1) > (r2, s2),
            },
            _ => false,
        }
    }

    /// `SAMSortOrderChecker.getSortKey`: `name:start` (`*:0` unmapped), or the read name.
    fn printed(&self, key: &Key) -> String {
        match key {
            Key::Queryname(name) => name.clone(),
            Key::Coordinate(index, start) => {
                let name = usize::try_from(*index)
                    .ok()
                    .and_then(|i| self.names.get(i))
                    .map_or("*", String::as_str);
                format!("{name}:{start}")
            }
        }
    }
}

impl Output {
    /// Opens `path` for `header`. `opened_first` says whether the tool's whole-file writer had
    /// already created the file when a record failed, so a failure leaves it empty, not absent.
    pub fn create(path: &str, header: &SamHeader, opened_first: bool) -> io::Result<Output> {
        Output::open(path, header, opened_first, path.ends_with(".sam"))
    }

    /// [`Self::create`] for a tool that writes a BAM whatever the name, a CRAM name aside.
    pub fn bam(path: &str, header: &SamHeader, opened_first: bool) -> io::Result<Output> {
        Output::open(path, header, opened_first, false)
    }

    fn open(path: &str, header: &SamHeader, opened_first: bool, sam: bool) -> io::Result<Output> {
        let sink = if crate::cram_output::is_cram(path) {
            if opened_first {
                File::create(path)?;
            }
            Sink::Cram(crate::output::bam_writer(Vec::new(), header)?)
        } else if sam {
            let mut out = BufWriter::new(File::create(path)?);
            let text = header.encode();
            out.write_all(text.as_bytes())?;
            Sink::Sam(out, header.clone(), text.len())
        } else {
            let out = BufWriter::new(File::create(path)?);
            Sink::Bam(crate::output::bam_writer(out, header)?)
        };
        if !matches!(sink, Sink::Cram(_)) {
            OPEN.lock()
                .expect("not poisoned")
                .push((path.to_string(), opened_first));
        }
        Ok(Output {
            path: path.to_string(),
            sink: Some(sink),
            checker: Checker::for_header(header),
        })
    }

    /// Writes `record`; a record the format cannot carry ends the output as a failure leaves it,
    /// with the message the whole-file writer returned.
    pub fn write(&mut self, record: &BamRecord) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(mut checker) = self.checker.take() {
            let key = checker.key(record);
            if let Some(previous) = &checker.previous {
                if Checker::after(previous, &key) {
                    // `isSorted` has already made this record the previous one when the message
                    // asks for it, so both keys printed are this record's.
                    let message =
                        format!(
                        "Alignments added out of order in SAMFileWriterImpl.addAlignment for {}. \
                         Sort order is {}. Offending records are at [{}] and [{}]",
                        self.writer_name(),
                        if checker.queryname { "queryname" } else { "coordinate" },
                        checker.printed(&key),
                        checker.printed(&key)
                    );
                    // `addAlignment` throws before the record is written, out of a writer that is
                    // never closed.
                    self.leave_unclosed()?;
                    eprintln!(
                        "Exception in thread \"main\" java.lang.IllegalArgumentException: {message}"
                    );
                    std::process::exit(1);
                }
            }
            checker.previous = Some(key);
            self.checker = Some(checker);
        }
        let failed = match self.sink.as_mut().expect("open until finished") {
            Sink::Sam(out, header, _) => match crate::bam_stream::sam_line(header, record) {
                Some(line) => {
                    out.write_all(line.as_bytes())?;
                    out.write_all(b"\n")?;
                    return Ok(());
                }
                None => "records failed to re-encode as SAM".to_string(),
            },
            Sink::Bam(writer) => match writer.write(record) {
                Ok(()) => return Ok(()),
                Err(e) => format!("{e:?}"),
            },
            Sink::Cram(writer) => return writer.write(record).map_err(|e| format!("{e:?}").into()),
        };
        self.abandon();
        Err(failed.into())
    }

    /// Ends a failed run: the file removed, or left empty when the tool had opened it first.
    fn abandon(&mut self) {
        // The writer goes first, or what it still buffers would land in the emptied file.
        if self.sink.take().is_some() {
            if let Some(empty) = close(&self.path) {
                leave(&self.path, empty);
            }
        }
    }

    /// `getFilename()`: a BAM writer's absolute path as a URI, a SAM writer's as a path, neither
    /// with its links resolved.
    /// shortcut: the URI is not percent-encoded; encode it once a path with a space matters.
    fn writer_name(&self) -> String {
        let absolute = std::path::absolute(&self.path)
            .map(|p| p.display().to_string())
            .unwrap_or_else(|_| self.path.clone());
        match self.sink {
            Some(Sink::Sam(..)) => absolute,
            _ => format!("file://{absolute}"),
        }
    }

    /// Ends a run that threw after the records written so far: the file as htsjdk's unclosed
    /// writer leaves it ([`java_writer`]), whatever the format.
    pub fn abandon_unclosed(mut self) -> io::Result<()> {
        self.leave_unclosed()
    }

    fn leave_unclosed(&mut self) -> io::Result<()> {
        close(&self.path);
        match self.sink.take().expect("open until finished") {
            Sink::Sam(out, _, header_len) => {
                let file = out.into_inner().map_err(|e| e.into_error())?;
                let records = file.metadata()?.len() as usize - header_len;
                let buffers = records / java_writer::BUFFER_SIZE * java_writer::BUFFER_SIZE;
                file.set_len((header_len + buffers) as u64)
            }
            Sink::Bam(writer) => {
                writer.finish()?.flush()?;
                let finished = std::fs::read(&self.path)?;
                std::fs::write(&self.path, java_writer::bam_on_disk(&finished))
            }
            // The whole-file writer wrote the partial BAM under the CRAM's name.
            Sink::Cram(writer) => {
                std::fs::write(&self.path, java_writer::bam_on_disk(&writer.finish()?))
            }
        }
    }

    /// Closes the output. A CRAM is indexed when `create_index` (`CREATE_INDEX`) is set, whatever
    /// its order, and a BAM when `index_bam` is: the tool decides whether its order asks for one.
    pub fn finish(
        mut self,
        create_index: bool,
        index_bam: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        close(&self.path);
        match self.sink.take().expect("open until finished") {
            Sink::Sam(mut out, ..) => out.flush()?,
            Sink::Bam(writer) => {
                writer.finish()?.flush()?;
                if index_bam {
                    let bam = std::fs::read(&self.path)?;
                    let index = build_bam_index(&bam).map_err(|e| format!("{e:?}"))?;
                    std::fs::write(index_path(&self.path), index)?;
                }
            }
            Sink::Cram(writer) => {
                crate::cram_output::write_from_bam(&self.path, &writer.finish()?, create_index)?
            }
        }
        Ok(())
    }
}

impl Drop for Output {
    fn drop(&mut self) {
        self.abandon();
    }
}

/// `BAMFileWriter.createBamIndex`: the `.bam` suffix is dropped, anything else is kept.
pub fn index_path(output: &str) -> String {
    match output.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.bai"),
        None => format!("{output}.bai"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::sam_file::{read_sam, write_sam};

    const SAM: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n\
                       r1\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
                       r2\t0\tchr1\t9\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

    fn scratch(name: &str) -> String {
        let dir = std::env::temp_dir().join(format!("stream_output_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name).display().to_string()
    }

    #[test]
    fn a_finished_sam_is_the_whole_file_writers() {
        let (header, records) = read_sam(SAM).unwrap();
        let path = scratch("whole.sam");
        let mut out = Output::create(&path, &header, false).unwrap();
        for record in &records {
            out.write(record).unwrap();
        }
        out.finish(false, false).unwrap();
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            write_sam(&header, &records).unwrap()
        );
    }

    #[test]
    fn a_dropped_output_is_removed_or_emptied() {
        let (header, records) = read_sam(SAM).unwrap();
        for (opened_first, name) in [(false, "dropped.bam"), (true, "emptied.bam")] {
            let path = scratch(name);
            let mut out = Output::create(&path, &header, opened_first).unwrap();
            out.write(&records[0]).unwrap();
            drop(out);
            match opened_first {
                false => assert!(!std::path::Path::new(&path).exists()),
                true => assert_eq!(std::fs::metadata(&path).unwrap().len(), 0),
            }
        }
    }

    #[test]
    fn an_unclosed_sam_keeps_its_header() {
        let (header, records) = read_sam(SAM).unwrap();
        let path = scratch("unclosed.sam");
        let mut out = Output::create(&path, &header, false).unwrap();
        out.write(&records[0]).unwrap();
        out.abandon_unclosed().unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), header.encode());
    }
}
