//! What `LiftoverVcf` shares with other tools: htsjdk's `VCFEncoder.formatVCFDouble`, and the
//! `SortingCollection` the tool writes its output through. The tool itself is
//! `bin/LiftoverVcf.rs`.

/// `VCFEncoder.formatVCFDouble`, which is how a computed number reaches the file.
pub fn format_vcf_double(value: f64) -> String {
    if value < 1.0 {
        if value < 0.01 {
            if value.abs() >= 1e-20 {
                return format!("{value:.3e}");
            }
            return "0.00".to_string();
        }
        return format!("{value:.3}");
    }
    format!("{value:.2}")
}

/// `SortingCollection` over encoded lines: each held as its target contig index, its start and
/// its bytes, a full buffer of `max_records_in_ram` sorted (stably) and spilled to an unlinked
/// file in `tmp_dir`, the runs merged back with ties going to the earlier run. Unsorted, every
/// key is equal and the lines come back in the order they were added.
pub struct LineSpool {
    sort: bool,
    max_records_in_ram: usize,
    tmp_dir: std::path::PathBuf,
    ram: Vec<(u32, i64, Box<[u8]>)>,
    runs: Vec<std::fs::File>,
}

impl LineSpool {
    pub fn new(sort: bool, max_records_in_ram: usize, tmp_dir: std::path::PathBuf) -> Self {
        LineSpool {
            sort,
            max_records_in_ram: max_records_in_ram.max(1),
            tmp_dir,
            ram: Vec::new(),
            runs: Vec::new(),
        }
    }

    pub fn add(&mut self, contig: u32, start: i64, line: &[u8]) -> std::io::Result<()> {
        if self.ram.len() == self.max_records_in_ram {
            self.spill()?;
        }
        let (contig, start) = if self.sort { (contig, start) } else { (0, 0) };
        self.ram.push((contig, start, line.into()));
        Ok(())
    }

    fn spill(&mut self) -> std::io::Result<()> {
        use std::io::{Seek, Write};
        self.ram.sort_by_key(|(contig, start, _)| (*contig, *start));
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let path = self
            .tmp_dir
            .join(format!("liftovervcf.{}.{n}.tmp", std::process::id()));
        let file = std::fs::File::options()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)?;
        // shortcut: unlinking an open file keeps it readable on Unix only; Windows is not a target.
        std::fs::remove_file(&path)?;
        let mut out = std::io::BufWriter::new(file);
        for (contig, start, line) in self.ram.drain(..) {
            out.write_all(&contig.to_le_bytes())?;
            out.write_all(&start.to_le_bytes())?;
            out.write_all(&(line.len() as u32).to_le_bytes())?;
            out.write_all(&line)?;
        }
        let mut file = out.into_inner().map_err(|e| e.into_error())?;
        file.seek(std::io::SeekFrom::Start(0))?;
        self.runs.push(file);
        Ok(())
    }

    /// `doneAdding` and `iterator`: the lines in order, from RAM when nothing spilled.
    pub fn into_lines(mut self) -> std::io::Result<SpoolLines> {
        if self.runs.is_empty() {
            self.ram.sort_by_key(|(contig, start, _)| (*contig, *start));
            return Ok(SpoolLines::Ram(self.ram.into_iter()));
        }
        if !self.ram.is_empty() {
            self.spill()?;
        }
        // shortcut: every run open at once, as htsjdk has them; a file-descriptor limit caps the
        // input at about that many times MAX_RECORDS_IN_RAM.
        let mut runs: Vec<_> = self.runs.into_iter().map(std::io::BufReader::new).collect();
        let mut heap = std::collections::BinaryHeap::new();
        for (index, run) in runs.iter_mut().enumerate() {
            if let Some((contig, start, line)) = read_spilled(run)? {
                heap.push(std::cmp::Reverse((contig, start, index, line)));
            }
        }
        Ok(SpoolLines::Merge { runs, heap })
    }
}

type Spilled = (u32, i64, Box<[u8]>);

/// A run's next line, keyed so the heap's top is the least, ties going to the earlier run.
type Head = std::cmp::Reverse<(u32, i64, usize, Box<[u8]>)>;

fn read_spilled(run: &mut impl std::io::BufRead) -> std::io::Result<Option<Spilled>> {
    if run.fill_buf()?.is_empty() {
        return Ok(None);
    }
    let mut key = [0u8; 16];
    run.read_exact(&mut key)?;
    let len = u32::from_le_bytes(key[12..].try_into().expect("four bytes")) as usize;
    let mut line = vec![0u8; len];
    run.read_exact(&mut line)?;
    Ok(Some((
        u32::from_le_bytes(key[..4].try_into().expect("four bytes")),
        i64::from_le_bytes(key[4..12].try_into().expect("eight bytes")),
        line.into(),
    )))
}

/// A [`LineSpool`]'s lines in order.
pub enum SpoolLines {
    Ram(std::vec::IntoIter<Spilled>),
    Merge {
        runs: Vec<std::io::BufReader<std::fs::File>>,
        heap: std::collections::BinaryHeap<Head>,
    },
}

impl Iterator for SpoolLines {
    type Item = std::io::Result<Box<[u8]>>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            SpoolLines::Ram(lines) => lines.next().map(|(_, _, line)| Ok(line)),
            SpoolLines::Merge { runs, heap } => {
                let std::cmp::Reverse((_, _, index, line)) = heap.pop()?;
                match read_spilled(&mut runs[index]) {
                    Ok(Some((contig, start, next))) => {
                        heap.push(std::cmp::Reverse((contig, start, index, next)))
                    }
                    Ok(None) => {}
                    Err(e) => return Some(Err(e)),
                }
                Some(Ok(line))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::LineSpool;

    fn drain(spool: LineSpool) -> Vec<String> {
        spool
            .into_lines()
            .unwrap()
            .map(|l| String::from_utf8(l.unwrap().into()).unwrap())
            .collect()
    }

    #[test]
    fn a_spilled_sort_is_the_stable_sort() {
        let keys = [(1, 5), (0, 9), (1, 5), (0, 2), (2, 1), (0, 9), (1, 1)];
        for max in [1, 2, 3, 100] {
            let mut sorted = LineSpool::new(true, max, std::env::temp_dir());
            let mut unsorted = LineSpool::new(false, max, std::env::temp_dir());
            for (i, (contig, start)) in keys.iter().enumerate() {
                let line = format!("{contig}:{start}#{i}");
                sorted.add(*contig, *start, line.as_bytes()).unwrap();
                unsorted.add(*contig, *start, line.as_bytes()).unwrap();
            }
            let mut expected: Vec<_> = keys.iter().enumerate().collect();
            expected.sort_by_key(|(_, k)| **k);
            let expected: Vec<String> = expected
                .iter()
                .map(|(i, (c, s))| format!("{c}:{s}#{i}"))
                .collect();
            assert_eq!(drain(sorted), expected, "max {max}");
            let in_order: Vec<String> = keys
                .iter()
                .enumerate()
                .map(|(i, (c, s))| format!("{c}:{s}#{i}"))
                .collect();
            assert_eq!(drain(unsorted), in_order, "max {max}");
        }
    }
}
