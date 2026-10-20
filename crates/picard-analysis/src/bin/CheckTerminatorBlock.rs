//! `CheckTerminatorBlock` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CheckTerminatorBlock.doWork` at tag 3.4.0. The classification is
//! `picard_analysis::check_terminator_block`; this prints its name on standard error, as
//! `System.err.println(term.name())` does, and exits 100 when the file is defective. The other
//! arguments the array varies (the stringency, the index and deflater switches, the record cap)
//! are read by nothing in `doWork`, so they change nothing here either.

use picard_analysis::check_terminator_block::check_terminator_block;
use std::io::{Read, Seek, SeekFrom};

/// The last block's worth of the file: `checkTermination` reads no further back than that.
fn tail(path: &str) -> std::io::Result<Vec<u8>> {
    let mut file = std::fs::File::open(path)?;
    let len = file.metadata()?.len();
    file.seek(SeekFrom::Start(
        len.saturating_sub(htsjdk_bgzf::MAX_COMPRESSED_BLOCK_SIZE as u64),
    ))?;
    let mut data = Vec::new();
    file.read_to_end(&mut data)?;
    Ok(data)
}

fn main() {
    picard_analysis::argument_check::check("CheckTerminatorBlock");
    picard_analysis::output::install();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let input = args
        .iter()
        .find_map(|a| a.strip_prefix("INPUT=").or_else(|| a.strip_prefix("I=")))
        .unwrap_or_else(|| {
            eprintln!("INPUT= is required");
            std::process::exit(1);
        });
    let data = tail(input).unwrap_or_else(|e| {
        eprintln!("Exception in thread \"main\" htsjdk.samtools.SAMException: {e}");
        std::process::exit(1);
    });
    let (name, code) = check_terminator_block(&data);
    eprintln!("{name}");
    std::process::exit(code);
}
