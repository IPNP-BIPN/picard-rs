//! `CheckTerminatorBlock` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.sam.CheckTerminatorBlock.doWork` at tag 3.4.0. The classification is
//! `picard_analysis::check_terminator_block`; this prints its name on standard error, as
//! `System.err.println(term.name())` does, and exits 100 when the file is defective. The other
//! arguments the array varies (the stringency, the index and deflater switches, the record cap)
//! are read by nothing in `doWork`, so they change nothing here either.

use picard_analysis::check_terminator_block::check_terminator_block;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let input = args
        .iter()
        .find_map(|a| a.strip_prefix("INPUT=").or_else(|| a.strip_prefix("I=")))
        .unwrap_or_else(|| {
            eprintln!("INPUT= is required");
            std::process::exit(1);
        });
    let data = std::fs::read(input).unwrap_or_else(|e| {
        eprintln!("Exception in thread \"main\" htsjdk.samtools.SAMException: {e}");
        std::process::exit(1);
    });
    let (name, code) = check_terminator_block(&data);
    eprintln!("{name}");
    std::process::exit(code);
}
