//! `FifoBuffer` as a runnable binary: the covering array's port side.
//!
//! Ports `picard.util.FifoBuffer.doWork` at tag 3.4.0: standard input is copied to standard output
//! through a circular buffer of `BUFFER_SIZE` bytes, `IO_SIZE` bytes at a time. The bytes that come
//! out are the bytes that went in whatever the shape of the buffer, so the arguments that matter
//! are the ones that could break that: a buffer smaller than one read, a read of a single byte.
//! `DEBUG_FREQUENCY` and `NAME` only drive a log thread on standard error, which a port that has
//! no second thread has nothing to say about.

use picard_analysis::fifo_buffer::pump;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .rev()
        .find_map(|a| a.strip_prefix(key).map(str::to_string))
}

fn number(args: &[String], key: &str, default: i64) -> i64 {
    match arg(args, key) {
        None => default,
        Some(value) => value.parse().unwrap_or_else(|_| {
            eprintln!(
                "Argument '{}' cannot be set to '{value}': it is not a number",
                &key[..key.len() - 1]
            );
            std::process::exit(1)
        }),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    picard_analysis::argument_check::check("FifoBuffer");
    let args: Vec<String> = std::env::args().skip(1).collect();
    let buffer = number(&args, "BUFFER_SIZE=", 512 * 1024 * 1024);
    let io = number(&args, "IO_SIZE=", 64 * 1024);
    let _debug_frequency = number(&args, "DEBUG_FREQUENCY=", 0);
    if buffer <= 0 || io <= 0 {
        // `new byte[n]` with a negative `n` is a NegativeArraySizeException on a Fifo thread, and
        // a zero is a loop that never ends; neither is an answer a comparison can hold.
        eprintln!(
            "Exception in thread \"main\" picard.PicardException: Exception on input thread."
        );
        std::process::exit(1);
    }

    // The ring is never larger than one read needs: the reference allocates all of BUFFER_SIZE
    // (512 MiB by default), and the bytes that come out do not depend on it.
    let size = (buffer as usize).min(io as usize);
    pump(
        &mut std::io::stdin().lock(),
        &mut std::io::BufWriter::new(std::io::stdout().lock()),
        size,
        io as usize,
    )?;
    Ok(())
}
