//! `FifoBuffer`: the smallest tool in Picard.
//!
//! It copies its input to its output through a CIRCULAR buffer so a pipeline can decouple a slow
//! writer from a slow reader. What there is to reproduce is that the bytes come out unchanged
//! whatever the buffer's shape, including a buffer smaller than the input: the reading and writing
//! halves take turns rather than needing room for the whole of it.
//!
//! Ported from `picard.util.FifoBuffer` and `picard.util.CircularByteBuffer` in Picard 3.5.0;
//! both files are unchanged since 3.4.0.
//!
//! The reference runs two threads around the buffer, synchronised with `wait`/`notify`. Here one
//! loop takes their turns in order. Only the output bytes are observable, and they are what the
//! golden (`tests/data/fifo_buffer.txt.gz`, eight runs) records; how many bytes each call moved
//! is not recorded and differs (see [`CircularBuffer::write`]).
//!
//! # Reading this file without reading Rust
//!
//! - `&[u8]` is a borrowed, read-only view of a byte array and `&mut [u8]` a writable one; the
//!   caller keeps ownership. `&data[..n]` is the view of its first `n` bytes.
//! - `usize` is an unsigned index type; a size here cannot be negative.
//! - `input: &mut impl std::io::Read` accepts any readable stream (a file, standard input, a
//!   slice), as a Java parameter typed `InputStream` would.
//! - `std::io::Result<()>` is "done, or an I/O error"; `?` hands an error to the caller at once,
//!   Java's `throws IOException` at a single call site.
//! - `0..n` is a half-open range of indices, `start` included and `end` excluded.

/// `CircularByteBuffer`: a fixed ring with a read cursor, a write cursor and a closed flag.
///
/// The reference's `nextWritePos` and `nextReadPos` are `write_at` and `read_at`; its
/// `bytesAvailableToRead` is `filled`, and `bytesAvailableToWrite` is `capacity - filled`.
#[derive(Debug)]
pub struct CircularBuffer {
    bytes: Vec<u8>,
    read_at: usize,
    write_at: usize,
    filled: usize,
    closed: bool,
}

impl CircularBuffer {
    /// A buffer of the given size, which the reference caps at nothing: one byte is legal.
    ///
    /// A size of zero is raised to one here. The reference's `new byte[size]` would make a
    /// zero-byte ring that `write` can never put a byte into.
    pub fn new(size: usize) -> Self {
        Self {
            bytes: vec![0; size.max(1)],
            read_at: 0,
            write_at: 0,
            filled: 0,
            closed: false,
        }
    }

    pub fn capacity(&self) -> usize {
        self.bytes.len()
    }

    pub fn is_closed(&self) -> bool {
        self.closed
    }

    /// `close()`, which the reading half does when its input runs out.
    pub fn close(&mut self) {
        self.closed = true;
    }

    /// `write(buffer, start, length)`: as much as fits, which may be none.
    ///
    /// The reference copies one contiguous run and waits while the ring is full. Its length is
    /// the smallest of `size`, `bytesAvailableToWrite` and the distance to the end of the array
    /// (`CircularByteBuffer.java:42-44`), which begins:
    ///
    /// ```java
    /// final int distanceToEnd = this.capacity - writePos;
    /// ```
    ///
    /// This one wraps past the end within the same call and returns zero instead of waiting. The
    /// same bytes enter the ring in the same order; only the number of calls differs.
    pub fn write(&mut self, data: &[u8]) -> usize {
        let room = self.capacity() - self.filled;
        let taken = room.min(data.len());
        // `%` wraps the cursor back to the start of the ring, as `% this.capacity` does.
        for byte in &data[..taken] {
            let at = self.write_at;
            self.bytes[at] = *byte;
            self.write_at = (at + 1) % self.capacity();
        }
        self.filled += taken;
        taken
    }

    /// `read(buffer, start, length)`: as much as is there, which may be none.
    ///
    /// Like [`CircularBuffer::write`], it wraps within one call where the reference stops at the
    /// end of its array (`CircularByteBuffer.java:66-68`), and returns zero where the reference
    /// waits for the writer.
    pub fn read(&mut self, into: &mut [u8]) -> usize {
        let taken = self.filled.min(into.len());
        // `iter_mut().take(taken)` visits the first `taken` slots of `into` writably; `*slot = ..`
        // stores through the reference.
        for slot in into.iter_mut().take(taken) {
            let at = self.read_at;
            *slot = self.bytes[at];
            self.read_at = (at + 1) % self.capacity();
        }
        self.filled -= taken;
        taken
    }
}

/// `doWork`: the whole tool, with the two threads' turns taken in order.
///
/// The reference runs the halves as two threads and lets the buffer block; here they alternate,
/// which produces the same bytes because the buffer is the only thing between them. A buffer of
/// ONE byte still copies everything, one byte at a time.
pub fn copy(input: &[u8], buffer_size: usize, io_size: usize) -> Vec<u8> {
    let mut output = Vec::with_capacity(input.len());
    let mut reader = input;
    // `expect` turns an error into a panic with this message; reading a slice and writing a
    // vector cannot fail, so there is nothing to hand back.
    pump(&mut reader, &mut output, buffer_size, io_size).expect("a slice and a vector do not fail");
    output
}

/// The same turns over streams: read up to `IO_SIZE` bytes, put what fits into the ring, and write
/// out what the ring holds, until the input ends and the ring is empty. Nothing is held beyond the
/// ring, so the tool keeps up with a pipeline rather than waiting for the end of its input.
///
/// The reference's two loops, Picard 3.5.0 `src/main/java/picard/util/FifoBuffer.java:116-120`
/// (input thread) and `:138-139` (output thread):
///
/// ```java
/// while ((read = inputStream.read(buffer)) > -1) {
///     int start = 0;
///     while (start < read) {
///         start += fifo.write(buffer, start, read-start);
///     }
/// ```
///
/// ```java
/// while ((read = fifo.read(buffer, 0, buffer.length)) > 0 || !fifo.isClosed()) {
///     outputStream.write(buffer, 0, read);
/// ```
///
/// `incoming` is the input thread's `buffer` and `pending` the part of it not yet in the ring
/// (`start..read`); `scratch` is the output thread's `buffer`. An `io_size` of zero is raised to
/// one.
pub fn pump(
    input: &mut impl std::io::Read,
    output: &mut impl std::io::Write,
    buffer_size: usize,
    io_size: usize,
) -> std::io::Result<()> {
    let mut buffer = CircularBuffer::new(buffer_size);
    let chunk = io_size.max(1);
    let mut incoming = vec![0u8; chunk];
    let mut pending = 0..0;
    let mut scratch = vec![0u8; chunk];
    let mut ended = false;
    loop {
        // The input half's turn: once the previous read is all in the ring, read again.
        if pending.is_empty() && !ended {
            // Retry a read the OS interrupted; `break other?` leaves the inner loop with the
            // byte count, or returns the error from `pump`. A Rust `read` returns 0 only at the
            // end of the input, where Java's returns -1.
            let n = loop {
                match input.read(&mut incoming) {
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                    other => break other?,
                }
            };
            ended = n == 0;
            pending = 0..n;
        }
        // `start += fifo.write(...)`. A range is not copied implicitly in Rust, so indexing with
        // it takes a `clone()` to keep `pending` usable afterwards.
        if !pending.is_empty() {
            pending.start += buffer.write(&incoming[pending.clone()]);
        }
        // The output half's turn.
        let taken = buffer.read(&mut scratch);
        output.write_all(&scratch[..taken])?;
        // The output thread's exit: the ring is drained and nothing more will be written to it.
        if ended && pending.is_empty() && buffer.filled == 0 {
            break;
        }
    }
    buffer.close();
    output.flush()
}
