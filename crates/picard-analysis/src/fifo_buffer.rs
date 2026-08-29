//! `FifoBuffer`: the smallest tool in Picard.
//!
//! It copies its input to its output through a CIRCULAR buffer so a pipeline can decouple a slow
//! writer from a slow reader. What there is to reproduce is that the bytes come out unchanged
//! whatever the buffer's shape, including a buffer smaller than the input: the reading and writing
//! halves take turns rather than needing room for the whole of it.
//!
//! Ported from `picard.util.FifoBuffer` and `htsjdk.samtools.util.CircularByteBuffer` in
//! Picard 3.4.0.

/// `CircularByteBuffer`: a fixed ring with a read cursor, a write cursor and a closed flag.
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
    pub fn write(&mut self, data: &[u8]) -> usize {
        let room = self.capacity() - self.filled;
        let taken = room.min(data.len());
        for byte in &data[..taken] {
            let at = self.write_at;
            self.bytes[at] = *byte;
            self.write_at = (at + 1) % self.capacity();
        }
        self.filled += taken;
        taken
    }

    /// `read(buffer, start, length)`: as much as is there, which may be none.
    pub fn read(&mut self, into: &mut [u8]) -> usize {
        let taken = self.filled.min(into.len());
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
    let mut buffer = CircularBuffer::new(buffer_size);
    let mut output = Vec::with_capacity(input.len());
    let chunk = io_size.max(1);
    let mut written = 0;
    let mut scratch = vec![0u8; chunk];
    while written < input.len() || buffer.filled > 0 {
        if written < input.len() {
            let end = (written + chunk).min(input.len());
            written += buffer.write(&input[written..end]);
        }
        let taken = buffer.read(&mut scratch);
        output.extend_from_slice(&scratch[..taken]);
        if taken == 0 && written >= input.len() {
            break;
        }
    }
    buffer.close();
    output
}
