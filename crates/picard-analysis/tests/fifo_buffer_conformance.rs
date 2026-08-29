//! Conformance for `FifoBuffer` against Picard 3.4.0.
//!
//! Golden from `tools/fifobuffer-conformance/FifoBufferDump.java`, eight runs whose outputs are
//! recorded as base64 because the point of the tool is that it does not care what it is carrying.
//!
//! # What this suite is for
//!
//!  * **the bytes coming out unchanged, whatever the buffer's shape**;
//!  * **a buffer smaller than the input still copying all of it**;
//!  * **bytes that are not text surviving**;
//!  * **and an empty input being an empty output.**

use std::io::Read;

use picard_analysis::fifo_buffer::{copy, CircularBuffer};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fifo_buffer.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn field(text: &str, kind: &str, case: &str) -> Option<String> {
    let prefix = format!("{kind}\t{case}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
}

/// The golden's base64, decoded.
fn output(text: &str, case: &str) -> Vec<u8> {
    let encoded = field(text, "out", case).unwrap_or_else(|| panic!("{case}"));
    let table: Vec<u8> =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/".to_vec();
    let mut bits = 0u32;
    let mut count = 0;
    let mut bytes = Vec::new();
    for character in encoded.bytes() {
        if character == b'=' {
            break;
        }
        let value = table
            .iter()
            .position(|entry| *entry == character)
            .unwrap_or_else(|| panic!("a base64 character, not {character}"));
        bits = (bits << 6) | value as u32;
        count += 6;
        if count >= 8 {
            count -= 8;
            bytes.push((bits >> count) as u8);
        }
    }
    bytes
}

/// The inputs the dump used, rebuilt.
fn long_input() -> Vec<u8> {
    (0..100_000).map(|index| (index % 251) as u8).collect()
}

/// Every case's bytes, through the port.
#[test]
fn the_bytes_come_out_unchanged() {
    let text = corpus();
    let line = b"the quick brown fox\n".to_vec();
    let binary = vec![0u8, 13, 200, 10, 9, 255];

    for (case, input, buffer, io) in [
        ("a-line-of-text", line.clone(), 512 * 1024 * 1024, 64 * 1024),
        ("nothing-at-all", Vec::new(), 512 * 1024 * 1024, 64 * 1024),
        ("a-small-buffer", line.clone(), 8, 64 * 1024),
        ("a-buffer-of-one-byte", line.clone(), 1, 64 * 1024),
        ("an-io-size-above-the-buffer", line.clone(), 8, 64),
        (
            "bytes-that-are-not-text",
            binary.clone(),
            512 * 1024 * 1024,
            64 * 1024,
        ),
        (
            "a-hundred-thousand-bytes",
            long_input(),
            512 * 1024 * 1024,
            64 * 1024,
        ),
        (
            "a-hundred-thousand-bytes-through-a-small-buffer",
            long_input(),
            1024,
            64 * 1024,
        ),
    ] {
        let produced = copy(&input, buffer, io);
        assert_eq!(produced, input, "{case}: the port changed the bytes");
        assert_eq!(produced, output(&text, case), "{case}");
    }
}

/// The buffer itself: a ring that takes what fits and gives back what is there.
#[test]
fn the_buffer_is_a_ring() {
    let mut buffer = CircularBuffer::new(4);
    assert_eq!(buffer.capacity(), 4);
    assert_eq!(buffer.write(b"abcdef"), 4);
    let mut into = [0u8; 2];
    assert_eq!(buffer.read(&mut into), 2);
    assert_eq!(&into, b"ab");
    // The two bytes read make room for two more, which wrap around the end.
    assert_eq!(buffer.write(b"gh"), 2);
    let mut rest = [0u8; 4];
    assert_eq!(buffer.read(&mut rest), 4);
    assert_eq!(&rest, b"cdgh");
    // A size of zero is a size of one: the reference has no empty buffer.
    assert_eq!(CircularBuffer::new(0).capacity(), 1);
    assert!(!buffer.is_closed());
    buffer.close();
    assert!(buffer.is_closed());
}

/// A buffer of one byte copies a hundred thousand, which is what circular means.
#[test]
fn a_buffer_of_one_byte_copies_everything() {
    let input = long_input();
    assert_eq!(copy(&input, 1, 1), input);
    assert_eq!(copy(&input, 1, 64 * 1024), input);
    // And the golden agrees about the twenty-byte case it measured.
    let text = corpus();
    assert_eq!(
        output(&text, "a-buffer-of-one-byte"),
        b"the quick brown fox\n".to_vec()
    );
}
