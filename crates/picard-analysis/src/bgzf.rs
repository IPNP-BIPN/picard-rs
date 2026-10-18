//! Whole-file BGZF decompression on every core, with the serial reader as the arbiter.
//!
//! The inflating is noodles-bgzf's multithreaded reader, on libdeflate. Inflating is
//! deterministic, so the bytes cannot depend on which inflater produced them. What could differ is
//! the answer to a damaged file: noodles and htsjdk do not check the same fields (noodles checks
//! the CRC, htsjdk by default does not) or reject with the same words. So noodles only ever answers
//! for a file whose every block passes the checks `BgzfReader::next_block` makes, and whose
//! inflated length is the sum of the declared ones. Anything else hands the whole file to
//! [`htsjdk_bgzf::decompress_all`], whose answer, success or error, is the one returned.
//!
//! The declared lengths also size the output exactly, so it is allocated once.

use htsjdk_bgzf::BgzfError;
use std::io::Read;
use std::num::NonZeroUsize;
use std::sync::Arc;

const HEADER: usize = 18;
const FOOTER: usize = 8;
const MAX_BLOCK: usize = 64 * 1024;

/// [`htsjdk_bgzf::decompress_all`], in parallel. Same bytes, same errors.
///
/// It takes the compressed bytes by value: the inflating threads must own what they read, and a
/// caller that is done with the file then holds only one copy of it, the inflated one.
pub fn decompress_owned(bytes: Vec<u8>) -> Result<Vec<u8>, BgzfError> {
    let Some(total) = scan(&bytes) else {
        return htsjdk_bgzf::decompress_all(&bytes);
    };
    let bytes = Arc::new(bytes);
    let workers = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    let mut reader = noodles_bgzf::io::MultithreadedReader::with_worker_count(
        workers,
        Shared {
            bytes: Arc::clone(&bytes),
            pos: 0,
        },
    );
    // One spare byte, so a stream longer than declared shows as one rather than as a realloc.
    let mut out = Vec::with_capacity(total + 1);
    let read = reader.read_to_end(&mut out);
    drop(reader);
    if read.is_ok() && out.len() == total {
        return Ok(out);
    }
    drop(out);
    htsjdk_bgzf::decompress_all(&bytes)
}

/// The compressed bytes as the inflating threads read them, still reachable from here when they
/// fail.
struct Shared {
    bytes: Arc<Vec<u8>>,
    pos: usize,
}

impl Read for Shared {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let rest = &self.bytes[self.pos.min(self.bytes.len())..];
        let n = rest.len().min(buf.len());
        buf[..n].copy_from_slice(&rest[..n]);
        self.pos += n;
        Ok(n)
    }
}

/// [`decompress_owned`] for a caller that still needs the compressed bytes afterwards.
pub fn decompress_all(bytes: &[u8]) -> Result<Vec<u8>, BgzfError> {
    decompress_owned(bytes.to_vec())
}

/// The inflated length of a well-formed file, or `None` at the first byte the serial reader would
/// have treated as anything other than a clean block.
fn scan(bytes: &[u8]) -> Option<usize> {
    let mut total = 0usize;
    let mut p = 0usize;
    while p < bytes.len() {
        let h = bytes.get(p..p + HEADER)?;
        // The checks `BgzfReader::next_block` makes, and no others.
        if h[0] != 0x1f || h[1] != 0x8b || h[2] != 8 || h[3] != 4 {
            return None;
        }
        if u16::from_le_bytes([h[10], h[11]]) != 6 {
            return None;
        }
        let block_length = u16::from_le_bytes([h[16], h[17]]) as usize + 1;
        if !(HEADER + FOOTER..=MAX_BLOCK).contains(&block_length) {
            return None;
        }
        let end = p + block_length;
        let footer = bytes.get(end - FOOTER..end)?;
        let out_len = u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]) as usize;
        total = total.checked_add(out_len)?;
        p = end;
    }
    Some(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bgzf(chunks: &[&[u8]]) -> Vec<u8> {
        let mut w = htsjdk_bgzf::BgzfWriter::new(Vec::new());
        for chunk in chunks {
            std::io::Write::write_all(&mut w, chunk).unwrap();
            std::io::Write::flush(&mut w).unwrap();
        }
        w.into_inner().unwrap()
    }

    #[test]
    fn agrees_with_the_serial_reader_on_a_clean_file() {
        let big: Vec<u8> = (0..300_000u32).map(|i| (i * 7 % 251) as u8).collect();
        let file = bgzf(&[b"BAM\x01", &big, b"", b"tail"]);
        let ours = decompress_all(&file).unwrap();
        assert_eq!(ours, htsjdk_bgzf::decompress_all(&file).unwrap());
        assert!(
            scan(&file).is_some(),
            "a clean file takes the parallel path"
        );
    }

    #[test]
    fn empty_input_is_empty_output() {
        assert_eq!(decompress_all(&[]).unwrap(), Vec::<u8>::new());
    }

    #[test]
    fn damage_is_answered_by_the_serial_reader() {
        let file = bgzf(&[b"some bytes worth compressing, some bytes worth compressing"]);
        // Truncated mid-block, a corrupt payload, and a lying ISIZE.
        let truncated = &file[..file.len() - 30];
        let mut corrupt = file.clone();
        corrupt[HEADER + 1] ^= 0xff;
        let mut lying = file.clone();
        let first_end = u16::from_le_bytes([file[16], file[17]]) as usize + 1;
        lying[first_end - 4] = lying[first_end - 4].wrapping_add(1);
        for damaged in [truncated, &corrupt[..], &lying[..]] {
            assert_eq!(
                format!("{:?}", decompress_all(damaged)),
                format!("{:?}", htsjdk_bgzf::decompress_all(damaged))
            );
        }
    }
}
