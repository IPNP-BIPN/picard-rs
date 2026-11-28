//! `htsjdk.samtools.util.Murmur3`, for the scoring strategy that is not random at all.
//!
//! `DUPLICATE_SCORING_STRATEGY=RANDOM` does not call a random number generator: it hashes the read
//! name with Murmur3 seeded at 1 and keeps the low fourteen bits, so the same file scores the same
//! way on every run and both ends of a pair get the same score. Reproducing the choice of which
//! read in a duplicate set is kept therefore means reproducing this hash exactly.
//!
//! It is Murmur3-32 over UTF-16 code units, two at a time, which is Guava's `hashUnencodedChars`
//! and not Murmur3 over the name's UTF-8 bytes: the length fed to the finalizer is `2 * length`
//! in CHARACTERS, and a name of odd length mixes its last character on its own.
//!
//! The caller, htsjdk 5.0.0's `DuplicateScoringStrategy`
//! (`src/main/java/htsjdk/samtools/DuplicateScoringStrategy.java:43`):
//!
//! ```java
//! private static final Murmur3 hasher = new Murmur3(1);
//! ```
//!
//! and, in its `RANDOM` case (`:102`):
//!
//! ```java
//! score += (short) (hasher.hashUnencodedChars(record.getReadName()) & 0b11_1111_1111_1111);
//! ```
//!
//! # Why not a Rust Murmur3 crate
//!
//! A crate's Murmur3_x86_32 hashes bytes. Fed the UTF-16LE bytes of the name it gives exactly this
//! hash (two characters are one little-endian four-byte block, and an odd last character is the
//! two-byte tail), but a Rust `&str` holds UTF-8, so the natural call hashes the name's UTF-8
//! bytes: for an ASCII name, half as many bytes, and a different hash. For `"a"` with seed 1 the
//! UTF-16 hash is -810024386 (the oracle value the test below checks) and the UTF-8 one is
//! 1485495528. Spelling the algorithm out over `encode_utf16` keeps the input Java's.
//!
//! # Reading this file without reading Rust
//!
//! - The state is `u32`, not Java's `int`: the bits are the same, and `>>` on a `u32` shifts in
//!   zeros, which is Java's `>>>`.
//! - `wrapping_mul` and `wrapping_add` are Java's `int` `*` and `+`, which wrap on overflow; plain
//!   Rust `*` would panic on overflow in a debug build. `rotate_left` is `Integer.rotateLeft`.
//! - `x as u32` and `x as i32` reinterpret the 32 bits between signed and unsigned; `as u32` on a
//!   `u16` zero-extends it, as Java's `char` widens to `int`.

/// `mixK1`'s `c1` and `c2` (`Murmur3.java:94-95`).
const C1: u32 = 0xcc9e_2d51;
const C2: u32 = 0x1b87_3593;

/// `Murmur3.mixK1` (`Murmur3.java:93-100`): multiply, rotate left 15, multiply.
fn mix_k1(mut k1: u32) -> u32 {
    k1 = k1.wrapping_mul(C1);
    k1 = k1.rotate_left(15);
    k1.wrapping_mul(C2)
}

/// `Murmur3.mixH1` (`Murmur3.java:102-107`): `h1 = h1 * 5 + 0xe6546b64` after the xor and rotate.
fn mix_h1(mut h1: u32, k1: u32) -> u32 {
    h1 ^= k1;
    h1 = h1.rotate_left(13);
    h1.wrapping_mul(5).wrapping_add(0xe654_6b64)
}

/// `Murmur3.fmix` (`Murmur3.java:110-118`), the finalization mix.
fn fmix(mut h1: u32, length: u32) -> u32 {
    h1 ^= length;
    h1 ^= h1 >> 16;
    h1 = h1.wrapping_mul(0x85eb_ca6b);
    h1 ^= h1 >> 13;
    h1 = h1.wrapping_mul(0xc2b2_ae35);
    h1 ^ (h1 >> 16)
}

/// `Murmur3(seed).hashUnencodedChars(input)`, as a Java `int`.
///
/// htsjdk 5.0.0, `src/main/java/htsjdk/samtools/util/Murmur3.java:55-68`:
///
/// ```java
/// for (int i = 1; i < length; i += 2) {
///     int k1 = input.charAt(i - 1) | (input.charAt(i) << 16);
///     k1 = mixK1(k1);
///     h1 = mixH1(h1, k1);
/// }
///
/// // deal with any remaining characters
/// if ((length & 1) == 1) {
///     int k1 = input.charAt(length - 1);
///     k1 = mixK1(k1);
///     h1 ^= k1;
/// }
///
/// return fmix(h1, 2 * length);
/// ```
///
/// `length` is `input.length()`, a count of UTF-16 code units, so the name is first re-encoded
/// into the units Java's `charAt` would return.
pub fn hash_unencoded_chars(input: &str, seed: i32) -> i32 {
    let chars: Vec<u16> = input.encode_utf16().collect();
    let mut h1 = seed as u32;
    let length = chars.len();

    // Two characters per block, low character in the low half.
    let mut i = 1;
    while i < length {
        let k1 = (chars[i - 1] as u32) | ((chars[i] as u32) << 16);
        h1 = mix_h1(h1, mix_k1(k1));
        i += 2;
    }
    // The odd last character is mixed and xored in without `mixH1`, as in the Java.
    if length & 1 == 1 {
        h1 ^= mix_k1(chars[length - 1] as u32);
    }
    fmix(h1, (2 * length) as u32) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Taken from the oracle container: `new Murmur3(1).hashUnencodedChars(name)`, run against
    /// htsjdk itself rather than derived from the algorithm.
    #[test]
    fn matches_htsjdks_hash() {
        assert_eq!(hash_unencoded_chars("", 1), 1_364_076_727);
        assert_eq!(hash_unencoded_chars("a", 1), -810_024_386);
        assert_eq!(hash_unencoded_chars("read0000", 1), 123_561_743);
        assert_eq!(hash_unencoded_chars("read0314", 1), -1_542_954_484);
        assert_eq!(
            hash_unencoded_chars("INST:1:FLOWCELL:1:1101:1000:2000", 1),
            -911_918_772
        );
    }

    /// The scoring strategy keeps the low fourteen bits, which is what makes the score small.
    #[test]
    fn the_low_fourteen_bits_are_what_the_score_uses() {
        assert_eq!(
            hash_unencoded_chars("read0000", 1) & 0b11_1111_1111_1111,
            9_999
        );
        assert_eq!(
            hash_unencoded_chars("read0314", 1) & 0b11_1111_1111_1111,
            8_716
        );
    }
}
