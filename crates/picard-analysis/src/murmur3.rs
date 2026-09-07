//! `htsjdk.samtools.util.Murmur3`, for the scoring strategy that is not random at all.
//!
//! `DUPLICATE_SCORING_STRATEGY=RANDOM` does not call a random number generator: it hashes the read
//! name with Murmur3 seeded at 1 and keeps the low fourteen bits, so the same file scores the same
//! way on every run and both ends of a pair get the same score. Reproducing the choice of which
//! read in a duplicate set is kept therefore means reproducing this hash exactly.
//!
//! It is Murmur3-32 over UTF-16 code units, two at a time, which is Guava's `hashUnencodedChars`
//! and not the byte-oriented Murmur3 of the same name: the length fed to the finalizer is `2 *
//! length` in CHARACTERS, and a name of odd length mixes its last character on its own.

const C1: u32 = 0xcc9e_2d51;
const C2: u32 = 0x1b87_3593;

fn mix_k1(mut k1: u32) -> u32 {
    k1 = k1.wrapping_mul(C1);
    k1 = k1.rotate_left(15);
    k1.wrapping_mul(C2)
}

fn mix_h1(mut h1: u32, k1: u32) -> u32 {
    h1 ^= k1;
    h1 = h1.rotate_left(13);
    h1.wrapping_mul(5).wrapping_add(0xe654_6b64)
}

fn fmix(mut h1: u32, length: u32) -> u32 {
    h1 ^= length;
    h1 ^= h1 >> 16;
    h1 = h1.wrapping_mul(0x85eb_ca6b);
    h1 ^= h1 >> 13;
    h1 = h1.wrapping_mul(0xc2b2_ae35);
    h1 ^ (h1 >> 16)
}

/// `Murmur3(seed).hashUnencodedChars(input)`, as a Java `int`.
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
