//! The pieces of Apache Commons Math the fingerprint tools lean on.
//!
//! `CalculateFingerprintMetrics` reports chi-squared p-values and a number drawn from a hundred
//! random permutations of the fingerprint, so reproducing its rows means reproducing the library's
//! own arithmetic: the regularized incomplete gamma function behind the chi-squared distribution,
//! and the Mersenne Twister the permutations are drawn from, seeded with the constant the tool
//! does not expose.
//!
//! Ported from `org.apache.commons.math3.special.Gamma`,
//! `org.apache.commons.math3.util.ContinuedFraction`,
//! `org.apache.commons.math3.stat.inference.ChiSquareTest`,
//! `org.apache.commons.math3.distribution.ChiSquaredDistribution`,
//! `org.apache.commons.math3.distribution.GammaDistribution`,
//! `org.apache.commons.math3.distribution.UniformIntegerDistribution`,
//! `org.apache.commons.math3.random.MersenneTwister`,
//! `org.apache.commons.math3.random.RandomDataGenerator` and
//! `org.apache.commons.math3.util.MathArrays` in Commons Math 3.5.

// The constants below are transcribed digit for digit from the reference, so that they can be
// checked against it by eye. Several carry more digits than a double holds, and the compiler
// rounds them the way Java's parser does.
#![allow(clippy::excessive_precision)]

/// `Gamma.LANCZOS_G`.
const LANCZOS_G: f64 = 607.0 / 128.0;

/// `Gamma.LANCZOS`.
const LANCZOS: [f64; 15] = [
    0.999_999_999_999_997_09,
    57.156_235_665_862_923_517,
    -59.597_960_355_475_491_248,
    14.136_097_974_741_747_174,
    -0.491_913_816_097_620_199_78,
    0.339_946_499_848_118_886_99e-4,
    0.465_236_289_270_485_756_65e-4,
    -0.983_744_753_048_795_646_77e-4,
    0.158_088_703_224_912_488_84e-3,
    -0.210_264_441_724_104_883_19e-3,
    0.217_439_618_115_212_643_20e-3,
    -0.164_318_106_536_763_890_22e-3,
    0.844_182_239_838_527_432_93e-4,
    -0.261_908_384_015_814_086_70e-4,
    0.368_991_826_595_316_227_04e-5,
];

const HALF_LOG_2_PI: f64 = 0.5 * 1.837_877_066_409_345_5;

const DEFAULT_EPSILON: f64 = 10e-15;
const MAX_ITERATIONS: usize = i32::MAX as usize;

// The constants of `DGAM1`, copied from the NSWC library the way the reference copied them.
const A0: f64 = 0.611_609_510_448_141_581_788e-08;
const A1: f64 = 0.624_730_830_116_465_516_210e-08;
const B1: f64 = 0.203_610_414_066_806_987_300e+00;
const B2: f64 = 0.266_205_348_428_949_217_746e-01;
const B3: f64 = 0.493_944_979_382_446_875_238e-03;
const B4: f64 = -0.851_419_432_440_314_906_588e-05;
const B5: f64 = -0.643_045_481_779_353_022_248e-05;
const B6: f64 = 0.992_641_840_672_773_722_196e-06;
const B7: f64 = -0.607_761_895_722_825_260_739e-07;
const B8: f64 = 0.195_755_836_614_639_731_882e-09;
const P0: f64 = 0.611_609_510_448_141_581_786_1e-08;
const P1: f64 = 0.687_167_411_306_719_873_615_2e-08;
const P2: f64 = 0.682_016_166_849_617_065_791_8e-09;
const P3: f64 = 0.468_684_332_294_884_803_108_0e-10;
const P4: f64 = 0.157_283_302_771_044_628_699_5e-11;
const P5: f64 = -0.124_944_157_227_636_621_322_2e-12;
const P6: f64 = 0.434_352_993_740_859_425_517_8e-14;
const Q1: f64 = 0.305_696_107_836_522_102_500_9e+00;
const Q2: f64 = 0.546_421_308_604_229_653_601_6e-01;
const Q3: f64 = 0.495_683_009_382_588_731_202_0e-02;
const Q4: f64 = 0.269_236_946_618_636_119_287_6e-03;
const C: f64 = -0.422_784_335_098_467_139_393_487_909_917_598e+00;
const C0: f64 = 0.577_215_664_901_532_860_606_512_090_082_402e+00;
const C1: f64 = -0.655_878_071_520_253_881_077_019_515_145_390e+00;
const C2: f64 = -0.420_026_350_340_952_355_290_039_348_754_298e-01;
const C3: f64 = 0.166_538_611_382_291_489_501_700_795_102_105e+00;
const C4: f64 = -0.421_977_345_555_443_367_482_083_012_891_874e-01;
const C5: f64 = -0.962_197_152_787_697_356_211_492_167_234_820e-02;
const C6: f64 = 0.721_894_324_666_309_954_239_501_034_044_657e-02;
const C7: f64 = -0.116_516_759_185_906_511_211_397_108_401_839e-02;
const C8: f64 = -0.215_241_674_114_950_972_815_729_963_053_648e-03;
const C9: f64 = 0.128_050_282_388_116_186_153_198_626_328_164e-03;
const C10: f64 = -0.201_348_547_807_882_386_556_893_914_210_218e-04;
const C11: f64 = -0.125_049_348_214_267_065_734_535_947_383_309e-05;
const C12: f64 = 0.113_302_723_198_169_588_237_412_962_033_074e-05;
const C13: f64 = -0.205_633_841_697_760_710_345_015_413_002_057e-06;

/// `Gamma.invGamma1pm1`: one over gamma of one plus x, less one, for x in `[-0.5, 1.5]`.
pub fn inv_gamma_1p_m1(x: f64) -> f64 {
    let t = if x <= 0.5 { x } else { (x - 0.5) - 0.5 };
    if t < 0.0 {
        let a = A0 + t * A1;
        let mut b = B8;
        for coefficient in [B7, B6, B5, B4, B3, B2, B1] {
            b = coefficient + t * b;
        }
        b = 1.0 + t * b;
        let mut c = C13 + t * (a / b);
        for coefficient in [C12, C11, C10, C9, C8, C7, C6, C5, C4, C3, C2, C1, C] {
            c = coefficient + t * c;
        }
        if x > 0.5 {
            t * c / x
        } else {
            x * ((c + 0.5) + 0.5)
        }
    } else {
        let mut p = P6;
        for coefficient in [P5, P4, P3, P2, P1, P0] {
            p = coefficient + t * p;
        }
        let mut q = Q4;
        for coefficient in [Q3, Q2, Q1] {
            q = coefficient + t * q;
        }
        q = 1.0 + t * q;
        let mut c = C13 + (p / q) * t;
        for coefficient in [C12, C11, C10, C9, C8, C7, C6, C5, C4, C3, C2, C1, C0] {
            c = coefficient + t * c;
        }
        if x > 0.5 {
            (t / x) * ((c - 0.5) - 0.5)
        } else {
            x * c
        }
    }
}

/// `Gamma.logGamma1p`.
pub fn log_gamma_1p(x: f64) -> f64 {
    -inv_gamma_1p_m1(x).ln_1p()
}

/// `Gamma.lanczos`.
fn lanczos(x: f64) -> f64 {
    let mut sum = 0.0;
    for index in (1..LANCZOS.len()).rev() {
        sum += LANCZOS[index] / (x + index as f64);
    }
    sum + LANCZOS[0]
}

/// `Gamma.logGamma`, by the same four branches the reference uses.
pub fn log_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 0.5 {
        return log_gamma_1p(x) - x.ln();
    }
    if x <= 2.5 {
        return log_gamma_1p((x - 0.5) - 0.5);
    }
    if x <= 8.0 {
        let n = (x - 1.5).floor() as i32;
        let mut product = 1.0;
        for index in 1..=n {
            product *= x - f64::from(index);
        }
        return log_gamma_1p(x - f64::from(n + 1)) + product.ln();
    }
    let sum = lanczos(x);
    let tmp = x + LANCZOS_G + 0.5;
    ((x + 0.5) * tmp.ln()) - tmp + HALF_LOG_2_PI + (sum / x).ln()
}

/// `ContinuedFraction.evaluate`, for the fraction the upper gamma uses.
fn continued_fraction(a: f64, x: f64, epsilon: f64, max_iterations: usize) -> f64 {
    let small = 1e-50;
    let term_a = |n: usize| ((2.0 * n as f64) + 1.0) - a + x;
    let term_b = |n: usize| n as f64 * (a - n as f64);

    let mut h_previous = term_a(0);
    if (h_previous - 0.0).abs() < small {
        h_previous = small;
    }
    let mut n = 1;
    let mut d_previous = 0.0;
    let mut c_previous = h_previous;
    let mut h = h_previous;
    while n < max_iterations {
        let mut d = term_a(n) + term_b(n) * d_previous;
        if (d - 0.0).abs() < small {
            d = small;
        }
        let mut c = term_a(n) + term_b(n) / c_previous;
        if (c - 0.0).abs() < small {
            c = small;
        }
        d = 1.0 / d;
        let delta = c * d;
        h = h_previous * delta;
        if (delta - 1.0).abs() < epsilon {
            break;
        }
        d_previous = d;
        c_previous = c;
        h_previous = h;
        n += 1;
    }
    h
}

/// `Gamma.regularizedGammaP`, the lower incomplete gamma as a fraction of the whole.
pub fn regularized_gamma_p(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x >= a + 1.0 {
        // The upper one converges faster there, and the reference hands the work over.
        return 1.0 - regularized_gamma_q(a, x);
    }
    let mut n = 0.0;
    let mut an = 1.0 / a;
    let mut sum = an;
    while (an / sum).abs() > DEFAULT_EPSILON && n < MAX_ITERATIONS as f64 && sum.is_finite() {
        n += 1.0;
        an *= x / (a + n);
        sum += an;
    }
    if sum.is_infinite() {
        return 1.0;
    }
    (-x + (a * x.ln()) - log_gamma(a)).exp() * sum
}

/// `Gamma.regularizedGammaQ`, the upper tail.
pub fn regularized_gamma_q(a: f64, x: f64) -> f64 {
    if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 1.0;
    }
    if x < a + 1.0 {
        return 1.0 - regularized_gamma_p(a, x);
    }
    let fraction = 1.0 / continued_fraction(a, x, DEFAULT_EPSILON, MAX_ITERATIONS);
    (-x + (a * x.ln()) - log_gamma(a)).exp() * fraction
}

/// `ChiSquaredDistribution.cumulativeProbability`, which is the gamma distribution's with a shape
/// of half the degrees of freedom and a scale of two.
pub fn chi_squared_cumulative_probability(x: f64, degrees_of_freedom: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    regularized_gamma_p(degrees_of_freedom / 2.0, x / 2.0)
}

/// `ChiSquareTest.chiSquare`: the statistic, with the expectations rescaled to the observations
/// where the two do not already sum to the same number.
pub fn chi_square(expected: &[f64], observed: &[i64]) -> f64 {
    let sum_expected: f64 = expected.iter().sum();
    let sum_observed: f64 = observed.iter().map(|count| *count as f64).sum();
    let rescale = (sum_expected - sum_observed).abs() > 10e-6;
    let ratio = if rescale {
        sum_observed / sum_expected
    } else {
        1.0
    };
    let mut sum_of_squares = 0.0;
    for (count, expectation) in observed.iter().zip(expected.iter()) {
        let expectation = if rescale {
            ratio * expectation
        } else {
            *expectation
        };
        let deviation = *count as f64 - expectation;
        sum_of_squares += deviation * deviation / expectation;
    }
    sum_of_squares
}

/// `ChiSquareTest.chiSquareTest`: the p-value, on one degree of freedom fewer than there are bins.
pub fn chi_square_test(expected: &[f64], observed: &[i64]) -> f64 {
    let statistic = chi_square(expected, observed);
    1.0 - chi_squared_cumulative_probability(statistic, expected.len() as f64 - 1.0)
}

/// `MersenneTwister`, which is what the fingerprint metrics' sampling is drawn from.
pub struct MersenneTwister {
    state: [i32; MersenneTwister::N],
    index: usize,
}

impl MersenneTwister {
    const N: usize = 624;
    const M: usize = 397;
    const MAG01: [i32; 2] = [0x0, -0x66f7_4f21];

    /// `setSeed(int)`, the seeding of the 2002 C version.
    pub fn new(seed: i32) -> MersenneTwister {
        let mut state = [0i32; MersenneTwister::N];
        let mut value = i64::from(seed);
        state[0] = value as i32;
        for (index, slot) in state.iter_mut().enumerate().skip(1) {
            value = (1_812_433_253i64 * (value ^ (value >> 30)) + index as i64) & 0xffff_ffff;
            *slot = value as i32;
        }
        MersenneTwister {
            state,
            index: MersenneTwister::N,
        }
    }

    /// `next(int bits)`: the twist, the tempering, and the top bits.
    pub fn next_bits(&mut self, bits: u32) -> i32 {
        if self.index >= MersenneTwister::N {
            let mut next = self.state[0];
            for k in 0..MersenneTwister::N - MersenneTwister::M {
                let current = next;
                next = self.state[k + 1];
                let y = (current & i32::MIN) | (next & i32::MAX);
                self.state[k] = self.state[k + MersenneTwister::M]
                    ^ ((y as u32) >> 1) as i32
                    ^ MersenneTwister::MAG01[(y & 0x1) as usize];
            }
            for k in MersenneTwister::N - MersenneTwister::M..MersenneTwister::N - 1 {
                let current = next;
                next = self.state[k + 1];
                let y = (current & i32::MIN) | (next & i32::MAX);
                self.state[k] = self.state[k + MersenneTwister::M - MersenneTwister::N]
                    ^ ((y as u32) >> 1) as i32
                    ^ MersenneTwister::MAG01[(y & 0x1) as usize];
            }
            let y = (next & i32::MIN) | (self.state[0] & i32::MAX);
            self.state[MersenneTwister::N - 1] = self.state[MersenneTwister::M - 1]
                ^ ((y as u32) >> 1) as i32
                ^ MersenneTwister::MAG01[(y & 0x1) as usize];
            self.index = 0;
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= ((y as u32) >> 11) as i32;
        y ^= (y << 7) & -0x62d3_a980;
        y ^= (y << 15) & -0x1039_a000;
        y ^= ((y as u32) >> 18) as i32;
        ((y as u32) >> (32 - bits)) as i32
    }

    /// `BitsStreamGenerator.nextInt(int n)`.
    pub fn next_int(&mut self, n: i32) -> i32 {
        if n & -n == n {
            return ((i64::from(n) * i64::from(self.next_bits(31))) >> 31) as i32;
        }
        loop {
            let bits = self.next_bits(31);
            let value = bits % n;
            if bits - value + (n - 1) >= 0 {
                return value;
            }
        }
    }
}

/// `RandomDataGenerator.nextPermutation`, which shuffles the naturals from the tail down.
///
/// Each step draws a target from the range that is still unplaced, and the draw is a uniform
/// integer, which is one number off the generator per step except the last, which takes none.
pub fn next_permutation(rng: &mut MersenneTwister, n: usize) -> Vec<usize> {
    let mut index: Vec<usize> = (0..n).collect();
    for position in (0..n).rev() {
        let target = if position == 0 {
            0
        } else {
            rng.next_int(position as i32 + 1) as usize
        };
        index.swap(target, position);
    }
    index
}
