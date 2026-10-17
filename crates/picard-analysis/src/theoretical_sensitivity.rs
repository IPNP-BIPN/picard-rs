//! `TheoreticalSensitivity.hetSNPSensitivity`, the `HET_SNP_SENSITIVITY` column of every
//! `WgsMetrics` row.
//!
//! The figure is a Monte-Carlo estimate, so reproducing it means reproducing the random stream:
//! `RouletteWheel` draws through a `java.util.Random` seeded with 51, and every draw it makes is
//! two `nextDouble` calls (one to pick a bin, one to accept it) until a bin is accepted or 600
//! attempts have failed. The draws are consumed in the order the Java consumes them -- sample by
//! sample, summand by summand -- because a different order is a different estimate.
//!
//! Ported from `picard.analysis.TheoreticalSensitivity` (Picard 3.4.0) and `java.util.Random`.

/// `java.util.Random`: the 48-bit linear congruential generator.
#[derive(Debug, Clone)]
pub struct JavaRandom {
    seed: i64,
    next_next_gaussian: Option<f64>,
}

const MULTIPLIER: i64 = 0x5DEECE66D;
const ADDEND: i64 = 0xB;
const MASK: i64 = (1 << 48) - 1;

impl JavaRandom {
    /// `new Random(seed)`, which scrambles the seed before use.
    pub fn new(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ MULTIPLIER) & MASK,
            next_next_gaussian: None,
        }
    }

    /// `Random.next(bits)`.
    fn next(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND)) & MASK;
        ((self.seed as u64) >> (48 - bits)) as i32
    }

    /// `Random.nextInt()`: the next 32 bits, signed.
    pub fn next_int_unbounded(&mut self) -> i32 {
        self.next(32)
    }

    /// `Random.nextInt(bound)`: a power of two takes the high bits, anything else rejects the
    /// draws that would bias the remainder.
    pub fn next_int(&mut self, bound: i32) -> i32 {
        assert!(bound > 0, "bound must be positive");
        if bound & -bound == bound {
            return ((bound as i64 * self.next(31) as i64) >> 31) as i32;
        }
        loop {
            let bits = self.next(31);
            let val = bits % bound;
            if bits.wrapping_sub(val).wrapping_add(bound - 1) >= 0 {
                return val;
            }
        }
    }

    /// `Random.nextGaussian`: the polar method, which draws PAIRS and hands out the second one on
    /// the next call.
    pub fn next_gaussian(&mut self) -> f64 {
        if let Some(cached) = self.next_next_gaussian.take() {
            return cached;
        }
        loop {
            let v1 = 2.0 * self.next_double() - 1.0;
            let v2 = 2.0 * self.next_double() - 1.0;
            let s = v1 * v1 + v2 * v2;
            if s < 1.0 && s != 0.0 {
                let multiplier = (-2.0 * jmath::strict_log::log(s) / s).sqrt();
                self.next_next_gaussian = Some(v2 * multiplier);
                return v1 * multiplier;
            }
        }
    }

    /// `Random.nextDouble`: 53 random bits scaled into [0, 1).
    pub fn next_double(&mut self) -> f64 {
        let high = (self.next(26) as i64) << 27;
        let low = self.next(27) as i64;
        (high + low) as f64 * (1.0f64 / (1u64 << 53) as f64)
    }
}

const SAMPLING_MAX: i32 = 600;
const MAX_CONSIDERED_DEPTH_HET_SENS: usize = 1000;
const RANDOM_SEED: i64 = 51;

/// `MathUtil.max`, which is `nums[indexOfMax(nums)]`: a strict `>` scan from the first element,
/// so a leading NaN is never replaced and is what comes back.
fn java_max(nums: &[f64]) -> f64 {
    let mut max = nums[0];
    for &n in &nums[1..] {
        if n > max {
            max = n;
        }
    }
    max
}

/// `TheoreticalSensitivity.RouletteWheel`.
struct RouletteWheel {
    probabilities: Vec<f64>,
    n: usize,
    count: i32,
    rng: JavaRandom,
}

impl RouletteWheel {
    /// Refuses an all-zero distribution. An all-NaN one (a histogram with no observations,
    /// normalized by a zero sum) is NOT refused: its maximum is NaN, which is not `== 0`.
    fn new(weights: &[f64]) -> Result<Self, String> {
        let w_max = java_max(weights);
        if w_max == 0.0 {
            return Err("picard.PicardException: Quality score distribution is empty.".to_string());
        }
        Ok(RouletteWheel {
            probabilities: weights.iter().map(|w| w / w_max).collect(),
            n: weights.len(),
            count: 0,
            rng: JavaRandom::new(RANDOM_SEED),
        })
    }

    fn draw(&mut self) -> i32 {
        loop {
            let n = (self.n as f64 * self.rng.next_double()) as i32;
            self.count += 1;
            if self.rng.next_double() < self.probabilities[n as usize] {
                self.count = 0;
                return n;
            } else if self.count >= SAMPLING_MAX {
                self.count = 0;
                return 0;
            }
        }
    }

    fn sample_cumulative_sums(&mut self, max_summands: usize, sample_size: i32) -> Vec<Vec<i32>> {
        let mut result: Vec<Vec<i32>> = (0..max_summands)
            .map(|_| Vec::with_capacity(sample_size.max(0) as usize))
            .collect();
        for _ in 0..sample_size {
            let mut cumulative: i32 = 0;
            for row in result.iter_mut() {
                row.push(cumulative);
                cumulative = cumulative.wrapping_add(self.draw());
            }
        }
        result
    }
}

/// `TheoreticalSensitivity.normalizeHistogram`, over a histogram whose bins are `0..len`.
pub fn normalize(values: &[f64]) -> Vec<f64> {
    let sum: f64 = values.iter().sum();
    values.iter().map(|v| v / sum).collect()
}

fn proportions_above_thresholds(lists: &mut [Vec<i32>], thresholds: &[f64]) -> Vec<Vec<f64>> {
    let mut result = Vec::with_capacity(lists.len());
    for list in lists.iter_mut() {
        let mut row = vec![0.0; thresholds.len()];
        list.sort_unstable();
        let mut n = 0;
        let mut j = 0;
        while n < thresholds.len() && j < list.len() {
            if thresholds[n] > list[j] as f64 {
                j += 1;
            } else {
                row[n] = (list.len() - j) as f64 / list.len() as f64;
                n += 1;
            }
        }
        result.push(row);
    }
    result
}

fn het_alt_depth_distribution(n_max: usize) -> Vec<Vec<f64>> {
    let mut table: Vec<Vec<f64>> = Vec::with_capacity(n_max);
    for n in 0..n_max {
        let mut row = Vec::with_capacity(n + 1);
        row.push(0.5f64.powi(n as i32));
        for m in 1..n {
            let value = (n as f64 * 0.5 / m as f64) * table[n - 1][m - 1];
            row.push(value);
        }
        if n > 0 {
            row.push(row[0]);
        }
        table.push(row);
    }
    table
}

/// `TheoreticalSensitivity.hetSNPSensitivity(depth, quality, sampleSize, logOddsThreshold)`.
pub fn het_snp_sensitivity(
    depth_distribution: &[f64],
    quality_distribution: &[f64],
    sample_size: i32,
    log_odds_threshold: f64,
) -> Result<f64, String> {
    let n_max = depth_distribution
        .len()
        .min(MAX_CONSIDERED_DEPTH_HET_SENS + 1);
    let mut sampler = RouletteWheel::new(quality_distribution)?;
    let mut quality_sums = sampler.sample_cumulative_sums(n_max, sample_size);
    let log_10 = 2f64.log10();
    let thresholds: Vec<f64> = (0..n_max)
        .map(|n| 10.0 * (n as f64 * log_10 + log_odds_threshold))
        .collect();
    let probability = proportions_above_thresholds(&mut quality_sums, &thresholds);
    let alt = het_alt_depth_distribution(n_max);
    let mut result = 0.0;
    for n in 0..n_max {
        for m in 0..=n {
            result += depth_distribution[n] * alt[n][m] * probability[m][n];
        }
    }
    Ok(result)
}

/// `QualityUtil.getPhredScoreFromErrorProbability`: `(int) Math.round(-10 * log10(p))`.
///
/// `Math.round` returns a `long` and the cast keeps its low 32 bits, so a probability of zero --
/// a sensitivity of exactly one -- gives `Long.MAX_VALUE`, which is `-1` as an `int`.
pub fn phred_from_error_probability(probability: f64) -> i32 {
    java_round(-10.0 * probability.log10()) as i32
}

/// `Math.round(double)`: `floor(x + 1/2)` computed exactly, NaN to zero, saturating.
pub fn java_round(x: f64) -> i64 {
    if x.is_nan() {
        return 0;
    }
    if x >= i64::MAX as f64 {
        return i64::MAX;
    }
    if x <= i64::MIN as f64 {
        return i64::MIN;
    }
    let floor = x.floor();
    if x - floor >= 0.5 {
        floor as i64 + 1
    } else {
        floor as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_random_matches_the_jdk_sequence() {
        // new Random(42).nextDouble() == 0.7275636800328681 on every JDK.
        let mut rng = JavaRandom::new(42);
        assert_eq!(rng.next_double(), 0.7275636800328681);
    }

    #[test]
    fn rounding_follows_math_round() {
        assert_eq!(java_round(2.5), 3);
        assert_eq!(java_round(-2.5), -2);
        assert_eq!(java_round(f64::NAN), 0);
        assert_eq!(phred_from_error_probability(0.0), -1);
    }
}

/// `org.apache.commons.math3.random.Well19937c` (Commons Math 3.5), as `new Well19937c(long)`
/// builds it: the seed split into two ints and spread over the pool by the Mersenne Twister's
/// initialisation recurrence.
pub struct Well19937c {
    v: Vec<i32>,
    index: usize,
    i_rm1: Vec<usize>,
    i_rm2: Vec<usize>,
    i1: Vec<usize>,
    i2: Vec<usize>,
    i3: Vec<usize>,
}

impl Well19937c {
    pub fn new(seed: i64) -> Self {
        const K: usize = 19937;
        const R: usize = K.div_ceil(32);
        let (m1, m2, m3) = (70, 179, 449);
        let mut v = vec![0i32; R];
        let init = [(seed as u64 >> 32) as i32, seed as i32];
        v[..2].copy_from_slice(&init);
        for i in 2..R {
            let l = v[i - 2] as i64;
            v[i] = (1_812_433_253i64
                .wrapping_mul(l ^ (l >> 30))
                .wrapping_add(i as i64)
                & 0xffff_ffff) as i32;
        }
        Well19937c {
            v,
            index: 0,
            i_rm1: (0..R).map(|j| (j + R - 1) % R).collect(),
            i_rm2: (0..R).map(|j| (j + R - 2) % R).collect(),
            i1: (0..R).map(|j| (j + m1) % R).collect(),
            i2: (0..R).map(|j| (j + m2) % R).collect(),
            i3: (0..R).map(|j| (j + m3) % R).collect(),
        }
    }

    /// `Well19937c.next(bits)`, the Matsumoto-Kurita tempering included.
    fn next(&mut self, bits: u32) -> i32 {
        let index_rm1 = self.i_rm1[self.index];
        let index_rm2 = self.i_rm2[self.index];
        let v0 = self.v[self.index];
        let vm1 = self.v[self.i1[self.index]];
        let vm2 = self.v[self.i2[self.index]];
        let vm3 = self.v[self.i3[self.index]];
        let ushr = |x: i32, n: u32| ((x as u32) >> n) as i32;
        let z0 = (0x8000_0000u32 as i32 & self.v[index_rm1]) ^ (0x7FFF_FFFF & self.v[index_rm2]);
        let z1 = (v0 ^ (v0 << 25)) ^ (vm1 ^ ushr(vm1, 27));
        let z2 = ushr(vm2, 9) ^ (vm3 ^ ushr(vm3, 1));
        let z3 = z1 ^ z2;
        let mut z4 = z0 ^ (z1 ^ (z1 << 9)) ^ (z2 ^ (z2 << 21)) ^ (z3 ^ ushr(z3, 21));
        self.v[self.index] = z3;
        self.v[index_rm1] = z4;
        self.v[index_rm2] &= 0x8000_0000u32 as i32;
        self.index = index_rm1;
        z4 ^= (z4 << 7) & 0xe46e_1700u32 as i32;
        z4 ^= (z4 << 15) & 0x9b86_8000u32 as i32;
        ushr(z4, 32 - bits)
    }

    /// `BitsStreamGenerator.nextDouble`: two 26-bit draws, 52 bits in all.
    pub fn next_double(&mut self) -> f64 {
        let high = (self.next(26) as i64) << 26;
        let low = self.next(26) as i64;
        (high | low) as f64 * f64::powi(2.0, -52)
    }
}

/// `TheoreticalSensitivity.trimDistribution`: the trailing zeros removed.
fn trim_distribution(distribution: &[f64]) -> &[f64] {
    let mut end = distribution.len();
    while end > 0 && distribution[end - 1] == 0.0 {
        end -= 1;
    }
    &distribution[..end]
}

/// `TheoreticalSensitivity.sensitivityAtConstantDepth` with its seed of 51: simulated variants at
/// a depth, their alt counts binomial (inverted from a `Well19937c`), their summed qualities drawn
/// from the quality distribution below ten alt reads and from a Gaussian above.
fn sensitivity_at_constant_depth(
    depth: i32,
    quality: &JavaHistogram,
    log_odds_threshold: f64,
    sample_size: i32,
    allele_fraction: f64,
) -> Result<f64, String> {
    if depth == 0 {
        return Ok(0.0);
    }
    let mut wheel = RouletteWheel::new(trim_distribution(&quality.normalized()))?;
    let mut gaussian = JavaRandom::new(RANDOM_SEED);
    let mut uniform = Well19937c::new(RANDOM_SEED);
    let (average, deviation) = (quality.mean, quality.deviation);
    let mut called = 0;
    for _ in 0..sample_size {
        let alt_depth = jmath::binomial::inverse_cumulative_probability(
            depth,
            allele_fraction,
            uniform.next_double(),
        )
        .map_err(|e| format!("{e:?}"))?;
        let sum_of_qualities: i32 = if alt_depth < 10 {
            (0..alt_depth).map(|_| wheel.draw()).sum()
        } else {
            (alt_depth as f64 * average
                + gaussian.next_gaussian() * (alt_depth as f64).sqrt() * deviation)
                as i32
        };
        let threshold = 10.0
            * (alt_depth as f64 * -jmath::math::log10(allele_fraction)
                + (depth - alt_depth) as f64 * -jmath::math::log10(1.0 - allele_fraction)
                + log_odds_threshold);
        if sum_of_qualities as f64 > threshold {
            called += 1;
        }
    }
    Ok(called as f64 / sample_size as f64)
}

/// `TheoreticalSensitivityMetrics`.
#[derive(Debug, Clone, PartialEq)]
pub struct SensitivityRow {
    pub allele_fraction: f64,
    pub theoretical_sensitivity: f64,
    pub theoretical_sensitivity_q: i32,
    pub log_odds_threshold: f64,
    pub sample_size: i32,
}

/// A `Histogram<Integer>` as `TheoreticalSensitivity` reads it: `normalizeHistogram`'s view, which
/// is `get(i) / getSumOfValues()` for `i` below the number of BINS (so a histogram with a gap loses
/// its top keys there but not from the sum), and the mean and deviation over every bin.
pub struct JavaHistogram {
    pub values: Vec<f64>,
    pub sum: f64,
    pub mean: f64,
    pub deviation: f64,
}

impl JavaHistogram {
    /// From the bins, in key order.
    pub fn from_bins(bins: &[(f64, f64)]) -> Self {
        let by_key: std::collections::HashMap<i64, f64> =
            bins.iter().map(|(k, v)| (*k as i64, *v)).collect();
        let values = (0..bins.len())
            .map(|i| by_key.get(&(i as i64)).copied().unwrap_or(0.0))
            .collect();
        let count: f64 = bins.iter().map(|(_, v)| v).sum();
        let mean = bins.iter().map(|(k, v)| k * v).sum::<f64>() / count;
        let squares: f64 = bins.iter().map(|(k, v)| v * (k - mean).powi(2)).sum();
        JavaHistogram {
            values,
            sum: count,
            mean,
            deviation: (squares / (count - 1.0)).sqrt(),
        }
    }

    /// From an array with a bin at every index.
    pub fn from_dense(values: &[i64]) -> Self {
        let bins: Vec<(f64, f64)> = values
            .iter()
            .enumerate()
            .map(|(k, &v)| (k as f64, v as f64))
            .collect();
        Self::from_bins(&bins)
    }

    fn normalized(&self) -> Vec<f64> {
        self.values.iter().map(|v| v / self.sum).collect()
    }
}

/// `TheoreticalSensitivity.calculateSensitivities` over a depth and a base-quality histogram with a
/// bin at every index, which is how the WGS collectors hand them over.
pub fn calculate_sensitivities(
    sample_size: i32,
    depth: &[i64],
    quality: &[i64],
    allele_fractions: &[f64],
) -> Result<Vec<SensitivityRow>, String> {
    calculate_sensitivities_with(
        sample_size,
        &JavaHistogram::from_dense(depth),
        &JavaHistogram::from_dense(quality),
        allele_fractions,
    )
}

/// The same over `normalizeHistogram`'s view of a depth histogram, for a collector whose
/// histograms have bins only where something was counted.
pub fn calculate_sensitivities_with(
    sample_size: i32,
    depth: &JavaHistogram,
    quality: &JavaHistogram,
    allele_fractions: &[f64],
) -> Result<Vec<SensitivityRow>, String> {
    const DEPTH_BIN_WIDTH: f64 = 0.01;
    let log_odds_threshold = 6.2;
    let distribution = depth.normalized();
    let mut rows = Vec::new();
    for &allele_fraction in allele_fractions {
        if !(0.0..=1.0).contains(&allele_fraction) {
            return Err(
                "java.lang.IllegalArgumentException: Allele fractions must be between 0 and 1."
                    .to_string(),
            );
        }
        let mut sensitivity = 0.0;
        let mut current = 0usize;
        let mut right = 0.0;
        while current < distribution.len() {
            let mut delta = 0.0;
            while (delta == 0.0 && current < distribution.len())
                || (delta < DEPTH_BIN_WIDTH
                    && current < distribution.len()
                    && distribution[current] < DEPTH_BIN_WIDTH / 2.0)
            {
                delta += distribution[current];
                current += 1;
            }
            let left = right;
            right = sensitivity_at_constant_depth(
                current as i32,
                quality,
                log_odds_threshold,
                sample_size,
                allele_fraction,
            )?;
            sensitivity += delta * (left + right) / 2.0;
        }
        rows.push(SensitivityRow {
            allele_fraction,
            theoretical_sensitivity: sensitivity,
            theoretical_sensitivity_q: phred_from_error_probability(1.0 - sensitivity),
            log_odds_threshold,
            sample_size,
        });
    }
    Ok(rows)
}

#[cfg(test)]
mod well_tests {
    use super::*;

    /// The first draws of `new Well19937c(51L)`, of a `BinomialDistribution(37, 0.1)` sampled from
    /// another, and of `new Random(51).nextGaussian()`, printed by Commons Math 3.5 and JDK 17 in
    /// the oracle container.
    #[test]
    fn the_generators_draw_what_the_reference_draws() {
        let mut well = Well19937c::new(51);
        let doubles: Vec<f64> = (0..5).map(|_| well.next_double()).collect();
        assert_eq!(
            doubles,
            [
                0.6236127482818385,
                0.31954229124778966,
                0.960406052784996,
                0.3707656267778141,
                0.541059975047363
            ]
        );
        let mut uniform = Well19937c::new(51);
        let samples: Vec<i32> = (0..20)
            .map(|_| {
                jmath::binomial::inverse_cumulative_probability(37, 0.1, uniform.next_double())
                    .unwrap()
            })
            .collect();
        assert_eq!(
            samples,
            [4, 3, 7, 3, 4, 8, 2, 5, 4, 5, 2, 6, 3, 3, 3, 4, 4, 6, 4, 5]
        );
        let mut random = JavaRandom::new(51);
        let gaussians: Vec<f64> = (0..3).map(|_| random.next_gaussian()).collect();
        assert_eq!(
            gaussians,
            [0.7059203462238288, 0.8721702467442611, -0.8754653493198823]
        );
    }
}
