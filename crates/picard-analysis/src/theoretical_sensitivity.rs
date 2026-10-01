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
}

const MULTIPLIER: i64 = 0x5DEECE66D;
const ADDEND: i64 = 0xB;
const MASK: i64 = (1 << 48) - 1;

impl JavaRandom {
    /// `new Random(seed)`, which scrambles the seed before use.
    pub fn new(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ MULTIPLIER) & MASK,
        }
    }

    /// `Random.next(bits)`.
    fn next(&mut self, bits: u32) -> i32 {
        self.seed = (self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND)) & MASK;
        ((self.seed as u64) >> (48 - bits)) as i32
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
