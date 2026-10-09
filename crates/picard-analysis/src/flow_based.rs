//! `picard.flow`: a read as the flows of a flow-based sequencer rather than as bases.
//!
//! A flow sequencer (Ultima, 454) dispenses one nucleotide at a time in a fixed cycle, the flow
//! order, and measures how many of that nucleotide were incorporated: the homopolymer length in
//! that flow. A read's KEY is that list of lengths, a zero wherever a flow added nothing. Its
//! FLOW MATRIX holds, for every flow, the probability of each length from zero to `maxHmer`, and
//! it is built from the base qualities and the `tp` tag (which says, per base, how far the true
//! length is from the called one) and, unless told to ignore it, the `t0` tag (the quality of a
//! one-to-zero error, which has no base to sit on).
//!
//! Ported from `picard.flow.FlowBasedRead`, `FlowBasedKeyCodec.baseArrayToKey`,
//! `FlowReadGroupInfo` and `FlowBasedArgumentCollection` in Picard 3.4.0.

use htsjdk_bam::header::ReadGroup;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

/// `FlowBasedRead.MAX_CLASS`, the longest homopolymer kept when the read group does not say.
pub const MAX_CLASS: i32 = 12;
/// `FlowBasedRead.MINIMAL_CALL_PROB`.
const MINIMAL_CALL_PROB: f64 = 0.1;

/// `FlowReadGroupInfo`.
#[derive(Debug, Clone)]
pub struct FlowReadGroupInfo {
    pub flow_order: Option<String>,
    pub max_class: i32,
    pub is_flow_platform: bool,
}

impl FlowReadGroupInfo {
    /// The constructor: a flow platform is Ultima or LS454, its `mc` attribute is the longest
    /// homopolymer, and an Ultima group without a flow order is malformed.
    pub fn new(group: &ReadGroup) -> Result<Self, String> {
        let platform = group.attributes.get("PL");
        let is_flow_platform = matches!(platform, Some("ULTIMA") | Some("LS454"));
        let max_class = if is_flow_platform {
            match group.attributes.get("mc") {
                None => MAX_CLASS,
                Some(text) => text.parse::<i32>().map_err(|_| {
                    format!("java.lang.NumberFormatException: For input string: \"{text}\"")
                })?,
            }
        } else {
            0
        };
        let flow_order = group.attributes.get("FO").map(str::to_string);
        if platform == Some("ULTIMA") && flow_order.is_none() {
            return Err(format!(
                "java.lang.RuntimeException: Malformed Ultima read group identified, aborting: {}",
                group_text(group)
            ));
        }
        Ok(FlowReadGroupInfo {
            flow_order,
            max_class,
            is_flow_platform,
        })
    }
}

/// `SAMReadGroupRecord.toString()`: `SAMReadGroupRecord[ID=...]` as htsjdk prints a record's
/// attributes, ID first.
fn group_text(group: &ReadGroup) -> String {
    let mut text = format!("SAMReadGroupRecord{{ID={}", group.id);
    for (key, value) in group.attributes.iter() {
        text.push_str(&format!(", {key}={value}"));
    }
    text.push('}');
    text
}

/// `FlowBasedKeyCodec.baseArrayToKey`.
pub fn base_array_to_key(bases: &[u8], flow_order: &str) -> Result<Vec<i32>, String> {
    let order = flow_order.as_bytes();
    let period = order.len();
    let mut result = Vec::new();
    let mut loc = 0usize;
    let mut flow_number = 0usize;
    let mut period_guard = 0usize;
    while loc < bases.len() {
        let flow_base = order[flow_number % period];
        if bases[loc] != flow_base && bases[loc] != b'N' {
            result.push(0);
            period_guard += 1;
            if period_guard > period {
                return Err(format!(
                    "picard.PicardException: baseArrayToKey periodGuard tripped, on {}, flowOrder: \
                     {flow_order} This probably indicates the presence of a base (value) in the \
                     sequence that is not included in the provided flow order",
                    String::from_utf8_lossy(bases)
                ));
            }
        } else {
            let mut count = 0;
            while loc < bases.len() && (bases[loc] == flow_base || bases[loc] == b'N') {
                loc += 1;
                count += 1;
            }
            result.push(count);
            period_guard = 0;
        }
        flow_number += 1;
    }
    Ok(result)
}

/// `FlowBasedArgumentCollection`.
#[derive(Debug, Clone, Copy)]
pub struct FlowArguments {
    /// `--flow-ignore-t0-tag`.
    pub ignore_t0_tag: bool,
    /// `--flow-fill-empty-bins-value`; zero means "estimate it from the read".
    pub filling_value: f64,
}

impl Default for FlowArguments {
    fn default() -> Self {
        FlowArguments {
            ignore_t0_tag: false,
            filling_value: 0.0,
        }
    }
}

/// `FlowBasedRead`: the key and the flow matrix.
#[derive(Debug, Clone)]
pub struct FlowBasedRead {
    key: Vec<i32>,
    max_hmer: i32,
    per_hmer_min_error_probability: f64,
    /// `[hmer length][flow]`.
    flow_matrix: Vec<Vec<f64>>,
}

fn pow10(exponent: f64) -> f64 {
    jmath::strict_math::pow(10.0, exponent)
}

impl FlowBasedRead {
    /// The constructor: matrix from the `tp` and `t0` tags, then the boundary flows of an
    /// unmapped or hard-clipped read spread across the lengths at least as long as the call.
    pub fn new(
        record: &BamRecord,
        flow_order: &str,
        max_hmer: i32,
        arguments: &FlowArguments,
    ) -> Result<Self, String> {
        let tp = match record.tags.get(Tag::new(b"tp")) {
            None => {
                return Err("picard.PicardException: read missing flow matrix attribute: tp".into())
            }
            Some(TagValue::ByteArray {
                values,
                unsigned: false,
            }) => values.clone(),
            Some(_) => {
                return Err(
                    "java.lang.IllegalArgumentException: tp is not a signed byte array".into(),
                )
            }
        };
        let mut read = FlowBasedRead {
            key: Vec::new(),
            max_hmer,
            per_hmer_min_error_probability: 0.0,
            flow_matrix: Vec::new(),
        };
        read.read_flow_matrix(record, flow_order, &tp, arguments)?;

        let cigar = &record.cigar;
        let unmapped = record.flags & 0x4 != 0;
        let first_is_hard_clip = cigar
            .elements
            .first()
            .is_some_and(|e| e.op == htsjdk_bam::cigar::Op::H && e.length > 0);
        let last_is_hard_clip = cigar
            .elements
            .last()
            .is_some_and(|e| e.op == htsjdk_bam::cigar::Op::H && e.length > 0);
        if unmapped || first_is_hard_clip {
            read.spread_flow_length_probs_across_counts_at_flow(first_non_zero(&read.key))?;
        }
        if unmapped || last_is_hard_clip {
            read.spread_flow_length_probs_across_counts_at_flow(last_non_zero(&read.key))?;
        }
        Ok(read)
    }

    fn spread_flow_length_probs_across_counts_at_flow(&mut self, flow: i64) -> Result<(), String> {
        if flow < 0 {
            return Ok(());
        }
        let flow = flow as usize;
        let call = self.key[flow];
        if call == 0 {
            return Err(
                "java.lang.IllegalStateException: Boundary key value should not be zero for the \
                 spreading"
                    .into(),
            );
        }
        let number_to_fill = self.max_hmer - call + 1;
        let mut total = 0.0;
        for i in call..=self.max_hmer {
            total += self.flow_matrix[i as usize][flow];
        }
        let fill_prob = (total / number_to_fill as f64).max(self.per_hmer_min_error_probability);
        for i in call..=self.max_hmer {
            self.flow_matrix[i as usize][flow] = fill_prob;
        }
        Ok(())
    }

    fn read_flow_matrix(
        &mut self,
        record: &BamRecord,
        flow_order: &str,
        tp: &[i8],
        arguments: &FlowArguments,
    ) -> Result<(), String> {
        let qualities = &record.base_qualities;
        let total_min_error_probability = if arguments.filling_value == 0.0 {
            estimate_filling_value(qualities)
        } else {
            arguments.filling_value
        };
        self.per_hmer_min_error_probability = total_min_error_probability / self.max_hmer as f64;

        self.key = base_array_to_key(&record.read_bases, flow_order)?;
        let max_hmer = self.max_hmer as usize;
        self.flow_matrix =
            vec![vec![self.per_hmer_min_error_probability; self.key.len()]; max_hmer + 1];

        // `SAMUtils.fastqToPhred(getStringAttribute("t0"))`: null when the tag is absent.
        let t0: Option<Vec<i32>> = match record.tags.get(Tag::new(b"t0")) {
            Some(TagValue::Str(text)) => {
                let mut scores = Vec::with_capacity(text.len());
                for c in text.chars() {
                    let code = c as u32;
                    if !(33..=126).contains(&code) {
                        return Err(format!(
                            "java.lang.IllegalArgumentException: Invalid fastq character: {c}"
                        ));
                    }
                    scores.push((code - 33) as i32);
                }
                Some(scores)
            }
            _ => None,
        };
        let special_treatment_for_zero_calls = t0.is_some() && !arguments.ignore_t0_tag;
        if special_treatment_for_zero_calls && t0.as_ref().map(|t| t.len()) != Some(tp.len()) {
            return Err(format!(
                "picard.PicardException: Illegal read len(t0)!=len(qual): {}",
                record.read_name
            ));
        }
        let mut probs = vec![0.0f64; qualities.len()];
        let mut t0probs = vec![0.0f64; qualities.len()];
        for i in 0..qualities.len() {
            probs[i] = pow10(-(qualities[i] as f64) / 10.0);
            if special_treatment_for_zero_calls {
                let qq = t0.as_ref().map(|t| t[i]).unwrap_or(0) as f64;
                t0probs[i] = pow10(-qq / 10.0);
            }
        }

        let mut qual_ofs = 0usize;
        for i in 0..self.key.len() {
            let run = self.key[i];
            if run > 0 {
                self.parse_single_hmer(&probs, tp, i, run, qual_ofs)?;
            }
            if run == 0 && special_treatment_for_zero_calls {
                self.parse_zero_quals(&t0probs, i, qual_ofs, total_min_error_probability);
            }
            let mut total_error_prob = 0.0;
            for k in 0..max_hmer {
                total_error_prob += self.flow_matrix[k][i];
            }
            let call_prob = MINIMAL_CALL_PROB.max(1.0 - total_error_prob);
            self.flow_matrix[(run.min(self.max_hmer)) as usize][i] = call_prob;
            qual_ofs += run as usize;
        }
        self.clip_probs();
        Ok(())
    }

    fn parse_single_hmer(
        &mut self,
        probs: &[f64],
        tp: &[i8],
        flow_idx: usize,
        flow_call: i32,
        qual_ofs: usize,
    ) -> Result<(), String> {
        for (i, &probability) in probs
            .iter()
            .enumerate()
            .skip(qual_ofs)
            .take(flow_call as usize)
        {
            let shift = *tp.get(i).ok_or_else(|| {
                format!(
                    "java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {}",
                    tp.len()
                )
            })?;
            if shift != 0 {
                let loc = (flow_call + shift as i32).min(self.max_hmer).max(0) as usize;
                if self.flow_matrix[loc][flow_idx] == self.per_hmer_min_error_probability {
                    self.flow_matrix[loc][flow_idx] = probability;
                } else {
                    self.flow_matrix[loc][flow_idx] += probability;
                }
            }
        }
        Ok(())
    }

    fn parse_zero_quals(
        &mut self,
        probs: &[f64],
        flow_idx: usize,
        qual_ofs: usize,
        total_min_error_probability: f64,
    ) {
        if qual_ofs == 0 || qual_ofs == probs.len() {
            return;
        }
        let neighbours = probs[qual_ofs - 1].min(probs[qual_ofs]);
        if neighbours <= total_min_error_probability {
            self.flow_matrix[1][flow_idx] =
                self.flow_matrix[1][flow_idx].max(self.per_hmer_min_error_probability);
        } else {
            self.flow_matrix[1][flow_idx] = self.flow_matrix[1][flow_idx].max(neighbours);
        }
    }

    fn clip_probs(&mut self) {
        for i in 0..self.max_hmer as usize {
            for j in 0..self.key.len() {
                if self.flow_matrix[i][j] <= self.per_hmer_min_error_probability
                    && self.key[j] != i as i32
                {
                    self.flow_matrix[i][j] = self.per_hmer_min_error_probability;
                }
            }
        }
    }

    pub fn max_hmer(&self) -> i32 {
        self.max_hmer
    }

    pub fn key(&self) -> &[i32] {
        &self.key
    }

    /// `getProb`: the cell, capped at one.
    pub fn prob(&self, flow: usize, hmer: i32) -> f64 {
        let row = if hmer < self.max_hmer {
            hmer
        } else {
            self.max_hmer
        };
        let prob = self.flow_matrix[row as usize][flow];
        if prob <= 1.0 {
            prob
        } else {
            1.0
        }
    }
}

/// `estimateFillingValue`: the error probability of the read's best base quality, or of 40 when
/// no base has one.
fn estimate_filling_value(qualities: &[u8]) -> f64 {
    let mut max_qual = 0.0f64;
    for &q in qualities {
        if q as f64 > max_qual {
            max_qual = q as f64;
        }
    }
    if max_qual == 0.0 {
        max_qual = 40.0;
    }
    pow10(-max_qual / 10.0)
}

fn first_non_zero(array: &[i32]) -> i64 {
    array
        .iter()
        .position(|&v| v != 0)
        .map(|i| i as i64)
        .unwrap_or(-1)
}

fn last_non_zero(array: &[i32]) -> i64 {
    array
        .iter()
        .rposition(|&v| v != 0)
        .map(|i| i as i64)
        .unwrap_or(-1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_key_counts_each_flows_homopolymer() {
        // TGCA flows over TTGAAA: two Ts, a G, no C, three As.
        assert_eq!(
            base_array_to_key(b"TTGAAA", "TGCA").unwrap(),
            vec![2, 1, 0, 3]
        );
    }

    #[test]
    fn an_n_joins_the_run_it_is_in() {
        assert_eq!(base_array_to_key(b"TNG", "TGCA").unwrap(), vec![2, 1]);
    }
}
