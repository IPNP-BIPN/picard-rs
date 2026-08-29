//! `CombineGenotypingArrayVcfs`: the lockstep an array VCF merge walks in.
//!
//! The tool puts single-sample array VCFs side by side into one multi-sample file. It is not a
//! merge by position: the nth variant of each input is merged with the nth of every other, so two
//! files holding the same loci in a different order are a refusal rather than a reordering.
//!
//! What is ported is the deciding: which pairs line up, which attributes have to agree, what the
//! merged header keeps, and the twelve refusals around all of that.
//!
//! Ported from `picard.arrays.CombineGenotypingArrayVcfs`.

/// The attributes the agreement check skips, which may differ freely between inputs.
///
/// Three of them are recalculated after the merge and the other four are documented as carrying
/// minor allowable differences.
pub const EXEMPT_ATTRIBUTES: [&str; 7] =
    ["AC", "AF", "AN", "devX_AB", "devY_AB", "SOURCE", "refSNP"];

/// The header lines the merged file drops, each of them one sample's own.
pub const SAMPLE_SPECIFIC_HEADER_LINES: [&str; 9] = [
    "pipelineVersion",
    "analysisVersionNumber",
    "autocallDate",
    "autocallGender",
    "chipWellBarcode",
    "expectedGender",
    "extendedIlluminaManifestVersion",
    "fingerprintGender",
    "gtcCallRate",
];

/// Why a merge stopped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Refusal {
    /// A sample name appeared in two inputs, named with the file that repeated it.
    RepeatedSample {
        file: String,
        sample: String,
    },
    /// One input ran out of variants before another.
    VariantCount,
    /// Two variants of one step sit at different loci.
    Locus,
    Id,
    ReferenceAllele,
    AlternateAlleleCount,
    /// The alternate allele, which is the one refusal that names where it happened.
    AlternateAllele {
        contig: String,
        start: i64,
    },
    /// An attribute one input has and another does not.
    AttributeMissing(String),
    /// An attribute whose values disagree.
    AttributeDisagrees(String),
    /// The index was asked for and no sequence dictionary was reachable.
    NoSequenceDictionary,
    /// The depth: the one attribute the merge adds up, and the one that makes the tool throw.
    ///
    /// The sum is written back into the map `VariantContext.getAttributes` returned, which is
    /// unmodifiable, so the run ends in an `UnsupportedOperationException` carrying no message at
    /// all. No depth ever reaches an output.
    DepthIsUnwritable,
}

impl Refusal {
    /// The message the reference writes, which is empty for the depth because the exception the
    /// reference throws there has none.
    pub fn message(&self) -> String {
        match self {
            Refusal::RepeatedSample { file, sample } => format!(
                "Input file {file} contains a sample entry ({sample}) that appears in another \
                 input file."
            ),
            Refusal::VariantCount => "Mismatch in number of variants among input VCFs".to_string(),
            Refusal::Locus => "Mismatch in loci among input VCFs".to_string(),
            Refusal::Id => "Mismatch in ID field among input VCFs".to_string(),
            Refusal::ReferenceAllele => "Mismatch in REF allele among input VCFs".to_string(),
            Refusal::AlternateAlleleCount => {
                "Mismatch in ALT allele count among input VCFs".to_string()
            }
            Refusal::AlternateAllele { contig, start } => {
                format!("Mismatch in ALT allele among input VCFs for {contig}.{start}")
            }
            Refusal::AttributeMissing(key) => format!("Attribute '{key}' not found in all VCFs"),
            Refusal::AttributeDisagrees(key) => {
                format!("Values for attribute '{key}' disagrees among input VCFs")
            }
            Refusal::NoSequenceDictionary => "A sequence dictionary must be available (either \
                 through the input file or by setting it explicitly) when creating indexed output."
                .to_string(),
            Refusal::DepthIsUnwritable => String::new(),
        }
    }

    /// The exception class the reference throws, which is not the same for all of them.
    pub fn class(&self) -> &'static str {
        match self {
            Refusal::RepeatedSample { .. } => "java.lang.IllegalArgumentException",
            Refusal::DepthIsUnwritable => "java.lang.UnsupportedOperationException",
            _ => "picard.PicardException",
        }
    }
}

/// `VCFConstants.DEPTH_KEY`.
pub const DEPTH_KEY: &str = "DP";

/// Whether an attribute has to agree between the inputs.
pub fn attribute_must_agree(key: &str) -> bool {
    !EXEMPT_ATTRIBUTES.contains(&key)
}

/// Whether the merged header keeps a header line.
pub fn header_line_is_kept(key: &str) -> bool {
    !SAMPLE_SPECIFIC_HEADER_LINES.contains(&key)
}

/// One variant of one input, as far as the lockstep looks at it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub contig: String,
    pub start: i64,
    pub id: String,
    pub reference: String,
    pub alternates: Vec<String>,
    pub attributes: Vec<(String, String)>,
}

/// `checkThatAllelesMatch` and the two checks around it, in the reference's own order.
pub fn check_step(first: &Site, other: &Site) -> Result<(), Refusal> {
    if first.contig != other.contig || first.start != other.start {
        return Err(Refusal::Locus);
    }
    if first.id != other.id {
        return Err(Refusal::Id);
    }
    if first.reference != other.reference {
        return Err(Refusal::ReferenceAllele);
    }
    if first.alternates.len() != other.alternates.len() {
        return Err(Refusal::AlternateAlleleCount);
    }
    if first.alternates != other.alternates {
        return Err(Refusal::AlternateAllele {
            contig: first.contig.clone(),
            start: first.start,
        });
    }
    Ok(())
}

/// The attribute loop, which runs over the OTHER inputs' attributes against the first's.
///
/// The direction is the whole of it: a key only the FIRST file has is never looked for, so it is
/// kept without comment and reaches the output, and a key a LATER file has is refused by name. The
/// check reads as symmetric and is not.
///
/// A depth is refused twice over: by this loop when the two disagree, and by the write-back when
/// they agree, which is why a run that carries one never produces an output either way.
pub fn check_attributes(first: &Site, other: &Site) -> Result<(), Refusal> {
    for (key, value) in &other.attributes {
        if !attribute_must_agree(key) {
            continue;
        }
        match first.attributes.iter().find(|(name, _)| name == key) {
            None => return Err(Refusal::AttributeMissing(key.clone())),
            Some((_, extant)) if extant != value => {
                return Err(Refusal::AttributeDisagrees(key.clone()))
            }
            Some(_) => {}
        }
    }
    if other.attributes.iter().any(|(key, _)| key == DEPTH_KEY)
        || first.attributes.iter().any(|(key, _)| key == DEPTH_KEY)
    {
        return Err(Refusal::DepthIsUnwritable);
    }
    Ok(())
}

/// Whether the samples of the inputs, in order, are a valid multi-sample set.
pub fn sample_list(inputs: &[(String, Vec<String>)]) -> Result<Vec<String>, Refusal> {
    let mut samples: Vec<String> = Vec::new();
    for (file, names) in inputs {
        for name in names {
            if samples.contains(name) {
                return Err(Refusal::RepeatedSample {
                    file: file.clone(),
                    sample: name.clone(),
                });
            }
            samples.push(name.clone());
        }
    }
    Ok(samples)
}
