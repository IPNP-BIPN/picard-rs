//! Illumina's binary array files as Picard reads them: the bead pool manifest (`.bpm`), the
//! cluster file (`.egt`) and the genotype call file (`.gtc`).
//!
//! Ported from `picard.arrays.illumina.InfiniumDataFile`, `IlluminaBPMFile`,
//! `IlluminaBPMLocusEntry`, `InfiniumEGTFile`, `InfiniumGTCFile` and `InfiniumFileTOC` at tag
//! 3.4.0. Everything is little-endian, strings carry a seven-bits-a-byte length, and a GTC is a
//! table of contents the reader seeks through. A read past the end is the `EOFException` the
//! `DataInputStream` throws, which the tools wrap as an `IOException`.

/// What a parse can fail with: an `IOException` (the stream ran out, or a version the format
/// refuses) or a `PicardException` with its message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadError {
    Io(String),
    Picard(String),
}

/// A `DataInputStream` over a byte array.
pub struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Reader<'a> {
    pub fn new(bytes: &'a [u8]) -> Self {
        Reader { bytes, at: 0 }
    }

    pub fn seek(&mut self, at: usize) {
        self.at = at;
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], ReadError> {
        let end = self.at.checked_add(n).filter(|e| *e <= self.bytes.len());
        match end {
            Some(end) => {
                let s = &self.bytes[self.at..end];
                self.at = end;
                Ok(s)
            }
            None => Err(ReadError::Io("java.io.EOFException".to_string())),
        }
    }

    pub fn byte(&mut self) -> Result<i8, ReadError> {
        Ok(self.take(1)?[0] as i8)
    }

    pub fn int(&mut self) -> Result<i32, ReadError> {
        let b = self.take(4)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    pub fn short(&mut self) -> Result<i32, ReadError> {
        let b = self.take(2)?;
        Ok(i32::from(b[0]) | (i32::from(b[1]) << 8))
    }

    pub fn float(&mut self) -> Result<f32, ReadError> {
        let b = self.take(4)?;
        Ok(f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    pub fn skip(&mut self, n: usize) {
        self.at = (self.at + n).min(self.bytes.len());
    }

    /// `parseString`: the length in seven-bit groups, the bytes as Latin-1 characters.
    pub fn string(&mut self) -> Result<String, ReadError> {
        let mut total: i64 = 0;
        let mut partial = self.byte()?;
        let mut groups = 0;
        while (i32::from(partial) & 0x80) > 0 {
            total += i64::from(i32::from(partial) & 0x7F) << (7 * groups);
            partial = self.byte()?;
            groups += 1;
        }
        total += i64::from(partial) << (7 * groups);
        if total == 0 {
            return Ok(String::new());
        }
        let available = self.bytes.len() - self.at;
        let want = total as usize;
        if available < want {
            return Err(ReadError::Io(format!(
                "java.io.IOException: Did not fully read string. Read {available} out of {want}."
            )));
        }
        let b = self.take(want)?;
        Ok(b.iter().map(|&c| char::from(c)).collect())
    }
}

/// `IlluminaManifestRecord.IlluminaStrand.valueOf`, by name.
pub fn illumina_strand(name: &str) -> Result<String, ReadError> {
    match name {
        "" => Ok("NONE".to_string()),
        "TOP" | "BOT" | "PLUS" | "MINUS" | "NONE" => Ok(name.to_string()),
        other => Err(ReadError::Picard(format!(
            "java.lang.IllegalArgumentException: No enum constant picard.arrays.illumina.IlluminaManifestRecord.IlluminaStrand.{other}"
        ))),
    }
}

/// One `IlluminaBPMLocusEntry`.
#[derive(Debug, Clone, PartialEq)]
pub struct BpmLocus {
    pub version: i32,
    pub ilmn_id: String,
    pub name: String,
    pub index: usize,
    pub ilmn_strand: String,
    pub snp: String,
    pub chrom: String,
    pub ploidy: String,
    pub species: String,
    pub map_info: i32,
    pub top_genomic_seq: String,
    pub customer_strand: String,
    pub address_a: i32,
    pub address_b: i32,
    pub allele_a_probe_seq: String,
    pub allele_b_probe_seq: String,
    pub genome_build: String,
    pub source: String,
    pub source_version: String,
    pub source_strand: String,
    pub source_seq: String,
    pub exp_clusters: i32,
    pub intensity_only: bool,
    pub assay_type: i32,
    pub frac: [f32; 4],
    pub ref_strand: Option<char>,
    pub normalization_id: i32,
}

/// `IlluminaBPMFile`.
#[derive(Debug, Clone, PartialEq)]
pub struct Bpm {
    pub manifest_name: String,
    pub control_config: String,
    pub loci: Vec<BpmLocus>,
    pub all_normalization_ids: Vec<i32>,
    /// A `TreeSet`, so sorted.
    pub unique_normalization_ids: Vec<i32>,
}

fn locus_entry(r: &mut Reader) -> Result<BpmLocus, ReadError> {
    let version = r.int()?;
    if !(6..=8).contains(&version) {
        return Err(ReadError::Picard(format!(
            "Unsupported Locus version: {version}"
        )));
    }
    let ilmn_id = r.string()?;
    let name = r.string()?;
    r.string()?;
    r.string()?;
    r.string()?;
    let index = r.int()? - 1;
    r.string()?;
    let ilmn_strand = illumina_strand(&r.string()?)?;
    let snp = r.string()?;
    let chrom = r.string()?;
    let ploidy = r.string()?;
    let species = r.string()?;
    let map_text = r.string()?;
    let map_info: i32 = map_text.parse().map_err(|_| {
        ReadError::Picard(format!(
            "java.lang.NumberFormatException: For input string: \"{map_text}\""
        ))
    })?;
    let top_genomic_seq = r.string()?;
    let customer_strand = r.string()?;
    let address_a = r.int()?;
    let address_b = r.int()?;
    let allele_a_probe_seq = r.string()?;
    let allele_b_probe_seq = r.string()?;
    let genome_build = r.string()?;
    let source = r.string()?;
    let source_version = r.string()?;
    let source_strand = illumina_strand(&r.string()?)?;
    let source_seq = r.string()?;
    r.byte()?;
    let exp_clusters = i32::from(r.byte()?);
    let intensity = i32::from(r.byte()?);
    if intensity != 0 && intensity != 1 {
        return Err(ReadError::Picard(format!(
            "Unexpected value ('{intensity}') for intensity_only field"
        )));
    }
    let assay_type = i32::from(r.byte()?);
    let mut frac = [0.0f32; 4];
    if version >= 7 {
        for f in &mut frac {
            *f = r.float()?;
        }
    }
    let mut ref_strand = None;
    if version == 8 {
        ref_strand = r.string()?.chars().next();
    }
    if !(0..=2).contains(&assay_type) {
        return Err(ReadError::Picard(format!(
            "Invalid assay_type '{assay_type}' in BPM file"
        )));
    }
    if (assay_type != 0 && address_b == 0) || (assay_type == 0 && address_b != 0) {
        return Err(ReadError::Picard(format!(
            "Invalid assay_type '{assay_type}' for address B '{address_b}' in BPM file"
        )));
    }
    Ok(BpmLocus {
        version,
        ilmn_id,
        name,
        index: index.max(0) as usize,
        ilmn_strand,
        snp,
        chrom,
        ploidy,
        species,
        map_info,
        top_genomic_seq,
        customer_strand,
        address_a,
        address_b,
        allele_a_probe_seq,
        allele_b_probe_seq,
        genome_build,
        source,
        source_version,
        source_strand,
        source_seq,
        exp_clusters,
        intensity_only: intensity == 1,
        assay_type,
        frac,
        ref_strand,
        normalization_id: 0,
    })
}

impl Bpm {
    /// `new IlluminaBPMFile(file)`.
    pub fn parse(bytes: &[u8]) -> Result<Bpm, ReadError> {
        let mut r = Reader::new(bytes);
        let mut id = String::new();
        for _ in 0..3 {
            id.push(char::from(r.byte()? as u8));
        }
        if id != "BPM" {
            return Err(ReadError::Picard(format!(
                "Invalid identifier '{id}' for BPM file"
            )));
        }
        let file_version = r.byte()?;
        if file_version != 1 {
            return Err(ReadError::Picard(format!(
                "Unknown BPM version ({file_version})"
            )));
        }
        let mut version = r.int()?;
        if version & 0x1000 == 0x1000 {
            version ^= 0x1000;
        }
        if !(3..=5).contains(&version) {
            return Err(ReadError::Picard(format!(
                "Unsupported BPM version ({version})"
            )));
        }
        let manifest_name = r.string()?;
        let control_config = r.string()?;
        let n = r.int()?.max(0) as usize;
        r.skip(4 * n);
        let mut names = Vec::with_capacity(n);
        for _ in 0..n {
            names.push(r.string()?);
        }
        let mut all: Vec<i32> = (0..n)
            .map(|_| r.byte().map(i32::from))
            .collect::<Result<_, _>>()?;
        let mut unique: Vec<i32> = Vec::new();
        let mut loci: Vec<Option<BpmLocus>> = vec![None; n];
        for _ in 0..n {
            let mut locus = locus_entry(&mut r)?;
            let norm = *all.get(locus.index).ok_or_else(|| {
                ReadError::Picard(format!(
                    "java.lang.ArrayIndexOutOfBoundsException: Index {} out of bounds for length {n}",
                    locus.index
                ))
            })?;
            if norm > 100 {
                return Err(ReadError::Picard(format!(
                    "Invalid normalization ID: {norm} for name: {}",
                    locus.name
                )));
            }
            locus.normalization_id = norm + 100 * locus.assay_type;
            all[locus.index] = locus.normalization_id;
            if !unique.contains(&locus.normalization_id) {
                unique.push(locus.normalization_id);
            }
            if loci[locus.index].is_some() {
                return Err(ReadError::Picard(format!(
                    "Duplicate locus entry for index: {} '{}",
                    locus.index, locus.name
                )));
            }
            if names[locus.index] != locus.name {
                return Err(ReadError::Picard(format!(
                    "Mismatch in names at index: {}",
                    locus.index
                )));
            }
            let at = locus.index;
            loci[at] = Some(locus);
        }
        unique.sort_unstable();
        Ok(Bpm {
            manifest_name,
            control_config,
            loci: loci.into_iter().flatten().collect(),
            all_normalization_ids: all,
            unique_normalization_ids: unique,
        })
    }
}

/// `InfiniumEGTFile`.
#[derive(Debug, Clone, PartialEq)]
pub struct Egt {
    pub manifest_name: String,
    pub n: Vec<[i32; 3]>,
    pub dev_r: Vec<[f32; 3]>,
    pub mean_r: Vec<[f32; 3]>,
    pub dev_theta: Vec<[f32; 3]>,
    pub mean_theta: Vec<[f32; 3]>,
    pub total_score: Vec<f32>,
    pub rs_names: Vec<String>,
}

impl Egt {
    /// `new InfiniumEGTFile(file)`.
    pub fn parse(bytes: &[u8]) -> Result<Egt, ReadError> {
        let mut r = Reader::new(bytes);
        let file_version = r.int()?;
        for _ in 0..5 {
            r.string()?;
        }
        r.byte()?;
        if file_version == 2 {
            return Err(ReadError::Io(
                "java.io.IOException: Version '2' unsupported".to_string(),
            ));
        }
        r.string()?;
        let version = r.int()?;
        if version > 9 {
            return Err(ReadError::Io(format!(
                "java.io.IOException: Error. Cannot read file - unknown version {version} for gentrain data type"
            )));
        }
        let manifest_name = r.string()?;
        let codes = r.int()?.max(0) as usize;
        let mut egt = Egt {
            manifest_name,
            n: Vec::with_capacity(codes),
            dev_r: Vec::with_capacity(codes),
            mean_r: Vec::with_capacity(codes),
            dev_theta: Vec::with_capacity(codes),
            mean_theta: Vec::with_capacity(codes),
            total_score: Vec::with_capacity(codes),
            rs_names: Vec::with_capacity(codes),
        };
        let triple_f = |r: &mut Reader| -> Result<[f32; 3], ReadError> {
            Ok([r.float()?, r.float()?, r.float()?])
        };
        for _ in 0..codes {
            egt.n.push([r.int()?, r.int()?, r.int()?]);
            egt.dev_r.push(triple_f(&mut r)?);
            egt.mean_r.push(triple_f(&mut r)?);
            egt.dev_theta.push(triple_f(&mut r)?);
            egt.mean_theta.push(triple_f(&mut r)?);
            r.skip(15 * 4);
        }
        for _ in 0..codes {
            r.skip(4);
            egt.total_score.push(r.float()?);
            r.skip(4);
            r.skip(1);
        }
        for _ in 0..codes {
            // `skipString`: one length byte, however long the string.
            let len = r.byte()?;
            r.skip(len.max(0) as usize);
        }
        for _ in 0..codes {
            let name = r.string()?;
            if egt.rs_names.contains(&name) {
                return Err(ReadError::Picard(format!(
                    "Non-unique rsName '{name}' found in cluster file"
                )));
            }
            egt.rs_names.push(name);
        }
        Ok(egt)
    }
}
