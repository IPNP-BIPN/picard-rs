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

/// `InfiniumGTCFile`: every field a getter returns, `None` where the table of contents had no
/// entry for it (the reference's null, or 0 for a primitive).
#[derive(Debug, Clone, Default)]
pub struct Gtc {
    pub identifier: String,
    pub file_version: i32,
    pub number_of_snps: i32,
    pub ploidy: i32,
    pub ploidy_type: i32,
    pub sample_name: Option<String>,
    pub sample_plate: Option<String>,
    pub sample_well: Option<String>,
    pub cluster_file: Option<String>,
    pub snp_manifest: Option<String>,
    pub imaging_date: Option<String>,
    pub auto_call_date: Option<String>,
    pub auto_call_version: Option<String>,
    pub transformations: Vec<crate::gtc_to_vcf::Transformation>,
    pub raw_control_x: Option<Vec<i32>>,
    pub raw_control_y: Option<Vec<i32>>,
    pub raw_x: Option<Vec<i32>>,
    pub raw_y: Option<Vec<i32>>,
    pub genotypes: Option<Vec<i8>>,
    pub base_calls: Option<Vec<[i8; 2]>>,
    pub genotype_scores: Option<Vec<f32>>,
    pub scanner_name: Option<String>,
    pub pmt_green: i32,
    pub pmt_red: i32,
    pub scanner_version: Option<String>,
    pub imaging_user: Option<String>,
    pub call_rate: f64,
    pub gender: Option<String>,
    pub log_r_dev: f32,
    pub p10_gc: f32,
    pub dx: i32,
    pub p50_gc: f32,
    pub num_calls: i32,
    pub num_no_calls: i32,
    pub num_intensity_only: i32,
    pub red_percentiles: Option<[i32; 3]>,
    pub green_percentiles: Option<[i32; 3]>,
    pub sentrix_barcode: Option<String>,
    pub b_allele_freqs: Option<Vec<f32>>,
    pub log_r_ratios: Option<Vec<f32>>,
    pub normalized_x: Vec<f32>,
    pub normalized_y: Vec<f32>,
    pub r_ilmn: Vec<f32>,
    pub theta_ilmn: Vec<f32>,
    pub aa_calls: i64,
    pub ab_calls: i32,
    pub bb_calls: i64,
}

impl Gtc {
    /// `new InfiniumGTCFile(gtc, bpm)`: the manifest's normalization ids are read first, then the
    /// table of contents, each entry from the start of the file.
    pub fn parse(bytes: &[u8], bpm: &Bpm) -> Result<Gtc, ReadError> {
        let mut g = Gtc::default();
        let mut r = Reader::new(bytes);
        for _ in 0..3 {
            g.identifier.push(char::from(r.byte()? as u8));
        }
        if g.identifier != "gtc" {
            return Err(ReadError::Picard(format!(
                "Invalid identifier '{}' for GTC file",
                g.identifier
            )));
        }
        g.file_version = i32::from(r.byte()?);
        let entries = r.int()?.max(0) as usize;
        let mut toc = Vec::with_capacity(entries);
        for _ in 0..entries {
            let id = r.short()? as i16;
            let offset = r.int()?;
            toc.push((id, offset));
        }
        for (id, offset) in toc {
            let mut s = Reader::new(bytes);
            s.seek(offset.max(0) as usize);
            let ushorts = |s: &mut Reader| -> Result<Vec<i32>, ReadError> {
                let n = s.int()?.max(0) as usize;
                (0..n).map(|_| s.short()).collect()
            };
            let floats = |s: &mut Reader| -> Result<Vec<f32>, ReadError> {
                let n = s.int()?.max(0) as usize;
                (0..n).map(|_| s.float()).collect()
            };
            match id {
                1 => g.number_of_snps = offset,
                2 => g.ploidy = offset,
                3 => g.ploidy_type = offset,
                10 => g.sample_name = Some(s.string()?),
                11 => g.sample_plate = Some(s.string()?),
                12 => g.sample_well = Some(s.string()?),
                100 => g.cluster_file = Some(s.string()?),
                101 => g.snp_manifest = Some(s.string()?),
                200 => g.imaging_date = Some(s.string()?),
                201 => g.auto_call_date = Some(s.string()?),
                300 => g.auto_call_version = Some(s.string()?),
                400 => {
                    let n = s.int()?.max(0) as usize;
                    g.transformations.clear();
                    for _ in 0..n {
                        s.int()?;
                        let t = crate::gtc_to_vcf::Transformation {
                            offset_x: s.float()?,
                            offset_y: s.float()?,
                            scale_x: s.float()?,
                            scale_y: s.float()?,
                            shear: s.float()?,
                            theta: s.float()?,
                        };
                        for _ in 0..6 {
                            s.float()?;
                        }
                        g.transformations.push(t);
                    }
                }
                500 => g.raw_control_x = Some(ushorts(&mut s)?),
                501 => g.raw_control_y = Some(ushorts(&mut s)?),
                1000 => g.raw_x = Some(ushorts(&mut s)?),
                1001 => g.raw_y = Some(ushorts(&mut s)?),
                1002 => {
                    let n = s.int()?.max(0) as usize;
                    let calls: Vec<i8> = (0..n).map(|_| s.byte()).collect::<Result<_, _>>()?;
                    for c in &calls {
                        match c {
                            1 => g.aa_calls += 1,
                            2 => g.ab_calls += 1,
                            3 => g.bb_calls += 1,
                            _ => {}
                        }
                    }
                    g.genotypes = Some(calls);
                }
                1003 => {
                    let n = s.int()?.max(0) as usize;
                    let mut calls = Vec::with_capacity(n);
                    for _ in 0..n {
                        let mut pair = [s.byte()?, s.byte()?];
                        for b in &mut pair {
                            if *b == 0 {
                                *b = b'-' as i8;
                            }
                        }
                        calls.push(pair);
                    }
                    g.base_calls = Some(calls);
                }
                1004 => g.genotype_scores = Some(floats(&mut s)?),
                1005 => {
                    g.scanner_name = Some(s.string()?);
                    g.pmt_green = s.int()?;
                    g.pmt_red = s.int()?;
                    g.scanner_version = Some(s.string()?);
                    g.imaging_user = Some(s.string()?);
                }
                1006 => g.call_rate = f64::from(s.float()?),
                1007 => {
                    g.gender = Some(match bytes.get(offset.max(0) as usize) {
                        Some(b) => char::from(*b).to_string(),
                        None => char::from_u32(0xFFFF).map(String::from).unwrap_or_default(),
                    })
                }
                1008 => g.log_r_dev = s.float()?,
                1009 => g.p10_gc = s.float()?,
                1010 => g.dx = s.int()?,
                1011 => {
                    g.p50_gc = s.float()?;
                    g.num_calls = s.int()?;
                    g.num_no_calls = s.int()?;
                    g.num_intensity_only = s.int()?;
                }
                1012 => g.b_allele_freqs = Some(floats(&mut s)?),
                1013 => g.log_r_ratios = Some(floats(&mut s)?),
                1014 => g.red_percentiles = Some([s.short()?, s.short()?, s.short()?]),
                1015 => g.green_percentiles = Some([s.short()?, s.short()?, s.short()?]),
                1016 => g.sentrix_barcode = Some(s.string()?),
                _ => {}
            }
        }
        if g.num_calls == 0 {
            g.num_calls = (g.aa_calls + i64::from(g.ab_calls) + g.bb_calls) as i32;
        }
        g.normalize(bpm)?;
        Ok(g)
    }

    /// `normalizeIntensities` and `calculateRandTheta`.
    fn normalize(&mut self, bpm: &Bpm) -> Result<(), ReadError> {
        let n = self.number_of_snps.max(0) as usize;
        let raw_x = self.raw_x.clone().unwrap_or_default();
        let raw_y = self.raw_y.clone().unwrap_or_default();
        self.normalized_x = vec![0.0; n];
        self.normalized_y = vec![0.0; n];
        for (i, &x) in raw_x.iter().enumerate() {
            let oob = |len: usize| {
                ReadError::Picard(format!(
                    "java.lang.ArrayIndexOutOfBoundsException: Index {i} out of bounds for length {len}"
                ))
            };
            let y = *raw_y.get(i).ok_or_else(|| oob(raw_y.len()))?;
            let index = bpm
                .all_normalization_ids
                .get(i)
                .and_then(|id| bpm.unique_normalization_ids.iter().position(|u| u == id));
            let (nx, ny) = match index {
                Some(at) => {
                    let t = self
                        .transformations
                        .get(at)
                        .ok_or_else(|| oob(self.transformations.len()))?;
                    crate::gtc_to_vcf::normalize(x, y, t)
                }
                None => (x as f32, y as f32),
            };
            if i >= n {
                return Err(oob(n));
            }
            self.normalized_x[i] = nx;
            self.normalized_y[i] = ny;
        }
        self.r_ilmn = Vec::with_capacity(n);
        self.theta_ilmn = Vec::with_capacity(n);
        for i in 0..n {
            let (r, theta) =
                crate::gtc_to_vcf::r_and_theta(self.normalized_x[i], self.normalized_y[i]);
            self.r_ilmn.push(r);
            self.theta_ilmn.push(theta);
        }
        Ok(())
    }
}
