/*
 * Builds the fixture corpus the covering arrays run against.
 *
 * Usage: MakeFixtures <output directory>
 *
 * A covering array cannot invent a file path: the value has to be a file that exists and holds
 * content the tool accepts (gatk-rs tools/coverage/domains.py excludes every path-typed argument
 * for exactly that reason, which is most of what it excludes). This produces that corpus, small
 * and deterministic, at fixed paths, so a row of the array can be turned into a command line.
 *
 * Three properties matter and are the reason this is a program rather than a directory of
 * committed files:
 *
 *   1. Deterministic. Fixed seed, fixed content, no timestamps, no temp directories. Two runs
 *      produce the same bytes, so a divergence between the oracle and the port is about the tool
 *      and not about the input.
 *   2. Small. Every row of a covering array runs the tool once; HaplotypeCaller's t=2 array is 62
 *      rows and the whole inventory is 19,437. The corpus is sized for that, not for realism.
 *   3. Branchy. Uniform perfect reads exercise one path. These reads carry unmapped mates,
 *      duplicates, secondary and supplementary alignments, soft clips, an indel, no-calls, both
 *      strands, two read groups and two libraries, so a row that flips a filtering argument
 *      actually changes the answer. A corpus where every argument produces the same output would
 *      make a covering array look green while testing nothing.
 *
 * The reference is two short contigs, which is enough for a sequence dictionary, an interval
 * list, and a tool that needs REFERENCE_SEQUENCE.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import java.io.File;
import java.io.PrintWriter;
import java.util.Random;

public class MakeFixtures {

    static final int CHR1 = 2_000;
    static final int CHR2 = 1_000;
    static final int READS = 400;
    static final int READ_LENGTH = 50;

    public static void main(String[] args) throws Exception {
        File dir = new File(args.length > 0 ? args[0] : "fixtures");
        dir.mkdirs();

        String chr1 = reference(CHR1, 20260729L);
        String chr2 = reference(CHR2, 20260730L);

        writeFasta(new File(dir, "ref.fasta"), chr1, chr2);
        writeFai(new File(dir, "ref.fasta.fai"), chr1, chr2);

        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate);
        writeBam(new File(dir, "small.bam"), header, reads(header, true), true);
        writeSam(new File(dir, "small.sam"), header, reads(header, true));

        SAMFileHeader queryname = header(SAMFileHeader.SortOrder.queryname);
        writeBam(new File(dir, "queryname.bam"), queryname, reads(queryname, false), false);

        SAMFileHeader oneGroup = header(SAMFileHeader.SortOrder.coordinate);
        oneGroup.setReadGroups(java.util.Collections.singletonList(oneGroup.getReadGroup("rg1")));
        java.util.List<SAMRecord> oneGroupReads = reads(oneGroup, true);
        for (SAMRecord r : oneGroupReads) r.setAttribute("RG", "rg1");
        writeBam(new File(dir, "one_read_group.bam"), oneGroup, oneGroupReads, false);

        // Illumina-style read names, for the tools that read a physical location out of one.
        // `PositionBasedDownsampleSam` keeps reads by where they sit on the flowcell, and the
        // names everywhere else in this corpus (`read0314`) carry no location at all: the parser
        // refuses them, every read lands on the same defaulted tile, and the mask becomes
        // all-or-nothing. These names give it tiles and coordinates to work on, in a file of its
        // own so that no existing fixture's bytes move.
        SAMFileHeader tiled = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> tiledReads = reads(tiled, true);
        for (SAMRecord r : tiledReads) {
            // The new name is derived from the old one, so the two ends of a pair still share it:
            // `read0314` gives 314, and both mates carry that number.
            int n = Integer.parseInt(r.getReadName().replaceAll("[^0-9]", ""));
            int tile = 1101 + (n % 4);
            int x = 1000 + ((n * 137) % 20000);
            int y = 1000 + ((n * 251) % 20000);
            r.setReadName(String.format("INST:1:FLOWCELL:1:%d:%d:%d", tile, x, y));
        }
        writeBam(new File(dir, "tiled.bam"), tiled, tiledReads, false);

        // Reads carrying `MC`, for the duplicate markers that read the mate's cigar instead of
        // waiting for the mate. Without the tag `SimpleMarkDuplicatesWithMateCigar` refuses every
        // file outright and `MarkDuplicatesWithMateCigar` skips every pair, so a corpus without it
        // measures the refusal and nothing else. `SamPairUtil.setMateInformation` writes the tag
        // the same way `FixMateInformation --ADD_MATE_CIGAR` does; the file is new, so no existing
        // fixture's bytes move.
        SAMFileHeader mateCigarHeader = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> mateCigarReads = reads(mateCigarHeader, true);
        java.util.Map<String, java.util.List<SAMRecord>> byName = new java.util.LinkedHashMap<>();
        for (SAMRecord r : mateCigarReads) {
            byName.computeIfAbsent(r.getReadName(), k -> new java.util.ArrayList<>()).add(r);
        }
        for (java.util.List<SAMRecord> template : byName.values()) {
            if (template.size() == 2) {
                htsjdk.samtools.SamPairUtil.setMateInfo(template.get(0), template.get(1), true);
            }
        }
        writeBam(new File(dir, "mate_cigar.bam"), mateCigarHeader, mateCigarReads, true);

        // Records of one query name that DISAGREE about their duplicate flag, for
        // `CheckDuplicateMarking`. Everywhere else in this corpus both ends of a pair carry the
        // same flag, so the tool finds nothing whatever `MODE` is asked for: the array covers the
        // argument and observes none of it.
        //
        // Which record of the pair is flipped decides which modes see the disagreement, and that
        // is what makes the four values four answers. A flipped SECONDARY or SUPPLEMENTARY record
        // is seen by `ALL` alone; a flipped UNMAPPED one is also seen by `PRIMARY_ONLY`; a flipped
        // record of a pair that is not proper is seen by those and by `PRIMARY_MAPPED_ONLY`; and a
        // flipped ordinary mate is seen by all four.
        //
        // Coordinate-sorted, so the tool has to sort it into query-name order itself: the order it
        // sorts into decides which record of a name is the one the others are compared against.
        SAMFileHeader inconsistentHeader = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> inconsistentReads = reads(inconsistentHeader, true);
        java.util.Map<String, java.util.List<SAMRecord>> inconsistentTemplates = new java.util.LinkedHashMap<>();
        for (SAMRecord r : inconsistentReads) {
            inconsistentTemplates.computeIfAbsent(r.getReadName(), k -> new java.util.ArrayList<>()).add(r);
        }
        int flipped = 0;
        for (java.util.List<SAMRecord> template : inconsistentTemplates.values()) {
            if (template.size() != 2) continue;
            SAMRecord one = template.get(0).getFirstOfPairFlag() ? template.get(0) : template.get(1);
            SAMRecord two = template.get(0).getFirstOfPairFlag() ? template.get(1) : template.get(0);
            int n = Integer.parseInt(one.getReadName().replaceAll("[^0-9]", ""));
            SAMRecord flip;
            if (n % 16 == 0) {
                flip = one;          // secondary
            } else if (n % 18 == 0) {
                flip = two;          // supplementary
            } else if (n % 20 == 0) {
                flip = two;          // unmapped
            } else if (n % 5 == 0) {
                flip = two;          // not a proper pair
            } else if (n % 6 == 0) {
                flip = two;          // an ordinary mate
            } else {
                continue;
            }
            flip.setDuplicateReadFlag(!flip.getDuplicateReadFlag());
            flipped++;
        }
        if (flipped == 0) throw new IllegalStateException("no duplicate flag was flipped");
        writeBam(new File(dir, "inconsistent_duplicates.bam"), inconsistentHeader, inconsistentReads, false);

        // Reads whose ends really are Illumina adapters, for `MarkIlluminaAdapters`. On the random
        // bases of the other fixtures no adapter is ever found, so every accepted row produces the
        // same output and the array covers the search without running it.
        //
        // Three families are planted, and which family a pair gets is what makes `--ADAPTERS`
        // observable. Picard's parser APPENDS to a list argument's default, so every row searches
        // INDEXED, DUAL_INDEXED and PAIRED_END whatever it asks for: a corpus carrying only a
        // PAIRED_END adapter answers the same thing for all nine values. NEXTERA_V2 and
        // TRUSEQ_SMALLRNA are not in that default, so the pairs carrying them are marked only by
        // the rows that name them.
        //
        // Two planting shapes, because the paired rule has two branches. A one-sided plant puts
        // the three prime adapter in read one and leaves read two alone, which is the branch that
        // re-checks the single match against twice the minimum and then marks BOTH reads. A
        // two-sided plant puts the three prime adapter in read one and the five prime adapter in
        // read-order in read two, at the same offset, which is the branch where the two indices
        // agree and the pair is marked immediately.
        SAMFileHeader adapterHeader = header(SAMFileHeader.SortOrder.queryname);
        java.util.List<SAMRecord> adapterReads = reads(adapterHeader, false);
        java.util.Map<String, java.util.List<SAMRecord>> adapterTemplates = new java.util.LinkedHashMap<>();
        for (SAMRecord r : adapterReads) {
            adapterTemplates.computeIfAbsent(r.getReadName(), k -> new java.util.ArrayList<>()).add(r);
        }
        int planted = 0;
        for (java.util.List<SAMRecord> template : adapterTemplates.values()) {
            if (template.size() != 2) continue;
            SAMRecord one = template.get(0).getFirstOfPairFlag() ? template.get(0) : template.get(1);
            SAMRecord two = template.get(0).getFirstOfPairFlag() ? template.get(1) : template.get(0);
            int n = Integer.parseInt(one.getReadName().replaceAll("[^0-9]", ""));
            String fivePrime, threePrime;
            boolean twoSided;
            if (n % 12 == 0) {
                fivePrime = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
                threePrime = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG";
                twoSided = false;
            } else if (n % 12 == 4) {
                fivePrime = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
                threePrime = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
                twoSided = true;
            } else if (n % 12 == 8) {
                fivePrime = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC";
                threePrime = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
                twoSided = false;
            } else {
                continue;
            }
            // 24 bases, which is well over the paired minimum and over the single-end one too.
            plantAdapter(one, threePrime, 24);
            if (twoSided) {
                plantAdapter(two, SequenceUtil.reverseComplement(fivePrime), 24);
            }
            planted++;
        }
        if (planted == 0) throw new IllegalStateException("no adapter was planted");
        writeBam(new File(dir, "adapters.bam"), adapterHeader, adapterReads, false);

        // Read pairs that really do repeat, for `EstimateLibraryComplexity`. The tool groups pairs
        // by the first MIN_IDENTICAL_BASES of both ends and counts how big each group of duplicates
        // is; on the random bases of the other fixtures every group holds one pair, every bin is
        // one, and MIN_GROUP_COUNT then drops the lot, so the metrics are zeros whatever the
        // arguments say.
        //
        // Twenty-four families of one to four identical pairs each, which gives bins of one, two,
        // three and four. Every fourth family mutates its later copies -- two bases, inside the
        // default MAX_DIFF_RATE's allowance, or eight, outside it -- so the rate decides whether
        // those copies join the group. Every sixth family is written at quality fifteen, below the
        // default MIN_MEAN_QUALITY, so the quality filter has something to drop.
        //
        // The names are Illumina-style because the optical duplicate finder reads a tile and a
        // position out of them: the first two copies of a family share a tile twenty pixels apart,
        // which is inside the default pixel distance and outside a smaller one, and the later
        // copies sit on tiles of their own.
        //
        // The two read groups alternate by family, never within one, so a family stays in one
        // library and the library split does not cut a group in half.
        SAMFileHeader complexityHeader = header(SAMFileHeader.SortOrder.queryname);
        java.util.List<SAMRecord> complexityReads = new java.util.ArrayList<>();
        Random complexityRng = new Random(20260908L);
        String complexityBases = "ACGT";
        for (int family = 0; family < 24; family++) {
            byte[] readOne = new byte[READ_LENGTH];
            byte[] readTwo = new byte[READ_LENGTH];
            for (int b = 0; b < READ_LENGTH; b++) {
                readOne[b] = (byte) complexityBases.charAt(complexityRng.nextInt(4));
                readTwo[b] = (byte) complexityBases.charAt(complexityRng.nextInt(4));
            }
            int copies = 1 + (family % 4);
            String group = (family % 2 == 0) ? "rg1" : "rg2";
            byte quality = (byte) (family % 6 == 5 ? 15 : 35);
            for (int copy = 0; copy < copies; copy++) {
                byte[] one = readOne.clone();
                if (copy > 0 && family % 4 == 3) {
                    int errors = (family % 8 == 3) ? 2 : 8;
                    for (int e = 0; e < errors; e++) {
                        one[20 + e] = mutateBase(one[20 + e]);
                    }
                }
                int tile = 1101 + (copy < 2 ? 0 : copy);
                int x = 1000 + family * 7 + (copy < 2 ? copy * 20 : copy * 5000);
                int y = 2000 + family * 11 + (copy < 2 ? copy * 20 : copy * 5000);
                String name = String.format("INST:1:FLOWCELL:1:%d:%d:%d", tile, x, y);

                SAMRecord first = new SAMRecord(complexityHeader);
                SAMRecord second = new SAMRecord(complexityHeader);
                for (SAMRecord r : new SAMRecord[] {first, second}) {
                    byte[] quals = new byte[READ_LENGTH];
                    java.util.Arrays.fill(quals, quality);
                    r.setReadName(name);
                    r.setBaseQualities(quals);
                    r.setReferenceIndex(0);
                    r.setCigarString(READ_LENGTH + "M");
                    r.setMappingQuality(60);
                    r.setAttribute("RG", group);
                    r.setReadPairedFlag(true);
                    r.setProperPairFlag(true);
                }
                first.setFirstOfPairFlag(true);
                second.setSecondOfPairFlag(true);
                first.setAlignmentStart(100 + family);
                second.setAlignmentStart(400 + family);
                first.setReadBases(one);
                // The second end is on the negative strand, so the file stores it reverse
                // complemented and the tool has to complement it back before it compares: a
                // corpus whose ends are all forward never runs that path.
                byte[] two = readTwo.clone();
                SequenceUtil.reverseComplement(two);
                second.setReadBases(two);
                second.setReadNegativeStrandFlag(true);
                first.setMateNegativeStrandFlag(true);
                SamPairUtil.setMateInfo(first, second, false);
                complexityReads.add(first);
                complexityReads.add(second);
            }
        }
        complexityReads.sort(new SAMRecordQueryNameComparator());
        writeBam(new File(dir, "complexity.bam"), complexityHeader, complexityReads, false);

        SAMFileHeader unmappedHeader = header(SAMFileHeader.SortOrder.unsorted);
        writeBam(new File(dir, "unmapped.bam"), unmappedHeader, unmapped(unmappedHeader), false);

        // A read-name list, for FilterSamReads' includeReadList / excludeReadList. Every fourth
        // pair, so both filters keep something and drop something.
        try (PrintWriter out = new PrintWriter(new File(dir, "read_names.txt"), "UTF-8")) {
            for (int i = 0; i < READS; i += 8) out.printf("read%04d%n", i);
        }

        // Three small VCFs, for the tools that read variants rather than reads. `variants.vcf` has
        // two samples and a mix of SNPs and an indel, half of its sites also in `dbsnp.vcf`, so a
        // tool that partitions by novelty sees both partitions; `single_sample.vcf` is the
        // one-sample file the tools that refuse more than one need. Written through htsjdk's own
        // writer, index and all, because a hand-written VCF is a fixture whose bugs become the
        // tool's answers.
        writeVcf(new File(dir, "variants.vcf"), chr1, chr2, true);
        writeVcf(new File(dir, "dbsnp.vcf"), chr1, chr2, false);
        writeVcf(new File(dir, "single_sample.vcf"), chr1, chr2, true, 1);

        // Two `CollectQualityYieldMetrics` outputs, for the tools that accumulate metrics files
        // rather than reads. The header comments are what that tool writes, command line and
        // timestamp included, because the accumulator reads past them to the table and a fixture
        // that dropped them would not be the file it is given in practice. The two differ in every
        // counter, so a row that reads one is a different answer from a row that reads the other.
        writeQualityYield(new File(dir, "quality_yield_one.metrics"),
                400, 400, 50, 20000, 20000, 10547, 10547, 5200, 5200, 20456, 20456);
        writeQualityYield(new File(dir, "quality_yield_two.metrics"),
                150, 120, 60, 9000, 7200, 4100, 3300, 2000, 1600, 8800, 7000);

        // Reads carrying cell and molecular barcodes, for `SamToFastqWithTags`. That tool writes
        // the ordinary read FASTQ and, beside it, one FASTQ per SEQUENCE_TAG_GROUP whose reads are
        // built from TAG VALUES rather than from bases, so a corpus needs records that carry the
        // tags and records that do not: a group naming a tag a read is missing is a refusal, and
        // that is a row of the array.
        //
        // Every read carries every tag the array's groups name: a read missing one is a refusal
        // ("does have a value for tag"), and a corpus that refused most rows would measure the
        // check rather than the writing.
        SAMFileHeader taggedHeader = header(SAMFileHeader.SortOrder.queryname);
        java.util.List<SAMRecord> taggedReads = reads(taggedHeader, false);
        for (SAMRecord r : taggedReads) {
            int n = Integer.parseInt(r.getReadName().replaceAll("[^0-9]", ""));
            r.setAttribute("CB", "ACGTAC" + (char) ('A' + (n % 4)));
            r.setAttribute("CY", "IIIIIII");
            r.setAttribute("UR", "TTGCA");
            r.setAttribute("UY", "IIIII");
        }
        writeBam(new File(dir, "tagged.bam"), taggedHeader, taggedReads, false);

        writeIntervals(new File(dir, "targets.interval_list"));
        writeBed(new File(dir, "targets.bed"));
        writeMixedBed(new File(dir, "targets_mixed.bed"));
        writeMixedIntervals(new File(dir, "targets_mixed.interval_list"));
        writeDescribedFasta(new File(dir, "described.fasta"), chr2);
        writeDict(new File(dir, "ref.dict"), chr1, chr2);
        writeFastq(new File(dir, "reads_1.fastq"), 1);
        writeFastq(new File(dir, "reads_2.fastq"), 2);

        System.out.println("fixtures written to " + dir.getAbsolutePath());
        for (File f : dir.listFiles()) {
            System.out.printf("%s\t%d%n", f.getName(), f.length());
        }
    }

    /** A reference with a fixed but non-uniform base composition, so GC-dependent tools vary. */
    static String reference(int length, long seed) {
        Random rng = new Random(seed);
        char[] bases = new char[length];
        String alphabet = "ACGT";
        for (int i = 0; i < length; i++) {
            // A GC-rich stretch in the middle, and a run of Ns, so the reference is not featureless.
            if (i > length / 2 && i < length / 2 + 100) bases[i] = (i % 2 == 0) ? 'G' : 'C';
            else if (i > length - 60 && i < length - 40) bases[i] = 'N';
            else bases[i] = alphabet.charAt(rng.nextInt(4));
        }
        return new String(bases);
    }

    static void writeFasta(File f, String chr1, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            writeContig(p, "chr1", chr1);
            writeContig(p, "chr2", chr2);
        }
    }

    static void writeContig(PrintWriter p, String name, String bases) {
        p.println(">" + name);
        for (int i = 0; i < bases.length(); i += 60) {
            p.println(bases.substring(i, Math.min(i + 60, bases.length())));
        }
    }

    /**
     * The FASTA index, whose offsets have to be the file's real ones.
     *
     * The byte length of a contig is its bases plus one newline per line, and the last line is
     * short: `ceil(len / 60) * 61` over-counts it by `60 - (len % 60)`. chr1 is 2000 bases, so the
     * old arithmetic put chr2's offset 40 bytes past where chr2 begins.
     *
     * Nothing noticed until a tool took the indexed path. `ReferenceSequenceFileFactory` opens
     * `IndexedFastaSequenceFile` only when the caller asks for names truncated at whitespace
     * ("Using faidx requires truncateNamesAtWhitespace"), so NormalizeFasta read this index only
     * with TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE=true, and then sliced chr2 from the wrong byte:
     * its output carried the file's own line terminators as if they were bases, one of them
     * producing an empty line. The covering array is what ran that combination.
     */
    static void writeFai(File f, String chr1, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            int lineWidth = 60, lineBytes = 61;
            long offset1 = ">chr1\n".length();
            long chr1Bytes = bytesOnDisk(chr1.length(), lineWidth);
            long offset2 = offset1 + chr1Bytes + ">chr2\n".length();
            p.printf("chr1\t%d\t%d\t%d\t%d%n", chr1.length(), offset1, lineWidth, lineBytes);
            p.printf("chr2\t%d\t%d\t%d\t%d%n", chr2.length(), offset2, lineWidth, lineBytes);
        }
    }

    /**
     * The same corpus under a single read group.
     *
     * Every other fixture carries the same two read groups, so a tool whose whole output is a
     * function of the header's `@RG` records answers identically on all of them:
     * CalculateReadGroupChecksum's array was nine rows and one digest, which covers its arguments
     * without testing them. This file differs in exactly the thing that tool reads.
     */
    static SAMFileHeader header(SAMFileHeader.SortOrder order) {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", CHR1));
        d.addSequence(new SAMSequenceRecord("chr2", CHR2));
        h.setSequenceDictionary(d);
        h.setSortOrder(order);
        // Two read groups in two libraries: the multi-level collectors have a LIBRARY and a
        // READ_GROUP accumulation level, and one read group would leave both untested.
        for (String[] rg : new String[][] {{"rg1", "lib1", "sample1"}, {"rg2", "lib2", "sample2"}}) {
            SAMReadGroupRecord r = new SAMReadGroupRecord(rg[0]);
            r.setLibrary(rg[1]);
            r.setSample(rg[2]);
            r.setPlatform("ILLUMINA");
            r.setPlatformUnit("unit-" + rg[0]);
            h.addReadGroup(r);
        }
        return h;
    }

    /** The next base in ACGT order, for planting a mismatch that is still a real base. */
    static byte mutateBase(byte base) {
        switch (base) {
            case 'A': return 'C';
            case 'C': return 'G';
            case 'G': return 'T';
            default: return 'A';
        }
    }

    /** One `QualityYieldMetrics` row, written the way `CollectQualityYieldMetrics` writes it. */
    static void writeQualityYield(File f, long totalReads, long pfReads, int readLength,
                                  long totalBases, long pfBases, long q20, long pfQ20, long q30,
                                  long pfQ30, long q20Yield, long pfQ20Yield) throws Exception {
        try (PrintWriter out = new PrintWriter(f, "UTF-8")) {
            out.print("## htsjdk.samtools.metrics.StringHeader\n");
            out.print("# CollectQualityYieldMetrics INPUT=/work/fixtures/small.bam OUTPUT=/work/out/output.txt\n");
            out.print("## htsjdk.samtools.metrics.StringHeader\n");
            out.print("# Started on: Mon Sep 07 00:00:00 UTC 2026\n");
            out.print("\n");
            out.print("## METRICS CLASS\tpicard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics\n");
            out.print("TOTAL_READS\tPF_READS\tREAD_LENGTH\tTOTAL_BASES\tPF_BASES\tQ20_BASES\tPF_Q20_BASES\tQ30_BASES\tPF_Q30_BASES\tQ20_EQUIVALENT_YIELD\tPF_Q20_EQUIVALENT_YIELD\n");
            out.printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n",
                    totalReads, pfReads, readLength, totalBases, pfBases, q20, pfQ20, q30, pfQ30,
                    q20Yield, pfQ20Yield);
            out.print("\n");
        }
    }

    /**
     * Put the first `length` bases of an adapter at the end of a read, IN READ ORDER.
     *
     * A record on the negative strand stores its bases reverse complemented, and the tool searches
     * what the sequencer read rather than what the file stores: it reverse complements a copy
     * before it looks. Planting into the stored bases would therefore put the adapter at the front
     * of half the reads, where the search never looks, so the plant is done on the read-order copy
     * and complemented back.
     */
    static void plantAdapter(SAMRecord read, String adapter, int length) throws Exception {
        byte[] bases = read.getReadBases();
        if (read.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);
        byte[] planted = adapter.substring(0, length).getBytes("UTF-8");
        System.arraycopy(planted, 0, bases, bases.length - planted.length, planted.length);
        if (read.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);
        read.setReadBases(bases);
    }

    static java.util.List<SAMRecord> reads(SAMFileHeader header, boolean coordinateSorted) {
        Random rng = new Random(20260729L);
        java.util.List<SAMRecord> out = new java.util.ArrayList<>();
        String alphabet = "ACGT";

        for (int i = 0; i < READS; i += 2) {
            String name = String.format("read%04d", i);
            boolean chr2 = i % 8 == 0;
            int contig = chr2 ? 1 : 0;
            int limit = (chr2 ? CHR2 : CHR1) - READ_LENGTH - 10;
            int start = 1 + rng.nextInt(limit);

            SAMRecord first = new SAMRecord(header);
            SAMRecord second = new SAMRecord(header);
            for (SAMRecord r : new SAMRecord[] {first, second}) {
                byte[] bases = new byte[READ_LENGTH];
                byte[] quals = new byte[READ_LENGTH];
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases[b] = (byte) alphabet.charAt(rng.nextInt(4));
                    // Qualities span the CollectQualityYieldMetrics thresholds (Q20, Q30) rather
                    // than sitting above both, which would make those counters constant.
                    quals[b] = (byte) (2 + rng.nextInt(38));
                }
                // A no-call every so often, so the base-distribution and N-handling paths run.
                if (i % 10 == 0) bases[3] = 'N';
                r.setReadName(name);
                r.setReadBases(bases);
                r.setBaseQualities(quals);
                r.setReferenceIndex(contig);
                r.setMappingQuality(i % 12 == 0 ? 0 : 20 + rng.nextInt(40));
                r.setAttribute("RG", i % 4 == 0 ? "rg2" : "rg1");
            }

            first.setAlignmentStart(start);
            second.setAlignmentStart(start + 100 <= limit ? start + 100 : start);
            first.setCigarString(cigarFor(i));
            second.setCigarString(READ_LENGTH + "M");
            first.setReadPairedFlag(true);
            second.setReadPairedFlag(true);
            first.setFirstOfPairFlag(true);
            second.setSecondOfPairFlag(true);
            first.setReadNegativeStrandFlag(i % 3 == 0);
            second.setReadNegativeStrandFlag(!first.getReadNegativeStrandFlag());
            first.setProperPairFlag(i % 5 != 0);
            second.setProperPairFlag(i % 5 != 0);
            first.setDuplicateReadFlag(i % 14 == 0);
            second.setDuplicateReadFlag(i % 14 == 0);
            // Secondary and supplementary are properties of an *alignment*, so htsjdk's validation
            // rejects them on an unmapped read ("Supplementary alignment flag should not be set for
            // unaligned read"). The first covering-array run hit that: nine of eleven rows failed,
            // and one class of failure was this fixture rather than the argument under test. A
            // corpus that is invalid under STRICT tests the validator, not the tool.
            boolean secondUnmapped = i % 20 == 0;
            if (i % 16 == 0) first.setNotPrimaryAlignmentFlag(true);
            if (i % 18 == 0 && !secondUnmapped) second.setSupplementaryAlignmentFlag(true);

            // One pair in twenty has an unmapped mate: the paired-metrics and mate-info paths
            // behave differently there, and a corpus of clean pairs never reaches them.
            if (secondUnmapped) {
                second.setReadUnmappedFlag(true);
                second.setAlignmentStart(first.getAlignmentStart());
                second.setCigarString("*");
                second.setMappingQuality(0);
                first.setMateUnmappedFlag(true);
            }
            SamPairUtil.setMateInfo(first, second, false);
            out.add(first);
            out.add(second);
        }

        if (coordinateSorted) {
            out.sort(new SAMRecordCoordinateComparator());
        } else {
            out.sort(new SAMRecordQueryNameComparator());
        }
        return out;
    }

    /** Soft clips, an insertion and a deletion, so cigar-walking tools take more than one branch. */
    static String cigarFor(int i) {
        switch (i % 6) {
            case 0: return READ_LENGTH + "M";
            case 1: return "5S" + (READ_LENGTH - 5) + "M";
            case 2: return (READ_LENGTH - 8) + "M8S";
            case 3: return "20M2I" + (READ_LENGTH - 22) + "M";
            case 4: return "20M3D" + (READ_LENGTH - 20) + "M";
            default: return "10M5N" + (READ_LENGTH - 10) + "M";
        }
    }

    static java.util.List<SAMRecord> unmapped(SAMFileHeader header) {
        java.util.List<SAMRecord> out = new java.util.ArrayList<>();
        Random rng = new Random(20260731L);
        for (int i = 0; i < 40; i += 2) {
            SAMRecord first = new SAMRecord(header);
            SAMRecord second = new SAMRecord(header);
            for (SAMRecord r : new SAMRecord[] {first, second}) {
                byte[] bases = new byte[READ_LENGTH];
                byte[] quals = new byte[READ_LENGTH];
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases[b] = (byte) "ACGT".charAt(rng.nextInt(4));
                    quals[b] = (byte) (2 + rng.nextInt(38));
                }
                r.setReadName(String.format("unmapped%04d", i));
                r.setReadBases(bases);
                r.setBaseQualities(quals);
                r.setReadUnmappedFlag(true);
                r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                r.setMappingQuality(0);
                r.setAttribute("RG", "rg1");
                r.setReadPairedFlag(true);
                r.setMateUnmappedFlag(true);
            }
            first.setFirstOfPairFlag(true);
            second.setSecondOfPairFlag(true);
            out.add(first);
            out.add(second);
        }
        return out;
    }

    /**
     * The caller passes -Dsamjdk.try_use_intel_deflater=false: the fixture must be
     * byte-reproducible, and the GKL deflater emits different bytes than zlib for the same input.
     * The oracle contract pins the JDK deflater for the same reason.
     *
     * Only the coordinate-sorted BAM is indexed; indexing a queryname-sorted or unsorted file is
     * an error, not an option.
     */
    static void writeBam(File f, SAMFileHeader header, java.util.List<SAMRecord> records,
                         boolean index) {
        SAMFileWriterFactory factory = new SAMFileWriterFactory().setUseAsyncIo(false);
        factory.setCreateIndex(index);
        try (SAMFileWriter w = factory.makeBAMWriter(header, true, f)) {
            for (SAMRecord r : records) w.addAlignment(r);
        }
    }

    static void writeSam(File f, SAMFileHeader header, java.util.List<SAMRecord> records) {
        try (SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(header, true, f)) {
            for (SAMRecord r : records) w.addAlignment(r);
        }
    }

    static void writeIntervals(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("@HD\tVN:1.6");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", CHR1);
            p.printf("@SQ\tSN:chr2\tLN:%d%n", CHR2);
            p.println("chr1\t100\t400\t+\ttarget1");
            p.println("chr1\t900\t1200\t+\ttarget2");
            p.println("chr2\t50\t200\t-\ttarget3");
        }
    }

    /**
     * The same three targets as the interval list, in BED coordinates.
     *
     * A BED start is 0-based and its end is exclusive, where an interval list is 1-based and
     * inclusive, so the same target is one lower on the left here. Writing both from one set of
     * numbers is the point: a tool that converts between them can then be checked against a
     * fixture that says what the answer is, rather than against its own output.
     */
    static void writeBed(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("chr1\t99\t400\ttarget1\t0\t+");
            p.println("chr1\t899\t1200\ttarget2\t0\t+");
            p.println("chr2\t49\t200\ttarget3\t0\t-");
        }
    }

    /**
     * A FASTA whose headers carry a description and whose lines are not the output length.
     *
     * ref.fasta has bare contig names and is already wrapped at the length NormalizeFasta writes,
     * so TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE has nothing to truncate and normalizing is the
     * identity: the array covers both arguments without observing either. Here each header is
     * "name description", so truncation changes the header line, and the bases are wrapped at 37
     * rather than 100, so normalizing rewraps them.
     *
     * It deliberately has no .fai beside it. With one, ReferenceSequenceFileFactory opens the
     * indexed reader, whose index would have to agree with the names; without one it opens
     * FastaSequenceFile, which is the path this port reproduces.
     */
    static void writeDescribedFasta(File f, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println(">seq1 first sequence, described");
            for (int i = 0; i < chr2.length(); i += 37) {
                p.println(chr2.substring(i, Math.min(i + 37, chr2.length())));
            }
            p.println(">seq2\ta tab-separated description");
            for (int i = 0; i < 120; i += 37) {
                p.println(chr2.substring(i, Math.min(i + 37, 120)));
            }
        }
    }

    /**
     * An interval list whose order is not the coordinate order.
     *
     * targets.interval_list is already sorted, so SORT produces the same file with it on or off
     * and an array over that argument covers it without testing it. Here chr2 leads, the chr1
     * entries are out of order, and both strands appear, so sorting moves lines and the
     * strand-then-name tiebreak of IntervalCoordinateComparator is reachable.
     */
    static void writeMixedIntervals(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("@HD\tVN:1.6");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", CHR1);
            p.printf("@SQ\tSN:chr2\tLN:%d%n", CHR2);
            p.println("chr2\t50\t200\t-\ttargetB");
            p.println("chr1\t300\t500\t+\ttargetC");
            p.println("chr1\t100\t400\t+\ttargetA");
            p.println("chr1\t600\t700\t-\ttargetD");
        }
    }

    /**
     * A BED the interval tools' arguments can actually be observed on.
     *
     * targets.bed is already sorted, disjoint and length-nonzero, so SORT, UNIQUE and
     * KEEP_LENGTH_ZERO_INTERVALS all produce the same file on it: the array covers those
     * arguments without testing them, which the runner says out loud. This one is built so that
     * each of the three changes the output.
     *
     * Out of coordinate order, so SORT moves lines. Two overlapping features and two abutting
     * ones, so UNIQUE merges and concatenates names. One feature whose BED start equals its end,
     * which becomes `start == end + 1` and is dropped unless KEEP_LENGTH_ZERO_INTERVALS is set.
     */
    static void writeMixedBed(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("chr2\t49\t200\ttargetB\t0\t-");
            p.println("chr1\t299\t500\ttargetC\t0\t+");
            p.println("chr1\t99\t400\ttargetA\t0\t+");
            p.println("chr1\t599\t700\ttargetD\t0\t+");
            p.println("chr1\t699\t800\ttargetE\t0\t+");
            p.println("chr1\t900\t900\tzeroLength\t0\t+");
        }
    }

    /**
     * The sequence dictionary, as its own file.
     *
     * `SAMSequenceDictionaryExtractor` reads a FASTA's dictionary through
     * `ReferenceSequenceFileFactory`, which does not derive one: it looks for the `.dict` beside
     * the reference and throws "Could not find dictionary next to reference file" when there is
     * none. Every tool taking a SEQUENCE_DICTIONARY therefore needed this file before it could be
     * given an array at all.
     */
    static void writeVcf(File f, String chr1, String chr2, boolean withGenotypes) throws Exception {
        writeVcf(f, chr1, chr2, withGenotypes, 2);
    }

    static void writeVcf(File f, String chr1, String chr2, boolean withGenotypes, int sampleCount)
            throws Exception {
        htsjdk.samtools.SAMSequenceDictionary dict = new htsjdk.samtools.SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", chr1.length()));
        dict.addSequence(new SAMSequenceRecord("chr2", chr2.length()));

        java.util.Set<htsjdk.variant.vcf.VCFHeaderLine> lines = new java.util.LinkedHashSet<>();
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "GT", 1, htsjdk.variant.vcf.VCFHeaderLineType.String, "Genotype"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "GQ", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Genotype quality"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "DP", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Depth"));
        lines.add(new htsjdk.variant.vcf.VCFInfoHeaderLine(
                "AC", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Allele count"));
        lines.add(new htsjdk.variant.vcf.VCFFilterHeaderLine("LowQual", "Low quality"));

        java.util.List<String> samples = new java.util.ArrayList<>();
        if (withGenotypes) {
            for (int i = 1; i <= sampleCount; i++) samples.add("sample" + i);
        }
        htsjdk.variant.vcf.VCFHeader header = new htsjdk.variant.vcf.VCFHeader(lines, samples);
        header.setSequenceDictionary(dict);

        try (htsjdk.variant.variantcontext.writer.VariantContextWriter writer =
                     new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder()
                             .setOutputFile(f)
                             .setReferenceDictionary(dict)
                             .setOption(htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY)
                             .build()) {
            writer.writeHeader(header);
            // The first four are on chr1 (2,000 bases), the last two on chr2 (1,000).
            int[] positions = {100, 300, 500, 700, 200, 600};
            for (int i = 0; i < positions.length; i++) {
                // The sites-only file keeps every other variant, so half of the full file is known.
                if (!withGenotypes && i % 2 == 1) continue;
                String contig = i < 4 ? "chr1" : "chr2";
                int position = positions[i];
                String reference = String.valueOf((i < 4 ? chr1 : chr2).charAt(position - 1));
                boolean indel = i == 3;
                htsjdk.variant.variantcontext.Allele ref = htsjdk.variant.variantcontext.Allele
                        .create(indel ? reference + "AT" : reference, true);
                htsjdk.variant.variantcontext.Allele alt = htsjdk.variant.variantcontext.Allele
                        .create(indel ? reference : (reference.equals("A") ? "G" : "A"), false);
                htsjdk.variant.variantcontext.VariantContextBuilder builder =
                        new htsjdk.variant.variantcontext.VariantContextBuilder()
                                .chr(contig)
                                .start(position)
                                .stop(position + ref.length() - 1)
                                .alleles(java.util.Arrays.asList(ref, alt))
                                .attribute("AC", 1 + (i % 2));
                if (i == 5) builder.filter("LowQual");
                if (withGenotypes) {
                    java.util.List<htsjdk.variant.variantcontext.Genotype> genotypes =
                            new java.util.ArrayList<>();
                    for (int g = 0; g < samples.size(); g++) {
                        java.util.List<htsjdk.variant.variantcontext.Allele> called = (i + g) % 3 == 0
                                ? java.util.Arrays.asList(ref, ref)
                                : ((i + g) % 3 == 1
                                        ? java.util.Arrays.asList(ref, alt)
                                        : java.util.Arrays.asList(alt, alt));
                        genotypes.add(new htsjdk.variant.variantcontext.GenotypeBuilder(
                                samples.get(g), called)
                                .GQ(20 + 7 * ((i + g) % 5))
                                .DP(10 + ((i + g) % 4))
                                .make());
                    }
                    builder.genotypes(genotypes);
                }
                writer.add(builder.make());
            }
        }
    }

    static void writeDict(File f, String chr1, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("@HD\tVN:1.6\tSO:unsorted");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", chr1.length());
            p.printf("@SQ\tSN:chr2\tLN:%d%n", chr2.length());
        }
    }

    /** A contig's bytes in the file: its bases, plus the newline that ends each line. */
    static long bytesOnDisk(int bases, int lineWidth) {
        long lines = (bases + lineWidth - 1) / lineWidth;
        return bases + lines;
    }

    static void writeFastq(File f, int end) throws Exception {
        Random rng = new Random(20260800L + end);
        try (PrintWriter p = new PrintWriter(f)) {
            for (int i = 0; i < 40; i++) {
                StringBuilder bases = new StringBuilder();
                StringBuilder quals = new StringBuilder();
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases.append("ACGT".charAt(rng.nextInt(4)));
                    quals.append((char) (33 + 2 + rng.nextInt(38)));
                }
                p.printf("@fq%04d/%d%n%s%n+%n%s%n", i, end, bases, quals);
            }
        }
    }
}
