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

        SAMFileHeader unmappedHeader = header(SAMFileHeader.SortOrder.unsorted);
        writeBam(new File(dir, "unmapped.bam"), unmappedHeader, unmapped(unmappedHeader), false);

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
