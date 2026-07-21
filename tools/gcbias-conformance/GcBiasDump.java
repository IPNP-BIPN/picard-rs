/*
 * Builds BAMs and reference FASTAs, runs Picard's CollectGcBiasMetrics, and emits all of it.
 *
 * Output, one line per case:
 *   bam     <TAB> <case> <TAB> <hex of the input BAM>
 *   fasta   <TAB> <case> <TAB> <the reference bases, one contig>
 *   detail  <TAB> <case> <TAB> <the detail metrics file, escaped>
 *   summary <TAB> <case> <TAB> <the summary metrics file, escaped>
 *
 * The cases are built around where the collector and the definitions disagree: which window a
 * read is charged to on each strand, what happens at the edges of the computed window range, and
 * the bins where a divide-by-zero is avoided by a guard rather than by the data.
 */

import htsjdk.samtools.*;
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.*;

public class GcBiasDump {

    static final int REF_LENGTH = 2000;

    static String hex(byte[] b) {
        StringBuilder sb = new StringBuilder();
        for (byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    static String esc(String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static SAMFileHeader header() {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", REF_LENGTH));
        h.setSequenceDictionary(d);
        h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("s"); rg.setLibrary("l"); rg.setPlatform("ILLUMINA");
        h.addReadGroup(rg);
        return h;
    }

    static SAMRecord read(SAMFileHeader h, String name, int start, String cigar, String bases,
                          int flags, int mapq) {
        SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setFlags(flags);
        r.setReferenceIndex(0);
        r.setAlignmentStart(start);
        r.setMappingQuality(mapq);
        r.setCigarString(cigar);
        r.setReadString(bases);
        byte[] q = new byte[bases.length()];
        Arrays.fill(q, (byte) 30);
        r.setBaseQualities(q);
        r.setAttribute("RG", "rg1");
        if ((flags & 0x4) != 0) {
            // A truly unmapped record: no reference, no position, so it sorts last under
            // coordinate order rather than colliding with a placed read.
            r.setReadUnmappedFlag(true);
            r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            r.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            r.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
        }
        return r;
    }

    static String repeat(char c, int n) {
        char[] a = new char[n];
        Arrays.fill(a, c);
        return new String(a);
    }

    /** A reference that sweeps through GC content so many bins are populated. */
    static String gradientReference() {
        StringBuilder sb = new StringBuilder();
        Random rng = new Random(20260721L);
        for (int block = 0; block < 20; block++) {
            int gcPercent = block * 5;
            for (int i = 0; i < 100; i++) {
                sb.append(rng.nextInt(100) < gcPercent ? (rng.nextBoolean() ? 'G' : 'C')
                                                       : (rng.nextBoolean() ? 'A' : 'T'));
            }
        }
        while (sb.length() < REF_LENGTH) sb.append('A');
        return sb.substring(0, REF_LENGTH);
    }

    static void emit(String name, List<SAMRecord> records, String refBases) throws Exception {
        SAMFileHeader h = header();
        File bam = File.createTempFile("gcbias", ".bam");
        try (SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(false)
                .makeBAMWriter(h, true, bam)) {
            for (SAMRecord r : records) w.addAlignment(r);
        }
        File fasta = File.createTempFile("gcbias", ".fasta");
        try (PrintWriter p = new PrintWriter(fasta)) {
            p.println(">chr1");
            for (int i = 0; i < refBases.length(); i += 60)
                p.println(refBases.substring(i, Math.min(i + 60, refBases.length())));
        }
        try (PrintWriter p = new PrintWriter(fasta.getPath() + ".fai")) {
            p.printf("chr1\t%d\t6\t60\t61%n", refBases.length());
        }
        try (PrintWriter p = new PrintWriter(fasta.getPath().replaceAll("\\.fasta$", "") + ".dict")) {
            p.println("@HD\tVN:1.6\tSO:unsorted");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", refBases.length());
        }
        File detail = File.createTempFile("gcbias", ".detail.txt");
        File summary = File.createTempFile("gcbias", ".summary.txt");
        File chart = File.createTempFile("gcbias", ".pdf");
        int rc = new picard.analysis.CollectGcBiasMetrics().instanceMain(new String[] {
                "INPUT=" + bam.getPath(),
                "OUTPUT=" + detail.getPath(),
                "SUMMARY_OUTPUT=" + summary.getPath(),
                "CHART_OUTPUT=" + chart.getPath(),
                "REFERENCE_SEQUENCE=" + fasta.getPath(),
                "ASSUME_SORTED=true",
        });
        if (rc != 0) { System.err.println("case " + name + " exited " + rc); return; }

        System.out.println("bam\t" + name + "\t" + hex(Files.readAllBytes(bam.toPath())));
        System.out.println("fasta\t" + name + "\t" + refBases);
        System.out.println("detail\t" + name + "\t" + esc(new String(Files.readAllBytes(detail.toPath()))));
        System.out.println("summary\t" + name + "\t" + esc(new String(Files.readAllBytes(summary.toPath()))));
        bam.delete(); fasta.delete(); detail.delete(); summary.delete(); chart.delete();
    }

    public static void main(String[] args) throws Exception {
        SAMFileHeader h = header();
        String ref = gradientReference();

        emit("one_forward_read",
             List.of(read(h, "r1", 500, "100M", repeat('A', 100), 0, 60)), ref);

        // The same bases on each strand, which the collector charges to different windows.
        emit("forward_and_reverse_same_bases", List.of(
                read(h, "f", 500, "100M", repeat('A', 100), 0, 60),
                read(h, "r", 500, "100M", repeat('A', 100), 0x10, 60)), ref);

        // A read whose window start is past lastWindowStart, charged to bin 0.
        emit("past_last_window_start",
             List.of(read(h, "r1", 1950, "40M", repeat('A', 40), 0, 60)), ref);

        // A reverse read near the contig start, dropped from the bins but counted as aligned.
        emit("reverse_near_contig_start",
             List.of(read(h, "r1", 1, "50M", repeat('A', 50), 0x10, 60)), ref);

        // Unmapped reads bump totalClusters without reaching a bin.
        emit("unmapped_counts_as_cluster", List.of(
                read(h, "m", 500, "100M", repeat('A', 100), 0, 60),
                read(h, "u", 1, "*", repeat('A', 100), 0x4, 0)), ref);

        // Indels and mismatches, which feed errorsByGc and so MEAN_BASE_QUALITY.
        emit("with_errors", List.of(
                read(h, "e1", 500, "40M5D40M", repeat('G', 80), 0, 60),
                read(h, "e2", 600, "40M5I45M", repeat('C', 90), 0, 60)), ref);

        // A spread of reads over the gradient, so many bins are populated and the dropout and
        // normalized-coverage arithmetic has something to work on.
        List<SAMRecord> spread = new ArrayList<>();
        Random rng = new Random(4242L);
        for (int i = 0; i < 400; i++) {
            int start = 1 + rng.nextInt(REF_LENGTH - 200);
            int flags = rng.nextInt(100) < 50 ? 0x10 : 0;
            if (rng.nextInt(100) < 50) flags |= 0x1 | 0x2 | (rng.nextBoolean() ? 0x40 : 0x80);
            SAMRecord r = read(h, "s" + i, start, "100M", repeat('A', 100), flags, 60);
            if ((flags & 0x1) != 0) {
                r.setMateReferenceIndex(0);
                r.setMateAlignmentStart(Math.max(1, start + 150));
                r.setInferredInsertSize(250);
            }
            spread.add(r);
        }
        spread.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        emit("spread_over_gradient", spread, ref);

        // A read binned into a GC bin the reference never produces. The reference below has no
        // window at 0% GC, and the read starts past lastWindowStart so it is charged to bin 0
        // anyway. calculateGcNormCoverage skips bins with no windows from *both* its numerator
        // and its window total, so this is the only case that distinguishes that skip from
        // including them - without it, a port that included them passes.
        // Exactly one G every 50 bases, so every 100-base window holds exactly two and every
        // window bins at 2. Bin 0 therefore has *no* windows while the 0-19 range as a whole
        // has many - which is what makes the zero-window skip in calculateGcNormCoverage
        // observable. A random reference does not: at high GC the whole 0-19 range is empty and
        // the function returns 0 either way, and at low GC bin 0 picks up windows of its own.
        StringBuilder noZero = new StringBuilder();
        for (int i = 0; i < REF_LENGTH; i++) noZero.append(i % 50 == 0 ? 'G' : 'A');
        emit("read_in_a_bin_the_reference_never_reaches",
             List.of(read(h, "r1", 1950, "40M", repeat('A', 40), 0, 60)), noZero.toString());

        // Paired reads, where only the first of pair bumps totalClusters.
        List<SAMRecord> pairs = new ArrayList<>();
        for (int i = 0; i < 20; i++) {
            SAMRecord a = read(h, "p" + i, 300 + i * 10, "100M", repeat('A', 100), 0x1 | 0x2 | 0x40, 60);
            SAMRecord b = read(h, "p" + i, 300 + i * 10, "100M", repeat('A', 100), 0x1 | 0x2 | 0x80, 60);
            for (SAMRecord r : List.of(a, b)) {
                r.setMateReferenceIndex(0);
                r.setMateAlignmentStart(300 + i * 10);
                r.setInferredInsertSize(100);
            }
            pairs.add(a); pairs.add(b);
        }
        emit("paired_clusters", pairs, ref);
    }
}
