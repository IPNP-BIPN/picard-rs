/*
 * Probe: does the iteration order of CollectRnaSeqMetrics' overlapping-gene set reach the output?
 *
 * geneOverlapDetector.getOverlaps(readInterval) returns a Set<Gene> (a HashSet). The collector
 * iterates it (RnaSeqMetricsCollector line 197) and, for each overlapping transcript, writes a
 * per-base LocusFunction via Gene.Transcript.assignLocusFunctionForRange.
 *
 * Reading that method (Gene.java line 161) shows it only overwrites a base with a *higher*
 * ordinal LocusFunction - a max reduction, which is commutative. And coverage is accumulated per
 * transcript into that transcript's own array, which is also commutative. If both readings are
 * right, the HashSet's iteration order cannot reach the output, and the 1273-line red-black
 * IntervalTree behind the OverlapDetector does not need a byte-exact port: any correct
 * overlap-set structure suffices.
 *
 * The claim is testable without reproducing the tree: run the tool twice on the same input and
 * compare the bytes. HashSet<Gene> iteration order is not stable across constructions (Gene does
 * not override hashCode, so it inherits identity hashCode, which varies run to run), so if the
 * order reached the output the two runs would differ. Determinism across runs is therefore
 * evidence that it does not.
 *
 * The input is built so that two genes overlap the same bases with *different* locus functions,
 * which is the only situation where last-writer-wins and max-reduction diverge - the case that
 * would expose an order dependency if one existed.
 *
 * Prints TWO_RUNS_IDENTICAL=true|false and, on false, the first differing line.
 */

import htsjdk.samtools.*;
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.*;

public class OverlapOrderProbe {

    static final int REF_LENGTH = 5000;

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

    static void writeRefFlat(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            // geneName transcript chrom strand txStart txEnd cdsStart cdsEnd exonCount starts ends
            // Two genes overlapping the same region, one coding, one with a UTR there, on opposite
            // strands, so their per-base locus functions differ.
            p.println("GENEA\ttxA\tchr1\t+\t1000\t2000\t1200\t1800\t1\t1000,\t2000,");
            p.println("GENEB\ttxB\tchr1\t-\t1000\t2000\t1000\t1000\t1\t1000,\t2000,");
            // A dozen more overlapping genes to make the HashSet non-trivial and its order matter
            // if it ever could.
            for (int i = 0; i < 12; i++) {
                p.printf("GENE%d\ttx%d\tchr1\t+\t%d\t%d\t%d\t%d\t1\t%d,\t%d,%n",
                        i, i, 1000 + i * 20, 2000 - i * 10, 1100 + i * 20, 1900 - i * 10,
                        1000 + i * 20, 2000 - i * 10);
            }
        }
    }

    static File makeBam(SAMFileHeader h) throws Exception {
        File bam = File.createTempFile("rnaprobe", ".bam");
        List<SAMRecord> reads = new ArrayList<>();
        Random rng = new Random(11L);
        for (int i = 0; i < 500; i++) {
            SAMRecord r = new SAMRecord(h);
            r.setReadName("r" + i);
            r.setReferenceIndex(0);
            int start = 1000 + rng.nextInt(1000);
            r.setAlignmentStart(start);
            r.setCigarString("100M");
            r.setMappingQuality(60);
            char[] b = new char[100]; Arrays.fill(b, 'A');
            r.setReadString(new String(b));
            byte[] q = new byte[100]; Arrays.fill(q, (byte) 30);
            r.setBaseQualities(q);
            r.setFlags(rng.nextInt(100) < 50 ? 0x10 : 0);
            r.setAttribute("RG", "rg1");
            reads.add(r);
        }
        reads.sort(Comparator.comparingInt(SAMRecord::getAlignmentStart));
        try (SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(false)
                .makeBAMWriter(h, true, bam)) {
            for (SAMRecord r : reads) w.addAlignment(r);
        }
        return bam;
    }

    static String run(File bam, File refFlat, File fasta) throws Exception {
        File out = File.createTempFile("rnaout", ".txt");
        int rc = new picard.analysis.CollectRnaSeqMetrics().instanceMain(new String[] {
                "INPUT=" + bam.getPath(),
                "OUTPUT=" + out.getPath(),
                "REF_FLAT=" + refFlat.getPath(),
                "STRAND_SPECIFICITY=NONE",
                "REFERENCE_SEQUENCE=" + fasta.getPath(),
                "ASSUME_SORTED=true",
        });
        if (rc != 0) return "EXIT_" + rc;
        String text = new String(Files.readAllBytes(out.toPath()));
        out.delete();
        // Drop the two run-time header lines.
        StringBuilder sb = new StringBuilder();
        for (String line : text.split("\n"))
            if (!line.startsWith("# CollectRnaSeqMetrics") && !line.startsWith("# Started on:"))
                sb.append(line).append("\n");
        return sb.toString();
    }

    public static void main(String[] args) throws Exception {
        SAMFileHeader h = header();

        char[] ref = new char[REF_LENGTH];
        Arrays.fill(ref, 'A');
        for (int i = 0; i < REF_LENGTH; i += 3) ref[i] = 'G';
        File fasta = File.createTempFile("rnaref", ".fasta");
        try (PrintWriter p = new PrintWriter(fasta)) {
            p.println(">chr1");
            for (int i = 0; i < REF_LENGTH; i += 60)
                p.println(new String(ref, i, Math.min(60, REF_LENGTH - i)));
        }
        try (PrintWriter p = new PrintWriter(fasta.getPath() + ".fai")) {
            p.printf("chr1\t%d\t6\t60\t61%n", REF_LENGTH);
        }
        try (PrintWriter p = new PrintWriter(fasta.getPath().replaceAll("\\.fasta$", "") + ".dict")) {
            p.println("@HD\tVN:1.6\tSO:unsorted");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", REF_LENGTH);
        }

        File refFlat = File.createTempFile("refflat", ".txt");
        writeRefFlat(refFlat);
        File bam = makeBam(h);

        // Two independent runs, each in a fresh JVM-object graph, so any identity-hash-driven
        // HashSet order would differ between them.
        String a = run(bam, refFlat, fasta);
        String b = run(bam, refFlat, fasta);

        System.out.println("TWO_RUNS_IDENTICAL=" + a.equals(b));
        if (!a.equals(b)) {
            String[] al = a.split("\n"), bl = b.split("\n");
            for (int i = 0; i < Math.min(al.length, bl.length); i++) {
                if (!al[i].equals(bl[i])) {
                    System.out.println("  run1: " + al[i]);
                    System.out.println("  run2: " + bl[i]);
                    break;
                }
            }
        }
        // Print the metric line so the probe also shows the tool did real work.
        for (String line : a.split("\n")) {
            if (line.startsWith("PF_BASES") || (line.length() > 0 && Character.isDigit(line.charAt(0)))) {
                System.out.println("  metrics: " + line.substring(0, Math.min(80, line.length())));
            }
        }
    }
}
