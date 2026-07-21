/*
 * Probe: which reference window does CollectGcBiasMetrics charge a read to?
 *
 * GcBiasMetricsCollector.addRead:
 *
 *     final int pos = rec.getReadNegativeStrandFlag()
 *         ? rec.getAlignmentEnd() - scanWindowSize
 *         : rec.getAlignmentStart();
 *     if (pos > 0) { final int windowGc = gc[pos]; ... }
 *
 * So a forward read is charged to the window starting at its alignment start, and a reverse read
 * to the window *ending* at its alignment end. riker's ERRATA describes "a forward-strand window
 * binning fix", so this is a place an independent reimplementation chose to differ.
 *
 * Two consequences are checked here, both against a reference built so the answer is legible:
 * the first 100 bases are pure AT (GC 0) and the next 100 pure GC (GC 100), so which window a
 * read lands in is readable straight off the GC bin it reports.
 *
 *   1. A forward and a reverse read covering the *same* bases bin differently.
 *   2. `gc[pos]` is read from an array whose entries at index 0 and at or past lastWindowStart
 *      were never computed and are still zero. A read whose window falls there is charged to
 *      GC bin 0 rather than skipped.
 */

import htsjdk.samtools.*;
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.*;

public class WindowBinningProbe {

    static final int REF_LENGTH = 400;

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

    static SAMRecord read(SAMFileHeader h, String name, int start, int len, boolean reverse) {
        SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setReferenceIndex(0);
        r.setAlignmentStart(start);
        r.setCigarString(len + "M");
        r.setMappingQuality(60);
        char[] b = new char[len];
        Arrays.fill(b, 'A');
        r.setReadString(new String(b));
        byte[] q = new byte[len];
        Arrays.fill(q, (byte) 30);
        r.setBaseQualities(q);
        r.setFlags(reverse ? 0x10 : 0);
        r.setAttribute("RG", "rg1");
        return r;
    }

    static void run(String label, List<SAMRecord> reads) throws Exception {
        SAMFileHeader h = header();

        // 100 AT, then 100 GC, then 200 AT. Window size is 100, so window 1 is ~0% GC and
        // window 101 is ~100% GC.
        StringBuilder ref = new StringBuilder();
        for (int i = 0; i < 100; i++) ref.append('A');
        for (int i = 0; i < 100; i++) ref.append('G');
        for (int i = 0; i < 200; i++) ref.append('A');

        File fasta = File.createTempFile("gcprobe", ".fasta");
        try (PrintWriter p = new PrintWriter(fasta)) {
            p.println(">chr1");
            for (int i = 0; i < REF_LENGTH; i += 60)
                p.println(ref.substring(i, Math.min(i + 60, REF_LENGTH)));
        }
        try (PrintWriter p = new PrintWriter(fasta.getPath() + ".fai")) {
            p.printf("chr1\t%d\t6\t60\t61%n", REF_LENGTH);
        }
        try (PrintWriter p = new PrintWriter(fasta.getPath().replaceAll("\\.fasta$", "") + ".dict")) {
            p.println("@HD\tVN:1.6\tSO:unsorted");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", REF_LENGTH);
        }

        File bam = File.createTempFile("gcprobe", ".bam");
        try (SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(false)
                .makeBAMWriter(h, true, bam)) {
            for (SAMRecord r : reads) w.addAlignment(r);
        }

        File detail = File.createTempFile("gcprobe", ".detail.txt");
        File summary = File.createTempFile("gcprobe", ".summary.txt");
        File chart = File.createTempFile("gcprobe", ".pdf");
        new picard.analysis.CollectGcBiasMetrics().instanceMain(new String[] {
                "INPUT=" + bam.getPath(),
                "OUTPUT=" + detail.getPath(),
                "SUMMARY_OUTPUT=" + summary.getPath(),
                "CHART_OUTPUT=" + chart.getPath(),
                "REFERENCE_SEQUENCE=" + fasta.getPath(),
                "ASSUME_SORTED=true",
        });

        // Report the GC bins that received read starts.
        String text = new String(Files.readAllBytes(detail.toPath()));
        StringBuilder bins = new StringBuilder();
        boolean inRows = false;
        for (String line : text.split("\n")) {
            if (line.startsWith("ACCUMULATION_LEVEL")) { inRows = true; continue; }
            if (!inRows || line.isEmpty() || line.startsWith("#")) continue;
            String[] f = line.split("\t");
            // ACCUMULATION_LEVEL READS_USED GC WINDOWS READ_STARTS ...
            if (f.length > 4 && !f[4].equals("0")) bins.append(" gc=").append(f[2])
                    .append(" starts=").append(f[4]);
        }
        System.out.println(label + ":" + bins);

        bam.delete(); fasta.delete(); detail.delete(); summary.delete(); chart.delete();
    }

    public static void main(String[] args) throws Exception {
        SAMFileHeader h = header();

        // Both cover reference bases 101..200, the pure-GC stretch. Forward is charged to the
        // window at its start (101, GC 100); reverse to the window at end-100 (100, GC ~0).
        run("forward_over_gc_block", List.of(read(h, "f", 101, 100, false)));
        run("reverse_over_gc_block", List.of(read(h, "r", 101, 100, true)));

        // A read whose window start is past lastWindowStart, where calculateAllGcs never wrote.
        run("past_last_window_start", List.of(read(h, "p", 320, 50, false)));

        // A reverse read near the contig start, where alignmentEnd - windowSize is <= 0 and the
        // read is dropped from the GC bins entirely while still counting as an aligned read.
        run("reverse_near_contig_start", List.of(read(h, "n", 1, 50, true)));
    }
}
