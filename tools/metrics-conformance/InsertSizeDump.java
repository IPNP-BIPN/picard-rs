/*
 * Builds BAMs of paired reads, runs Picard's CollectInsertSizeMetrics, and emits both.
 *
 * Output, one line per case:
 *   bam     <TAB> <case> <TAB> <hex of the input BAM>
 *   metrics <TAB> <case> <TAB> <metrics file, \n and \t escaped>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class InsertSizeDump {

    static String hex(byte[] b) {
        StringBuilder sb = new StringBuilder();
        for (byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    static SAMFileHeader header() {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", 250000000));
        h.setSequenceDictionary(d);
        h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return h;
    }

    /** One pair. `orientation` is 0 for FR, 1 for RF, 2 for TANDEM. */
    static void addPair(SAMFileHeader h, List<SAMRecord> out, String name,
                        int start, int insert, int orientation, int extraFlags) {
        final int len = 50;
        for (int mate = 0; mate < 2; mate++) {
            SAMRecord r = new SAMRecord(h);
            r.setReadName(name);
            r.setReferenceIndex(0);
            r.setMateReferenceIndex(0);
            r.setMappingQuality(60);
            r.setCigarString(len + "M");
            byte[] bases = new byte[len];
            java.util.Arrays.fill(bases, (byte) 'A');
            r.setReadBases(bases);
            byte[] q = new byte[len];
            java.util.Arrays.fill(q, (byte) 30);
            r.setBaseQualities(q);

            boolean first = (mate == 0);
            int flags = 0x1 | 0x2 | extraFlags | (first ? 0x40 : 0x80);
            // FR: first forward, second reverse. RF: first reverse, second forward.
            // TANDEM: both on the same strand.
            boolean firstReverse = (orientation == 1);
            boolean secondReverse = (orientation == 0);
            if (orientation == 2) { firstReverse = false; secondReverse = false; }
            boolean thisReverse = first ? firstReverse : secondReverse;
            boolean mateReverse = first ? secondReverse : firstReverse;
            if (thisReverse) flags |= 0x10;
            if (mateReverse) flags |= 0x20;
            r.setFlags(flags);

            int leftStart = start;
            int rightStart = start + insert - len;
            r.setAlignmentStart(first ? leftStart : rightStart);
            r.setMateAlignmentStart(first ? rightStart : leftStart);
            r.setInferredInsertSize(first ? insert : -insert);
            out.add(r);
        }
    }

    static void emit(String name, List<SAMRecord> records, String... extra) throws Exception {
        SAMFileHeader h = header();
        File bam = File.createTempFile("is-", ".bam");
        bam.deleteOnExit();
        SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false).setCreateMd5File(false).setUseAsyncIo(false)
                .makeBAMWriter(h, false, bam);
        for (SAMRecord r : records) w.addAlignment(r);
        w.close();

        File out = File.createTempFile("is-", ".txt");
        out.deleteOnExit();
        List<String> args = new ArrayList<>();
        args.add("INPUT=" + bam.getPath());
        args.add("OUTPUT=" + out.getPath());
        args.add("USE_JDK_DEFLATER=true");
        args.add("USE_JDK_INFLATER=true");
        args.add("ASSUME_SORTED=true");
        File chart = File.createTempFile("is-", ".pdf");
        chart.deleteOnExit();
        args.add("HISTOGRAM_FILE=" + chart.getPath());
        for (String a : extra) args.add(a);

        int rc = new picard.analysis.CollectInsertSizeMetrics().instanceMain(args.toArray(new String[0]));
        if (rc != 0) throw new IllegalStateException(name + ": exit " + rc);

        System.out.println("bam\t" + name + "\t" + hex(Files.readAllBytes(bam.toPath())));
        System.out.println("metrics\t" + name + "\t"
                + new String(Files.readAllBytes(out.toPath()))
                        .replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
    }

    public static void main(String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        SAMFileHeader h = header();

        List<SAMRecord> one = new ArrayList<>();
        addPair(h, one, "p1", 1000, 300, 0, 0);
        emit("one_pair", one);

        // A unimodal distribution, the normal case.
        List<SAMRecord> normal = new ArrayList<>();
        Random rng = new Random(20260721L);
        for (int i = 0; i < 2000; i++) {
            int insert = (int) Math.max(60, Math.round(350 + rng.nextGaussian() * 60));
            addPair(h, normal, "n" + i, 1000 + i * 20, insert, 0, 0);
        }
        emit("normal", normal);

        // All three orientations present.
        List<SAMRecord> mixed = new ArrayList<>();
        for (int i = 0; i < 600; i++) {
            int insert = 200 + (i % 100);
            addPair(h, mixed, "m" + i, 1000 + i * 20, insert, i % 3, 0);
        }
        emit("mixed_orientations", mixed);

        // A minority orientation below the default MINIMUM_PCT, which must be dropped.
        List<SAMRecord> rare = new ArrayList<>();
        for (int i = 0; i < 1000; i++) addPair(h, rare, "a" + i, 1000 + i * 20, 300 + (i % 50), 0, 0);
        for (int i = 0; i < 3; i++)    addPair(h, rare, "b" + i, 900000 + i * 20, 500, 1, 0);
        emit("rare_orientation", rare);
        emit("rare_orientation_pct0", rare, "MINIMUM_PCT=0");

        // Duplicates, excluded by default.
        List<SAMRecord> dups = new ArrayList<>();
        for (int i = 0; i < 500; i++) addPair(h, dups, "d" + i, 1000 + i * 20, 300, 0, 0);
        for (int i = 0; i < 500; i++) addPair(h, dups, "x" + i, 1000 + i * 20, 900, 0, 0x400);
        emit("duplicates", dups);
        emit("duplicates_included", dups, "INCLUDE_DUPLICATES=true");

        // A long tail, so the histogram trim actually removes bins.
        List<SAMRecord> tail = new ArrayList<>();
        for (int i = 0; i < 1000; i++) addPair(h, tail, "t" + i, 1000 + i * 20, 300 + (i % 20), 0, 0);
        for (int i = 0; i < 20; i++)   addPair(h, tail, "o" + i, 800000 + i * 20, 5000 + i * 500, 0, 0);
        emit("long_tail", tail);
        emit("long_tail_fixed_width", tail, "HISTOGRAM_WIDTH=1000");
        emit("long_tail_deviations", tail, "DEVIATIONS=2");
    }
}
