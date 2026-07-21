/*
 * Builds BAMs, runs Picard's CollectQualityYieldMetrics on them, and emits both.
 *
 * Output, one line per case:
 *
 *     bam     <TAB> <case> <TAB> <hex of the input BAM>
 *     metrics <TAB> <case> <TAB> <metrics file, newlines as \n>
 *
 * The input BAM is emitted so the Rust side can confirm it is measuring the same bytes.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class QualityYieldDump {

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

    static void emit(String name, List<SAMRecord> records, String... extraArgs) throws Exception {
        SAMFileHeader h = header();
        File bam = File.createTempFile("qy-", ".bam");
        bam.deleteOnExit();
        SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false).setCreateMd5File(false).setUseAsyncIo(false)
                .makeBAMWriter(h, true, bam);
        for (SAMRecord r : records) w.addAlignment(r);
        w.close();

        File out = File.createTempFile("qy-", ".txt");
        out.deleteOnExit();
        List<String> args = new ArrayList<>();
        args.add("INPUT=" + bam.getPath());
        args.add("OUTPUT=" + out.getPath());
        args.add("USE_JDK_DEFLATER=true");
        args.add("USE_JDK_INFLATER=true");
        for (String a : extraArgs) args.add(a);

        int rc = new picard.analysis.CollectQualityYieldMetrics()
                .instanceMain(args.toArray(new String[0]));
        if (rc != 0) throw new IllegalStateException(name + ": tool exited " + rc);

        System.out.println("bam\t" + name + "\t" + hex(Files.readAllBytes(bam.toPath())));
        System.out.println("metrics\t" + name + "\t"
                + new String(Files.readAllBytes(out.toPath()))
                        .replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
    }

    static SAMRecord read(SAMFileHeader h, String name, int start, byte[] quals, int flags) {
        SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setReferenceIndex(0);
        r.setAlignmentStart(start);
        r.setMappingQuality(60);
        r.setCigarString(quals.length + "M");
        byte[] bases = new byte[quals.length];
        for (int i = 0; i < bases.length; i++) bases[i] = (byte) "ACGT".charAt(i % 4);
        r.setReadBases(bases);
        r.setBaseQualities(quals);
        r.setFlags(flags);
        return r;
    }

    static byte[] quals(int n, int base) {
        byte[] q = new byte[n];
        for (int i = 0; i < n; i++) q[i] = (byte) ((base + i) % 45);
        return q;
    }

    public static void main(String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        SAMFileHeader h = header();

        emit("empty", new ArrayList<>());

        emit("one_read", List.of(read(h, "r1", 100, quals(50, 20), 0)));

        // A spread of qualities across the Q20 and Q30 thresholds.
        List<SAMRecord> spread = new ArrayList<>();
        for (int i = 0; i < 200; i++) spread.add(read(h, "r" + i, 1 + i * 7, quals(50, i), 0));
        emit("quality_spread", spread);

        // Vendor-failed reads, which count toward TOTAL but not PF.
        List<SAMRecord> pf = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            pf.add(read(h, "ok" + i, 1 + i * 7, quals(50, i), 0));
            pf.add(read(h, "fail" + i, 1 + i * 7, quals(50, i), 0x200));
        }
        emit("vendor_failed", pf);

        // Secondary and supplementary alignments, excluded by default.
        List<SAMRecord> flagged = new ArrayList<>();
        for (int i = 0; i < 50; i++) {
            flagged.add(read(h, "p" + i, 1 + i * 7, quals(50, i), 0));
            flagged.add(read(h, "s" + i, 1 + i * 7, quals(50, i), 0x100));
            flagged.add(read(h, "u" + i, 1 + i * 7, quals(50, i), 0x800));
        }
        emit("secondary_supplementary", flagged);
        emit("secondary_included", flagged, "INCLUDE_SECONDARY_ALIGNMENTS=true");
        emit("supplemental_included", flagged, "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true");

        // Varying read lengths, so READ_LENGTH is an integer division with a remainder.
        List<SAMRecord> lengths = new ArrayList<>();
        for (int i = 0; i < 30; i++) lengths.add(read(h, "L" + i, 1 + i * 7, quals(10 + i * 3, i), 0));
        emit("varied_lengths", lengths);

        // OQ tags, which USE_ORIGINAL_QUALITIES prefers by default.
        List<SAMRecord> oq = new ArrayList<>();
        for (int i = 0; i < 50; i++) {
            SAMRecord r = read(h, "o" + i, 1 + i * 7, quals(50, i), 0);
            r.setOriginalBaseQualities(quals(50, i + 10));
            oq.add(r);
        }
        emit("original_qualities", oq);
        emit("original_qualities_off", oq, "USE_ORIGINAL_QUALITIES=false");
    }
}
