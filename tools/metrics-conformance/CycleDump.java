/*
 * Builds BAMs and runs both cycle-based collectors on them: MeanQualityByCycle and
 * CollectBaseDistributionByCycle.
 *
 * Two tools in one harness on purpose. They are stratum-mates in the calibration gate, and
 * running them over identical inputs is what makes the second one's cost a delta rather than a
 * separate measurement.
 *
 * Output: <tool> <TAB> <case> <TAB> <metrics file, escaped>, plus one `bam` line per case.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class CycleDump {

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
        d.addSequence(new SAMSequenceRecord("chr1", 250000000));
        h.setSequenceDictionary(d);
        h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return h;
    }

    static SAMRecord read(SAMFileHeader h, String name, int start, byte[] bases, byte[] quals,
                          int flags, String oq) {
        SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setReferenceIndex(0);
        r.setAlignmentStart(start);
        r.setMappingQuality(60);
        r.setCigarString(bases.length + "M");
        r.setReadBases(bases);
        r.setBaseQualities(quals);
        r.setFlags(flags);
        if (oq != null) r.setAttribute("OQ", oq);
        if ((flags & 0x1) != 0) {
            // A paired record needs a resolvable mate, or htsjdk refuses to write it.
            r.setMateReferenceIndex(0);
            r.setMateAlignmentStart(start + 100);
            r.setInferredInsertSize(((flags & 0x40) != 0) ? 250 : -250);
        }
        return r;
    }

    static void emit(String name, List<SAMRecord> records) throws Exception {
        SAMFileHeader h = header();
        File bam = File.createTempFile("cy-", ".bam");
        bam.deleteOnExit();
        SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false).setCreateMd5File(false).setUseAsyncIo(false)
                .makeBAMWriter(h, false, bam);
        for (SAMRecord r : records) w.addAlignment(r);
        w.close();
        System.out.println("bam\t" + name + "\t" + hex(Files.readAllBytes(bam.toPath())));

        for (String tool : new String[]{"MeanQualityByCycle", "CollectBaseDistributionByCycle"}) {
            File out = File.createTempFile("cy-", ".txt");
            out.deleteOnExit();
            File chart = File.createTempFile("cy-", ".pdf");
            chart.deleteOnExit();
            List<String> args = new ArrayList<>();
            args.add("INPUT=" + bam.getPath());
            args.add("OUTPUT=" + out.getPath());
            args.add("CHART=" + chart.getPath());
            args.add("USE_JDK_DEFLATER=true");
            args.add("USE_JDK_INFLATER=true");
            args.add("ASSUME_SORTED=true");

            picard.cmdline.CommandLineProgram p = tool.equals("MeanQualityByCycle")
                    ? new picard.analysis.MeanQualityByCycle()
                    : new picard.analysis.CollectBaseDistributionByCycle();
            int rc = p.instanceMain(args.toArray(new String[0]));
            if (rc != 0) throw new IllegalStateException(name + "/" + tool + ": exit " + rc);
            System.out.println(tool + "\t" + name + "\t"
                    + esc(new String(Files.readAllBytes(out.toPath()))));
        }
    }

    static byte[] bases(int n, int seed) {
        byte[] b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = (byte) "ACGTN".charAt((seed + i) % 5);
        return b;
    }

    /* Lower case and IUPAC ambiguity codes. baseToInt folds case and drops everything that is
     * not one of the four into PCT_N, and a corpus of only upper-case ACGTN cannot tell a
     * case-folding implementation from a case-sensitive one. Found by sabotage. */
    static byte[] mixedBases(int n, int seed) {
        final String alphabet = "AcGtNaCgTnMRSVWYHKDB=";
        byte[] b = new byte[n];
        for (int i = 0; i < n; i++) b[i] = (byte) alphabet.charAt((seed + i) % alphabet.length());
        return b;
    }

    static byte[] quals(int n, int base) {
        byte[] q = new byte[n];
        for (int i = 0; i < n; i++) q[i] = (byte) ((base + i) % 45);
        return q;
    }

    public static void main(String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        SAMFileHeader h = header();

        emit("single", List.of(read(h, "r1", 100, bases(50, 0), quals(50, 20), 0, null)));

        // Unpaired reads on both strands: the cycle is reversed for a reverse-strand read.
        List<SAMRecord> strands = new ArrayList<>();
        for (int i = 0; i < 200; i++) {
            strands.add(read(h, "s" + i, 100 + i * 60, bases(50, i), quals(50, i),
                    (i % 2 == 0) ? 0 : 0x10, null));
        }
        emit("both_strands", strands);

        // Paired reads: the second end's cycles are offset past the first end's length.
        List<SAMRecord> paired = new ArrayList<>();
        for (int i = 0; i < 200; i++) {
            paired.add(read(h, "p" + i, 100 + i * 60, bases(50, i), quals(50, i), 0x1 | 0x2 | 0x40, null));
            paired.add(read(h, "p" + i, 200 + i * 60, bases(50, i + 1), quals(50, i + 3), 0x1 | 0x2 | 0x80, null));
        }
        emit("paired", paired);

        // Reads of differing lengths, which grow the per-cycle arrays.
        List<SAMRecord> lengths = new ArrayList<>();
        Random rng = new Random(20260721L);
        for (int i = 0; i < 300; i++) {
            int n = 20 + rng.nextInt(120);
            lengths.add(read(h, "L" + i, 100 + i * 200, bases(n, i), quals(n, i), 0, null));
        }
        emit("varied_lengths", lengths);

        // OQ tags, which MeanQualityByCycle reports as a second histogram column.
        List<SAMRecord> oq = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < 50; j++) sb.append((char) (33 + ((i + j) % 40)));
            oq.add(read(h, "o" + i, 100 + i * 60, bases(50, i), quals(50, i), 0, sb.toString()));
        }
        emit("original_qualities", oq);

        // Vendor-failed, unmapped, secondary and supplementary records.
        List<SAMRecord> flagged = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            flagged.add(read(h, "a" + i, 100 + i * 60, bases(50, i), quals(50, i), 0, null));
            flagged.add(read(h, "b" + i, 100 + i * 60, bases(50, i), quals(50, i), 0x200, null));
            flagged.add(read(h, "c" + i, 100 + i * 60, bases(50, i), quals(50, i), 0x100, null));
            flagged.add(read(h, "d" + i, 100 + i * 60, bases(50, i), quals(50, i), 0x800, null));
        }
        emit("flagged", flagged);

        // Lower case and ambiguity codes, so a case-sensitive baseToInt is distinguishable.
        List<SAMRecord> mixed = new ArrayList<>();
        for (int i = 0; i < 150; i++) {
            mixed.add(read(h, "x" + i, 100 + i * 60, mixedBases(50, i), quals(50, i),
                    (i % 3 == 0) ? 0x10 : 0, null));
        }
        emit("mixed_case_bases", mixed);
    }
}
