/*
 * Builds BAMs, runs Picard's QualityScoreDistribution, and emits both.
 *
 * Output, one line per case:
 *   bam     <TAB> <case> <TAB> <hex of the input BAM>
 *   metrics <TAB> <case> <TAB> <metrics file, escaped>
 *
 * The tool writes no metric rows at all, only histograms, so this exercises a shape none of the
 * other suites do: a MetricsFile whose entire body is one or two histogram tables keyed on
 * java.lang.Byte. The OQ histogram is emitted only when it is non-empty, so the same tool
 * produces a one-column table for most files and a two-column one when OQ tags are present -
 * and the two histograms' key sets differ, which is exactly the union the writer has to sort.
 */

import htsjdk.samtools.*;
import java.io.File;
import java.nio.file.Files;
import java.util.*;

public class QualityScoreDistributionDump {

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
        d.addSequence(new SAMSequenceRecord("chr1", 100000));
        h.setSequenceDictionary(d);
        h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("s"); rg.setLibrary("lib1"); rg.setPlatform("ILLUMINA");
        h.addReadGroup(rg);
        return h;
    }

    static SAMRecord read(SAMFileHeader h, String name, int start, String bases, byte[] quals,
                          int flags, String oq) {
        SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setFlags(flags);
        r.setReferenceIndex(0);
        r.setAlignmentStart(start);
        r.setMappingQuality(60);
        r.setCigarString(bases.length() + "M");
        r.setReadString(bases);
        r.setBaseQualities(quals);
        r.setAttribute("RG", "rg1");
        if (oq != null) r.setAttribute("OQ", oq);
        if ((flags & 0x4) != 0) {
            r.setReadUnmappedFlag(true);
            r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            r.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            r.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
        }
        return r;
    }

    static byte[] quals(int... v) {
        byte[] q = new byte[v.length];
        for (int i = 0; i < v.length; i++) q[i] = (byte) v[i];
        return q;
    }

    static String oqOf(int... v) {
        StringBuilder sb = new StringBuilder();
        for (int x : v) sb.append((char) (x + 33));
        return sb.toString();
    }

    static void emit(String name, List<SAMRecord> records, String... extraArgs) throws Exception {
        SAMFileHeader h = header();
        File bam = File.createTempFile("qsdist", ".bam");
        try (SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(false)
                .makeBAMWriter(h, true, bam)) {
            for (SAMRecord r : records) w.addAlignment(r);
        }
        File out = File.createTempFile("qsdist", ".txt");
        File chart = File.createTempFile("qsdist", ".pdf");
        List<String> args = new ArrayList<>(List.of(
                "INPUT=" + bam.getPath(), "OUTPUT=" + out.getPath(),
                "CHART_OUTPUT=" + chart.getPath(), "ASSUME_SORTED=true"));
        args.addAll(Arrays.asList(extraArgs));
        int rc = new picard.analysis.QualityScoreDistribution()
                .instanceMain(args.toArray(new String[0]));
        if (rc != 0) { System.err.println("case " + name + " exited " + rc); return; }
        System.out.println("bam\t" + name + "\t" + hex(Files.readAllBytes(bam.toPath())));
        System.out.println("metrics\t" + name + "\t" + esc(new String(Files.readAllBytes(out.toPath()))));
        bam.delete(); out.delete(); chart.delete();
    }

    public static void main(String[] args) throws Exception {
        SAMFileHeader h = header();

        emit("simple", List.of(
                read(h, "r1", 1, "ACGT", quals(30, 30, 20, 10), 0, null)));

        // A spread of qualities, so the histogram has many bins in ascending byte order.
        List<SAMRecord> spread = new ArrayList<>();
        Random rng = new Random(7L);
        for (int i = 0; i < 50; i++) {
            char[] b = new char[20];
            byte[] q = new byte[20];
            for (int j = 0; j < 20; j++) {
                b[j] = "ACGT".charAt(rng.nextInt(4));
                q[j] = (byte) rng.nextInt(41);
            }
            spread.add(read(h, "s" + i, 1 + i, new String(b), q, 0, null));
        }
        emit("spread", spread);

        // No-call bases are skipped by default, so the qualities at those positions vanish.
        emit("no_calls_excluded", List.of(
                read(h, "r1", 1, "ANNT", quals(30, 5, 6, 30), 0, null)));
        emit("no_calls_included", List.of(
                read(h, "r1", 1, "ANNT", quals(30, 5, 6, 30), 0, null)),
             "INCLUDE_NO_CALLS=true");

        // With OQ, a second histogram appears - and its key set differs from the first's, which
        // is the union the writer has to sort rather than concatenate.
        emit("with_oq", List.of(
                read(h, "r1", 1, "ACGT", quals(30, 30, 30, 30), 0, oqOf(2, 3, 40, 41))));

        // OQ on one read only: the OQ histogram is sparse next to the Q one.
        emit("oq_on_one_read_only", List.of(
                read(h, "r1", 1, "ACGT", quals(30, 30, 30, 30), 0, oqOf(11, 12, 13, 14)),
                read(h, "r2", 2, "ACGT", quals(31, 31, 31, 31), 0, null)));

        // Secondary and supplementary records are filtered out entirely.
        emit("secondary_and_supplementary_filtered", List.of(
                read(h, "r1", 1, "ACGT", quals(30, 30, 30, 30), 0, null),
                read(h, "r2", 2, "ACGT", quals(1, 1, 1, 1), 0x100, null),
                read(h, "r3", 3, "ACGT", quals(2, 2, 2, 2), 0x800, null)));

        // Vendor-failed and unmapped reads are included by default and excluded on request.
        List<SAMRecord> mixed = List.of(
                read(h, "pass", 1, "ACGT", quals(30, 30, 30, 30), 0, null),
                read(h, "fail", 2, "ACGT", quals(11, 11, 11, 11), 0x200, null),
                read(h, "unmap", 1, "ACGT", quals(12, 12, 12, 12), 0x4, null));
        emit("all_reads_by_default", mixed);
        emit("pf_reads_only", mixed, "PF_READS_ONLY=true");
        emit("aligned_reads_only", mixed, "ALIGNED_READS_ONLY=true");

        // Quality zero, the lowest bin.
        emit("quality_zero", List.of(
                read(h, "r1", 1, "ACGT", quals(0, 0, 30, 30), 0, null)));
    }
}
