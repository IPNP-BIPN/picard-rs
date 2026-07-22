/*
 * Oracle dump harness for SamToFastq (default options) conformance in picard-rs.
 *
 * Emits an escaped TSV to stdout, two cases:
 *   unpaired: a `bam` row and a `fastq` row (single FASTQ output).
 *   paired:   a `bam` row and `fastq_r1`/`fastq_r2` rows (FASTQ + SECOND_END_FASTQ).
 *
 * The committed corpus is `java ... SamToFastqDump | gzip > ...`, regenerated and compared in CI.
 * FASTQ output has no header, so every row is compared raw.
 *
 *   java -cp picard-fat.jar:. SamToFastqDump
 */
import htsjdk.samtools.*;

import java.io.File;
import java.nio.file.Files;

public class SamToFastqDump {
    public static void main(final String[] args) throws Exception {
        emitUnpaired();
        emitPaired();
    }

    /** Unpaired reads to a single FASTQ, exercising RE_REVERSE and the default filters. */
    private static void emitUnpaired() throws Exception {
        final SAMFileHeader header = header();
        final File bam = File.createTempFile("s2f-u-", ".bam");
        bam.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory().makeBAMWriter(header, false, bam);
        w.addAlignment(read(header, "fwd", 0, "ACGTACGT", new int[] {40, 40, 30, 30, 20, 20, 41, 45}));
        w.addAlignment(read(header, "rev", SAMFlag.READ_REVERSE_STRAND.intValue(), "ACGTN",
                new int[] {40, 30, 20, 10, 5}));
        w.addAlignment(read(header, "withN", 0, "ANNGT", new int[] {40, 2, 2, 30, 20}));
        w.addAlignment(read(header, "sec", SAMFlag.SECONDARY_ALIGNMENT.intValue(), "AC", new int[] {40, 40}));
        w.addAlignment(read(header, "sup", SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue(), "AC", new int[] {40, 40}));
        w.addAlignment(read(header, "qcfail", SAMFlag.READ_FAILS_VENDOR_QUALITY_CHECK.intValue(), "AC", new int[] {40, 40}));
        w.close();

        final File out = File.createTempFile("s2f-u-", ".fastq");
        out.deleteOnExit();
        run(new String[] {"INPUT=" + bam.getAbsolutePath(), "FASTQ=" + out.getAbsolutePath()});

        emit("bam", "unpaired", hex(Files.readAllBytes(bam.toPath())));
        emit("fastq", "unpaired", esc(new String(Files.readAllBytes(out.toPath()))));
    }

    /** Paired reads to two FASTQs, exercising the /1 /2 split and RE_REVERSE on a mate. */
    private static void emitPaired() throws Exception {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.queryname);
        builder.setRandomSeed(0);
        final int idx = builder.getHeader().getSequenceIndex("chr1");
        builder.addPair("pairA", idx, 100, 300);
        builder.addPair("pairB", idx, 500, 700);

        final File bam = File.createTempFile("s2f-p-", ".bam");
        bam.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory().makeBAMWriter(builder.getHeader(), false, bam);
        for (final SAMRecord r : builder.getRecords()) w.addAlignment(r);
        w.close();

        final File r1 = File.createTempFile("s2f-p1-", ".fastq");
        final File r2 = File.createTempFile("s2f-p2-", ".fastq");
        r1.deleteOnExit();
        r2.deleteOnExit();
        run(new String[] {
                "INPUT=" + bam.getAbsolutePath(),
                "FASTQ=" + r1.getAbsolutePath(),
                "SECOND_END_FASTQ=" + r2.getAbsolutePath(),
        });

        emit("bam", "paired", hex(Files.readAllBytes(bam.toPath())));
        emit("fastq_r1", "paired", esc(new String(Files.readAllBytes(r1.toPath()))));
        emit("fastq_r2", "paired", esc(new String(Files.readAllBytes(r2.toPath()))));
    }

    private static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("s1");
        header.addReadGroup(rg);
        return header;
    }

    private static SAMRecord read(final SAMFileHeader header, final String name, final int extraFlags,
                                  final String bases, final int[] quals) {
        final SAMRecord r = new SAMRecord(header);
        r.setReadName(name);
        r.setReferenceName("chr1");
        r.setAlignmentStart(100);
        r.setCigarString(bases.length() + "M");
        r.setReadBases(bases.getBytes());
        final byte[] q = new byte[quals.length];
        for (int i = 0; i < quals.length; i++) q[i] = (byte) quals[i];
        r.setBaseQualities(q);
        r.setAttribute("RG", "rg1");
        r.setFlags(r.getFlags() | extraFlags);
        return r;
    }

    private static void run(final String[] toolArgs) {
        final int rc = new picard.sam.SamToFastq().instanceMain(toolArgs);
        if (rc != 0) {
            System.err.println("SamToFastq exited " + rc);
            System.exit(rc);
        }
    }

    private static String hex(final byte[] b) {
        final StringBuilder sb = new StringBuilder(b.length * 2);
        for (final byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    private static String esc(final String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    private static void emit(final String kind, final String kase, final String payload) {
        System.out.println(kind + "\t" + kase + "\t" + payload);
    }
}
