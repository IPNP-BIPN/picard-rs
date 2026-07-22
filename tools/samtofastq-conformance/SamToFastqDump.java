/*
 * Oracle dump harness for SamToFastq (unpaired, default options) conformance in picard-rs.
 *
 * Emits an escaped TSV to stdout: a `bam` row (hex of the input BAM) and a `fastq` row (the FASTQ
 * file SamToFastq produced). The committed corpus is `java ... SamToFastqDump | gzip > ...`,
 * regenerated and compared in CI. The Rust test decodes the BAM, runs the unpaired path, and must
 * reproduce the FASTQ bytes.
 *
 * The reads are all unpaired and mapped (the SAM validator forbids secondary/supplementary flags on
 * unmapped reads). They exercise a forward read, a negative-strand read (reverse-complemented with
 * its qualities reversed), an N base, and secondary / supplementary / vendor-fail reads that are
 * dropped by default.
 *
 *   java -cp picard-fat.jar:. SamToFastqDump
 */
import htsjdk.samtools.*;

import java.io.File;
import java.nio.file.Files;

public class SamToFastqDump {
    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("s1");
        header.addReadGroup(rg);

        final File bam = File.createTempFile("s2f-", ".bam");
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

        final File out = File.createTempFile("s2f-", ".fastq");
        out.deleteOnExit();
        final int rc = new picard.sam.SamToFastq().instanceMain(new String[] {
                "INPUT=" + bam.getAbsolutePath(),
                "FASTQ=" + out.getAbsolutePath(),
        });
        if (rc != 0) {
            System.err.println("SamToFastq exited " + rc);
            System.exit(rc);
        }

        emit("bam", "unpaired", hex(Files.readAllBytes(bam.toPath())));
        emit("fastq", "unpaired", esc(new String(Files.readAllBytes(out.toPath()))));
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
