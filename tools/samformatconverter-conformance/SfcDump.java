/*
 * Oracle dump harness for SamFormatConverter (BAM -> SAM) conformance in picard-rs.
 *
 * Builds a small coordinate-sorted BAM with varied reads (mixed flags, cigars, tags), converts it to
 * SAM with SamFormatConverter, and emits a `bam` row (the BAM bytes, hex) and a `sam` row (the SAM
 * text). SamFormatConverter adds no @PG and no timestamp, so the SAM is compared raw; the BAM is the
 * input the port decodes.
 *
 * The JDK deflater is forced before writing the BAM so the committed hex matches what CI regenerates
 * (GKL's igzip emits different BGZF bytes than zlib).
 *
 *   java -cp picard-fat.jar:. SfcDump
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.file.Files;

public class SfcDump {
    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dict);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("s");
        header.addReadGroup(rg);

        final File bam = File.createTempFile("sfc-", ".bam");
        bam.deleteOnExit();
        try (SAMFileWriter w = new SAMFileWriterFactory().makeBAMWriter(header, true, bam)) {
            w.addAlignment(rec(header, "r0", 0, 100, "100M", 60, false));
            w.addAlignment(rec(header, "r1", 147, 200, "50M5D50M", 60, true));
            w.addAlignment(rec(header, "r2", 99, 300, "40M10I50M", 30, true));
            w.addAlignment(rec(header, "r3", 4, 0, "*", 0, false)); // unmapped
        }

        final File sam = File.createTempFile("sfc-", ".sam");
        sam.deleteOnExit();
        final int rc = new picard.sam.SamFormatConverter().instanceMain(new String[] {
                "INPUT=" + bam.getAbsolutePath(),
                "OUTPUT=" + sam.getAbsolutePath(),
                "VALIDATION_STRINGENCY=SILENT",
        });
        if (rc != 0) { System.err.println("SamFormatConverter exited " + rc); System.exit(rc); }

        emit("bam", "case", hex(Files.readAllBytes(bam.toPath())));
        emit("sam", "case", esc(new String(Files.readAllBytes(sam.toPath()))));
    }

    static SAMRecord rec(final SAMFileHeader h, final String name, final int flags, final int pos,
                         final String cigar, final int mapq, final boolean tag) {
        final SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setFlags(flags);
        if (pos > 0) {
            r.setReferenceIndex(0);
            r.setAlignmentStart(pos);
            r.setCigarString(cigar);
            r.setMappingQuality(mapq);
        } else {
            r.setReferenceIndex(-1);
            r.setAlignmentStart(0);
        }
        r.setReadBases("ACGTACGTAC".getBytes());
        r.setBaseQualities(new byte[] {40, 40, 30, 30, 20, 20, 41, 45, 10, 15});
        r.setAttribute("RG", "rg1");
        if (tag) r.setAttribute("MQ", 42);
        return r;
    }

    static String hex(final byte[] b) {
        final StringBuilder sb = new StringBuilder(b.length * 2);
        for (final byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    static String esc(final String p) {
        return p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static void emit(final String kind, final String kase, final String payload) {
        System.out.println(kind + "\t" + kase + "\t" + payload);
    }
}
