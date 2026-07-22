/*
 * Oracle dump harness for AddOrReplaceReadGroups (SAM I/O, default sort order) conformance in
 * picard-rs.
 *
 * Emits an escaped TSV to stdout: an `input` row (the SAM, already carrying an @RG ID:1 and RG:Z:1
 * tags) and an `output` row (after replacing the read group with ID:2). The committed corpus is
 * `java ... AddReplaceRgDump | gzip > ...`, regenerated and compared in CI. The tool adds no @PG and
 * no timestamp and does not re-sort with the default SORT_ORDER, so both SAMs are compared raw.
 *
 *   java -cp picard-fat.jar:. AddReplaceRgDump
 */
import htsjdk.samtools.*;

import java.io.File;
import java.nio.file.Files;

public class AddReplaceRgDump {
    public static void main(final String[] args) throws Exception {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        builder.setRandomSeed(0);
        final int chr1 = builder.getHeader().getSequenceIndex("chr1");
        builder.addFrag("r1", chr1, 100, false);
        builder.addFrag("r2", chr1, 200, false);
        builder.addFrag("r3", chr1, 300, true);

        final File input = File.createTempFile("arrg-in-", ".sam");
        input.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, input);
        for (final SAMRecord rec : builder.getRecords()) w.addAlignment(rec);
        w.close();

        final File out = File.createTempFile("arrg-out-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.AddOrReplaceReadGroups().instanceMain(new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "RGID=2",
                "RGLB=lib1",
                "RGPL=ILLUMINA",
                "RGPU=unit1",
                "RGSM=sample1",
        });
        if (rc != 0) {
            System.err.println("AddOrReplaceReadGroups exited " + rc);
            System.exit(rc);
        }

        emit("input", "case", new String(Files.readAllBytes(input.toPath())));
        emit("output", "case", new String(Files.readAllBytes(out.toPath())));
    }

    private static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
