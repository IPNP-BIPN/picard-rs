/*
 * Oracle dump harness for FastqToSam (unpaired, default options, SAM output) conformance in
 * picard-rs.
 *
 * Emits an escaped TSV to stdout: a `fastq` row (the input) and a `sam` row (the SAM file
 * FastqToSam produced). The committed corpus is `java ... FastqToSamDump | gzip > ...`, regenerated
 * and compared in CI. FastqToSam writes no @PG and no timestamp, so the whole SAM is compared raw.
 *
 * The input reads are in a deliberately non-queryname order (r2, r10, r1, r3) so the output must be
 * re-sorted, and one carries a /1 suffix to exercise the read-name cleanup.
 *
 *   java -cp picard-fat.jar:. FastqToSamDump
 */
import htsjdk.samtools.SAMFileHeader;

import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class FastqToSamDump {
    public static void main(final String[] args) throws Exception {
        final File fastq = File.createTempFile("f2s-", ".fastq");
        fastq.deleteOnExit();
        try (PrintStream ps = new PrintStream(fastq)) {
            ps.print("@r2 a comment\nACGT\n+\nIIII\n");
            ps.print("@r10\nAACC\n+\n5555\n");
            ps.print("@r1\nGGTTA\n+\n#$%&'\n");
            ps.print("@r3/1\nACGTAC\n+\nIIIIII\n");
        }

        final File out = File.createTempFile("f2s-", ".sam");
        out.deleteOnExit();
        final int rc = new picard.sam.FastqToSam().instanceMain(new String[] {
                "FASTQ=" + fastq.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "SAMPLE_NAME=s1",
                "QUALITY_FORMAT=Standard",
        });
        if (rc != 0) {
            System.err.println("FastqToSam exited " + rc);
            System.exit(rc);
        }

        emit("fastq", "unpaired", new String(Files.readAllBytes(fastq.toPath())));
        emit("sam", "unpaired", new String(Files.readAllBytes(out.toPath())));
    }

    private static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
