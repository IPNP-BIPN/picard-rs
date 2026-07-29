/*
 * The rejections: rows the reference refuses, and the message it refuses them with.
 *
 * A covering array is generated over combinations the tool accepts, because rows spent being
 * rejected are rows not spent on the tool (picard-rs decision 0009). That leaves the rejections
 * themselves untested, and they are behaviour: a port that happily processes a queryname-sorted
 * BAM where Picard refuses is not a byte-identical port, it is a different tool with the same
 * name. This harness holds them.
 *
 * Each case runs the tool through PicardCommandLine, the same entry point the arrays use, and
 * records the exit code together with the exception class and message. The message is compared,
 * because it is what a user and a pipeline see.
 *
 * Usage: RejectionDump <fixture directory>
 */

import picard.cmdline.PicardCommandLine;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;

public class RejectionDump {

    /**
     * `PicardCommandLine.instanceMain` is protected and `main` calls System.exit, which would end
     * the dump after the first case. A subclass reaches the protected entry point without the
     * exit, so every case runs and each one's outcome is recorded.
     */
    static class Runner extends PicardCommandLine {
        int run(String[] argv) {
            Object code = instanceMain(argv);
            return code instanceof Integer ? (Integer) code : (code == null ? 0 : 1);
        }
    }


    /**
     * The output path is deliberately outside the fixture directory: the runner mounts fixtures
     * read-only, so a tool asked to write there fails on the mount instead of reaching the
     * behaviour under test. The fourth case caught that, because it asserts its message by name.
     */
    static final String OUT = "/tmp/rejection-out.txt";

    /** case name, then the argument list. */
    static String[][] cases(String dir) {
        return new String[][] {
            {
                "qualityyield_flow_mode_obsolete",
                "CollectQualityYieldMetrics",
                "--INPUT", dir + "/small.bam",
                "--OUTPUT", OUT,
                "--FLOW_MODE", "true",
            },
            {
                "qualityyield_queryname_not_assumed_sorted",
                "CollectQualityYieldMetrics",
                "--INPUT", dir + "/queryname.bam",
                "--OUTPUT", OUT,
                "--ASSUME_SORTED", "false",
            },
            {
                "alignmentsummary_queryname_not_assumed_sorted",
                "CollectAlignmentSummaryMetrics",
                "--INPUT", dir + "/queryname.bam",
                "--OUTPUT", OUT,
                "--REFERENCE_SEQUENCE", dir + "/ref.fasta",
                "--ASSUME_SORTED", "false",
            },
            {
                // ASSUME_SORTED=true gets past the header check, and the reference walk then fails
                // on the first record that goes backwards. A different message, a different code
                // path, and the one a user is most likely to hit by accident.
                "alignmentsummary_queryname_assumed_sorted",
                "CollectAlignmentSummaryMetrics",
                "--INPUT", dir + "/queryname.bam",
                "--OUTPUT", OUT,
                "--REFERENCE_SEQUENCE", dir + "/ref.fasta",
                "--ASSUME_SORTED", "true",
            },
        };
    }

    public static void main(String[] args) {
        String dir = args.length > 0 ? args[0] : "/work/fixtures";
        PrintStream stdout = System.out;

        stdout.println("# RejectionDump: rows the reference refuses, and how");
        for (String[] one : cases(dir)) {
            String name = one[0];
            String[] argv = new String[one.length - 1];
            System.arraycopy(one, 1, argv, 0, argv.length);

            // Picard writes progress to stdout and stderr; both are captured so the dump carries
            // only the rows, and a tool that starts logging differently does not become a
            // divergence in a file that is about exit codes and messages.
            PrintStream quiet = new PrintStream(new ByteArrayOutputStream());
            PrintStream realErr = System.err;
            String outcome;
            try {
                System.setOut(quiet);
                System.setErr(quiet);
                outcome = "EXIT=" + new Runner().run(argv);
            } catch (Throwable t) {
                Throwable root = t;
                while (root.getCause() != null) root = root.getCause();
                outcome = "THROWN=" + root.getClass().getName() + ": " + root.getMessage();
            } finally {
                System.setOut(stdout);
                System.setErr(realErr);
            }
            // The fixture directory is a per-run temporary path, so it is replaced by a token: the
            // message is the behaviour, the path is not.
            outcome = outcome.replace(dir, "<FIXTURES>");
            stdout.printf("reject\t%s\t%s%n", name, outcome.replace("\n", "\\n"));
            new File(OUT).delete();
        }
    }
}
