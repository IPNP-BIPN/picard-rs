/*
 * Oracle dump harness for SortSam spilling to TMP_DIR in picard-rs.
 *
 * SortSam sorts through htsjdk's SortingCollection, which holds MAX_RECORDS_IN_RAM records, spills
 * each full batch to TMP_DIR as a sorted run and merges the runs. This writes an unsorted SAM of
 * 600 records, most of which the comparators call equal, and sorts it with MAX_RECORDS_IN_RAM=50
 * (twelve runs): coordinate to an indexed BAM, queryname to a BAM and to a SAM, all with
 * USE_JDK_DEFLATER=true. It refuses to emit anything unless the same sorts without spilling write
 * the same bytes, and unless the small threshold really spills (a TMP_DIR that is a file fails).
 *
 *   java -cp picard-fat.jar:. SortSamSpillDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.Arrays;

public class SortSamSpillDump {
    public static void main(final String[] args) throws Exception {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:unsorted\n");
        sam.append("@SQ\tSN:chr1\tLN:100000\n");
        sam.append("@SQ\tSN:chr2\tLN:100000\n");
        sam.append("@RG\tID:rg1\tSM:s\n");
        final int[] flags = {99, 147, 0, 16};
        for (int i = 0; i < 600; i++) {
            final String name = "q" + (i * 7 % 13);
            final String tags = "RG:Z:rg1\tXI:i:" + i + "\tXF:f:" + (i % 3) + ".5\tXB:B:s," + (i % 7) + ",-" + (i % 5);
            if (i % 50 == 49) {
                sam.append(name).append("\t4\t*\t0\t0\t*\t*\t0\t0\tACGTAC\tIIIIII\t").append(tags).append('\n');
                continue;
            }
            final int flag = flags[i % 4];
            final String chr = "chr" + (1 + i % 2);
            final int pos = 100 + (i * 31 % 5) * 50;
            final boolean paired = (flag & 1) != 0;
            sam.append(name).append('\t').append(flag).append('\t').append(chr).append('\t').append(pos)
                .append("\t60\t6M\t").append(paired ? "=" : "*").append('\t').append(paired ? pos + 100 : 0)
                .append('\t').append(paired ? ((flag & 16) == 0 ? 106 : -106) : 0)
                .append("\tACGTAC\tIIIIII\t").append(tags).append('\n');
        }

        final File dir = Files.createTempDirectory("ss-spill-").toFile();
        final File in = new File(dir, "in.sam");
        try (PrintStream ps = new PrintStream(in)) { ps.print(sam.toString()); }
        final File tmp = new File(dir, "tmp");
        tmp.mkdir();

        emit("sam", "input", esc(sam.toString()));
        final String[][] runs = {{"coordinate", "bam"}, {"queryname", "bam"}, {"queryname", "sam"}};
        for (final String[] run : runs) {
            final String kase = run[0] + "-" + run[1];
            final File spilled = sort(in, new File(dir, kase + ".spilled." + run[1]), run[0], "50", tmp);
            final File whole = sort(in, new File(dir, kase + ".whole." + run[1]), run[0], "500000", tmp);
            if (!Arrays.equals(Files.readAllBytes(spilled.toPath()), Files.readAllBytes(whole.toPath()))) {
                System.err.println(kase + ": the spilled sort differs from the in-memory one");
                System.exit(2);
            }
            final byte[] out = Files.readAllBytes(spilled.toPath());
            emit(run[1], kase, run[1].equals("sam") ? esc(new String(out, "UTF-8")) : hex(out));
            if (run[0].equals("coordinate")) {
                final File bai = new File(spilled.getPath().replaceAll("\\.bam$", ".bai"));
                emit("bai", kase, hex(Files.readAllBytes(bai.toPath())));
            }
        }

        // A TMP_DIR that is a file cannot hold a run, so this fails only if the sort spills.
        final File notADirectory = new File(dir, "not-a-directory");
        notADirectory.createNewFile();
        boolean spilledToDisk;
        try {
            spilledToDisk = new picard.sam.SortSam().instanceMain(command(in, new File(dir, "x.bam"),
                    "coordinate", "50", notADirectory)) != 0;
        } catch (final RuntimeException e) {
            spilledToDisk = true;
        }
        if (!spilledToDisk) {
            System.err.println("MAX_RECORDS_IN_RAM=50 did not spill");
            System.exit(3);
        }
    }

    static String[] command(final File in, final File out, final String order, final String max, final File tmp) {
        return new String[] {
                "INPUT=" + in.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
                "SORT_ORDER=" + order,
                "MAX_RECORDS_IN_RAM=" + max,
                "TMP_DIR=" + tmp.getAbsolutePath(),
                "CREATE_INDEX=" + order.equals("coordinate"),
                "USE_JDK_DEFLATER=true",
                "VALIDATION_STRINGENCY=SILENT",
        };
    }

    static File sort(final File in, final File out, final String order, final String max, final File tmp) {
        final int rc = new picard.sam.SortSam().instanceMain(command(in, out, order, max, tmp));
        if (rc != 0) { System.err.println("SortSam exited " + rc); System.exit(rc); }
        return out;
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
