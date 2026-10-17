/*
 * Oracle dump for the input shapes the tools' own suites never carried, and that the port refused
 * or answered differently until they were measured here: secondary and supplementary records, mate
 * cigars (MC), missing and duplicated mates, unsorted input to a presorted writer, malformed tags,
 * and a run that throws after the writer has flushed part of its output.
 *
 * Each case runs one tool in-process (instanceMain, so nothing exits) on a SAM and emits the input,
 * the exit code or the exception, and what the output file holds afterwards.
 *
 *   java -cp picard.jar:. InputShapesDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class InputShapesDump {
    static final String H2 = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n@SQ\tSN:chr2\tLN:10000\n";
    static final String H150 = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:150\n@SQ\tSN:chr2\tLN:1000\n";

    public static void main(final String[] args) throws Exception {
        revertOriginal();
        cleanSam();
        samErrorMetrics();
    }

    static final String[] STRATIFIERS = {"ALL", "GC_CONTENT", "READ_ORDINALITY", "READ_BASE", "READ_DIRECTION",
        "PAIR_ORIENTATION", "REFERENCE_BASE", "PRE_DINUC", "POST_DINUC", "HOMOPOLYMER_LENGTH", "HOMOPOLYMER",
        "BINNED_HOMOPOLYMER", "FLOWCELL_TILE", "FLOWCELL_Y", "FLOWCELL_X", "READ_GROUP", "CYCLE", "BINNED_CYCLE",
        "SOFT_CLIPS", "INSERT_LENGTH", "BASE_QUALITY", "MAPPING_QUALITY", "MISMATCHES_IN_READ",
        "ONE_BASE_PADDED_CONTEXT", "TWO_BASE_PADDED_CONTEXT", "CONSENSUS", "NS_IN_READ", "INSERTIONS_IN_READ",
        "DELETIONS_IN_READ", "INDELS_IN_READ", "INDEL_LENGTH"};

    /**
     * CollectSamErrorMetrics on every stratifier, its default list (which names seven that were not
     * ported), and the two exceptions a stratifier or the loop throws: PAIR_PROPERNESS on a read
     * that is not paired, and a read stored without base qualities.
     */
    static void samErrorMetrics() throws Exception {
        long seed = 7;
        final StringBuilder ref = new StringBuilder();
        for (int i = 0; i < 3000; i++) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            ref.append(i % 400 < 6 ? 'N' : (i % 97 < 8 ? 'A' : "ACGT".charAt((int) ((seed >>> 33) & 3))));
        }
        final File fasta = File.createTempFile("shapes-ref-", ".fasta");
        fasta.deleteOnExit();
        final StringBuilder fa = new StringBuilder(">chr1\n");
        for (int i = 0; i < ref.length(); i += 60) fa.append(ref, i, Math.min(ref.length(), i + 60)).append('\n');
        Files.write(fasta.toPath(), fa.toString().getBytes());
        Files.write(new File(fasta.getPath() + ".fai").toPath(), ("chr1\t3000\t6\t60\t61\n").getBytes());
        final File dict = new File(fasta.getPath().replaceAll("\\.fasta$", ".dict"));
        Files.write(dict.toPath(), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:3000\n".getBytes());
        dict.deleteOnExit();
        new File(fasta.getPath() + ".fai").deleteOnExit();

        final String header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:3000\n@RG\tID:rg1\tSM:s\n@RG\tID:rg2\tSM:s\n";
        final java.util.List<String[]> rows = new java.util.ArrayList<>();
        for (int i = 0; i < 400; i++) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            final int start = 1 + (int) ((seed >>> 20) % 2800);
            final int kind = i % 6;
            final String cigar = kind == 0 ? "60M" : kind == 1 ? "20M2I38M" : kind == 2 ? "25M3D35M" : kind == 3 ? "5S55M" : kind == 4 ? "55M5S" : "60M";
            final StringBuilder bases = new StringBuilder(), quals = new StringBuilder();
            int r = start - 1;
            // Read bases follow the reference through the cigar, with a mismatch every 17th base.
            final java.util.regex.Matcher m = java.util.regex.Pattern.compile("(\\d+)([MIDS])").matcher(cigar);
            int n = 0;
            while (m.find()) {
                final int len = Integer.parseInt(m.group(1));
                final char op = m.group(2).charAt(0);
                for (int k = 0; k < len; k++) {
                    if (op == 'D') { r++; continue; }
                    char b = op == 'M' && r < ref.length() ? ref.charAt(r) : "ACGT".charAt((n + i) & 3);
                    if (++n % 17 == 0) b = b == 'A' ? 'C' : 'A';
                    bases.append(b);
                    seed = seed * 6364136223846793005L + 1442695040888963407L;
                    quals.append((char) (33 + 10 + (int) ((seed >>> 40) % 31)));
                    if (op == 'M') r++;
                }
            }
            final int flag = (i % 4 == 0) ? (i % 8 == 0 ? 0 : 16) : (1 | 2 | (i % 2 == 0 ? 64 : 128) | (i % 3 == 0 ? 16 : 32));
            final boolean paired = (flag & 1) != 0;
            final String name = "M01:7:FC1:1:" + (1101 + i % 3) + ":" + ((i * 7919) % 30000) + ":" + ((i * 104729) % 30000);
            String tags = "\tRG:Z:rg" + (1 + i % 2);
            if (i % 5 == 0) tags += "\taD:i:" + (i % 4) + "\tbD:i:" + (i % 3);
            if (i % 3 == 0) tags += "\tcD:i:" + (i % 5);
            if (i % 7 == 0) tags += "\tSA:Z:chr1,100,+,30M30S,60,0;";
            if (i % 2 == 0) tags += "\tNM:i:" + (i % 9);
            rows.add(new String[]{String.format("%010d", start), name + "\t" + flag + "\tchr1\t" + start + "\t" + (20 + i % 41) + "\t" + cigar
                + "\t" + (paired ? "=\t" + (start + 100) + "\t" + (i % 2 == 0 ? 160 : -160) : "*\t0\t0") + "\t" + bases + "\t" + quals + tags + "\n"});
        }
        rows.sort((a, b) -> a[0].compareTo(b[0]));
        final StringBuilder sam = new StringBuilder(header), pairedSam = new StringBuilder(header);
        for (final String[] row : rows) {
            sam.append(row[1]);
            if ((Integer.parseInt(row[1].split("\t")[1]) & 1) != 0) pairedSam.append(row[1]);
        }
        final String noQualities = header + "q1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n";

        final String tool = "SamErrorMetric.CollectSamErrorMetrics";
        final String r = "REFERENCE_SEQUENCE=" + fasta.getAbsolutePath();
        runMetrics(tool, "sem_defaults", sam.toString(), r);
        // LONG_HOMOPOLYMER is not varied here: it is read into a static Lazy the first time a
        // BINNED_HOMOPOLYMER stratifier is built, so in one JVM the first case would fix it for all.
        final String[] every = new String[STRATIFIERS.length + 5];
        every[0] = r;
        every[1] = "ERROR_METRICS=null";
        for (int i = 0; i < STRATIFIERS.length; i++) every[2 + i] = "ERROR_METRICS=ERROR:" + STRATIFIERS[i];
        every[2 + STRATIFIERS.length] = "ERROR_METRICS=INDEL_ERROR:INDEL_LENGTH";
        every[3 + STRATIFIERS.length] = "ERROR_METRICS=OVERLAPPING_ERROR:PAIR_ORIENTATION:BINNED_CYCLE";
        every[4 + STRATIFIERS.length] = "ERROR_METRICS=ERROR:READ_ORDINALITY:HOMOPOLYMER:CONSENSUS";
        runMetrics(tool, "sem_every_stratifier", sam.toString(), every);
        runMetrics(tool, "sem_pair_properness", pairedSam.toString(), r, "ERROR_METRICS=null",
            "ERROR_METRICS=ERROR:PAIR_PROPERNESS", "ERROR_METRICS=OVERLAPPING_ERROR:PAIR_PROPERNESS:READ_ORDINALITY");
        runMetrics(tool, "sem_pair_properness_unpaired", sam.toString(), r, "ERROR_METRICS=null", "ERROR_METRICS=ERROR:PAIR_PROPERNESS");
        runMetrics(tool, "sem_no_qualities", noQualities, r, "ERROR_METRICS=null", "ERROR_METRICS=ERROR:CYCLE");
    }

    /** A tool that writes one file per suffix next to OUTPUT: every file it leaves is a row. */
    static void runMetrics(final String tool, final String kase, final String sam, final String... extra) throws Exception {
        final File input = File.createTempFile("shapes-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam); }
        final File dir = Files.createTempDirectory("shapes-out-").toFile();
        final String prefix = new File(dir, "out").getAbsolutePath();
        final String[] argv = new String[2 + extra.length];
        argv[0] = "INPUT=" + input.getAbsolutePath();
        argv[1] = "OUTPUT=" + prefix;
        System.arraycopy(extra, 0, argv, 2, extra.length);
        emit("tool", kase, tool.substring(tool.lastIndexOf('.') + 1));
        final java.util.List<String> shown = new java.util.ArrayList<>();
        for (final String a : extra) shown.add(a.startsWith("REFERENCE_SEQUENCE=") ? "REFERENCE_SEQUENCE=<ref>" : a);
        emit("args", kase, String.join(" ", shown));
        emit("input", kase, sam);
        for (final String a : extra) {
            if (a.startsWith("REFERENCE_SEQUENCE=")) {
                emit("reference", kase, new String(Files.readAllBytes(new File(a.substring(19)).toPath())));
            }
        }
        try {
            final Object program = Class.forName("picard.sam." + tool).getDeclaredConstructor().newInstance();
            emit("exit", kase, Integer.toString(((picard.cmdline.CommandLineProgram) program).instanceMain(argv)));
        } catch (final Throwable t) {
            emit("thrown", kase, t.toString());
        }
        final File[] files = dir.listFiles();
        java.util.Arrays.sort(files);
        for (final File f : files) {
            // The metrics header carries the command line and a timestamp: only the table is compared.
            final StringBuilder body = new StringBuilder();
            for (final String line : Files.readAllLines(f.toPath())) if (!line.startsWith("# ")) body.append(line).append('\n');
            emit("file:" + f.getName().substring(4), kase, body.toString());
            f.delete();
        }
        dir.delete();
    }

    static void revertOriginal() throws Exception {
        final String tool = "RevertOriginalBaseQualitiesAndAddMateCigar";
        final String supplementary = H2
            + "p1\t97\tchr1\t100\t60\t50M\t=\t300\t250\t*\t*\tOQ:Z:5\n"
            + "p1\t2145\tchr1\t150\t20\t20M30S\tchr1\t300\t0\t*\t*\tSA:Z:chr1,100,+,50M,60,0;\n"
            + "p1\t353\tchr2\t500\t0\t50M\tchr1\t300\t0\t*\t*\n"
            + "p1\t145\tchr1\t300\t50\t40M10S\t=\t100\t-250\t*\t*\n"
            + "p1\t2193\tchr2\t900\t10\t10S40M\tchr1\t300\t0\t*\t*\tMC:Z:3M\tMQ:i:7\n"
            + "p2\t73\tchr1\t1000\t60\t50M\t=\t1000\t0\t*\t*\n"
            + "p2\t133\tchr1\t1000\t0\t*\t=\t1000\t0\t*\t*\n"
            + "p2\t2121\tchr1\t2000\t30\t25M25S\t*\t0\t0\t*\t*\n"
            + "u1\t0\tchr1\t3000\t60\t50M\t*\t0\t0\t*\t*\n"
            + "u1\t2048\tchr2\t3000\t60\t25M25S\t*\t0\t0\t*\t*\n";
        final String missingMate = H2
            + "m1\t99\tchr1\t100\t60\t50M\t=\t300\t250\t*\t*\n"
            + "m2\t0\tchr1\t200\t60\t50M\t*\t0\t0\t*\t*\n";
        final String twoFirst = H2
            + "t1\t99\tchr1\t100\t60\t50M\t=\t300\t250\t*\t*\n"
            + "t1\t99\tchr1\t150\t60\t50M\t=\t300\t200\t*\t*\n"
            + "t1\t147\tchr1\t300\t60\t50M\t=\t100\t-250\t*\t*\n";
        run(tool, "revert_supplementary", supplementary, ".sam", "VALIDATION_STRINGENCY=SILENT");
        run(tool, "revert_supplementary_queryname", supplementary, ".sam", "VALIDATION_STRINGENCY=SILENT", "SORT_ORDER=queryname");
        run(tool, "revert_supplementary_unsorted", supplementary, ".sam", "VALIDATION_STRINGENCY=SILENT", "SORT_ORDER=unsorted");
        run(tool, "revert_missing_mate", missingMate, ".sam", "VALIDATION_STRINGENCY=SILENT");
        run(tool, "revert_missing_mate_bam", missingMate, ".bam", "VALIDATION_STRINGENCY=SILENT");
        run(tool, "revert_two_first", twoFirst, ".sam", "VALIDATION_STRINGENCY=SILENT");
    }

    static void cleanSam() throws Exception {
        final String tool = "CleanSam";
        final String mateCigars = H150
            + "d\t0\tchr1\t10\t60\t10M\t*\t0\t0\t*\t*\tMC:Z:300M\n"
            + "e\t65\tchr1\t20\t60\t10M\tchr1\t145\t0\t*\t*\tMC:Z:2S20M3I2D4S\n"
            + "f\t65\tchr1\t20\t60\t10M\tchr1\t145\t0\t*\t*\tMC:Z:*\n"
            + "a\t99\tchr1\t60\t60\t36M\t=\t130\t106\t*\t*\tMC:Z:36M\n"
            + "b\t97\tchr1\t100\t60\t20M\tchr2\t900\t0\t*\t*\tMC:Z:10S100M5S\n"
            + "a\t147\tchr1\t130\t60\t36M\t=\t60\t-106\t*\t*\tMC:Z:36M\n"
            + "c\t73\tchr1\t140\t60\t10M5S\t=\t140\t0\t*\t*\tMC:Z:30M\n"
            + "g\t137\tchr1\t140\t60\t10M\t=\t140\t0\t*\t*\tMC:Z:30M\n"
            + "b\t145\tchr2\t900\t60\t10S100M5S\tchr1\t100\t0\t*\t*\tMC:Z:20M\n";
        run(tool, "clean_mate_cigars", mateCigars, ".sam");
        run(tool, "clean_unsorted", H150
            + "x\t0\tchr1\t50\t60\t10M\t*\t0\t0\t*\t*\n"
            + "y\t0\tchr1\t20\t60\t10M\t*\t0\t0\t*\t*\n", ".sam");
        run(tool, "clean_mc_not_a_string", H150
            + "g\t65\tchr1\t20\t60\t10M\tchr1\t100\t0\t*\t*\tMC:i:5\n", ".sam");
        run(tool, "clean_mc_unknown_operator", H150
            + "h0\t0\tchr1\t10\t60\t10M\t*\t0\t0\t*\t*\n"
            + "g\t65\tchr1\t20\t60\t10M\tchr1\t100\t0\t*\t*\tMC:Z:5Q\n", ".sam");
        run(tool, "clean_mc_no_operator", H150
            + "g\t65\tchr1\t20\t60\t10M\tchr1\t100\t0\t*\t*\tMC:Z:5\n", ".sam");
        run(tool, "clean_mc_no_length", H150
            + "g\t65\tchr1\t20\t60\t10M\tchr1\t100\t0\t*\t*\tMC:Z:M5\n", ".sam");
        run(tool, "clean_no_reference", H150
            + "h0\t0\tchr1\t10\t60\t10M\t*\t0\t0\t*\t*\n"
            + "n\t0\t*\t20\t60\t10M\t*\t0\t0\t*\t*\n", ".sam");
        // Enough record text before the throw to fill the SAM writer's 128 KiB buffer once.
        final StringBuilder big = new StringBuilder("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n");
        long seed = 1;
        for (int i = 0; i < 1000; i++) {
            final StringBuilder seq = new StringBuilder(), qual = new StringBuilder();
            for (int j = 0; j < 100; j++) {
                seed = (seed * 6364136223846793005L + 1442695040888963407L);
                seq.append("ACGT".charAt((int) ((seed >>> 33) & 3)));
                qual.append((char) (35 + (int) ((seed >>> 40) % 38)));
            }
            final int pos = 100 + i * 10;
            big.append("r").append(i).append("\t65\tchr1\t").append(pos).append("\t60\t100M\tchr1\t").append(pos + 50)
               .append("\t0\t").append(seq).append('\t').append(qual).append(i == 900 ? "\tMC:i:5\n" : "\tMC:Z:100M\n");
        }
        run(tool, "clean_throws_after_a_flush", big.toString(), ".sam");
    }

    static void run(final String tool, final String kase, final String sam, final String ext, final String... extra) throws Exception {
        final File input = File.createTempFile("shapes-in-", ".sam");
        input.deleteOnExit();
        try (PrintStream ps = new PrintStream(input)) { ps.print(sam); }
        final File out = File.createTempFile("shapes-out-", ext);
        out.delete();
        out.deleteOnExit();
        final String[] argv = new String[2 + extra.length];
        argv[0] = "INPUT=" + input.getAbsolutePath();
        argv[1] = "OUTPUT=" + out.getAbsolutePath();
        System.arraycopy(extra, 0, argv, 2, extra.length);
        emit("tool", kase, tool);
        emit("args", kase, String.join(" ", extra));
        emit("input", kase, sam);
        try {
            final Object program = Class.forName("picard.sam." + tool).getDeclaredConstructor().newInstance();
            final int rc = ((picard.cmdline.CommandLineProgram) program).instanceMain(argv);
            emit("exit", kase, Integer.toString(rc));
        } catch (final Throwable t) {
            emit("thrown", kase, t.toString());
        }
        if (!out.exists()) {
            emit("output_absent", kase, "");
        } else if (ext.equals(".sam")) {
            emit("output", kase, new String(Files.readAllBytes(out.toPath())));
        } else {
            emit("output_bytes", kase, Long.toString(out.length()));
        }
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
