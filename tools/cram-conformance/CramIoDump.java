/*
 * Oracle dump harness for CRAM input and output in picard-rs.
 *
 * Picard 3.5.0 bundles htsjdk 5.0.0, which writes CRAM 3.1 (rANS Nx16, fqzcomp, the name
 * tokeniser, and GZIP blocks deflated by libdeflate) and regenerates MD and NM when it reads one.
 * This builds a reference and a coordinate-sorted SAM with every read feature, tag type and pairing
 * the CRAM codec distinguishes, then runs Picard on them under /work/cram, a fixed directory,
 * because the file definition of every CRAM carries its output's URI:
 *
 *   sam            the input SAM
 *   cram           SamFormatConverter SAM -> CRAM, the bytes
 *   cram-sorted    SortSam CRAM -> CRAM with CREATE_INDEX, the bytes
 *   bai            its index, the bytes
 *   sam-from-cram  SamFormatConverter CRAM -> SAM, the text
 *   metrics        CollectAlignmentSummaryMetrics on the CRAM, without its `# ` lines
 *   shards         SplitSamByNumberOfReads on the CRAM: the names it writes
 *   error          each refusal, as its class and message: reading without a reference, writing
 *                  without one, reading with the wrong one, and indexing a CRAM that is not
 *                  coordinate-sorted. The wrong reference is tried twice: on the file above,
 *                  whose one slice spans both contigs and so is never checked against an MD5,
 *                  and on a file of chr1 reads alone, whose slice is
 *
 * The suite runs without -Dsamjdk.use_libdeflate=false: Picard's command line deflates CRAM GZIP
 * blocks with libdeflate, and that is what the port writes.
 *
 *   java -cp picard.jar:. CramIoDump
 */
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.Arrays;

public class CramIoDump {
    static final File DIR = new File("/work/cram");

    public static void main(final String[] args) throws Exception {
        DIR.mkdirs();
        final String chr1 = contig(5000, 11L);
        final String chr2 = contig(3000, 23L);
        writeFasta(new File(DIR, "ref.fasta"), chr1, chr2);
        // The same dictionary over other bases: the slice MD5s no longer match.
        writeFasta(new File(DIR, "wrong.fasta"), contig(5000, 12L), contig(3000, 24L));

        final String sam = buildSam(chr1, chr2);
        final File in = new File(DIR, "in.sam");
        try (PrintStream ps = new PrintStream(in)) { ps.print(sam); }
        emit("sam", "input", esc(sam));

        final String ref = "REFERENCE_SEQUENCE=" + new File(DIR, "ref.fasta").getPath();
        run(new picard.sam.SamFormatConverter(), "INPUT=" + in, "OUTPUT=/work/cram/out.cram", ref);
        emit("cram", "sam-to-cram", hex(new File(DIR, "out.cram")));

        run(new picard.sam.SortSam(), "INPUT=/work/cram/out.cram", "OUTPUT=/work/cram/sorted.cram",
                "SORT_ORDER=coordinate", "CREATE_INDEX=true", ref);
        emit("cram-sorted", "sortsam", hex(new File(DIR, "sorted.cram")));
        emit("bai", "sortsam", hex(new File(DIR, "sorted.cram.bai")));

        run(new picard.sam.SamFormatConverter(), "INPUT=/work/cram/out.cram", "OUTPUT=/work/cram/back.sam", ref);
        emit("sam-from-cram", "converter", esc(new String(Files.readAllBytes(new File(DIR, "back.sam").toPath()))));

        run(new picard.analysis.CollectAlignmentSummaryMetrics(), "INPUT=/work/cram/out.cram",
                "OUTPUT=/work/cram/summary.txt", ref);
        final StringBuilder metrics = new StringBuilder();
        for (final String line : Files.readAllLines(new File(DIR, "summary.txt").toPath())) {
            if (!line.startsWith("# ")) metrics.append(line).append('\n');
        }
        emit("metrics", "alignment-summary", esc(metrics.toString()));

        final File shards = new File(DIR, "shards");
        shards.mkdirs();
        run(new picard.sam.SplitSamByNumberOfReads(), "INPUT=/work/cram/out.cram",
                "OUTPUT=" + shards, "SPLIT_TO_N_FILES=3", ref);
        final String[] names = shards.list();
        Arrays.sort(names);
        emit("shards", "split-by-number", String.join(",", names));

        error("read-without-reference", new picard.sam.SamFormatConverter(),
                "INPUT=/work/cram/out.cram", "OUTPUT=/work/cram/noref.sam");
        error("write-without-reference", new picard.sam.SamFormatConverter(),
                "INPUT=" + in, "OUTPUT=/work/cram/noref.cram");
        final StringBuilder chr1Only = new StringBuilder();
        for (final String line : sam.split("\n")) {
            if (line.startsWith("@") || line.split("\t")[2].equals("chr1")) chr1Only.append(line).append('\n');
        }
        final File chr1In = new File(DIR, "chr1.sam");
        try (PrintStream ps = new PrintStream(chr1In)) { ps.print(chr1Only); }
        run(new picard.sam.SamFormatConverter(), "INPUT=" + chr1In, "OUTPUT=/work/cram/chr1.cram", ref);
        error("wrong-reference-single", new picard.sam.SamFormatConverter(),
                "INPUT=/work/cram/chr1.cram", "OUTPUT=/work/cram/wrong1.sam",
                "REFERENCE_SEQUENCE=" + new File(DIR, "wrong.fasta").getPath());
        error("wrong-reference", new picard.sam.SamFormatConverter(),
                "INPUT=/work/cram/out.cram", "OUTPUT=/work/cram/wrong.sam",
                "REFERENCE_SEQUENCE=" + new File(DIR, "wrong.fasta").getPath());
        error("index-not-coordinate", new picard.sam.SortSam(), "INPUT=/work/cram/out.cram",
                "OUTPUT=/work/cram/byname.cram", "SORT_ORDER=queryname", "CREATE_INDEX=true", ref);
    }

    /** A deterministic contig: a linear congruential walk over ACGT. */
    static String contig(final int length, final long seed) {
        final StringBuilder sb = new StringBuilder(length);
        long x = seed;
        for (int i = 0; i < length; i++) {
            x = x * 6364136223846793005L + 1442695040888963407L;
            sb.append("ACGT".charAt((int) ((x >>> 33) % 4)));
        }
        return sb.toString();
    }

    static void writeFasta(final File fasta, final String chr1, final String chr2) throws Exception {
        final StringBuilder text = new StringBuilder();
        final StringBuilder fai = new StringBuilder();
        for (final String[] c : new String[][] {{"chr1", chr1}, {"chr2", chr2}}) {
            text.append('>').append(c[0]).append('\n');
            final long offset = text.length();
            for (int i = 0; i < c[1].length(); i += 60) {
                text.append(c[1], i, Math.min(c[1].length(), i + 60)).append('\n');
            }
            fai.append(c[0]).append('\t').append(c[1].length()).append('\t').append(offset)
                    .append("\t60\t61\n");
        }
        Files.write(fasta.toPath(), text.toString().getBytes());
        Files.write(new File(fasta.getPath() + ".fai").toPath(), fai.toString().getBytes());
    }

    static String read(final String contig, final int start, final int length) {
        return contig.substring(start - 1, start - 1 + length);
    }

    static String quals(final int length, final int seed) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < length; i++) sb.append((char) ('!' + 2 + (i * 7 + seed) % 39));
        return sb.toString();
    }

    static String buildSam(final String chr1, final String chr2) {
        final StringBuilder s = new StringBuilder();
        s.append("@HD\tVN:1.6\tSO:coordinate\n");
        s.append("@SQ\tSN:chr1\tLN:").append(chr1.length()).append('\n');
        s.append("@SQ\tSN:chr2\tLN:").append(chr2.length()).append('\n');
        s.append("@RG\tID:rg1\tSM:s1\tLB:lib1\tPL:ILLUMINA\n");
        s.append("@RG\tID:rg2\tSM:s2\tLB:lib2\tPL:ILLUMINA\n");
        s.append("@CO\tcram conformance\n");
        // An exact match, with tags of every type the dictionary keys on.
        line(s, "exact", 0, "chr1", 101, 60, "50M", "*", 0, 0, read(chr1, 101, 50), quals(50, 1),
                "RG:Z:rg1\tXa:A:q\tXc:i:-5\tXC:i:200\tXs:i:-3000\tXS:i:40000\tXi:i:-100000\tXI:i:3000000000"
                        + "\tXf:f:0.25\tXZ:Z:text\tXH:H:1AE301\tXB:B:c,-1,2\tXb:B:C,1,250\tXt:B:s,-300,4"
                        + "\tXT:B:S,60000\tXn:B:i,-70000,1\tXN:B:I,4000000\tXF:B:f,0.5,1.5");
        // Substitutions, an N and an IUPAC base, with NM and MD as a reader would compute them,
        // which the writer drops, and a wrong NM, which it keeps.
        final StringBuilder sub = new StringBuilder(read(chr1, 201, 50));
        sub.setCharAt(5, sub.charAt(5) == 'A' ? 'C' : 'A');
        sub.setCharAt(17, 'N');
        sub.setCharAt(30, 'R');
        line(s, "substitutions", 0, "chr1", 201, 33, "50M", "*", 0, 0, sub.toString(), quals(50, 2),
                "RG:Z:rg2\tNM:i:7");
        // An insertion, a deletion, soft and hard clips, a skip, and no qualities.
        final String indel = read(chr1, 301, 20) + "TTT" + read(chr1, 324, 20);
        line(s, "indels", 16, "chr1", 301, 12, "20M3I3D20M", "*", 0, 0, indel, "*", "RG:Z:rg1");
        line(s, "clipped", 0, "chr1", 401, 50, "2H5S40M5S", "*", 0, 0,
                "GGGGG" + read(chr1, 401, 40) + "CCCCC", quals(50, 3), "RG:Z:rg1");
        line(s, "spliced", 0, "chr1", 501, 50, "25M200N25M", "*", 0, 0,
                read(chr1, 501, 25) + read(chr1, 726, 25), quals(50, 4), "RG:Z:rg1");
        // A pair whose template length is htslib's, which the writer attaches, and one whose is
        // not, which it leaves detached; a secondary and a supplementary; a duplicate.
        line(s, "pair-a", 99, "chr1", 1001, 60, "50M", "=", 1101, 150, read(chr1, 1001, 50), quals(50, 5), "RG:Z:rg1");
        line(s, "pair-b", 99, "chr1", 1051, 60, "50M", "=", 1151, 149, read(chr1, 1051, 50), quals(50, 6), "RG:Z:rg2");
        line(s, "pair-a", 147, "chr1", 1101, 60, "50M", "=", 1001, -150, read(chr1, 1101, 50), quals(50, 7), "RG:Z:rg1");
        line(s, "pair-b", 147, "chr1", 1151, 60, "50M", "=", 1051, -149, read(chr1, 1151, 50), quals(50, 8), "RG:Z:rg2");
        line(s, "secondary", 256, "chr1", 2001, 0, "50M", "*", 0, 0, read(chr1, 2001, 50), quals(50, 9), "RG:Z:rg1");
        line(s, "supplementary", 2048, "chr1", 2101, 3, "30M20H", "*", 0, 0, read(chr1, 2101, 30), quals(30, 10), "RG:Z:rg1");
        line(s, "duplicate", 1024, "chr1", 2201, 60, "50M", "*", 0, 0, read(chr1, 2201, 50), quals(50, 11), "RG:Z:rg2");
        // A mate on another contig, then chr2, a placed unmapped read with its mapped mate, and
        // the unplaced unmapped reads last, one with no sequence.
        line(s, "far", 65, "chr1", 4901, 40, "50M", "chr2", 101, 0, read(chr1, 4901, 50), quals(50, 12), "RG:Z:rg1");
        line(s, "far", 129, "chr2", 101, 40, "50M", "chr1", 4901, 0, read(chr2, 101, 50), quals(50, 13), "RG:Z:rg1");
        line(s, "placed", 73, "chr2", 501, 60, "50M", "=", 501, 0, read(chr2, 501, 50), quals(50, 14), "RG:Z:rg2");
        line(s, "placed", 133, "chr2", 501, 0, "*", "=", 501, 0, read(chr2, 520, 50), quals(50, 15), "RG:Z:rg2");
        line(s, "unplaced", 4, "*", 0, 0, "*", "*", 0, 0, "ACGTNACGTN", quals(10, 16), "RG:Z:rg1");
        line(s, "noseq", 4, "*", 0, 0, "*", "*", 0, 0, "*", "*", "RG:Z:rg2");
        return s.toString();
    }

    static void line(final StringBuilder s, final String name, final int flag, final String rname,
                     final int pos, final int mapq, final String cigar, final String rnext,
                     final int pnext, final int tlen, final String seq, final String qual,
                     final String tags) {
        s.append(name).append('\t').append(flag).append('\t').append(rname).append('\t').append(pos)
                .append('\t').append(mapq).append('\t').append(cigar).append('\t').append(rnext)
                .append('\t').append(pnext).append('\t').append(tlen).append('\t').append(seq)
                .append('\t').append(qual).append('\t').append(tags).append('\n');
    }

    static void run(final picard.cmdline.CommandLineProgram program, final String... args) {
        final int rc = program.instanceMain(args);
        if (rc != 0) {
            System.err.println(program.getClass().getSimpleName() + " exited " + rc);
            System.exit(rc);
        }
    }

    /** The exception a refusal escapes `instanceMain` with, as an uncaught one prints it. */
    static void error(final String kase, final picard.cmdline.CommandLineProgram program, final String... args) {
        try {
            program.instanceMain(args);
            emit("error", kase, "no exception");
        } catch (final Throwable t) {
            emit("error", kase, esc(t.getClass().getName() + ": " + t.getMessage()));
        }
    }

    static String hex(final File f) throws Exception {
        final byte[] b = Files.readAllBytes(f.toPath());
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
