/*
 * What a VCF writer leaves on disk when it refuses a record partway through, taken from the
 * reference.
 *
 * VariantContextWriterBuilder opens the output as IOUtil.maybeBufferOutputStream(
 * Files.newOutputStream(path), Defaults.BUFFER_SIZE), a 128 KiB BufferedOutputStream, and neither
 * BCF2Writer nor VCFWriter flushes it before close(). VcfFormatConverter never reaches close() when
 * a record throws, so the file holds exactly what the buffer had flushed: nothing on a small file,
 * whole buffers' worth on a large one, cut where a write call did not fit. No .idx is written.
 *
 * Each case is a VCF generated from a recipe, so the port builds the same bytes without the
 * golden carrying them: three samples, `n` records on chr1, `pad` extra contig lines to make the
 * header as large as the buffer, and record `bad` replaced by one the writer refuses. Each is
 * converted to .bcf and to .vcf without an index, and the undeclared case with one as well.
 *
 * Output:
 *
 *     input\t<case>\t<length> <md5>
 *     <step>\t<case>\t<length on disk> <md5 on disk> idx=<yes|no> <exception, or ok>
 *
 * Usage: PartialOutputDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;

public class PartialOutputDump {

    static final StringBuilder buf = new StringBuilder();
    static Path dir;

    static void emit(final String kind, final String key, final String value) {
        buf.append(kind).append('\t').append(key).append('\t').append(value).append('\n');
    }

    static String digest(final byte[] bytes) throws Exception {
        final StringBuilder b = new StringBuilder();
        for (final byte x : MessageDigest.getInstance("MD5").digest(bytes)) {
            b.append(Character.forDigit((x >> 4) & 0xf, 16)).append(Character.forDigit(x & 0xf, 16));
        }
        return bytes.length + " " + b;
    }

    /** The recipe; `PartialOutput` in the port's test builds the same text. */
    static String vcf(final int n, final int bad, final String kind, final int pad) {
        final StringBuilder b = new StringBuilder();
        b.append("##fileformat=VCFv4.2\n");
        b.append("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n");
        b.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
        b.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        b.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n");
        b.append("##contig=<ID=chr1,length=100000000>\n");
        for (int c = 0; c < pad; c++) {
            b.append("##contig=<ID=pad").append(c).append(",length=1000>\n");
        }
        b.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n");
        for (int i = 0; i < n; i++) {
            String chrom = "chr1";
            String alt = "G";
            String info = "DP=" + (i % 97) + ";AC=1";
            if (i == bad) {
                switch (kind) {
                    case "undeclared" -> info = "DP=3;XX=1";
                    case "count" -> info = "DP=3;AC=1,2";
                    case "contig" -> chrom = "chrZ";
                    case "missing" -> {
                        alt = "G,T";
                        info = "DP=3;AC=1,.";
                    }
                    default -> throw new IllegalArgumentException(kind);
                }
            }
            b.append(chrom).append('\t').append(i * 10 + 1).append("\trs").append(i).append("\tA\t").append(alt)
                    .append('\t').append(i % 60).append(".5\tPASS\t").append(info).append("\tGT:DP\t0/1:")
                    .append(i % 30).append("\t1/1:7\t0/0:12\n");
        }
        return b.toString();
    }

    static void run(final String name, final String step, final Path in, final Path out, final boolean index)
            throws Exception {
        final PrintStream realErr = System.err;
        final PrintStream realOut = System.out;
        String outcome;
        try {
            System.setErr(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final int code = new picard.vcf.VcfFormatConverter().instanceMain(new String[] {
                "I=" + in, "O=" + out, "REQUIRE_INDEX=false", "CREATE_INDEX=" + index});
            outcome = code == 0 ? "ok" : "exit " + code;
        } catch (final Throwable e) {
            outcome = e.getClass().getName() + ": "
                    + String.valueOf(e.getMessage()).replace(dir.toAbsolutePath().toString(), "<dir>")
                    .replace('\n', ' ').replace('\t', ' ');
        } finally {
            System.setErr(realErr);
            System.setOut(realOut);
        }
        final boolean idx = Files.exists(Path.of(out + ".idx"));
        emit(step, name, digest(Files.readAllBytes(out)) + " idx=" + (idx ? "yes" : "no") + " " + outcome);
    }

    public static void main(final String[] args) throws Exception {
        dir = Files.createTempDirectory("partial-dump");
        // name, n, bad, kind, pad
        final Object[][] cases = {
            {"late-undeclared", 20000, 18000, "undeclared", 0},
            {"late-count", 20000, 18000, "count", 0},
            {"late-contig", 20000, 18000, "contig", 0},
            {"late-missing", 20000, 18000, "missing", 0},
            {"big-header", 10, 5, "undeclared", 4000},
            {"small", 10, 5, "undeclared", 0},
        };
        for (final Object[] c : cases) {
            final String name = (String) c[0];
            final byte[] text = vcf((int) c[1], (int) c[2], (String) c[3], (int) c[4])
                    .getBytes(StandardCharsets.UTF_8);
            final Path in = dir.resolve(name + ".vcf");
            Files.write(in, text);
            emit("input", name, digest(text));
            run(name, "bcf", in, dir.resolve(name + ".bcf"), false);
            run(name, "vcf", in, dir.resolve(name + ".out.vcf"), false);
            if (c[3].equals("undeclared")) {
                run(name, "bcf-indexed", in, dir.resolve(name + ".indexed.bcf"), true);
                run(name, "vcf-indexed", in, dir.resolve(name + ".indexed.vcf"), true);
            }
        }
        System.out.print(buf);
    }
}
