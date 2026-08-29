/*
 * VcfToAdpc's binary output, taken from the reference.
 *
 * An `adpc.bin` is a sixteen-byte header and then one fixed-width record per sample per locus:
 * two unsigned shorts of raw intensity, three floats, and one more unsigned short for the
 * genotype. Everything is little-endian and nothing is compressed, so the file is byte-comparable
 * and what is interesting is which value each field takes.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE HEADER IS THE SIXTEEN CHARACTERS `1234567890123456`, which is a literal in the writer
 *     and not a magic number;
 *   - A RECORD IS TWENTY BYTES: two unsigned shorts, three floats, one unsigned short, in that
 *     order, and the floats are IEEE-754 little-endian;
 *   - THE GENOTYPE IS AN ILLUMINA ONE AND NOT THE VCF'S: AA, AB, BB and NN are 0, 1, 2 and 3, and
 *     which of them a call maps to is decided by the ALLELE_A and ALLELE_B header fields rather
 *     than by the reference and alternate alleles;
 *   - AN ALLELE THAT IS THE REFERENCE CARRIES A TRAILING `*` IN THOSE FIELDS, which the tool
 *     strips before matching, so a VCF that writes `A*` and one that writes `A` agree;
 *   - A NO-CALL IS `NN`, and it is the one genotype that needs no allele fields to be matched;
 *   - AN INTENSITY OVER 65535 IS TRUNCATED rather than refused, and the tool warns;
 *   - A NEGATIVE INTENSITY IS REFUSED, which is the other side of the same check;
 *   - A MISSING NORMALIZED INTENSITY IS WRITTEN AS `NaN`, the field being optional where the raw
 *     intensities are not;
 *   - A MISSING REQUIRED FIELD IS A REFUSAL, and so is a VCF with no records
 *     (`Found no records in VCF`) and a set of VCFs of differing lengths (`VCFs have differing
 *     number of loci`). Those messages are NOT in the golden: the tool catches its own exception
 *     and logs it, and the log does not reach a stream this dump can capture, so what a refusal
 *     leaves behind is an exit code of one;
 *   - THE SAMPLES ARE WRITTEN ONE PER LINE with no trailing newline, and the marker count is a
 *     bare number;
 *   - AND SEVERAL VCFS MUST AGREE ON THEIR NUMBER OF LOCI, a mismatch being refused after the
 *     records of the first have already been written.
 *
 * Output:
 *
 *     vcf\t<case>\t<the variant lines, escaped>
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     adpc\t<case>\t<the binary output, hex>
 *     text\t<case>/<file>\t<that text file's contents, escaped>
 *     error\t<case>\t<the reason, as the tool logged it>
 *
 * Usage: VcfToAdpcDump
 */

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class VcfToAdpcDump {

    static final StringBuilder buf = new StringBuilder();

    /** Everything the run wrote to the error stream, from before the first class was loaded. */
    static final java.io.ByteArrayOutputStream CAPTURED = new java.io.ByteArrayOutputStream();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static String hex(final byte[] bytes) {
        final StringBuilder text = new StringBuilder();
        for (final byte value : bytes) {
            text.append(String.format("%02x", value));
        }
        return text.toString();
    }

    /** One site of the fixture: its alleles, the Illumina A and B, and one sample's fields. */
    record Site(int position, String reference, String alternate, String alleleA, String alleleB,
                String gcScore, String genotype, String x, String y, String normX, String normY) {}

    static Site site(final int position, final String genotype) {
        return new Site(position, "A", "C", "A", "C", "0.75", genotype, "1000", "2000",
                "0.5", "1.5");
    }

    static String vcf(final List<Site> sites, final List<String> samples) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=1000>\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"Raw X\">\n");
        text.append("##FORMAT=<ID=Y,Number=1,Type=Integer,Description=\"Raw Y\">\n");
        text.append("##FORMAT=<ID=NORMX,Number=1,Type=Float,Description=\"Normalized X\">\n");
        text.append("##FORMAT=<ID=NORMY,Number=1,Type=Float,Description=\"Normalized Y\">\n");
        text.append("##INFO=<ID=GC_SCORE,Number=1,Type=Float,Description=\"GenTrain score\">\n");
        text.append("##INFO=<ID=ALLELE_A,Number=1,Type=String,Description=\"Illumina A allele\">\n");
        text.append("##INFO=<ID=ALLELE_B,Number=1,Type=String,Description=\"Illumina B allele\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (final String sample : samples) {
            text.append('\t').append(sample);
        }
        text.append('\n');
        for (final Site s : sites) {
            text.append("chr1\t").append(s.position()).append("\trs").append(s.position())
                    .append('\t').append(s.reference()).append('\t').append(s.alternate())
                    .append("\t100\tPASS\tGC_SCORE=").append(s.gcScore())
                    .append(";ALLELE_A=").append(s.alleleA())
                    .append(";ALLELE_B=").append(s.alleleB())
                    .append("\tGT:X:Y");
            if (s.normX() != null) {
                text.append(":NORMX:NORMY");
            }
            for (int i = 0; i < samples.size(); i++) {
                text.append('\t').append(s.genotype()).append(':').append(s.x()).append(':')
                        .append(s.y());
                if (s.normX() != null) {
                    text.append(':').append(s.normX()).append(':').append(s.normY());
                }
            }
            text.append('\n');
        }
        return text.toString();
    }

    static void run(final String name, final List<List<Site>> inputs, final List<String> samples)
            throws Exception {
        final Path dir = Files.createTempDirectory("vcftoadpc");
        final List<String> argv = new ArrayList<>();
        final StringBuilder written = new StringBuilder();
        for (int i = 0; i < inputs.size(); i++) {
            final Path file = dir.resolve("in" + i + ".vcf");
            final String text = vcf(inputs.get(i), samples);
            Files.writeString(file, text, StandardCharsets.UTF_8);
            if (i == 0) {
                written.append(text.lines().filter(line -> !line.startsWith("##"))
                        .reduce((a, b) -> a + "\n" + b).orElse(""));
            }
            argv.add("VCF=" + file);
        }
        emit("vcf", name, written.toString());
        final Path out = Files.createDirectory(dir.resolve("out"));
        argv.add("O=" + out.resolve("out.adpc.bin"));
        argv.add("SF=" + out.resolve("samples.txt"));
        argv.add("NMF=" + out.resolve("markers.txt"));

        // The tool LOGS its exception and returns 1 rather than throwing, and htsjdk's Log
        // writes to stdout, so both streams are captured: a refusal that reached neither would
        // leave the golden holding a bare exit code.
        final java.io.ByteArrayOutputStream errBytes = new java.io.ByteArrayOutputStream();
        final java.io.ByteArrayOutputStream outBytes = new java.io.ByteArrayOutputStream();
        final java.io.PrintStream realErr = System.err;
        final java.io.PrintStream realOut = System.out;
        final int code;
        try {
            System.setErr(new java.io.PrintStream(errBytes, true, StandardCharsets.UTF_8));
            System.setOut(new java.io.PrintStream(outBytes, true, StandardCharsets.UTF_8));
            code = new picard.arrays.VcfToAdpc().instanceMain(argv.toArray(new String[0]));
        } finally {
            System.err.flush();
            System.out.flush();
            System.setErr(realErr);
            System.setOut(realOut);
        }
        if (code != 0) {
            // The tool logs its exception and returns 1 rather than throwing, so the reason is on
            // the error stream and the run looks like an ordinary failure from outside.
            emit("error", name, "exit " + code + " "
                    + reason(errBytes.toString(StandardCharsets.UTF_8)
                            + "\n" + outBytes.toString(StandardCharsets.UTF_8))
                            .replace(dir.toString(), "<dir>"));
            return;
        }
        final List<String> names = new ArrayList<>();
        for (final File file : out.toFile().listFiles()) {
            names.add(file.getName());
        }
        Collections.sort(names);
        emit("files", name, String.join(" ", names));
        emit("adpc", name, hex(Files.readAllBytes(out.resolve("out.adpc.bin"))));
        for (final String file : List.of("samples.txt", "markers.txt")) {
            emit("text", name + "/" + file,
                    Files.readString(out.resolve(file), StandardCharsets.UTF_8));
        }
    }

    /** The line the tool logged, which is where its refusal ends up. */
    static String reason(final String stderr) {
        String last = "";
        for (final String line : stderr.split("\n", -1)) {
            final String trimmed = line.trim();
            if (trimmed.contains("Exception") || trimmed.startsWith("ERROR")) {
                last = trimmed;
            }
        }
        return last;
    }

    public static void main(final String[] args) throws Exception {
        final java.io.PrintStream realErr = System.err;
        System.setErr(new java.io.PrintStream(CAPTURED, true, StandardCharsets.UTF_8));
        try {
            dump();
        } finally {
            System.err.flush();
            System.setErr(realErr);
        }
        System.out.print(buf);
    }

    static void dump() throws Exception {
        final List<String> one = List.of("sample1");

        // The three called genotypes and the no-call, which are the four codes.
        run("homozygous-a", List.of(List.of(site(100, "0/0"))), one);
        run("heterozygous", List.of(List.of(site(100, "0/1"))), one);
        run("homozygous-b", List.of(List.of(site(100, "1/1"))), one);
        run("a-no-call", List.of(List.of(site(100, "./."))), one);

        // The allele fields, which decide the mapping rather than the reference and alternate.
        run("the-alleles-reversed", List.of(List.of(
                new Site(100, "A", "C", "C", "A", "0.75", "0/0", "1000", "2000", "0.5", "1.5"))),
                one);
        run("a-reference-allele-with-a-star", List.of(List.of(
                new Site(100, "A", "C", "A*", "C", "0.75", "0/0", "1000", "2000", "0.5", "1.5"))),
                one);

        // The intensities, on both sides of an unsigned short.
        run("an-intensity-at-the-limit", List.of(List.of(
                new Site(100, "A", "C", "A", "C", "0.75", "0/1", "65535", "0", "0.5", "1.5"))),
                one);
        run("an-intensity-over-the-limit", List.of(List.of(
                new Site(100, "A", "C", "A", "C", "0.75", "0/1", "70000", "0", "0.5", "1.5"))),
                one);

        // A missing normalized intensity, which is optional where the raw ones are not.
        run("without-the-normalized-intensities", List.of(List.of(
                new Site(100, "A", "C", "A", "C", "0.75", "0/1", "1000", "2000", null, null))),
                one);

        // Two loci and two samples, which is what the record order is read from.
        run("two-loci", List.of(List.of(site(100, "0/0"), site(200, "0/1"))), one);
        run("two-samples", List.of(List.of(site(100, "0/0"))), List.of("sample1", "sample2"));

        // Two VCFs, agreeing and not.
        run("two-vcfs", List.of(
                List.of(site(100, "0/0")), List.of(site(100, "1/1"))), one);
        run("two-vcfs-of-different-lengths", List.of(
                List.of(site(100, "0/0")), List.of(site(100, "1/1"), site(200, "0/1"))), one);

        // A VCF with no records at all.
        run("no-records", List.of(List.of()), one);
    }
}
