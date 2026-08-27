/*
 * MakeVcfSampleNameMap's map, taken from the reference.
 *
 * A TSV from VCF path to sample name, one line per input. What is measured is what goes on a
 * line, in what order the lines come out, and which inputs are refused.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE LINE IS THE PATH FIRST AND THE NAME SECOND, which is the opposite way round from the
 *     tool's name;
 *   - THE PATH ON THE LINE IS THE STRING AS GIVEN, never normalised or canonicalised, so a
 *     `<dir>/./a.vcf` argument comes back with its dot in it;
 *   - THE LINES ARE NOT IN INPUT ORDER: they come out of a HashMap keyed by the path string, so
 *     the order is that map's, which is the paths' hashes and not the argument list;
 *   - THE SAME PATH GIVEN TWICE IS ONE LINE, because that map is keyed by the path;
 *   - TWO PATHS THAT NAME THE SAME SAMPLE ARE BOTH KEPT, with only a warning, so the file may
 *     hold the same sample name twice;
 *   - A VCF WITH NO SAMPLE AT ALL IS REFUSED, by a message that names the input and the count;
 *   - SO IS ONE WITH TWO, by the same message;
 *   - A FILE THAT IS NOT A VCF IS REFUSED while it is being read rather than counted;
 *   - AND THE OUTPUT ALWAYS ENDS ON A NEWLINE, `Files.write` writing one after every line.
 *
 * Output:
 *
 *     vcf\t<name>\t<that input file, escaped>
 *     out\t<case>\t<the map file, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class MakeVcfSampleNameMapDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A minimal VCF whose header names the samples given, with one site each has a call at. */
    static String vcf(final List<String> samples) {
        final List<String> lines = new ArrayList<>();
        lines.add("##fileformat=VCFv4.2");
        lines.add("##contig=<ID=chr1,length=1000>");
        lines.add("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        final List<String> columns = new ArrayList<>(List.of(
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"));
        if (!samples.isEmpty()) {
            columns.add("FORMAT");
            columns.addAll(samples);
        }
        lines.add(String.join("\t", columns));
        final List<String> record = new ArrayList<>(List.of(
                "chr1", "100", ".", "A", "C", "50", "PASS", "."));
        if (!samples.isEmpty()) {
            record.add("GT");
            for (int i = 0; i < samples.size(); i++) {
                record.add("0/1");
            }
        }
        lines.add(String.join("\t", record));
        return String.join("\n", lines) + "\n";
    }

    /** Runs the tool over the inputs, which are named RELATIVE to the working directory. */
    static void run(final Path dir, final String name, final List<String> inputs) {
        final Path out = dir.resolve("out-" + name + ".sample_map");
        final List<String> argv = new ArrayList<>();
        for (final String input : inputs) {
            argv.add("I=" + input);
        }
        argv.add("O=" + out);
        try {
            final int code = new picard.vcf.MakeVcfSampleNameMap()
                    .instanceMain(argv.toArray(new String[0]));
            if (code != 0) {
                emit("error", name, "exit " + code);
                return;
            }
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        }
        try {
            emit("out", name, Files.readString(out).replace(dir.toString(), "<dir>"));
        } catch (final Exception e) {
            emit("error", name, "unreadable output");
        }
    }

    static Path write(final Path dir, final String name, final String content) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, content, StandardCharsets.UTF_8);
        return path;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        final Path dir = Files.createTempDirectory("vcfsamplenamemap");

        // Four single-sample VCFs whose file names are short enough to read in the output and
        // whose hashes do not order the way the argument list does.
        final String one = vcf(List.of("SAMPLE_ONE"));
        final String two = vcf(List.of("SAMPLE_TWO"));
        final String three = vcf(List.of("SAMPLE_THREE"));
        final String alsoOne = vcf(List.of("SAMPLE_ONE"));
        emit("vcf", "one", one);
        emit("vcf", "none", vcf(List.of()));
        emit("vcf", "two-samples", vcf(List.of("SAMPLE_ONE", "SAMPLE_TWO")));
        write(dir, "a.vcf", one);
        write(dir, "b.vcf", two);
        write(dir, "c.vcf", three);
        write(dir, "d.vcf", alsoOne);
        write(dir, "empty.vcf", vcf(List.of()));
        write(dir, "pair.vcf", vcf(List.of("SAMPLE_ONE", "SAMPLE_TWO")));
        write(dir, "not-a-vcf.vcf", "this is not a VCF at all\n");

        run(dir, "one-input", List.of(dir.resolve("a.vcf").toString()));
        run(dir, "three-inputs", List.of(
                dir.resolve("a.vcf").toString(),
                dir.resolve("b.vcf").toString(),
                dir.resolve("c.vcf").toString()));
        // The same three, named in the opposite order: the output order does not follow.
        run(dir, "three-inputs-reversed", List.of(
                dir.resolve("c.vcf").toString(),
                dir.resolve("b.vcf").toString(),
                dir.resolve("a.vcf").toString()));
        run(dir, "same-path-twice", List.of(
                dir.resolve("a.vcf").toString(),
                dir.resolve("a.vcf").toString()));
        run(dir, "same-sample-two-paths", List.of(
                dir.resolve("a.vcf").toString(),
                dir.resolve("d.vcf").toString()));
        run(dir, "no-sample", List.of(dir.resolve("empty.vcf").toString()));
        run(dir, "two-samples", List.of(dir.resolve("pair.vcf").toString()));
        run(dir, "not-a-vcf", List.of(dir.resolve("not-a-vcf.vcf").toString()));
        // A good input beside a bad one: the run stops on the bad one and writes nothing.
        run(dir, "good-then-bad", List.of(
                dir.resolve("a.vcf").toString(),
                dir.resolve("pair.vcf").toString()));
        run(dir, "missing-input", List.of(dir.resolve("gone.vcf").toString()));
        // The same file named twice, once plainly and once with a dot in the middle: the map is
        // keyed by the STRING, so these are two entries pointing at one file.
        run(dir, "unnormalised-path", List.of(
                dir.resolve("a.vcf").toString(),
                dir.resolve(".").resolve("a.vcf").toString()));

        System.out.print(buf);
    }
}
