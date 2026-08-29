/*
 * CombineGenotypingArrayVcfs' merge, taken from the reference.
 *
 * The tool puts several single-sample array VCFs side by side into one multi-sample VCF. It is not
 * a general merge: it walks the inputs IN LOCKSTEP, one variant from each per step, and refuses
 * anything that does not line up. What is measured is that alignment and the refusals around it.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE INPUTS ARE WALKED IN LOCKSTEP AND NOT MERGED BY POSITION: the nth variant of each file
 *     is merged with the nth of every other, so two files that hold the same loci in a different
 *     order are a refusal and not a reordering;
 *   - A DIFFERENT NUMBER OF VARIANTS IS REFUSED, and the message says so rather than naming a
 *     locus, because the loop finds out only when one iterator runs dry;
 *   - THE LOCUS, THE ID, THE REF AND THE ALT MUST ALL AGREE, each with its own message, and the
 *     ALT one names the locus where the others do not;
 *   - AN ALT COUNT THAT DIFFERS IS ITS OWN REFUSAL, separate from an ALT that differs;
 *   - A SAMPLE THAT APPEARS TWICE IS REFUSED BY NAME, with the file that repeated it;
 *   - THE SAMPLE ORDER IS THE INPUT ORDER, file by file and then within a file;
 *   - AN ATTRIBUTE PRESENT IN ONE FILE AND NOT ANOTHER IS A REFUSAL, and so is one whose values
 *     disagree, both naming the key;
 *   - SEVEN ATTRIBUTES ARE EXEMPT FROM THAT CHECK, `AC`, `AF`, `AN`, `devX_AB`, `devY_AB`,
 *     `SOURCE` and `refSNP`, which may differ freely;
 *   - `DP` IS THE ONE ATTRIBUTE THE MERGE MEANS TO ADD UP, AND IT IS THE ONE THAT MAKES THE TOOL
 *     THROW: the sum is written back with `firstAttributes.put`, and `getAttributes` returns an
 *     UNMODIFIABLE map, so any input carrying a depth ends in an UnsupportedOperationException
 *     with no message. The dump records the frame it was thrown from, because the message says
 *     nothing on its own;
 *   - THE SAMPLE-SPECIFIC HEADER LINES ARE DROPPED, so the merged header carries no
 *     `chipWellBarcode` and no `autocallDate`, which belonged to one sample each;
 *   - THE FILTERS OF EVERY INPUT ARE UNIONED onto the merged variant;
 *   - AN ATTRIBUTE IN THE FIRST FILE ONLY IS KEPT WITHOUT COMMENT, where one in a later file is
 *     refused: the loop runs over the other files' attributes and looks each up in the first's,
 *     so the direction decides whether there is a refusal at all;
 *   - AND THE OUTPUT IS INDEXED WHETHER OR NOT `CREATE_INDEX` ASKED FOR IT, the writer's builder
 *     indexing any file it has a sequence dictionary for. The argument only adds an option the
 *     builder already had.
 *
 * Output:
 *
 *     input\t<case>/<file>\t<the variant lines, escaped>
 *     merged\t<case>\t<the output's variant lines, escaped>
 *     header\t<case>\t<the output's header lines, escaped>
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CombineGenotypingArrayVcfsDump
 */

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CombineGenotypingArrayVcfsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 1000;

    /** One variant line of one input. */
    record Variant(int position, String id, String reference, String alternate, String info,
                   String filter, String genotype) {}

    static Variant variant(final int position, final String genotype) {
        return new Variant(position, "rs" + position, "A", "C", "BEADSET=7", "PASS", genotype);
    }

    /** One input file: its sample, and its variants. */
    record Input(String sample, List<Variant> variants) {}

    static String vcf(final Input input) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=").append(CONTIG_LENGTH).append(">\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
        text.append("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n");
        // The merge recalculates AC, AF and AN, and the writer refuses a key the header does not
        // define, so all three have to be there before the output can be written at all.
        text.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n");
        text.append("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n");
        text.append("##INFO=<ID=BEADSET,Number=1,Type=Integer,Description=\"Bead set\">\n");
        text.append("##INFO=<ID=EXTRA,Number=1,Type=Integer,Description=\"An extra key\">\n");
        text.append("##FILTER=<ID=LOW,Description=\"Low quality\">\n");
        // Two of the sample-specific header lines, which the merge is supposed to drop.
        text.append("##chipWellBarcode=barcode-").append(input.sample()).append('\n');
        text.append("##autocallDate=09/21/2016 20:40\n");
        text.append("##arrayType=TestArray-24v1-0_A1\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
                .append(input.sample()).append('\n');
        for (final Variant v : input.variants()) {
            text.append("chr1\t").append(v.position()).append('\t').append(v.id()).append('\t')
                    .append(v.reference()).append('\t').append(v.alternate()).append("\t100\t")
                    .append(v.filter()).append('\t').append(v.info()).append("\tGT\t")
                    .append(v.genotype()).append('\n');
        }
        return text.toString();
    }

    /** A file's lines, split into its header and its variants. */
    static String[] split(final Path file) throws Exception {
        final List<String> header = new ArrayList<>();
        final List<String> variants = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.isEmpty()) {
                continue;
            }
            (line.startsWith("#") ? header : variants).add(line);
        }
        return new String[]{String.join("\n", header), String.join("\n", variants)};
    }

    static void run(final String name, final List<Input> inputs, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("combinearrays");
        final List<String> argv = new ArrayList<>();
        for (int i = 0; i < inputs.size(); i++) {
            final Path file = dir.resolve("in" + i + ".vcf");
            final String text = vcf(inputs.get(i));
            Files.writeString(file, text, StandardCharsets.UTF_8);
            emit("input", name + "/in" + i, text.lines()
                    .filter(line -> !line.startsWith("##"))
                    .reduce((a, b) -> a + "\n" + b).orElse(""));
            argv.add("I=" + file);
        }
        final Path out = Files.createDirectory(dir.resolve("out"));
        argv.add("O=" + out.resolve("merged.vcf"));
        // The index is off unless a case asks for it: a plain text VCF written with the index on
        // the fly is what the `an-index` case is for.
        if (Arrays.stream(extra).noneMatch(a -> a.startsWith("CREATE_INDEX"))) {
            argv.add("CREATE_INDEX=false");
        }
        argv.addAll(Arrays.asList(extra));
        final java.io.ByteArrayOutputStream errBytes = new java.io.ByteArrayOutputStream();
        final java.io.PrintStream realErr = System.err;
        try {
            final int code;
            try {
                System.setErr(new java.io.PrintStream(errBytes, true, StandardCharsets.UTF_8));
                code = new picard.arrays.CombineGenotypingArrayVcfs()
                        .instanceMain(argv.toArray(new String[0]));
            } finally {
                System.err.flush();
                System.setErr(realErr);
            }
            if (code != 0) {
                emit("error", name, "exit " + code);
                return;
            }
        } catch (final Exception e) {
            System.setErr(realErr);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            // An exception with no message says nothing on its own, so the frame it was thrown
            // from is recorded beside it: that is what tells a reader WHICH operation the
            // reference refused.
            final StackTraceElement[] frames = cause.getStackTrace();
            final String where = frames.length == 0 ? ""
                    : " at " + frames[0].getClassName() + "." + frames[0].getMethodName();
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>") + where);
            return;
        }
        final List<String> written = new ArrayList<>();
        for (final File file : out.toFile().listFiles()) {
            written.add(file.getName());
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        final String[] parts = split(out.resolve("merged.vcf"));
        emit("header", name, parts[0]);
        emit("merged", name, parts[1]);
    }

    public static void main(final String[] args) throws Exception {
        final List<Variant> first = List.of(variant(100, "0/0"), variant(200, "0/1"));
        final List<Variant> second = List.of(variant(100, "1/1"), variant(200, "0/0"));

        // Two files, two samples, the same loci.
        run("two-samples", List.of(
                new Input("sampleA", first), new Input("sampleB", second)));
        // Three, to show the sample order is the input's.
        run("three-samples", List.of(
                new Input("sampleC", first), new Input("sampleA", second),
                new Input("sampleB", first)));
        // One file on its own, which is a merge of one.
        run("one-sample", List.of(new Input("sampleA", first)));

        // A sample that appears twice.
        run("a-repeated-sample", List.of(
                new Input("sampleA", first), new Input("sampleA", second)));

        // The lockstep: the same loci in the other order is a refusal and not a reordering.
        run("loci-in-another-order", List.of(
                new Input("sampleA", first),
                new Input("sampleB", List.of(variant(200, "0/1"), variant(100, "0/0")))));
        run("a-different-number-of-variants", List.of(
                new Input("sampleA", first),
                new Input("sampleB", List.of(variant(100, "0/0")))));
        run("a-different-locus", List.of(
                new Input("sampleA", first),
                new Input("sampleB", List.of(variant(100, "0/0"), variant(300, "0/1")))));
        run("a-different-id", List.of(
                new Input("sampleA", first),
                new Input("sampleB", List.of(
                        new Variant(100, "other", "A", "C", "DP=10", "PASS", "0/0"),
                        variant(200, "0/1")))));
        run("a-different-reference-allele", List.of(
                new Input("sampleA", first),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "T", "C", "DP=10", "PASS", "0/0"),
                        variant(200, "0/1")))));
        run("a-different-alternate-allele", List.of(
                new Input("sampleA", first),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "A", "G", "DP=10", "PASS", "0/0"),
                        variant(200, "0/1")))));
        run("a-different-alternate-count", List.of(
                new Input("sampleA", first),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "A", "C,G", "DP=10", "PASS", "0/0"),
                        variant(200, "0/1")))));

        // The attributes: one missing, one disagreeing, one exempt, and the one that is summed.
        // The loop runs over the OTHER files' attributes and looks each up in the FIRST file's,
        // so which file carries the extra key decides whether it is a refusal at all: one only the
        // first file has is kept without comment, and one a later file has is refused.
        run("an-attribute-in-the-first-file-only", List.of(
                new Input("sampleA", List.of(
                        new Variant(100, "rs100", "A", "C", "BEADSET=7;EXTRA=1", "PASS", "0/0"))),
                new Input("sampleB", List.of(variant(100, "0/1")))));
        run("an-attribute-in-a-later-file-only", List.of(
                new Input("sampleA", List.of(variant(100, "0/0"))),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "A", "C", "BEADSET=7;EXTRA=1", "PASS", "0/1")))));
        run("an-attribute-that-disagrees", List.of(
                new Input("sampleA", List.of(
                        new Variant(100, "rs100", "A", "C", "BEADSET=7", "PASS", "0/0"))),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "A", "C", "BEADSET=8", "PASS", "0/1")))));
        run("an-exempt-attribute-that-disagrees", List.of(
                new Input("sampleA", List.of(
                        new Variant(100, "rs100", "A", "C", "BEADSET=7;AC=1", "PASS", "0/0"))),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "A", "C", "BEADSET=7;AC=2", "PASS", "0/1")))));
        // `DP` is the one attribute the merge means to ADD UP, and it is the one input that makes
        // the tool throw: the sum is written back with `firstAttributes.put`, and
        // `VariantContext.getAttributes` returns an UNMODIFIABLE map. So any input carrying a
        // depth is refused by an UnsupportedOperationException with no message, from
        // `Collections$UnmodifiableMap.put`, and the depths never reach the output.
        run("a-depth-in-both-files", List.of(
                new Input("sampleA", List.of(
                        new Variant(100, "rs100", "A", "C", "DP=10", "PASS", "0/0"))),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "A", "C", "DP=10", "PASS", "0/1")))));
        run("depths-that-differ", List.of(
                new Input("sampleA", List.of(
                        new Variant(100, "rs100", "A", "C", "DP=10", "PASS", "0/0"))),
                new Input("sampleB", List.of(
                        new Variant(100, "rs100", "A", "C", "DP=30", "PASS", "0/1")))));
        run("a-depth-in-one-file-only", List.of(
                new Input("sampleA", List.of(
                        new Variant(100, "rs100", "A", "C", "DP=10", "PASS", "0/0"))),
                new Input("sampleB", List.of(variant(100, "0/1")))));

        // The index is written either way: `VariantContextWriterBuilder` indexes a file it has a
        // dictionary for, and `CREATE_INDEX` only adds the option the builder already had. Both
        // cases are here because the argument reads as if it decided something.
        run("an-index", List.of(
                new Input("sampleA", first), new Input("sampleB", second)),
                "CREATE_INDEX=true");

        // The filters, which are unioned.
        run("filters-that-differ", List.of(
                new Input("sampleA", List.of(
                        new Variant(100, "rs100", "A", "C", "BEADSET=7", "LOW", "0/0"))),
                new Input("sampleB", List.of(variant(100, "0/1")))));

        System.out.print(buf);
    }
}
