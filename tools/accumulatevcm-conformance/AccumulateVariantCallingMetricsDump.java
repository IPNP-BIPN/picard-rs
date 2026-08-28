/*
 * AccumulateVariantCallingMetrics' merged metrics, taken from the reference.
 *
 * The tool reads several pairs of variant-calling metrics files and writes one pair back. There is
 * no VCF involved at all: it is arithmetic over metrics tables, and the arithmetic is not a plain
 * sum.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE ARGUMENTS ARE PREFIXES AND NOT FILES, two fixed extensions being appended to each;
 *   - THE MERGE IS PER SAMPLE_ALIAS, so two files naming one sample give one row and two files
 *     naming two give two;
 *   - THE COUNTS ARE ADDED, so a sample seen in two files carries the sum of its SNPs;
 *   - THE RATIOS ARE NOT: they are RECONSTRUCTED from the printed columns before merging and
 *     recomputed after, so a round trip through the file is LOSSY. A single file of 301 SNPs at a
 *     TI/TV of 2.0 comes back at 2.01, having been merged with nothing at all;
 *   - THE RECONSTRUCTION ROUNDS, `invertFromRatio` being `Math.round(sum / (ratio + 1))`, so a
 *     ratio that does not divide its count loses a variant to rounding;
 *   - A RATIO OF NaN RECONSTRUCTS AS NOUGHT and comes back out as NOUGHT, so the round trip
 *     turns "no ratio" into "a ratio of zero";
 *   - THE SUMMARY'S REFERENCE BIAS IS RECONSTRUCTED AGAINST THE DETAIL FILE BESIDE IT, the
 *     total het depth being summed over that file's detail rows and handed to its summary, so a
 *     summary read without its detail would reconstruct differently;
 *   - A SUMMARY FILE OF MORE THAN ONE ROW IS REFUSED, by a message counting them;
 *   - A MISSING INPUT FILE IS REFUSED BY htsjdk AND NOT BY THE TOOL, so the message names the
 *     FILE it could not read rather than the prefix the tool's own catch would have named;
 *   - THE OUTPUT ROWS COME OUT OF A HashMap, which is why this dump sorts them;
 *   - AND ONE INPUT IS ENOUGH: merging a single pair is a round trip, and shows the loss by itself.
 *
 * Output:
 *
 *     in\t<name>\t<that metrics file, escaped>
 *     detail\t<case>\t<the merged detail table, its rows sorted, escaped>
 *     summary\t<case>\t<the merged summary table, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class AccumulateVariantCallingMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final String SUMMARY_HEADER = String.join("\t",
            "TOTAL_SNPS", "NUM_IN_DB_SNP", "NOVEL_SNPS", "FILTERED_SNPS", "PCT_DBSNP",
            "DBSNP_TITV", "NOVEL_TITV", "TOTAL_INDELS", "NOVEL_INDELS", "FILTERED_INDELS",
            "PCT_DBSNP_INDELS", "NUM_IN_DB_SNP_INDELS", "DBSNP_INS_DEL_RATIO",
            "NOVEL_INS_DEL_RATIO", "TOTAL_MULTIALLELIC_SNPS", "NUM_IN_DB_SNP_MULTIALLELIC",
            "TOTAL_COMPLEX_INDELS", "NUM_IN_DB_SNP_COMPLEX_INDELS", "SNP_REFERENCE_BIAS",
            "NUM_SINGLETONS");

    static final String DETAIL_HEADER = "SAMPLE_ALIAS\t" + SUMMARY_HEADER
            + "\tHET_HOMVAR_RATIO\tPCT_GQ0_VARIANTS\tTOTAL_GQ0_VARIANTS\tTOTAL_HET_DEPTH";

    /** One summary row: the twenty columns, most of them nought. */
    static String summaryRow(final long totalSnps, final long inDbSnp, final long novelSnps,
                             final String dbsnpTiTv, final String novelTiTv,
                             final String referenceBias) {
        return String.join("\t",
                Long.toString(totalSnps), Long.toString(inDbSnp), Long.toString(novelSnps), "0",
                "0", dbsnpTiTv, novelTiTv, "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0",
                referenceBias, "0");
    }

    static String detailRow(final String sample, final long totalSnps, final long inDbSnp,
                            final long novelSnps, final String dbsnpTiTv, final String novelTiTv,
                            final String referenceBias, final String hetHomVar,
                            final long hetDepth) {
        return sample + "\t"
                + summaryRow(totalSnps, inDbSnp, novelSnps, dbsnpTiTv, novelTiTv, referenceBias)
                + "\t" + hetHomVar + "\t0\t0\t" + hetDepth;
    }

    static String metricsFile(final String beanClass, final String header, final List<String> rows) {
        final List<String> lines = new ArrayList<>();
        lines.add("## htsjdk.samtools.metrics.StringHeader");
        lines.add("# a fixture");
        lines.add("");
        lines.add("## METRICS CLASS\t" + beanClass);
        lines.add(header);
        lines.addAll(rows);
        lines.add("");
        return String.join("\n", lines);
    }

    static final String DETAIL_CLASS =
            "picard.vcf.CollectVariantCallingMetrics$VariantCallingDetailMetrics";
    static final String SUMMARY_CLASS =
            "picard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics";

    /** A table without its comments, its data rows SORTED. */
    static String table(final String text) {
        final List<String> head = new ArrayList<>();
        final List<String> rows = new ArrayList<>();
        boolean seenHeader = false;
        for (final String line : text.split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            if (!seenHeader) {
                head.add(line);
                seenHeader = true;
            } else {
                rows.add(line);
            }
        }
        Collections.sort(rows);
        head.addAll(rows);
        return String.join("\n", head);
    }

    /** One input: a prefix, its detail rows and its summary row. */
    record Input(String prefix, List<String> detail, String summary) {}

    static void run(final String name, final List<Input> inputs, final boolean writeFiles)
            throws Exception {
        final Path dir = Files.createTempDirectory("accumulatevcm");
        final List<String> argv = new ArrayList<>();
        for (final Input input : inputs) {
            final Path base = dir.resolve(input.prefix());
            if (writeFiles) {
                Files.writeString(Path.of(base + ".variant_calling_detail_metrics"),
                        metricsFile(DETAIL_CLASS, DETAIL_HEADER, input.detail()),
                        StandardCharsets.UTF_8);
                Files.writeString(Path.of(base + ".variant_calling_summary_metrics"),
                        metricsFile(SUMMARY_CLASS, SUMMARY_HEADER,
                                input.summary() == null ? List.of()
                                        : List.of(input.summary().split("\n", -1))),
                        StandardCharsets.UTF_8);
            }
            argv.add("I=" + base);
        }
        final Path out = dir.resolve("out");
        argv.add("O=" + out);
        try {
            final int code = new picard.vcf.AccumulateVariantCallingMetrics()
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
        emit("detail", name, table(Files.readString(
                Path.of(out + ".variant_calling_detail_metrics"), StandardCharsets.UTF_8)));
        emit("summary", name, table(Files.readString(
                Path.of(out + ".variant_calling_summary_metrics"), StandardCharsets.UTF_8)));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // One sample, TI/TV of exactly 2 over a count that divides evenly.
        final Input even = new Input("even",
                List.of(detailRow("s1", 300, 300, 0, "2.0", "NaN", "0.5", "1.0", 100)),
                summaryRow(300, 300, 0, "2.0", "NaN", "0.5"));
        emit("in", "even-detail", metricsFile(DETAIL_CLASS, DETAIL_HEADER, even.detail()));
        emit("in", "even-summary", metricsFile(SUMMARY_CLASS, SUMMARY_HEADER,
                List.of(even.summary())));

        // A round trip of one input, which is where the loss shows by itself.
        run("one-input-even", List.of(even), true);

        // A count that does NOT divide evenly by the ratio.
        final Input odd = new Input("odd",
                List.of(detailRow("s1", 301, 301, 0, "2.0", "NaN", "0.5", "1.0", 100)),
                summaryRow(301, 301, 0, "2.0", "NaN", "0.5"));
        run("one-input-odd", List.of(odd), true);

        // Two inputs of the same sample, whose counts add.
        run("two-inputs-one-sample", List.of(even, new Input("even2", even.detail(),
                even.summary())), true);

        // Two inputs of two different samples, which stay apart.
        final Input other = new Input("other",
                List.of(detailRow("s2", 100, 100, 0, "1.0", "NaN", "0.25", "2.0", 40)),
                summaryRow(100, 100, 0, "1.0", "NaN", "0.25"));
        run("two-samples", List.of(even, other), true);

        // One input naming two samples in its detail file.
        run("two-samples-one-file", List.of(new Input("both",
                List.of(detailRow("s1", 300, 300, 0, "2.0", "NaN", "0.5", "1.0", 100),
                        detailRow("s2", 100, 100, 0, "1.0", "NaN", "0.25", "2.0", 40)),
                summaryRow(400, 400, 0, "1.75", "NaN", "0.4"))), true);

        // A NaN ratio, which reconstructs as nought.
        run("nan-ratio", List.of(new Input("nan",
                List.of(detailRow("s1", 300, 300, 0, "NaN", "NaN", "NaN", "1.0", 100)),
                summaryRow(300, 300, 0, "NaN", "NaN", "NaN"))), true);

        // A summary file of two rows.
        run("two-summary-rows", List.of(new Input("twosum",
                List.of(detailRow("s1", 300, 300, 0, "2.0", "NaN", "0.5", "1.0", 100)),
                summaryRow(300, 300, 0, "2.0", "NaN", "0.5") + "\n"
                        + summaryRow(1, 1, 0, "2.0", "NaN", "0.5"))), true);

        // A prefix whose files are not there.
        run("missing-input", List.of(new Input("gone", List.of(), null)), false);

        System.out.print(buf);
    }
}
