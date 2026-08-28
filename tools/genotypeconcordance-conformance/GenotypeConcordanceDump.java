/*
 * GenotypeConcordance's three tables, taken from the reference.
 *
 * The tool puts one sample's calls beside another's truth and counts the pairs of genotype states.
 * What is measured is which pair of states a site is filed under, what the three files carry, and
 * which arguments move any of it.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE OUTPUT ARGUMENT IS A BASENAME FOR THREE FILES, whose extensions are the tool's own;
 *   - THE DETAIL FILE IS ONE ROW PER PAIR OF STATES, so a truth het called hom-var is a row of its
 *     own and not a mismatch counter;
 *   - A SITE MISSING FROM ONE SIDE IS A STATE, `MISSING`, and --MISSING_SITES_HOM_REF turns it
 *     into HOM_REF instead, which moves the row and the summary with it;
 *   - THAT FLAG IS REFUSED WITHOUT --INTERVALS, because a site is only confidently hom-ref inside
 *     a confident region, and the interval list in turn requires both VCFs to be indexed;
 *   - --MIN_GQ AND --MIN_DP TURN A CALL INTO `LOW_GQ` AND `LOW_DP`, which are states and not
 *     exclusions: the site is still counted, under another name;
 *   - A FILTERED SITE IS ITS OWN STATE, `FILTERED`, and --IGNORE_FILTER_STATUS reads it as if it
 *     had passed;
 *   - THE SUMMARY IS PER VARIANT TYPE, so a file with a SNP and an indel gives two rows;
 *   - SENSITIVITY IS TRUE POSITIVES OVER TRUTH POSITIVES and specificity is left empty for SNPs,
 *     since a variant caller emits no true negatives;
 *   - THE CONTINGENCY FILE IS THE SAME DATA AS FOUR COUNTERS, TP, TN, FP and FN, and a state pair
 *     may contribute to more than one;
 *   - --OUTPUT_ALL_ROWS KEEPS THE PAIRS NOTHING WAS SEEN FOR, which is a hundred-odd rows rather
 *     than a handful;
 *   - --TRUTH_SAMPLE AND --CALL_SAMPLE NAME THE COLUMNS THE TWO FILES ARE READ FROM, so a file of
 *     two samples is read once per pairing;
 *   - AND A SITE WITH NO CALL ON EITHER SIDE IS `NO_CALL`, which is not the same as missing.
 *
 * Output:
 *
 *     vcf\t<name>\t<that input file, escaped>
 *     files\t<case>\t<the file names the basename produced, comma separated>
 *     summary\t<case>\t<the summary table, escaped>
 *     detail\t<case>\t<the detail rows that counted anything, escaped>
 *     contingency\t<case>\t<the contingency table, escaped>
 *     rows\t<case>\t<how many detail rows the file holds>
 *     error\t<case>\t<exception class>:<message>
 */

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenotypeConcordanceDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 100000;

    /** One VCF, with as many samples as the genotypes give. */
    static String vcf(final List<String> samples, final List<String> sites) {
        final List<String> lines = new ArrayList<>(List.of(
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">",
                "##FILTER=<ID=LowQual,Description=\"Low quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">"));
        lines.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + String.join("\t", samples));
        lines.addAll(sites);
        lines.add("");
        return String.join("\n", lines);
    }

    /** One site: its position, its alleles, its filter and one genotype per sample. */
    static String site(final int position, final String reference, final String alternate,
                       final String filter, final String... genotypes) {
        final List<String> columns = new ArrayList<>(List.of(
                "chr1", Integer.toString(position), ".", reference, alternate, "100.00",
                filter, ".", "GT:GQ:DP"));
        columns.addAll(Arrays.asList(genotypes));
        return String.join("\t", columns);
    }

    /** A genotype with a quality and a depth the two floors can reach. */
    static String genotype(final String alleles, final int quality, final int depth) {
        return alleles + ":" + quality + ":" + depth;
    }

    public static void main(final String[] args) throws Exception {
        final String plain = genotype("0/1", 60, 30);

        // The truth: a het SNP, a hom-var SNP, a het indel, a filtered site and a no-call.
        final List<String> truthSites = List.of(
                site(1000, "A", "G", "PASS", plain),
                site(2000, "C", "T", "PASS", genotype("1/1", 60, 30)),
                site(3000, "A", "ACGT", "PASS", plain),
                site(4000, "G", "T", "PASS", plain),
                site(5000, "A", "C", "PASS", plain),
                site(6000, "T", "A", "PASS", plain));
        final String truth = vcf(List.of("truth"), truthSites);
        emit("vcf", "truth", truth);

        // The calls: one agreeing, one called hom-var where the truth is het, one filtered, one
        // with a low quality, one with a low depth, and one missing altogether.
        final String calls = vcf(List.of("call"), List.of(
                site(1000, "A", "G", "PASS", plain),
                site(2000, "C", "T", "PASS", genotype("1/1", 60, 30)),
                site(3000, "A", "ACGT", "PASS", genotype("1/1", 60, 30)),
                site(4000, "G", "T", "LowQual", plain),
                site(5000, "A", "C", "PASS", genotype("0/1", 5, 30)),
                site(6000, "T", "A", "PASS", genotype("0/1", 60, 2))));
        emit("vcf", "calls", calls);

        run("plain", truth, calls, List.of());
        run("all-rows", truth, calls, List.of("OUTPUT_ALL_ROWS=true"));
        run("min-gq", truth, calls, List.of("MIN_GQ=20"));
        run("min-dp", truth, calls, List.of("MIN_DP=10"));
        run("ignore-filters", truth, calls, List.of("IGNORE_FILTER_STATUS=true"));

        // A call set that is missing a site the truth has.
        final String shortCalls = vcf(List.of("call"), List.of(
                site(1000, "A", "G", "PASS", plain)));
        emit("vcf", "short-calls", shortCalls);
        run("a-missing-site", truth, shortCalls, List.of());
        // The flag needs an interval list, and the interval list needs the VCFs indexed.
        run("missing-as-hom-ref-without-intervals", truth, shortCalls,
                List.of("MISSING_SITES_HOM_REF=true"));
        run("missing-as-hom-ref", truth, shortCalls, List.of("MISSING_SITES_HOM_REF=true"), true);
        run("intervals-alone", truth, shortCalls, List.of(), true);

        // A site neither side called.
        final String noCall = vcf(List.of("call"), List.of(
                site(1000, "A", "G", "PASS", genotype("./.", 60, 30))));
        emit("vcf", "no-call", noCall);
        run("a-no-call", vcf(List.of("truth"), List.of(site(1000, "A", "G", "PASS", plain))),
                noCall, List.of());

        // Two samples in one file, read once per pairing.
        final String pair = vcf(List.of("truth", "other"), List.of(
                "chr1\t1000\t.\tA\tG\t100.00\tPASS\t.\tGT:GQ:DP\t" + plain + "\t"
                        + genotype("1/1", 60, 30)));
        emit("vcf", "two-samples", pair);
        run("second-sample", pair, pair, List.of("TRUTH_SAMPLE=truth", "CALL_SAMPLE=other"));

        System.out.print(buf);
    }

    static void run(final String name, final String truth, final String calls,
                    final List<String> extra) throws Exception {
        run(name, truth, calls, extra, false);
    }

    /** One run, optionally with an interval list and the indexes that then become required. */
    static void run(final String name, final String truth, final String calls,
                    final List<String> extra, final boolean withIntervals) throws Exception {
        final Path dir = Files.createTempDirectory("concordance");
        final Path truthPath = write(dir, "truth.vcf", truth);
        final Path callPath = write(dir, "calls.vcf", calls);
        final Path prefix = dir.resolve("out");
        final List<String> argv = new ArrayList<>(List.of(
                "TRUTH_VCF=" + truthPath, "CALL_VCF=" + callPath, "O=" + prefix));
        if (extra.stream().noneMatch(value -> value.startsWith("TRUTH_SAMPLE"))) {
            argv.add("TRUTH_SAMPLE=truth");
            argv.add("CALL_SAMPLE=call");
        }
        argv.addAll(extra);
        if (withIntervals) {
            final Path intervals = write(dir, "confident.interval_list",
                    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH + "\n"
                            + "chr1\t500\t7000\t+\tconfident\n");
            argv.add("INTERVALS=" + intervals);
            index(truthPath);
            index(callPath);
        }
        final java.io.ByteArrayOutputStream errBytes = new java.io.ByteArrayOutputStream();
        final java.io.PrintStream realErr = System.err;
        try {
            final int code;
            try {
                System.setErr(new java.io.PrintStream(errBytes, true, StandardCharsets.UTF_8));
                code = new picard.vcf.GenotypeConcordance()
                        .instanceMain(argv.toArray(new String[0]));
            } finally {
                System.err.flush();
                System.setErr(realErr);
            }
            if (code != 0) {
                // The refusal is one line of a usage the golden has no reason to hold.
                final List<String> refusal = new ArrayList<>();
                for (final String line : errBytes.toString(StandardCharsets.UTF_8)
                        .split("\n", -1)) {
                    if (line.startsWith("You cannot use the MISSING_HOM option")
                            || line.startsWith("The index file was not found")) {
                        refusal.add(line);
                    }
                }
                emit("error", name, "exit " + code);
                emit("refusal", name, String.join("\n", refusal));
                return;
            }
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        }
        final List<String> written = new ArrayList<>();
        try (final var stream = Files.list(dir)) {
            stream.map(path -> path.getFileName().toString())
                    .filter(file -> file.startsWith("out."))
                    .sorted()
                    .forEach(written::add);
        }
        emit("files", name, String.join(",", written));
        emit("summary", name, table(dir.resolve(
                "out" + picard.vcf.GenotypeConcordance.SUMMARY_METRICS_FILE_EXTENSION)));
        emit("contingency", name, table(dir.resolve(
                "out" + picard.vcf.GenotypeConcordance.CONTINGENCY_METRICS_FILE_EXTENSION)));
        final List<String> detail = new ArrayList<>(List.of(table(dir.resolve(
                "out" + picard.vcf.GenotypeConcordance.DETAILED_METRICS_FILE_EXTENSION))
                .split("\n", -1)));
        emit("rows", name, Integer.toString(Math.max(0, detail.size() - 1)));
        // A hundred-odd rows of zeroes say nothing: the golden keeps the header and the rows that
        // counted something.
        final List<String> counted = new ArrayList<>(List.of(detail.get(0)));
        for (final String row : detail.subList(1, detail.size())) {
            if (!row.endsWith("\t0") && !row.isEmpty()) {
                counted.add(row);
            }
        }
        emit("detail", name, String.join("\n", counted));
    }

    /** A tribble index, which an interval list makes required. */
    static void index(final Path vcf) throws Exception {
        htsjdk.tribble.index.IndexFactory
                .createLinearIndex(vcf.toFile(), new htsjdk.variant.vcf.VCFCodec())
                .writeBasedOnFeatureFile(vcf.toFile());
    }

    static Path write(final Path dir, final String name, final String text) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, text, StandardCharsets.UTF_8);
        return path;
    }

    /** A metrics file without its comment lines. */
    static String table(final Path file) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (!line.startsWith("#") && !line.isEmpty()) {
                kept.add(line);
            }
        }
        return String.join("\n", kept);
    }
}
