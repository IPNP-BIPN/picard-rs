/*
 * CollectQualityYieldMetricsFlow's metrics, taken from the reference.
 *
 * The flow-based cousin of CollectQualityYieldMetrics: it counts flows rather than bases, and the
 * yield it reports is an equivalent one. What is measured is which reads are counted, how a flow
 * reaches a threshold, and what each of the three arguments changes.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE UNIT IS THE FLOW AND NOT THE BASE, so a read of eight bases in four homopolymers is
 *     four flows and one of eight singles is eight;
 *   - A READ THAT FAILS VENDOR QUALITY IS COUNTED IN TOTAL_READS AND NOT IN PF_READS, and its
 *     flows are counted nowhere;
 *   - A SKIPPED SECONDARY OR SUPPLEMENTARY READ IS COUNTED NOWHERE AT ALL, TOTAL_READS included:
 *     the two exclusions are not the same as the vendor one;
 *   - THE TWO INCLUDE ARGUMENTS ARE INDEPENDENT, so naming one leaves the other's records out;
 *   - THE Q20 AND Q30 COUNTS ARE OVER FLOWS THAT REACH THE THRESHOLD, and a run whose reads sit
 *     at quality 20 and 40 reports every flow under Q20 and half of them under Q30;
 *   - A QUALITY BETWEEN THE TWO THRESHOLDS COUNTS FOR ONE AND NOT THE OTHER;
 *   - PF_Q20_EQUIVALENT_YIELD IS NOT A COUNT OF ANYTHING: four flows at quality 40 yield 8 and
 *     eight flows at quality 25 yield 10, so it moves with the qualities and not only with the
 *     flows;
 *   - --INCLUDE_BQ_HISTOGRAM ADDS A HISTOGRAM TO THE FILE and changes none of the metrics;
 *   - A READ WITH NO FLOW MATRIX AT ALL IS REFUSED, by a message naming the `tp` attribute;
 *   - AND AN EMPTY FILE PRODUCES A METRICS FILE OF ZEROS rather than none.
 *
 * Output:
 *
 *     sam\t<case>\t<that bam as sam, without its header, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     histogram\t<case>\t<the histogram section, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class QualityYieldFlowDump {

    static final StringBuilder buf = new StringBuilder();
    static final String FLOW_ORDER = "TGCA";

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** name, bases, quality string, secondary, supplementary, vendor-failed, flow matrix. */
    record Read(String name, String bases, String quals, boolean secondary,
                boolean supplementary, boolean vendorFailed, boolean flowMatrix) { }

    static Read plain(final String name, final String bases, final String quals) {
        return new Read(name, bases, quals, false, false, false, true);
    }

    static SAMRecord record(final SAMFileHeader header, final Read spec, final int start) {
        final SAMRecord read = new SAMRecord(header);
        read.setReadName(spec.name());
        read.setReferenceName("chr1");
        read.setAlignmentStart(start);
        read.setCigarString(spec.bases().length() + "M");
        read.setReadString(spec.bases());
        read.setBaseQualityString(spec.quals());
        read.setMappingQuality(60);
        read.setAttribute("RG", "rg1");
        read.setSecondaryAlignment(spec.secondary());
        read.setSupplementaryAlignmentFlag(spec.supplementary());
        read.setReadFailsVendorQualityCheckFlag(spec.vendorFailed());
        if (spec.flowMatrix()) {
            read.setAttribute("tp", new byte[spec.bases().length()]);
        }
        return read;
    }

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("qyf");
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setPlatform("ULTIMA");
        group.setFlowOrder(FLOW_ORDER);
        group.setSample("sample");
        header.addReadGroup(group);
        final File bam = new File(dir.toFile(), "in.bam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setUseAsyncIo(false)
                .makeBAMWriter(header, false, bam)) {
            int start = 100;
            for (final Read spec : reads) {
                writer.addAlignment(record(header, spec, start));
                start += 100;
            }
        }
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam)) {
            for (final SAMRecord read : reader) {
                sam.append(read.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final File out = new File(dir.toFile(), "metrics.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + bam.getAbsolutePath(), "O=" + out.getAbsolutePath()));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.CollectQualityYieldMetricsFlow()
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
            emit("error", name, cause.getClass().getName() + ":" + cause.getMessage());
            return;
        }
        final StringBuilder table = new StringBuilder();
        final StringBuilder histogram = new StringBuilder();
        boolean inHistogram = false;
        for (final String line : Files.readString(out.toPath()).split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
                continue;
            }
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            if (inHistogram) {
                histogram.append(line).append('\n');
            } else {
                table.append(line).append('\n');
            }
        }
        emit("metrics", name, table.toString());
        emit("histogram", name, histogram.toString());
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Eight bases in four homopolymers, so four flows.
        final Read pairs = plain("pairs", "TTGGCCAA", "IIIIIIII");
        // Eight bases, eight flows.
        final Read singles = plain("singles", "TGCATGCA", "IIIIIIII");
        // The same bases at a lower quality.
        final Read low = plain("low", "TGCATGCA", "55555555");
        // A quality that sits between the two thresholds.
        final Read middle = plain("middle", "TGCATGCA", "::::::::");

        run("one-read", List.of(pairs));
        run("homopolymers-and-singles", List.of(pairs, singles));
        run("low-quality", List.of(singles, low));
        run("between-thresholds", List.of(middle));

        // A vendor-failed read.
        run("vendor-failed", List.of(singles,
                new Read("failed", "TGCATGCA", "IIIIIIII", false, false, true, true)));

        // Secondary and supplementary reads, in and out.
        final Read secondary = new Read("secondary", "TGCATGCA", "IIIIIIII", true, false, false,
                true);
        final Read supplementary = new Read("supplementary", "TGCATGCA", "IIIIIIII", false, true,
                false, true);
        run("secondary-out", List.of(singles, secondary));
        run("secondary-in", List.of(singles, secondary), "INCLUDE_SECONDARY_ALIGNMENTS=true");
        run("supplementary-out", List.of(singles, supplementary));
        run("supplementary-in", List.of(singles, supplementary),
                "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true");
        // Naming one leaves the other's records out.
        run("secondary-in-supplementary-out", List.of(singles, secondary, supplementary),
                "INCLUDE_SECONDARY_ALIGNMENTS=true");

        // The histogram.
        run("with-histogram", List.of(pairs, singles), "INCLUDE_BQ_HISTOGRAM=true");

        // A read with no flow matrix at all.
        run("no-flow-matrix", List.of(
                new Read("plain", "TGCATGCA", "IIIIIIII", false, false, false, false)));

        // No reads at all.
        run("empty", List.of());

        System.out.print(buf);
    }
}
