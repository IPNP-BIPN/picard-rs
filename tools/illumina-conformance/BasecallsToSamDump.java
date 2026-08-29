/*
 * `IlluminaBasecallsToSam`, taken from the reference.
 *
 * The same walk as `IlluminaBasecallsToFastq` with a different destination: an unaligned BAM,
 * whose records carry the run's identity in a read group rather than in a read name.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A CLUSTER BECOMES A RECORD PER TEMPLATE READ, so `2T2T` is two records a cluster and they
 *     are a PAIR: the flags say first and second of pair, which the FASTQ shape has no way to say;
 *   - `4T` IS ONE RECORD and it is NOT paired;
 *   - THE FILTER FILE BECOMES THE VENDOR CHECK FLAG rather than a dropped record, so the failing
 *     cluster is `0x200` and is still there;
 *   - `--INCLUDE_NON_PF_READS=false` IS WHAT DROPS IT;
 *   - A BARCODE SEGMENT NEEDS THE PARAMS FILE: the deprecated single-`--OUTPUT` form is declined
 *     for a read structure carrying one, because a barcode is only meaningful once something says
 *     which barcode goes where;
 *   - `--LIBRARY_PARAMS` SPLITS THE OUTPUT BY BARCODE, one file per row, and names the sample and
 *     the library of each;
 *   - AND THE READ GROUP CARRIES THE RUN BARCODE, THE LANE, THE PLATFORM AND THE CENTRE, so what
 *     the FASTQ put in every read name is here written once.
 *
 * Output:
 *
 *     files\t<case>\t<the files written, sorted, space separated>
 *     header\t<case>.<file>\t<the @RG lines, escaped>
 *     sam\t<case>.<file>\t<the records as sam text, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: BasecallsToSamDump
 */

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;

public class BasecallsToSamDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static void run(final String name, final String libraryParams, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("basecallstosam");
        final Path basecalls = MakeIlluminaFixtures.write(dir.resolve("run"));
        final Path out = Files.createDirectories(dir.resolve("out"));

        final List<String> argv = new ArrayList<>(List.of(
                "BASECALLS_DIR=" + basecalls, "LANE=1", "RUN_BARCODE=run17",
                "SEQUENCING_CENTER=centre", "PLATFORM=ILLUMINA", "SORT=false"));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("READ_STRUCTURE="))) {
            tail.add("READ_STRUCTURE=4T");
        }
        if (libraryParams == null) {
            tail.add("OUTPUT=" + out.resolve("reads.bam"));
            tail.add("SAMPLE_ALIAS=sample1");
            tail.add("LIBRARY_NAME=libraryA");
        } else {
            // The barcode each cluster was assigned is `ExtractIlluminaBarcodes`' answer, so a
            // multiplexed run is the reference's own pipeline: extract, then convert.
            final Path barcodes = dir.resolve("barcodes.tsv");
            Files.writeString(barcodes,
                    "barcode_sequence_1\tbarcode_name\tlibrary_name\nAG\tfirst\tlibraryA\n"
                            + "CT\tsecond\tlibraryB\n", StandardCharsets.UTF_8);
            new picard.illumina.ExtractIlluminaBarcodes().instanceMain(new String[]{
                    "BASECALLS_DIR=" + basecalls, "LANE=1", "READ_STRUCTURE=2T2B",
                    "BARCODE_FILE=" + barcodes,
                    "METRICS_FILE=" + dir.resolve("barcode-metrics.txt"),
                    "OUTPUT_DIR=" + basecalls});
            final Path params = dir.resolve("library.params");
            Files.writeString(params, libraryParams.replace("<out>", out.toString()),
                    StandardCharsets.UTF_8);
            tail.add("LIBRARY_PARAMS=" + params);
            tail.add("BARCODES_DIR=" + basecalls);
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        // A command line the tool declines is declined with a status and a printed reason rather
        // than an exception, so both streams are kept and the `ERROR:` line is read out of them.
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.illumina.IlluminaBasecallsToSam()
                    .instanceMain(argv.toArray(new String[0]));
            if (code != 0) {
                System.setOut(original);
                System.setErr(originalError);
                final List<String> reasons = new ArrayList<>();
                for (final String line : said.toString(StandardCharsets.UTF_8).split("\n", -1)) {
                    if (line.startsWith("ERROR:")) {
                        reasons.add(line);
                    }
                }
                emit("error", name, "exit " + code
                        + (reasons.isEmpty() ? "" : " " + String.join(" ", reasons)));
                return;
            }
        } catch (final Exception e) {
            System.setOut(original);
            System.setErr(originalError);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        } finally {
            System.setOut(original);
            System.setErr(originalError);
        }

        final List<String> written = new ArrayList<>();
        try (final Stream<Path> walk = Files.walk(out)) {
            walk.filter(Files::isRegularFile)
                    .forEach(path -> written.add(path.getFileName().toString()));
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            try (final SamReader reader = SamReaderFactory.makeDefault()
                    .open(out.resolve(file).toFile())) {
                final StringBuilder groups = new StringBuilder();
                for (final SAMReadGroupRecord group : reader.getFileHeader().getReadGroups()) {
                    groups.append(group.getSAMString());
                }
                emit("header", name + "." + file, groups.toString());
                final StringBuilder sam = new StringBuilder();
                for (final SAMRecord record : reader) {
                    sam.append(record.getSAMString());
                }
                emit("sam", name + "." + file, sam.toString());
            }
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // One record a cluster, then two that are a pair.
        run("one-read-of-four", null);
        run("two-reads-of-two", null, "READ_STRUCTURE=2T2T");
        // A barcode segment with a single OUTPUT is declined: a barcode read is only meaningful
        // once the params file says which barcode goes where, so the deprecated single-output form
        // has nothing to do with it.
        run("a-barcode-segment-without-params", null, "READ_STRUCTURE=2T2B");
        run("a-skipped-segment", null, "READ_STRUCTURE=2T2S");
        // The cluster that failed the filter: flagged, then dropped.
        run("without-the-non-pf-reads", null, "INCLUDE_NON_PF_READS=false");
        // The barcode in the read group's own identifier, which needs the params file for the
        // same reason.
        run("the-barcode-in-the-read-group",
                "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\n"
                        + "AG\t<out>/first.bam\tsampleA\tlibraryA\n"
                        + "CT\t<out>/second.bam\tsampleB\tlibraryB\n",
                "READ_STRUCTURE=2T2B", "INCLUDE_BC_IN_RG_TAG=true");

        // Split by barcode, with a sample and a library apiece.
        run("split-by-barcode",
                "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\n"
                        + "AG\t<out>/first.bam\tsampleA\tlibraryA\n"
                        + "CT\t<out>/second.bam\tsampleB\tlibraryB\n",
                "READ_STRUCTURE=2T2B");
        // One barcode declared, so the other clusters have nowhere to go.
        run("an-undeclared-barcode",
                "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\n"
                        + "AG\t<out>/first.bam\tsampleA\tlibraryA\n",
                "READ_STRUCTURE=2T2B");
        run("an-undeclared-barcode-ignored",
                "BARCODE_1\tOUTPUT\tSAMPLE_ALIAS\tLIBRARY_NAME\n"
                        + "AG\t<out>/first.bam\tsampleA\tlibraryA\n",
                "READ_STRUCTURE=2T2B", "IGNORE_UNEXPECTED_BARCODES=true");

        System.out.print(buf);
    }
}
