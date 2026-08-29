/*
 * `UmiAwareMarkDuplicatesWithMateCigar`, taken from the reference.
 *
 * A UMI is a barcode on the molecule rather than on the library, so two reads at one position that
 * carry different UMIs came from different molecules and are not duplicates of each other. This
 * tool is `SimpleMarkDuplicatesWithMateCigar` with that split applied inside each duplicate set,
 * and with a metrics file of its own about the UMIs it saw.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE UMI SPLITS A SET: two pairs at one position with different UMIs are both kept;
 *   - BUT ONLY BEYOND AN EDIT DISTANCE: `MAX_EDIT_DISTANCE_TO_JOIN` is one by default, so two UMIs
 *     that differ in a single base are the SAME molecule and one of the two pairs is marked;
 *   - THE ASSIGNED UMI IS WRITTEN BACK: `MOLECULAR_IDENTIFIER_TAG` names a tag that carries the
 *     UMI the set was assigned, which is not always the read's own;
 *   - AN N IS A BASE for this purpose, so a UMI with one joins the set its distance allows;
 *   - A MISSING UMI IS A REFUSAL unless `ALLOW_MISSING_UMIS` says otherwise, and the wording is
 *     the tool's own;
 *   - `UMI_METRICS` IS REQUIRED, so a command line without it is refused by the parser;
 *   - AND THE UMI METRICS THEMSELVES are a second table: observed and inferred counts, an entropy
 *     apiece, and an estimated error rate.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     marked\t<case>\t<the output as sam, without its header, escaped>
 *     metrics\t<case>\t<the duplication metrics table, escaped>
 *     umi\t<case>\t<the umi metrics table, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: UmiDuplicatesDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class UmiDuplicatesDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 4000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        group.setPlatformUnit("unit1");
        header.addReadGroup(group);
        return header;
    }

    /** One pair, both ends carrying the same UMI and each other's cigar. */
    static void pair(final List<SAMRecord> into, final String name, final int first,
                     final int second, final String umi) {
        final SAMFileHeader header = header();
        for (final boolean isFirst : new boolean[]{true, false}) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(name);
            record.setFlags(0x1 | 0x2 | (isFirst ? 0x40 | 0x20 : 0x80 | 0x10));
            record.setReferenceName("chr1");
            record.setAlignmentStart(isFirst ? first : second);
            record.setCigarString("20M");
            record.setMappingQuality(60);
            record.setReadString("A".repeat(20));
            record.setBaseQualityString("I".repeat(20));
            record.setAttribute("RG", "rg1");
            record.setMateReferenceName("chr1");
            record.setMateAlignmentStart(isFirst ? second : first);
            record.setAttribute("MC", "20M");
            if (umi != null) {
                record.setAttribute("RX", umi);
            }
            into.add(record);
        }
    }

    static void writeBam(final Path bam, final List<SAMRecord> records) {
        records.sort((a, b) -> {
            final int byStart = Integer.compare(a.getAlignmentStart(), b.getAlignmentStart());
            return byStart != 0 ? byStart : a.getReadName().compareTo(b.getReadName());
        });
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header(), true, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    /** A metrics file's table, without its comments and without its histogram. */
    static String table(final Path file) throws Exception {
        final List<String> rows = new ArrayList<>();
        boolean inHistogram = false;
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
            }
            if (inHistogram || line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            rows.add(line);
        }
        return String.join("\n", rows);
    }

    static void run(final String name, final List<SAMRecord> records, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("umidup");
        final Path in = dir.resolve("in.bam");
        writeBam(in, records);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final Path out = dir.resolve("out.bam");
        final Path metrics = dir.resolve("metrics.txt");
        final Path umi = dir.resolve("umi.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + out, "M=" + metrics));
        final List<String> tail = Arrays.asList(extra);
        if (!tail.contains("NO_UMI_METRICS")) {
            argv.add("UMI_METRICS=" + umi);
        }
        for (final String argument : tail) {
            if (!argument.equals("NO_UMI_METRICS")) {
                argv.add(argument);
            }
        }
        try {
            final int code = new picard.sam.markduplicates.UmiAwareMarkDuplicatesWithMateCigar()
                    .instanceMain(argv.toArray(new String[0]));
            if (code != 0) {
                emit("error", name, "exit " + code);
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
        final StringBuilder marked = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(out.toFile())) {
            for (final SAMRecord record : reader) {
                marked.append(record.getSAMString());
            }
        }
        emit("marked", name, marked.toString());
        emit("metrics", name, table(metrics));
        if (Files.exists(umi)) {
            emit("umi", name, table(umi));
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Two pairs at one position with UMIs four bases apart: two molecules, nothing marked.
        final List<SAMRecord> distinct = new ArrayList<>();
        pair(distinct, "HWI:1:FC:1:1:1:1", 100, 300, "AAAA");
        pair(distinct, "HWI:1:FC:1:1:9000:9000", 100, 300, "CCCC");
        run("two-umis", distinct);

        // The same two, one base apart: inside the default edit distance, so one is marked.
        final List<SAMRecord> close = new ArrayList<>();
        pair(close, "HWI:1:FC:1:1:1:1", 100, 300, "AAAA");
        pair(close, "HWI:1:FC:1:1:9000:9000", 100, 300, "AAAC");
        run("umis-one-base-apart", close);
        run("umis-one-base-apart-with-no-joining", close, "MAX_EDIT_DISTANCE_TO_JOIN=0");
        run("two-umis-joined-at-four", distinct, "MAX_EDIT_DISTANCE_TO_JOIN=4");

        // The same UMI twice, which is one molecule sequenced twice.
        final List<SAMRecord> same = new ArrayList<>();
        pair(same, "HWI:1:FC:1:1:1:1", 100, 300, "AAAA");
        pair(same, "HWI:1:FC:1:1:9000:9000", 100, 300, "AAAA");
        run("one-umi", same);
        run("one-umi-with-an-assigned-tag", same, "MOLECULAR_IDENTIFIER_TAG=MI");

        // A UMI with an N in it, which is a base like any other for the distance.
        final List<SAMRecord> ambiguous = new ArrayList<>();
        pair(ambiguous, "HWI:1:FC:1:1:1:1", 100, 300, "AAAA");
        pair(ambiguous, "HWI:1:FC:1:1:9000:9000", 100, 300, "AANA");
        run("a-umi-with-an-n", ambiguous);

        // No UMI at all, with and without the argument that allows it.
        final List<SAMRecord> missing = new ArrayList<>();
        pair(missing, "HWI:1:FC:1:1:1:1", 100, 300, null);
        pair(missing, "HWI:1:FC:1:1:9000:9000", 100, 300, null);
        run("no-umis", missing);
        run("no-umis-allowed", missing, "ALLOW_MISSING_UMIS=true");

        // A tag that is not the default one, and a command line with no UMI metrics file.
        final List<SAMRecord> tagged = new ArrayList<>();
        pair(tagged, "HWI:1:FC:1:1:1:1", 100, 300, "AAAA");
        pair(tagged, "HWI:1:FC:1:1:9000:9000", 100, 300, "CCCC");
        run("a-different-umi-tag", tagged, "UMI_TAG_NAME=RX");
        run("no-umi-metrics-file", distinct, "NO_UMI_METRICS");

        System.out.print(buf);
    }
}
