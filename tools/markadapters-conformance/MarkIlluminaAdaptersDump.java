/*
 * MarkIlluminaAdapters' XT tag and its histogram, taken from the reference.
 *
 * The tool does not clip: it MARKS, writing an `XT` tag whose value is the ONE-BASED position at
 * which the adapter starts, and a histogram of how many bases each read would lose. Everything
 * else it writes is the input's own record.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE TAG IS ONE-BASED, `index + 1`, so a read whose adapter starts at its first base carries
 *     `XT:i:1` and not `XT:i:0`;
 *   - THE SEARCH RUNS FROM THE END OF THE READ TOWARDS ITS FRONT, taking the LAST start that
 *     matches rather than the first, which is a different position when an adapter's prefix
 *     repeats;
 *   - A MATCH IS ONLY CONSIDERED FROM `read.length - minMatch`, so a read shorter than the minimum
 *     is never marked whatever it contains;
 *   - THE ERROR ALLOWANCE IS `(int)(length * maxErrorRate)`, TRUNCATED, so ten bases at a tenth
 *     allow ONE mismatch and nine bases allow none;
 *   - THE LENGTH IT IS COMPUTED FROM IS THE OVERLAP and not the adapter's length, so a match near
 *     the end of a read allows fewer mismatches than one further in;
 *   - AN `N` IN THE ADAPTER MATCHES ANYTHING, which is how the indexed adapters' eight Ns behave;
 *   - THE DEFAULT ADAPTER LIST IS THREE OF THE NINE the enum declares, tried in order, and the
 *     first that matches wins;
 *   - A PAIRED RUN USES A DIFFERENT MINIMUM AND A DIFFERENT ERROR RATE from a single one, and the
 *     two reads of a pair are marked together or not at all;
 *   - `--FIVE_PRIME_ADAPTER` AND `--THREE_PRIME_ADAPTER` REPLACE THE LIST rather than adding to
 *     it, and they are only read when the adapter list names `CUSTOM`;
 *   - THE HISTOGRAM COUNTS BASES CLIPPED per read, so a read marked at position one contributes
 *     its whole length;
 *   - THE OUTPUT RECORD IS THE INPUT'S with the tag added, so nothing else about it moves;
 *   - AND A READ THAT IS ALREADY MARKED LOSES ITS TAG when this run finds no adapter: the tag is
 *     SET from the search's answer rather than merged with what was there, so a file marked twice
 *     with different adapter lists carries the second run's answer and not the union.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     marked\t<case>\t<the output as sam, without its header, escaped>
 *     metrics\t<case>\t<the histogram without its comments, escaped>
 *     error\t<case>\t<the reason, as the run reported it>
 *
 * Usage: MarkIlluminaAdaptersDump
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

public class MarkIlluminaAdaptersDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** The three-prime adapter of the default list's first entry, which is what a read carries. */
    static final String INDEXED_THREE_PRIME =
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";

    /** One read: its name, its bases, and whether it is the first of a pair. */
    record Read(String name, String bases, Integer existingTag, int flags) {}

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 1000));
        header.setSequenceDictionary(dictionary);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setReadString(spec.bases());
                record.setBaseQualityString("I".repeat(spec.bases().length()));
                record.setFlags(spec.flags() | 0x4);
                record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                if (spec.existingTag() != null) {
                    record.setAttribute("XT", spec.existingTag());
                }
                writer.addAlignment(record);
            }
        }
    }

    /** A read of `length` bases whose adapter begins at the given ZERO-based offset. */
    static String withAdapterAt(final int length, final int offset, final String adapter) {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < length; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        for (int i = 0; i < adapter.length() && offset + i < length; i++) {
            final char base = adapter.charAt(i);
            // An `N` in the adapter matches anything, so the read carries a real base there.
            bases.setCharAt(offset + i, base == 'N' ? 'A' : base);
        }
        return bases.toString();
    }

    static String table(final Path file) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final List<Read> reads, final boolean paired,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("markadapters");
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final Path out = dir.resolve("out.bam");
        final Path metrics = dir.resolve("metrics.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + out, "M=" + metrics, "PAIRED_RUN=" + paired));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.illumina.MarkIlluminaAdapters()
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
    }

    static Read single(final String name, final String bases) {
        return new Read(name, bases, null, 0);
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // An adapter at the very start of the read, which is the one-based tag's clearest case.
        run("adapter-at-the-first-base", List.of(
                single("a", withAdapterAt(60, 0, INDEXED_THREE_PRIME))), false);
        // And one further in.
        run("adapter-halfway", List.of(
                single("a", withAdapterAt(60, 30, INDEXED_THREE_PRIME))), false);
        // A read with no adapter at all.
        run("no-adapter", List.of(single("a", "ACGT".repeat(15))), false);

        // The minimum match: an adapter shorter than it is never found.
        run("eleven-bases-of-adapter", List.of(
                single("a", withAdapterAt(60, 49, INDEXED_THREE_PRIME))), false);
        run("twelve-bases-of-adapter", List.of(
                single("a", withAdapterAt(60, 48, INDEXED_THREE_PRIME))), false);
        run("eleven-bases-with-a-lower-minimum", List.of(
                single("a", withAdapterAt(60, 49, INDEXED_THREE_PRIME))), false,
                "MIN_MATCH_BASES_SE=11");

        // The error allowance, which is truncated: a twelve-base overlap allows one mismatch.
        final String twelve = withAdapterAt(60, 48, INDEXED_THREE_PRIME);
        final StringBuilder oneOff = new StringBuilder(twelve);
        oneOff.setCharAt(50, oneOff.charAt(50) == 'A' ? 'C' : 'A');
        run("one-mismatch-in-twelve", List.of(single("a", oneOff.toString())), false);
        final StringBuilder twoOff = new StringBuilder(oneOff);
        twoOff.setCharAt(52, twoOff.charAt(52) == 'A' ? 'C' : 'A');
        run("two-mismatches-in-twelve", List.of(single("a", twoOff.toString())), false);
        run("two-mismatches-with-a-wider-rate", List.of(single("a", twoOff.toString())), false,
                "MAX_ERROR_RATE_SE=0.2");

        // A pair, whose two reads are marked together, and whose minimum is a different one.
        run("a-pair", List.of(
                new Read("p", withAdapterAt(60, 30, INDEXED_THREE_PRIME), null, 0x1 | 0x40 | 0x8),
                new Read("p", withAdapterAt(60, 30, INDEXED_THREE_PRIME), null, 0x1 | 0x80 | 0x8)),
                true);
        run("a-pair-with-one-adapter", List.of(
                new Read("p", withAdapterAt(60, 30, INDEXED_THREE_PRIME), null, 0x1 | 0x40 | 0x8),
                new Read("p", "ACGT".repeat(15), null, 0x1 | 0x80 | 0x8)), true);

        // The adapter list, and the custom pair that replaces it.
        run("one-adapter-named", List.of(
                single("a", withAdapterAt(60, 30, INDEXED_THREE_PRIME))), false,
                "ADAPTERS=null", "ADAPTERS=PAIRED_END");
        run("a-custom-adapter", List.of(
                single("a", withAdapterAt(60, 30, "TTTTTTTTTTTTTTTTTTTT"))), false,
                "ADAPTERS=null",
                "FIVE_PRIME_ADAPTER=AAAAAAAAAAAAAAAAAAAA",
                "THREE_PRIME_ADAPTER=TTTTTTTTTTTTTTTTTTTT");

        // A read that already carries a tag.
        run("an-existing-tag", List.of(
                new Read("a", "ACGT".repeat(15), 7, 0)), false);
        run("an-existing-tag-and-an-adapter", List.of(
                new Read("a", withAdapterAt(60, 30, INDEXED_THREE_PRIME), 7, 0)), false);

        System.out.print(buf);
    }
}
