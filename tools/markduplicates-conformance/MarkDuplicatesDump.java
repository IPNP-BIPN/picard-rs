/*
 * MarkDuplicates' flags and its metrics, taken from the reference.
 *
 * The tool decides which records are duplicates of which, marks them with the 0x400 flag, and
 * writes a metrics file counting what it found. It does not change a record otherwise, so the
 * output is the input with one bit moved, and the interesting questions are which record keeps the
 * bit clear and what the counters say.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - DUPLICATES ARE DECIDED BY THE 5' POSITION AND THE ORIENTATION, not by the bases: two reads
 *     with different sequences at the same unclipped start are duplicates of each other;
 *   - SOFT CLIPPING MOVES THAT POSITION: a read clipped at its front starts, for this purpose,
 *     where it would have started unclipped, so a clipped read and an unclipped one can be
 *     duplicates;
 *   - ONE RECORD OF EACH SET KEEPS THE BIT CLEAR, and which one is decided by the SCORING
 *     STRATEGY: the default sums the mapped reference length, and `SUM_OF_BASE_QUALITIES` picks a
 *     different read on the same set;
 *   - A PAIR IS A UNIT: both ends of the losing pair are marked, and a pair is never a duplicate
 *     of a single read;
 *   - AN UNPAIRED READ IS ITS OWN SET, so a file of singles marks duplicates the same way;
 *   - AN UNMAPPED READ IS NEVER A DUPLICATE, and neither is a secondary or supplementary one;
 *   - OPTICAL DUPLICATES ARE A SUBSET OF DUPLICATES, told apart by the tile coordinates in the
 *     read NAME and a pixel distance, so two reads at the same coordinates are optical and two far
 *     apart are not;
 *   - `--READ_NAME_REGEX null` TURNS THAT OFF, and the metrics then report no optical duplicates
 *     at all rather than refusing the file;
 *   - `--REMOVE_DUPLICATES` DROPS THEM instead of marking them, and
 *     `--REMOVE_SEQUENCING_DUPLICATES` drops only the optical ones;
 *   - `--TAGGING_POLICY` WRITES A `DT` TAG saying which kind a duplicate is, and `--CLEAR_DT`
 *     decides whether an incoming one survives;
 *   - `--BARCODE_TAG` MAKES TWO READS AT ONE POSITION TWO SETS, which is what a UMI is for;
 *   - AND THE METRICS COUNT PAIRS AND SINGLES APART, with the library as the row's key.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     marked\t<case>\t<the output as sam, without its header, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     histogram\t<case>\t<the histogram section, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MarkDuplicatesDump
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

public class MarkDuplicatesDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 2000;

    /** One record of the fixture, named the way the optical finder expects. */
    record Read(String name, int start, String cigar, String bases, String qualities, int flags,
                int mateStart, String barcode, String existingDt) {}

    static Read single(final String name, final int start, final String bases) {
        return new Read(name, start, bases.length() + "M", bases, "I".repeat(bases.length()), 0, 0,
                null, null);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        group.setPlatformUnit("unit1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = header();
        final List<SAMRecord> records = new ArrayList<>();
        for (final Read spec : reads) {
            final SAMRecord record = new SAMRecord(header);
            record.setReadName(spec.name());
            record.setFlags(spec.flags());
            if ((spec.flags() & 0x4) == 0) {
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setCigarString(spec.cigar());
                record.setMappingQuality(60);
            } else {
                record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                record.setMappingQuality(0);
            }
            record.setReadString(spec.bases());
            record.setBaseQualityString(spec.qualities());
            record.setAttribute("RG", "rg1");
            if ((spec.flags() & 0x1) != 0 && (spec.flags() & 0x8) == 0) {
                record.setMateReferenceName("chr1");
                record.setMateAlignmentStart(spec.mateStart());
                record.setAttribute("MC", spec.cigar());
            }
            if (spec.barcode() != null) {
                record.setAttribute("RX", spec.barcode());
            }
            if (spec.existingDt() != null) {
                record.setAttribute("DT", spec.existingDt());
            }
            records.add(record);
        }
        records.sort((a, b) -> {
            final int byStart = Integer.compare(
                    a.getReadUnmappedFlag() ? Integer.MAX_VALUE : a.getAlignmentStart(),
                    b.getReadUnmappedFlag() ? Integer.MAX_VALUE : b.getAlignmentStart());
            return byStart != 0 ? byStart : a.getReadName().compareTo(b.getReadName());
        });
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header(), true, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    /** A pair whose two ends sit where the case asks. */
    static List<Read> pair(final String name, final int first, final int second, final int length,
                           final String barcode) {
        final String bases = "A".repeat(length);
        return List.of(
                new Read(name, first, length + "M", bases, "I".repeat(length),
                        0x1 | 0x2 | 0x40 | 0x20, second, barcode, null),
                new Read(name, second, length + "M", bases, "I".repeat(length),
                        0x1 | 0x2 | 0x80 | 0x10, first, barcode, null));
    }

    /** The metrics table and the histogram under it. */
    static String[] split(final Path file) throws Exception {
        final List<String> table = new ArrayList<>();
        final List<String> histogram = new ArrayList<>();
        boolean inHistogram = false;
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
                continue;
            }
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            (inHistogram ? histogram : table).add(line);
        }
        return new String[]{String.join("\n", table), String.join("\n", histogram)};
    }

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("markduplicates");
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
                "I=" + in, "O=" + out, "M=" + metrics));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.sam.markduplicates.MarkDuplicates()
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
        final String[] parts = split(metrics);
        emit("metrics", name, parts[0]);
        emit("histogram", name, parts[1]);
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Two pairs at the same position: one of them keeps the bit clear.
        final List<Read> twoPairs = new ArrayList<>();
        twoPairs.addAll(pair("HWI:1:FC:1:1:1:1", 100, 300, 20, null));
        twoPairs.addAll(pair("HWI:1:FC:1:1:9000:9000", 100, 300, 20, null));
        run("two-pairs", twoPairs);
        // The same two with their tile coordinates close together, which makes them optical.
        final List<Read> optical = new ArrayList<>();
        optical.addAll(pair("HWI:1:FC:1:1:1:1", 100, 300, 20, null));
        optical.addAll(pair("HWI:1:FC:1:1:5:5", 100, 300, 20, null));
        run("optical-duplicates", optical);
        run("optical-duplicates-without-the-regex", optical, "READ_NAME_REGEX=null");
        run("optical-duplicates-with-a-smaller-distance", optical,
                "OPTICAL_DUPLICATE_PIXEL_DISTANCE=1");

        // A pair whose ends differ in length, which changes the default score.
        final List<Read> different = new ArrayList<>();
        different.addAll(pair("HWI:1:FC:1:1:1:1", 100, 300, 20, null));
        different.addAll(pair("HWI:1:FC:1:1:9000:9000", 100, 300, 30, null));
        run("pairs-of-different-lengths", different);
        run("pairs-of-different-lengths-by-quality", different,
                "DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES");

        // Unpaired reads at one position, and a soft-clipped one that starts where they do.
        run("three-singles", List.of(
                single("HWI:1:FC:1:1:1:1", 100, "ACGTACGTAC"),
                single("HWI:1:FC:1:1:2:2", 100, "TTTTTTTTTT"),
                single("HWI:1:FC:1:1:9000:9000", 100, "GGGGGGGGGG")));
        run("a-soft-clipped-read", List.of(
                single("HWI:1:FC:1:1:1:1", 105, "ACGTACGTAC"),
                new Read("HWI:1:FC:1:1:2:2", 105, "5S10M", "TTTTTTTTTTTTTTT",
                        "I".repeat(15), 0, 0, null, null)));

        // An unmapped read, a secondary one and a supplementary one, none of which is ever marked.
        run("an-unmapped-read", List.of(
                single("HWI:1:FC:1:1:1:1", 100, "ACGTACGTAC"),
                single("HWI:1:FC:1:1:2:2", 100, "ACGTACGTAC"),
                new Read("HWI:1:FC:1:1:3:3", 0, "*", "ACGTACGTAC", "IIIIIIIIII", 0x4, 0, null,
                        null)));
        run("a-secondary-read", List.of(
                single("HWI:1:FC:1:1:1:1", 100, "ACGTACGTAC"),
                new Read("HWI:1:FC:1:1:2:2", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", 0x100, 0,
                        null, null)));

        // Removing rather than marking.
        run("remove-duplicates", twoPairs, "REMOVE_DUPLICATES=true");
        run("remove-sequencing-duplicates", optical, "REMOVE_SEQUENCING_DUPLICATES=true");

        // The DT tag, and what happens to one that was already there.
        run("tagging-policy-all", twoPairs, "TAGGING_POLICY=All");
        run("tagging-policy-optical", optical, "TAGGING_POLICY=OpticalOnly");
        run("an-existing-dt-tag", List.of(
                new Read("HWI:1:FC:1:1:1:1", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", 0, 0, null,
                        "SQ"),
                single("HWI:1:FC:1:1:2:2", 100, "ACGTACGTAC")));
        run("an-existing-dt-tag-kept", List.of(
                new Read("HWI:1:FC:1:1:1:1", 100, "10M", "ACGTACGTAC", "IIIIIIIIII", 0, 0, null,
                        "SQ"),
                single("HWI:1:FC:1:1:2:2", 100, "ACGTACGTAC")), "CLEAR_DT=false");

        // A barcode, which splits one position into two sets.
        final List<Read> barcoded = new ArrayList<>();
        barcoded.addAll(pair("HWI:1:FC:1:1:1:1", 100, 300, 20, "AAAA"));
        barcoded.addAll(pair("HWI:1:FC:1:1:9000:9000", 100, 300, 20, "CCCC"));
        run("two-barcodes", barcoded, "BARCODE_TAG=RX");
        run("two-barcodes-ignored", barcoded);

        System.out.print(buf);
    }
}
