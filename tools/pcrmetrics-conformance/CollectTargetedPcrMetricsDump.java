/*
 * CollectTargetedPcrMetrics' three files, taken from the reference.
 *
 * The tool is CollectHsMetrics' collector under another set of column names: the baits are
 * AMPLICONS and the metrics say so, while the arithmetic underneath is the same. What is measured
 * is which columns the two tools differ in, and that the numbers agree wherever they do not.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE COLUMNS ARE AMPLICON'S AND NOT BAIT'S: ON_AMPLICON_BASES, NEAR_AMPLICON_BASES,
 *     OFF_AMPLICON_BASES, PCT_AMPLIFIED_BASES and ON_AMPLICON_VS_SELECTED, with no BAIT column at
 *     all;
 *   - THE AMPLICON ARITHMETIC IS THE SAME: the three amplicon columns partition the aligned bases
 *     exactly as the bait ones do, and the same fixture answers the same sixty, hundred and twenty
 *     and sixty under the other names;
 *   - AND THE TARGET ARITHMETIC IS NOT, because of one line in a constructor: CollectHsMetrics
 *     sets CLIP_OVERLAPPING_READS to true and this tool leaves the shared default of false, so on
 *     the same fixture the overlap of a pair is counted twice here and once there: a hundred and
 *     twenty on-target bases against a hundred, and a per-target mean of one against 0.833333;
 *   - --AMPLICON_INTERVALS AND --TARGET_INTERVALS ARE COUNTED SEPARATELY, so a target outside every
 *     amplicon is counted in the target columns and in none of the amplicon ones;
 *   - --NEAR_DISTANCE MOVES THE WINDOW, and at nought the near bases become off-amplicon ones;
 *   - --CUSTOM_AMPLICON_SET_NAME NAMES THE SET, and without it the name is the file's;
 *   - --MINIMUM_MAPPING_QUALITY AND --MINIMUM_BASE_QUALITY EMPTY THE COVERAGE and leave the
 *     amplicon counts alone;
 *   - --PER_TARGET_COVERAGE AND --PER_BASE_COVERAGE WRITE THE SAME TWO FILES the other tool writes;
 *   - A FILE WITH NO READS STILL REPORTS ITS TERRITORIES;
 *   - AND THE MAPPING-QUALITY DEFAULT IS ONE HERE TOO, so a pair at quality nought is out of the
 *     coverage before any argument is given.
 *
 * Output:
 *
 *     intervals\t<name>\t<that interval list, escaped>
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     per-target\t<case>\t<the per-target file, escaped, or `absent`>
 *     per-base\t<case>\t<the per-base file's first rows and its length, escaped, or `absent`>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CollectTargetedPcrMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 1000;

    /** The reference: a repeating pattern, so the per-target GC is not all one number. */
    static String fasta() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < CONTIG_LENGTH; i++) {
            bases.append("ACGTGGCCATAT".charAt(i % 12));
        }
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /** An interval list, in Picard's own format. */
    static String intervals(final List<int[]> spans, final List<String> names) {
        final List<String> lines = new ArrayList<>(List.of(
                "@HD\tVN:1.6\tSO:coordinate",
                "@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH));
        for (int i = 0; i < spans.size(); i++) {
            lines.add("chr1\t" + spans.get(i)[0] + "\t" + spans.get(i)[1] + "\t+\t" + names.get(i));
        }
        lines.add("");
        return String.join("\n", lines);
    }

    /** One read: where it sits, which end it is, and what its flags say. */
    record Read(String name, int start, int length, boolean first, int mappingQuality,
                int insertSize, String qualities) {}

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader head = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(head, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(head);
                record.setReadName(spec.name());
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setCigarString(spec.length() + "M");
                record.setMappingQuality(spec.mappingQuality());
                record.setReadPairedFlag(true);
                record.setProperPairFlag(true);
                record.setFirstOfPairFlag(spec.first());
                record.setSecondOfPairFlag(!spec.first());
                record.setReadNegativeStrandFlag(!spec.first());
                record.setMateNegativeStrandFlag(spec.first());
                record.setMateReferenceName("chr1");
                record.setMateAlignmentStart(spec.first()
                        ? spec.start() + spec.insertSize() - spec.length()
                        : spec.start() - spec.insertSize() + spec.length());
                record.setInferredInsertSize(spec.first() ? spec.insertSize() : -spec.insertSize());
                final StringBuilder bases = new StringBuilder();
                final StringBuilder quals = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    bases.append("ACGTGGCCATAT".charAt((spec.start() - 1 + i) % 12));
                    quals.append(spec.qualities() == null ? 'I' : spec.qualities().charAt(i));
                }
                record.setReadString(bases.toString());
                record.setBaseQualityString(quals.toString());
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /** A pair whose two ends are placed explicitly, so their overlap can be chosen. */
    static List<Read> pair(final String name, final int first, final int second, final int length) {
        final int insert = second + length - first;
        return List.of(
                new Read(name, first, length, true, 60, insert, null),
                new Read(name, second, length, false, 60, insert, null));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Three baits and three targets: two that agree, and one target outside every bait.
        final String baits = intervals(
                List.of(new int[]{101, 200}, new int[]{301, 400}),
                List.of("bait-a", "bait-b"));
        final String targets = intervals(
                List.of(new int[]{121, 180}, new int[]{321, 380}, new int[]{701, 760}),
                List.of("target-a", "target-b", "target-orphan"));
        emit("intervals", "amplicons", baits);
        emit("intervals", "targets", targets);

        // A pair over the first target, a pair over the orphan, and a pair over neither.
        final List<Read> reads = new ArrayList<>();
        reads.addAll(pair("on-target", 130, 150, 30));
        reads.addAll(pair("on-orphan", 705, 725, 30));
        reads.addAll(pair("off-bait", 500, 540, 30));
        // A pair just outside a bait, which the near window catches.
        reads.addAll(pair("near-bait", 210, 230, 30));

        run("plain", baits, targets, reads, false, false);
        run("per-target", baits, targets, reads, true, false);
        run("per-base", baits, targets, reads, false, true);
        run("near-distance-zero", baits, targets, reads, false, false, "NEAR_DISTANCE=0");
        run("amplicon-set-name", baits, targets, reads, false, false, "N=my-amplicon-set");
        run("mapping-quality-floor", baits, targets, reads, false, false, "MQ=61");
        run("base-quality-floor", baits, targets, reads, false, false, "Q=41");
        

        // A pair whose two ends overlap entirely, which is what the clipping argument is for.
        final List<Read> overlapping = new ArrayList<>(pair("overlap", 130, 135, 30));
        run("overlapping-pair", baits, targets, overlapping, true, false);
        run("overlapping-pair-clipped", baits, targets, overlapping, true, false,
                "CLIP_OVERLAPPING_READS=true");

        // No reads at all, and a read below the mapping-quality default of one.
        run("no-reads", baits, targets, List.of(), true, false);
        run("mapping-quality-zero", baits, targets, List.of(
                new Read("q0", 130, 30, true, 0, 60, null),
                new Read("q0", 160, 30, false, 0, 60, null)), true, false);

        System.out.print(buf);
    }

    static void run(final String name, final String baits, final String targets,
                    final List<Read> reads, final boolean perTarget, final boolean perBase,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("hsmetrics");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path baitPath = dir.resolve("amplicons.interval_list");
        final Path targetPath = dir.resolve("targets.interval_list");
        Files.writeString(baitPath, baits, StandardCharsets.UTF_8);
        Files.writeString(targetPath, targets, StandardCharsets.UTF_8);
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final Path metrics = dir.resolve("out.txt");
        final Path targetCoverage = dir.resolve("per-target.txt");
        final Path baseCoverage = dir.resolve("per-base.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + metrics, "R=" + reference,
                "AI=" + baitPath, "TI=" + targetPath));
        if (perTarget) {
            argv.add("PER_TARGET_COVERAGE=" + targetCoverage);
        }
        if (perBase) {
            argv.add("PER_BASE_COVERAGE=" + baseCoverage);
        }
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.directed.CollectTargetedPcrMetrics()
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
        emit("metrics", name, table(metrics, dir));
        emit("per-target", name, Files.exists(targetCoverage) ? table(targetCoverage, dir)
                : "absent");
        if (Files.exists(baseCoverage)) {
            // A row per target base: the golden holds the first few and the count.
            final List<String> rows = new ArrayList<>();
            for (final String line : Files.readString(baseCoverage, StandardCharsets.UTF_8)
                    .split("\n", -1)) {
                if (!line.startsWith("#") && !line.isEmpty()) {
                    rows.add(line);
                }
            }
            final List<String> kept = new ArrayList<>(rows.subList(0, Math.min(6, rows.size())));
            kept.add("# " + (rows.size() - 1) + " rows in all");
            emit("per-base", name, String.join("\n", kept));
        } else {
            emit("per-base", name, "absent");
        }
    }

    /** A metrics file without its comment lines, which carry the command line and the clock. */
    static String table(final Path file, final Path dir) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (!line.startsWith("#") && !line.isEmpty()) {
                kept.add(line.replace(dir.toString(), "<dir>"));
            }
        }
        return String.join("\n", kept);
    }
}
