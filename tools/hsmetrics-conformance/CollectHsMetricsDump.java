/*
 * CollectHsMetrics' three files, taken from the reference.
 *
 * The tool counts a hybrid-selection experiment against two interval lists: the BAITS that were
 * fished with and the TARGETS that were wanted. What is measured is which read reaches which
 * counter, what the per-target and per-base files carry, and which arguments move any of it.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE BAITS AND THE TARGETS ARE COUNTED SEPARATELY, and a read may be on one and not the
 *     other: ON_BAIT_BASES and ON_TARGET_BASES are different numbers on the same file;
 *   - NEAR_BAIT_BASES IS A WINDOW AROUND THE BAITS and --NEAR_DISTANCE moves it, so the same read
 *     is near bait at the default and off bait at nought;
 *   - THE THREE BAIT COLUMNS PARTITION THE ALIGNED BASES, so ON plus NEAR plus OFF is
 *     PF_BASES_ALIGNED;
 *   - --MINIMUM_MAPPING_QUALITY DROPS WHOLE READS FROM COVERAGE while leaving them in the bait
 *     counts: a floor of sixty-one takes ON_TARGET_BASES from a hundred to nought and leaves
 *     ON_BAIT_BASES at sixty. Its default is ONE, so a read at mapping quality nought is already
 *     out of the coverage before any argument is given;
 *   - --MINIMUM_BASE_QUALITY DROPS SINGLE BASES the same way, a floor of forty-one emptying the
 *     target coverage of reads whose qualities are all forty, and its default is NOUGHT, which is
 *     what makes this tool count bases CollectWgsMetrics would not;
 *   - --CLIP_OVERLAPPING_READS CHANGES NOTHING HERE, and the golden says so rather than assuming
 *     it: a pair whose ends overlap by twenty-five bases reports the same thirty-five on-target
 *     bases and the same per-target file with the argument and without it, because the coverage is
 *     counted per LOCUS and the overlap was never counted twice to begin with;
 *   - --PER_TARGET_COVERAGE WRITES A ROW PER TARGET with its own mean and its GC, and the rows are
 *     the target list's own order;
 *   - --PER_BASE_COVERAGE WRITES A ROW PER TARGET BASE, so a fifty-base target is fifty rows;
 *   - THE BAIT SET NAME IS THE BAIT FILE'S NAME WITHOUT ITS EXTENSION unless --BAIT_SET_NAME says
 *     otherwise;
 *   - A TARGET WITH NO READ OVER IT IS STILL A ROW, at zero, and it still counts towards
 *     TARGET_TERRITORY and the ZERO_CVG columns;
 *   - A TARGET FILE ROW CARRIES ITS OWN GC AND ITS NORMALISED COVERAGE, so the two covered
 *     targets read 0.833333 and 1.5 while the uncovered one reads nought and a `pct_0x` of one;
 *   - AND THE TWO INTERVAL LISTS MAY OVERLAP OR NOT: a target outside every bait is counted in the
 *     target columns and in none of the bait ones.
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

public class CollectHsMetricsDump {

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
        emit("intervals", "baits", baits);
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
        run("bait-set-name", baits, targets, reads, false, false, "N=my-bait-set");
        run("mapping-quality-floor", baits, targets, reads, false, false, "MQ=61");
        run("base-quality-floor", baits, targets, reads, false, false, "Q=41");
        run("clip-overlapping", baits, targets, reads, true, false,
                "CLIP_OVERLAPPING_READS=true");

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
        final Path baitPath = dir.resolve("baits.interval_list");
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
                "BI=" + baitPath, "TI=" + targetPath));
        if (perTarget) {
            argv.add("PER_TARGET_COVERAGE=" + targetCoverage);
        }
        if (perBase) {
            argv.add("PER_BASE_COVERAGE=" + baseCoverage);
        }
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.directed.CollectHsMetrics()
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
