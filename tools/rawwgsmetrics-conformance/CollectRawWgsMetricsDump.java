/*
 * CollectRawWgsMetrics' metrics, taken from the reference.
 *
 * The tool is CollectWgsMetrics with four defaults changed, and nothing else. What is measured is
 * exactly that: the same fixture through both tools, so the difference between the two files is
 * the difference between the two sets of defaults and no more.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE MAPPING-QUALITY FLOOR IS ZERO, so a read at quality five that CollectWgsMetrics excludes
 *     whole is counted here;
 *   - THE BASE-QUALITY FLOOR IS THREE AND NOT TWENTY, so bases at quality two are still excluded
 *     and bases at quality five are not;
 *   - AN `N` BASE IS STILL EXCLUDED, quality zero being under three as it is under twenty;
 *   - THE COVERAGE CAP IS A HUNDRED THOUSAND, so nothing this fixture can reach is capped. Only
 *     the `default-coverage-cap` case leaves the cap alone, and its histogram trailer counts the
 *     zero bins the default leaves behind; every other case names a cap of two hundred and fifty,
 *     because the reference writes one line per bin and a hundred thousand of them per case was
 *     most of what this suite cost to run;
 *   - THE DUPLICATE AND UNPAIRED RULES ARE UNCHANGED, so a duplicate is still excluded whole and
 *     an unpaired read still needs --COUNT_UNPAIRED;
 *   - THE OVERLAP RULE IS UNCHANGED, a pair's overlap still counting once;
 *   - AND THE COLUMNS ARE THE SAME, so the two tools' files differ only in their numbers.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     histogram\t<case>\t<the histogram section, escaped>
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

public class CollectRawWgsMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 200;

    /** The reference: a repeating pattern, with ten Ns at the end so the territory is smaller. */
    static String fasta() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < CONTIG_LENGTH - 10; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        bases.append("NNNNNNNNNN");
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /** One read: where it sits, how long it is, and what its flags and qualities say. */
    record Read(String name, int start, int length, int flags, int mappingQuality,
                String qualities, int mateStart, String bases) {}

    static Read read(final String name, final int start, final int length) {
        return new Read(name, start, length, 0, 60, null, 0, null);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads) {
        final SAMFileHeader header = header();
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(spec.flags());
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setMappingQuality(spec.mappingQuality());
                record.setCigarString(spec.length() + "M");
                final StringBuilder bases = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    bases.append(spec.bases() == null ? 'A' : spec.bases().charAt(i));
                }
                record.setReadString(bases.toString());
                final StringBuilder quals = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    quals.append(spec.qualities() == null ? 'I' : spec.qualities().charAt(i));
                }
                record.setBaseQualityString(quals.toString());
                if ((spec.flags() & 0x1) != 0) {
                    if ((spec.flags() & 0x8) != 0) {
                        record.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                        record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                    } else {
                        record.setMateReferenceName("chr1");
                        record.setMateAlignmentStart(spec.mateStart());
                    }
                }
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /**
     * The metrics table and the histogram, the histogram TRIMMED past its last non-zero bin.
     *
     * The coverage cap here is a hundred thousand, so the tool writes a hundred thousand and one
     * bins of which a handful are non-zero. Keeping them all made a golden of two and a half
     * megabytes and the slowest shard in CI, and every bin past the last non-zero one is zero by
     * construction: the cap itself is still visible in the metrics table.
     */
    static String[] split(final String text) {
        final List<String> table = new ArrayList<>();
        final List<String> histogram = new ArrayList<>();
        boolean inHistogram = false;
        for (final String line : text.split("\n", -1)) {
            if (line.startsWith("## HISTOGRAM")) {
                inHistogram = true;
                continue;
            }
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            (inHistogram ? histogram : table).add(line);
        }
        int last = 0;
        for (int i = 1; i < histogram.size(); i++) {
            final String[] columns = histogram.get(i).split("\t");
            boolean nonZero = false;
            for (int c = 1; c < columns.length; c++) {
                if (!columns[c].equals("0")) {
                    nonZero = true;
                }
            }
            if (nonZero) {
                last = i;
            }
        }
        final List<String> trimmed = new ArrayList<>(histogram.subList(0, last + 1));
        trimmed.add("# " + (histogram.size() - 1 - last) + " further bins, every one of them zero");
        return new String[]{String.join("\n", table), String.join("\n", trimmed)};
    }

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("rawwgsmetrics");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
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
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + metrics, "R=" + reference));
        // The default cap is a hundred thousand, and the reference writes a histogram line for
        // every bin up to it: a hundred thousand lines per case, which is most of what this suite
        // costs. One case leaves the cap alone, so the golden still holds the default in its
        // trailer of zero bins; the rest name a cap of two hundred and fifty, which changes no
        // number here because nothing this fixture can reach is capped either way.
        if (!Arrays.asList(extra).contains("DEFAULT_CAP")) {
            argv.add("COVERAGE_CAP=250");
        }
        for (final String value : extra) {
            if (!value.equals("DEFAULT_CAP")) {
                argv.add(value);
            }
        }
        try {
            final int code = new picard.analysis.CollectRawWgsMetrics()
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
        final String[] parts = split(Files.readString(metrics, StandardCharsets.UTF_8));
        emit("metrics", name, parts[0]);
        emit("histogram", name, parts[1]);
    }

    /** A pair whose two ends are given explicitly, so their overlap can be chosen. */
    static List<Read> pair(final String name, final int first, final int second, final int length) {
        return List.of(
                new Read(name, first, length, 0x1 | 0x2 | 0x40, 60, null, second, null),
                new Read(name, second, length, 0x1 | 0x2 | 0x80, 60, null, first, null));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // A read at mapping quality five, which CollectWgsMetrics excludes whole and this does
        // not: the floor here is zero.
        run("low-mapping-quality", List.of(
                new Read("a", 1, 20, 0, 5, null, 0, null)), "COUNT_UNPAIRED=true");

        // Bases at quality two and at quality five, which straddle the floor of three.
        run("quality-two", List.of(
                new Read("a", 1, 20, 0, 60, "IIIIIIIIII##########", 0, null)),
                "COUNT_UNPAIRED=true");
        run("quality-five", List.of(
                new Read("a", 1, 20, 0, 60, "IIIIIIIIII&&&&&&&&&&", 0, null)),
                "COUNT_UNPAIRED=true");

        // N bases, which quality zero still excludes.
        run("n-bases", List.of(
                new Read("a", 1, 20, 0, 60, null, 0, "AAAAANNNNNAAAAAAAAAA")),
                "COUNT_UNPAIRED=true");

        // The rules that did not change.
        run("one-unpaired-read", List.of(read("a", 1, 20)));
        run("one-unpaired-read-counted", List.of(read("a", 1, 20)), "COUNT_UNPAIRED=true");
        run("duplicate", List.of(
                new Read("a", 1, 20, 0x400, 60, null, 0, null)), "COUNT_UNPAIRED=true");
        run("pair-overlapping", pair("a", 1, 1, 20));
        run("pair-disjoint", pair("a", 1, 50, 20));

        // Ten reads over one span, which the cap of a hundred thousand does not reach.
        final List<Read> deep = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            deep.add(new Read("d" + i, 1, 20, 0, 60, null, 0, null));
        }
        run("deep", deep, "COUNT_UNPAIRED=true");
        // The same reads with the cap left alone, which is what shows the default.
        run("default-coverage-cap", deep, "COUNT_UNPAIRED=true", "DEFAULT_CAP");

        // A file with no reads at all.
        run("empty", List.of());

        System.out.print(buf);
    }
}
