/*
 * CollectSequencingArtifactMetrics' five files, taken from the reference.
 *
 * The tool counts, for every reference base and every context, how the reads that cover it
 * disagree with it, and splits those disagreements by the read's end and its strand: a pre-adapter
 * artifact is one whose direction follows the READ, a bait-bias artifact one whose direction
 * follows the REFERENCE STRAND. What is measured is which read reaches which counter, what the
 * five files carry, and which arguments change any of it.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE OUTPUT ARGUMENT IS A PREFIX FOR FIVE FILES, whose extensions are the metric class's own
 *     and not the tool's, and --FILE_EXTENSION appends to all five rather than replacing anything;
 *   - THE DETAIL FILES HOLD ONE ROW PER SUBSTITUTION PER CONTEXT, so a context size of one gives
 *     ninety-six rows per library and a context size of nought gives six;
 *   - --CONTEXTS_TO_PRINT CUTS THE TWO DETAIL FILES AND THE ERROR SUMMARY, and not the counting:
 *     naming one context takes the detail files from a hundred and ninety-two rows to three and
 *     the error summary from six to three, while the twelve summary rows come out unchanged;
 *   - THE SUMMARY FILES ARE THE WORST CONTEXT PER SUBSTITUTION, so twelve rows survive whatever
 *     the context size is;
 *   - A PRE-ADAPTER ARTIFACT FOLLOWS THE READ AND A BAIT-BIAS ONE FOLLOWS THE STRAND, which is
 *     what the same substitution in read one and in read two tells apart;
 *   - THE ERROR RATE IS FLOORED SO THE Q SCORE IS FINITE, and a substitution nothing was seen for
 *     still gets a row;
 *   - --MINIMUM_QUALITY_SCORE DROPS SINGLE BASES AND --USE_OQ MAKES IT READ THE `OQ` TAG;
 *   - --MINIMUM_MAPPING_QUALITY, THE INSERT-SIZE WINDOW, A DUPLICATE AND A SECONDARY ALIGNMENT
 *     EACH DROP A WHOLE READ, and --INCLUDE_DUPLICATES and --INCLUDE_UNPAIRED put two of them back;
 *   - AN UNPAIRED READ IS DROPPED BY THE INSERT-SIZE FILTER RATHER THAN BY A RULE OF ITS OWN,
 *     which is why --INCLUDE_UNPAIRED is what admits it;
 *   - --TANDEM_READS SWAPS READ TWO'S HALF OF EVERY SUM, which is invisible on a pair whose ends
 *     are on opposite strands (both conventions put such a pair on one side) and plain on a
 *     substitution carried by read two on the forward strand;
 *   - --DB_SNP AND --INTERVALS EACH REMOVE SITES BEFORE ANY COUNTING;
 *   - AND THE FIFTH FILE, THE ERROR SUMMARY, IS WRITTEN WHATEVER ELSE HAPPENS.
 *
 * Output:
 *
 *     ref\t<name>\t<the reference bases>
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     files\t<case>\t<the file names the prefix produced, comma separated>
 *     summary\t<case>\t<kind>=<that summary file's rows, escaped>
 *     detail\t<case>\t<kind>=<the detail rows that counted anything, escaped>
 *     rows\t<case>\t<kind>=<how many detail rows the file holds>
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
import java.util.Map;
import java.util.TreeMap;

public class CollectSequencingArtifactMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** The reference: every three-base context occurs, and the bases are not all one letter. */
    static final String BASES = "CACCGTTACGGTTAACCTGACGTACGTTGCAGCTAGCTAACGTGCATCGATCGTAGCTAGCTA"
            + "GGCCTTAACCGGTTAACCGGATCGATCGTTAAGGCCTTAAGGCCATCGATCGGCTA";

    static String fasta() {
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < BASES.length(); i += 60) {
            fasta.append(BASES, i, Math.min(i + 60, BASES.length())).append('\n');
        }
        return fasta.toString();
    }

    /** One read: which end, which strand, and where it disagrees with the reference. */
    record Read(String name, int start, int length, boolean first, boolean negative, boolean paired,
                int mappingQuality, int insertSize, Map<Integer, Character> substitutions,
                String qualities, String originalQualities, boolean duplicate, boolean secondary,
                boolean failsVendor) {}

    static Map<Integer, Character> none() {
        return new TreeMap<>();
    }

    static Map<Integer, Character> at(final int position, final char base) {
        final Map<Integer, Character> substitutions = new TreeMap<>();
        substitutions.put(position, base);
        return substitutions;
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", BASES.length()));
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
                record.setReadNegativeStrandFlag(spec.negative());
                if (spec.paired()) {
                    record.setReadPairedFlag(true);
                    record.setProperPairFlag(true);
                    record.setFirstOfPairFlag(spec.first());
                    record.setSecondOfPairFlag(!spec.first());
                    record.setMateNegativeStrandFlag(!spec.negative());
                    record.setMateReferenceName("chr1");
                    record.setMateAlignmentStart(spec.first() ? spec.start() + 90 : spec.start() - 90);
                    record.setInferredInsertSize(spec.insertSize());
                }
                record.setDuplicateReadFlag(spec.duplicate());
                record.setSecondaryAlignment(spec.secondary());
                record.setReadFailsVendorQualityCheckFlag(spec.failsVendor());
                final StringBuilder bases = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    final int position = spec.start() + i;
                    final Character substitution = spec.substitutions().get(position);
                    bases.append(substitution != null ? substitution : BASES.charAt(position - 1));
                }
                record.setReadString(bases.toString());
                final StringBuilder quals = new StringBuilder();
                for (int i = 0; i < spec.length(); i++) {
                    quals.append(spec.qualities() == null ? 'I' : spec.qualities().charAt(i));
                }
                record.setBaseQualityString(quals.toString());
                if (spec.originalQualities() != null) {
                    record.setAttribute("OQ", spec.originalQualities());
                }
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /** A pair, one end of which carries the substitutions on the strand asked for. */
    static List<Read> pair(final Map<Integer, Character> substitutions, final boolean inFirst,
                           final boolean negative) {
        return List.of(
                new Read("a", 20, 30, true, inFirst && negative, true, 60, 120,
                        inFirst ? substitutions : none(), null, null, false, false, false),
                new Read("a", 20, 30, false, !inFirst && negative, true, 60, -120,
                        inFirst ? none() : substitutions, null, null, false, false, false));
    }

    /** The clean pair every filter case carries alongside the read being filtered. */
    static final List<Read> CLEAN = List.of(
            new Read("b", 20, 30, true, false, true, 60, 120, none(), null, null, false, false, false),
            new Read("b", 20, 30, false, true, true, 60, -120, none(), null, null, false, false, false));

    static List<Read> with(final Read one) {
        final List<Read> reads = new ArrayList<>(CLEAN);
        reads.add(one);
        return reads;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        emit("ref", "chr1", BASES);

        // Position 22 is a G in the fixture, so a T there is a G>T substitution.
        run("plain", pair(none(), true, false));
        run("alt-read-one-forward", pair(at(22, 'T'), true, false));
        run("alt-read-two-forward", pair(at(22, 'T'), false, false));
        run("alt-read-one-reverse", pair(at(22, 'T'), true, true));
        run("tandem-reads", pair(at(22, 'T'), true, false), "TANDEM_READS=true");
        // A pair whose ends are on OPPOSITE strands lands on the same side under either
        // convention, so the swap only shows when the substitution is in read two on the forward
        // strand: that end changes side and the other does not.
        run("tandem-read-two-forward", pair(at(22, 'T'), false, false), "TANDEM_READS=true");

        // The detail files' size, and what cuts it.
        run("context-size-zero", pair(at(22, 'T'), true, false), "CONTEXT_SIZE=0");
        run("contexts-to-print", pair(at(22, 'T'), true, false), "CONTEXTS_TO_PRINT=ACG");
        run("file-extension", pair(at(22, 'T'), true, false), "FILE_EXTENSION=.txt");

        // The base-quality floor and the tag it can read instead.
        final String lowAt = "I".repeat(2) + "#" + "I".repeat(27);
        run("low-base-quality", List.of(
                new Read("a", 20, 30, true, false, true, 60, 120, at(22, 'T'), lowAt, null,
                        false, false, false),
                new Read("a", 20, 30, false, true, true, 60, -120, none(), null, null,
                        false, false, false)));
        run("original-qualities", List.of(
                new Read("a", 20, 30, true, false, true, 60, 120, at(22, 'T'), null, lowAt,
                        false, false, false),
                new Read("a", 20, 30, false, true, true, 60, -120, none(), null, null,
                        false, false, false)));
        run("original-qualities-ignored", List.of(
                new Read("a", 20, 30, true, false, true, 60, 120, at(22, 'T'), null, lowAt,
                        false, false, false),
                new Read("a", 20, 30, false, true, true, 60, -120, none(), null, null,
                        false, false, false)), "USE_OQ=false");

        // The whole-read filters, each alongside the clean pair.
        run("low-mapping-quality", with(new Read("a", 20, 30, true, false, true, 20, 120,
                at(22, 'T'), null, null, false, false, false)));
        run("insert-too-small", with(new Read("a", 20, 30, true, false, true, 60, 40,
                at(22, 'T'), null, null, false, false, false)));
        run("insert-too-large", with(new Read("a", 20, 30, true, false, true, 60, 900,
                at(22, 'T'), null, null, false, false, false)));
        final Read duplicate = new Read("a", 20, 30, true, false, true, 60, 120, at(22, 'T'),
                null, null, true, false, false);
        run("duplicate", with(duplicate));
        run("duplicate-included", with(duplicate), "INCLUDE_DUPLICATES=true");
        run("secondary", with(new Read("a", 20, 30, true, false, true, 60, 120, at(22, 'T'),
                null, null, false, true, false)));
        final Read vendor = new Read("a", 20, 30, true, false, true, 60, 120, at(22, 'T'),
                null, null, false, false, true);
        run("fails-vendor", with(vendor));
        run("fails-vendor-included", with(vendor), "INCLUDE_NON_PF_READS=true");
        final Read unpaired = new Read("a", 20, 30, true, false, false, 60, 0, at(22, 'T'),
                null, null, false, false, false);
        run("unpaired", with(unpaired));
        run("unpaired-included", with(unpaired), "INCLUDE_UNPAIRED=true");

        System.out.print(buf);
    }

    /** The five extensions the prefix produces, in the order the tool declares them. */
    static final List<String> KINDS = List.of(
            "pre_adapter_summary_metrics", "pre_adapter_detail_metrics",
            "bait_bias_summary_metrics", "bait_bias_detail_metrics", "error_summary_metrics");

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("artifacts");
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

        final Path prefix = dir.resolve("out");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + prefix, "R=" + reference));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.artifacts.CollectSequencingArtifactMetrics()
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
        final List<String> written = new ArrayList<>();
        try (final var stream = Files.list(dir)) {
            stream.map(path -> path.getFileName().toString())
                    .filter(file -> file.startsWith("out."))
                    .sorted()
                    .forEach(written::add);
        }
        emit("files", name, String.join(",", written));
        for (final String kind : KINDS) {
            for (final String suffix : List.of("", ".txt")) {
                final Path file = dir.resolve("out." + kind + suffix);
                if (!Files.exists(file)) {
                    continue;
                }
                final List<String> table = new ArrayList<>();
                for (final String line : Files.readString(file, StandardCharsets.UTF_8)
                        .split("\n", -1)) {
                    if (!line.startsWith("#") && !line.isEmpty()) {
                        table.add(line);
                    }
                }
                if (table.isEmpty()) {
                    emit("rows", name, kind + "=0");
                    continue;
                }
                emit("rows", name, kind + "=" + (table.size() - 1));
                if (kind.endsWith("detail_metrics")) {
                    // Ninety-six rows a library, of which a handful counted anything: the golden
                    // holds the header and the rows whose counts are not all nought.
                    final List<String> counted = new ArrayList<>(List.of(table.get(0)));
                    for (final String row : table.subList(1, table.size())) {
                        if (!allZero(row)) {
                            counted.add(row);
                        }
                    }
                    emit("detail", name, kind + "=" + String.join("\n", counted));
                } else {
                    emit("summary", name, kind + "=" + String.join("\n", table));
                }
            }
        }
    }

    /** Whether every integer column of a detail row is nought. */
    static boolean allZero(final String row) {
        for (final String column : row.split("\t", -1)) {
            if (column.matches("\\d+") && !column.equals("0")) {
                return false;
            }
        }
        return true;
    }
}
