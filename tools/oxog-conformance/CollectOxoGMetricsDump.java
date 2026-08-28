/*
 * CollectOxoGMetrics' counters, taken from the reference.
 *
 * The tool walks every C and G of the reference and asks, for each read covering one, whether the
 * base it carries is consistent with the 8-oxo-G damage that happens during shearing. What is
 * measured is which read reaches which of the eight counters, and what the derived rates make of
 * them.
 *
 * Thirteen behaviours this is built to catch.
 *
 *   - THE ROWS ARE ONE PER CONTEXT PER LIBRARY, and their order is a HashMap's over the context
 *     strings rather than a sorted one, which is the order that reaches the file;
 *   - A G SITE IS FOLDED INTO THE REVERSE COMPLEMENT OF ITS CONTEXT, so the C and the G halves of
 *     one context share a row and are told apart only by the C_REF and G_REF columns;
 *   - WHICH END OF THE PAIR CARRIES THE ALTERNATE DECIDES THE COUNTER: `G>T` in read one and
 *     `C>A` in read two are the oxidised state, and the same substitutions in the other end are
 *     the control;
 *   - THE READ'S OWN ORIENTATION IS UNDONE FIRST, a base on the negative strand being complemented
 *     back before the question is asked;
 *   - THE OXIDATION ERROR RATE IS FLOORED AT ONE BASE and not at zero, so a library with no
 *     alternate at all reports one over its total rather than nothing, and its Q is finite;
 *   - THE TWO REFERENCE-BIAS RATES ARE FLOORED AT 1e-10, which caps their Q at a hundred;
 *   - --MINIMUM_QUALITY_SCORE DROPS SINGLE BASES, and --USE_OQ MAKES IT READ THE `OQ` TAG when the
 *     read carries one, so the same read is counted or not according to a tag the alignment did
 *     not have to keep;
 *   - --MINIMUM_MAPPING_QUALITY DROPS WHOLE READS;
 *   - THE INSERT-SIZE WINDOW DROPS WHOLE PAIRS, at both ends of it;
 *   - A DUPLICATE AND A SECONDARY ALIGNMENT ARE DROPPED;
 *   - A SITE WITHIN THE CONTEXT OF EITHER END OF THE CONTIG IS SKIPPED, so a C at position one is
 *     never assayed;
 *   - --CONTEXT_SIZE 0 LEAVES ONE CONTEXT, `C`, and --CONTEXTS RESTRICTS THE ROWS to what it
 *     names, a context that never occurs still getting a row of zeroes;
 *   - AND A FILE WITH NO READ GROUP IS REFUSED, the analysis needing a library name to key on.
 *
 * Output:
 *
 *     ref\t<name>\t<the reference bases>
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     order\t<case>\t<the context of each row, in the file's own order>
 *     rows\t<case>\t<the rows with a site in them, escaped>
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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class CollectOxoGMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /**
     * The reference, whose C and G sites are the ones the cases assay.
     *
     * `ACCGTTACGGTTAAC...` repeated, so that several contexts occur and a C sits at position one,
     * which the end-of-contig rule then skips.
     */
    static final String BASES = "CACCGTTACGGTTAACCTGACGTACGTTGCAGCTAGCTAACGTGCATCGATCGTAGCTAGCTA"
            + "GGCCTTAACCGGTTAACCGGATCGATCGTTAAGGCCTTAAGGCCATCGATCGGCTA";

    static String fasta() {
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < BASES.length(); i += 60) {
            fasta.append(BASES, i, Math.min(i + 60, BASES.length())).append('\n');
        }
        return fasta.toString();
    }

    /**
     * One read of a pair.
     *
     * `substitutions` names the reference positions this read disagrees with and the base it
     * carries there, so a case can put a `T` in read one at a G site and the same `T` in read two.
     */
    record Read(String name, int start, int length, boolean first, boolean negative,
                int mappingQuality, int insertSize, Map<Integer, Character> substitutions,
                String qualities, String originalQualities, boolean duplicate, boolean secondary,
                String library) {}

    static Read read(final String name, final int start, final boolean first,
                     final Map<Integer, Character> substitutions) {
        return new Read(name, start, 30, first, !first, 60, first ? 120 : -120, substitutions,
                null, null, false, false, "lib1");
    }

    static SAMFileHeader header(final List<String> libraries) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", BASES.length()));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        for (final String library : libraries) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg-" + library);
            group.setSample("sample1");
            group.setLibrary(library);
            header.addReadGroup(group);
        }
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads, final List<String> libraries) {
        final SAMFileHeader head = header(libraries);
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
                record.setReadNegativeStrandFlag(spec.negative());
                record.setMateNegativeStrandFlag(!spec.negative());
                record.setMateReferenceName("chr1");
                record.setMateAlignmentStart(spec.first() ? spec.start() + 90 : spec.start() - 90);
                record.setInferredInsertSize(spec.insertSize());
                record.setDuplicateReadFlag(spec.duplicate());
                record.setSecondaryAlignment(spec.secondary());
                // A file with no read group at all is one of the cases, and a record may not
                // name a group the header does not declare.
                if (!libraries.isEmpty()) {
                    record.setAttribute("RG", "rg-" + spec.library());
                }
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
                writer.addAlignment(record);
            }
        }
    }

    /** The rows the metrics file holds, in its own order, keyed by context and library. */
    static void run(final String name, final List<Read> reads, final List<String> libraries,
                    final Map<String, String> files, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("oxog");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads, libraries);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + dir.resolve("out.txt"), "R=" + reference));
        for (final Map.Entry<String, String> file : files.entrySet()) {
            final Path path = dir.resolve(file.getKey());
            Files.writeString(path, file.getValue(), StandardCharsets.UTF_8);
            argv.add(file.getKey().endsWith(".vcf") ? "DB_SNP=" + path : "INTERVALS=" + path);
        }
        argv.addAll(Arrays.asList(extra));
        final java.io.ByteArrayOutputStream errBytes = new java.io.ByteArrayOutputStream();
        final java.io.PrintStream realErr = System.err;
        try {
            final int code;
            try {
                System.setErr(new java.io.PrintStream(errBytes, true, StandardCharsets.UTF_8));
                code = new picard.analysis.CollectOxoGMetrics()
                        .instanceMain(argv.toArray(new String[0]));
            } finally {
                System.err.flush();
                System.setErr(realErr);
            }
            if (code != 0) {
                emit("error", name, "exit " + code);
                // The refusal itself is one line of a usage the golden has no reason to hold.
                final List<String> refusal = new ArrayList<>();
                for (final String line : errBytes.toString(StandardCharsets.UTF_8).split("\n", -1)) {
                    if (line.startsWith("Middle base of context sequence")
                            || line.startsWith("Context ")) {
                        refusal.add(line);
                    }
                }
                emit("refusal", name, String.join("\n", refusal));
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
        final List<String> table = new ArrayList<>();
        for (final String line : Files.readString(dir.resolve("out.txt"),
                StandardCharsets.UTF_8).split("\n", -1)) {
            if (!line.startsWith("#") && !line.isEmpty()) {
                table.add(line);
            }
        }
        final String[] header = table.get(0).split("\t", -1);
        int contextColumn = 0;
        int sitesColumn = 0;
        for (int i = 0; i < header.length; i++) {
            if (header[i].equals("CONTEXT")) {
                contextColumn = i;
            }
            if (header[i].equals("TOTAL_SITES")) {
                sitesColumn = i;
            }
        }
        final List<String> order = new ArrayList<>();
        final List<String> populated = new ArrayList<>();
        populated.add(table.get(0));
        for (final String row : table.subList(1, table.size())) {
            final String[] columns = row.split("\t", -1);
            order.add(columns[contextColumn]);
            if (!columns[sitesColumn].equals("0")) {
                populated.add(row);
            }
        }
        emit("order", name, String.join(",", order));
        emit("rows", name, String.join("\n", populated));
    }

    /**
     * A pair whose named end carries the substitutions, on the strand asked for.
     *
     * Which end and which strand are the whole of the question the counters answer, so they are
     * given separately here rather than tied together as a real library's would be.
     */
    static List<Read> alt(final boolean inFirst, final boolean negative,
                          final Map<Integer, Character> substitutions) {
        return List.of(
                new Read("a", 20, 30, true, inFirst && negative, 60, 120,
                        inFirst ? substitutions : none(), null, null, false, false, "lib1"),
                new Read("a", 20, 30, false, !inFirst && negative, 60, -120,
                        inFirst ? none() : substitutions, null, null, false, false, "lib1"));
    }

    /** The clean pair, plus the one read a case is about. */
    static List<Read> with(final List<Read> clean, final Read one) {
        final List<Read> reads = new ArrayList<>(clean);
        reads.add(one);
        return reads;
    }

    static Map<Integer, Character> none() {
        return new LinkedHashMap<>();
    }

    static Map<Integer, Character> at(final int position, final char base) {
        final Map<Integer, Character> substitutions = new LinkedHashMap<>();
        substitutions.put(position, base);
        return substitutions;
    }

    /** A pair whose two ends cover the same window, so one site is seen by both. */
    static List<Read> pair(final String name, final int start,
                           final Map<Integer, Character> inFirst,
                           final Map<Integer, Character> inSecond) {
        return List.of(
                new Read(name, start, 30, true, false, 60, 120, inFirst, null, null, false, false,
                        "lib1"),
                new Read(name, start, 30, false, true, 60, -120, inSecond, null, null, false,
                        false, "lib1"));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        emit("ref", "chr1", BASES);
        final List<String> one = List.of("lib1");

        // A clean pair, which every filter case carries alongside the one being filtered, so a
        // dropped read shows as a missing base and not as an empty file.
        final List<Read> clean = List.of(
                new Read("b", 20, 30, true, false, 60, 120, none(), null, null, false, false,
                        "lib1"),
                new Read("b", 20, 30, false, true, 60, -120, none(), null, null, false, false,
                        "lib1"));

        // A pair carrying the reference everywhere, which fills the REF columns alone.
        run("reference-only", pair("a", 20, none(), none()), one, Map.of());

        // The alternate at a G site, in each end and on each strand. Position 22 is a G and
        // position 21 is a C in the fixture.
        run("alt-read-one-forward", alt(true, false, at(22, 'T')), one, Map.of());
        run("alt-read-two-reverse", alt(false, true, at(22, 'T')), one, Map.of());
        run("alt-read-one-reverse", alt(true, true, at(22, 'T')), one, Map.of());
        run("alt-read-two-forward", alt(false, false, at(22, 'T')), one, Map.of());

        // And at a C site, whose alternate is an A rather than a T.
        run("c-site-read-one-forward", alt(true, false, at(21, 'A')), one, Map.of());
        run("c-site-read-two-forward", alt(false, false, at(21, 'A')), one, Map.of());

        // The same alternate under the quality floor, and the OQ tag disagreeing with it.
        final String lowAt = "I".repeat(2) + "#" + "I".repeat(27);
        run("low-base-quality", with(clean, new Read("a", 20, 30, true, false, 60, 120,
                at(22, 'T'), lowAt, null, false, false, "lib1")), one, Map.of());
        run("original-qualities", with(clean, new Read("a", 20, 30, true, false, 60, 120,
                at(22, 'T'), null, lowAt, false, false, "lib1")), one, Map.of());
        run("original-qualities-ignored", with(clean, new Read("a", 20, 30, true, false, 60, 120,
                at(22, 'T'), null, lowAt, false, false, "lib1")), one, Map.of(), "USE_OQ=false");

        // The whole-read filters, each alongside the clean pair.
        run("low-mapping-quality", with(clean, new Read("a", 20, 30, true, false, 20, 120,
                at(22, 'T'), null, null, false, false, "lib1")), one, Map.of());
        run("insert-too-small", with(clean, new Read("a", 20, 30, true, false, 60, 40,
                at(22, 'T'), null, null, false, false, "lib1")), one, Map.of());
        run("insert-too-large", with(clean, new Read("a", 20, 30, true, false, 60, 900,
                at(22, 'T'), null, null, false, false, "lib1")), one, Map.of());
        run("duplicate", with(clean, new Read("a", 20, 30, true, false, 60, 120,
                at(22, 'T'), null, null, true, false, "lib1")), one, Map.of());
        run("secondary", with(clean, new Read("a", 20, 30, true, false, 60, 120,
                at(22, 'T'), null, null, false, true, "lib1")), one, Map.of());

        // The first base of the contig, which the context rule skips whatever covers it.
        run("contig-start", pair("a", 1, none(), none()), one, Map.of());

        // The context arguments, and one the validation refuses.
        run("context-size-zero", alt(true, false, at(22, 'T')), one, Map.of(),
                "CONTEXT_SIZE=0");
        run("restricted-contexts", alt(true, false, at(22, 'T')), one, Map.of(),
                "CONTEXTS=ACG", "CONTEXTS=TCG");
        run("context-without-a-c", alt(true, false, at(22, 'T')), one, Map.of(),
                "CONTEXTS=TTT");
        run("context-of-the-wrong-length", alt(true, false, at(22, 'T')), one, Map.of(),
                "CONTEXTS=ACGTA");

        // Two libraries, whose reads are told apart by the read group alone.
        run("two-libraries", List.of(
                new Read("a", 20, 30, true, false, 60, 120, at(22, 'T'), null, null, false, false,
                        "lib1"),
                new Read("a", 20, 30, false, true, 60, -120, none(), null, null, false, false,
                        "lib1"),
                new Read("b", 20, 30, true, false, 60, 120, none(), null, null, false, false,
                        "lib2"),
                new Read("b", 20, 30, false, true, 60, -120, none(), null, null, false, false,
                        "lib2")), List.of("lib1", "lib2"), Map.of());

        // The debugging stop, which cuts the walk short.
        run("stop-after", pair("a", 20, none(), none()), one, Map.of(), "STOP_AFTER=2");

        // A known site, which is skipped whatever covers it.
        run("dbsnp", alt(true, false, at(22, 'T')), one, Map.of("known.vcf",
                "##fileformat=VCFv4.2\n"
                        + "##contig=<ID=chr1,length=" + BASES.length() + ">\n"
                        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        + "chr1\t22\trs1\tG\tT\t.\t.\t.\n"));

        // An interval list that leaves the assayed site outside it.
        run("intervals", alt(true, false, at(22, 'T')), one,
                Map.of("targets.interval_list",
                        "@HD\tVN:1.6\tSO:coordinate\n"
                                + "@SQ\tSN:chr1\tLN:" + BASES.length() + "\n"
                                + "chr1\t40\t60\t+\ttarget\n"));

        // A file with no read group at all.
        run("no-read-group", pair("a", 20, none(), none()), List.of(), Map.of());

        System.out.print(buf);
    }
}
