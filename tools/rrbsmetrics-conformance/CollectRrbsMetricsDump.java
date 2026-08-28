/*
 * CollectRrbsMetrics' two metrics files, taken from the reference.
 *
 * The tool walks a bisulfite-converted alignment against its reference and counts cytosines twice
 * over: the ones in a CpG context, per site, and the ones outside it, in bulk. What is measured
 * here is which bases are counted at all, which is a set of rules that do not compose the way the
 * code around them reads.
 *
 * Ten behaviours this is built to catch.
 *
 *   - A CpG IS A REFERENCE C FOLLOWED BY A REFERENCE G, and it is looked for in the reference
 *     rather than in the read, so a read that reads TG over a reference CG is a CONVERTED CpG and
 *     not a mismatch;
 *   - THE LAST BASE OF AN ALIGNMENT BLOCK IS NEVER A CpG, the loop stopping one short so that the
 *     pair has a second base to be;
 *   - A VALID CpG NEEDS THREE THINGS: the read's C matching under bisulfite rules, the read's G
 *     matching exactly, and both qualities over their own thresholds, which are two different
 *     thresholds;
 *   - THE NON-CpG CYTOSINE BRANCH READS THE QUALITY OFF THE WHOLE READ, with the block's own
 *     index: `isAboveCytoQcThreshold(readQualities, i)` where the CpG branch passes the block's
 *     `readQualityFragment`. On a read whose alignment starts at its first base the two agree, and
 *     on any other read they do not;
 *   - THE COMMENT ABOVE THAT BRANCH SAYS THE NON-CpG COUNTS ARE HELD BACK until the read is known
 *     to carry a CpG, AND THE CODE DOES NOT: `nCytoTotal` is incremented straight into the totals,
 *     so a read with no CpG at all still contributes non-CpG bases;
 *   - A NEGATIVE-STRAND READ IS REVERSE COMPLEMENTED AND THEN TREATED AS POSITIVE, and the site it
 *     reports is `refStart + (blockLength - 1) - i - 1`, which is one further left than the
 *     mirrored index;
 *   - THE POSITION IN THE DETAIL FILE IS THE ZERO-BASED ONE, the block's reference start being
 *     taken as `getReferenceStart() - 1`;
 *   - A READ SHORTER THAN MINIMUM_READ_LENGTH AND ONE OVER MAX_MISMATCH_RATE ARE COUNTED AND
 *     DROPPED, the mismatch bound being `Math.round(length * rate)` and therefore a STRICTLY
 *     greater test against a rounded integer;
 *   - THE PREFIX GAINS A DOT IF IT HAS NONE, so the three files are named after it rather than
 *     concatenated to it;
 *   - AND THE PDF IS R'S, so its bytes are not something a golden can hold.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     summary\t<case>\t<the summary table without its comments, escaped>
 *     detail\t<case>\t<the detail table without its comments, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CollectRrbsMetricsDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CollectRrbsMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /**
     * The reference: `ACGTTCAACGTA` repeated, which is three cytosines per twelve bases of two
     * different kinds.
     *
     * `ACGT` repeated would have been simpler and would have measured half of the tool: every C in
     * it is followed by a G, so the non-CpG cytosine branch never runs and `NON_CPG_BASES` is zero
     * in every case. This pattern has two CpG sites (at offsets 1 and 8) and one isolated C (at
     * offset 5, followed by an A), so both branches are exercised by the same read.
     *
     * It is also not its own reverse complement, which `ACGT` repeated is. That matters for the
     * negative-strand case: on a palindrome the site the tool reports and the site a reader would
     * expect coincide by accident.
     */
    static final int CONTIG_LENGTH = 120;

    static final String MOTIF = "ACGTTCAACGTA";

    static String bases() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < CONTIG_LENGTH; i++) {
            bases.append(MOTIF.charAt(i % MOTIF.length()));
        }
        return bases.toString();
    }

    static String fasta() {
        final String bases = bases();
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /** One read: where it starts, what it reads, what its qualities are and which strand it is. */
    record Read(String name, int start, String bases, String qualities, boolean negative,
                String cigar) {}

    static Read read(final String name, final int start, final String bases) {
        return new Read(name, start, bases, null, false, null);
    }

    /** The reference over a span with the given ONE-BASED positions turned into a G.
     *
     * `A` is the base that appears at positions the motif leaves free, and a `G` read over an `A`
     * is a mismatch bisulfite cannot explain, where a `T` read over a `C` is a conversion.
     */
    static String mutated(final int length, final int... positions) {
        final char[] bases = reference(1, length).toCharArray();
        for (final int position : positions) {
            bases[position - 1] = 'G';
        }
        return new String(bases);
    }

    /** The reference's own bases over a span, which is what an unconverted read reads. */
    static String reference(final int start, final int length) {
        return bases().substring(start - 1, start - 1 + length);
    }

    static SAMFileHeader header(final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        header.addReadGroup(group);
        return header;
    }

    static void writeBam(final Path bam, final List<Read> reads,
                         final SAMFileHeader.SortOrder order) {
        final SAMFileHeader header = header(order);
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .makeBAMWriter(header, false, bam.toFile())) {
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setMappingQuality(60);
                record.setReadString(spec.bases());
                record.setCigarString(spec.cigar() == null
                        ? spec.bases().length() + "M" : spec.cigar());
                final StringBuilder quals = new StringBuilder();
                for (int i = 0; i < spec.bases().length(); i++) {
                    quals.append(spec.qualities() == null ? 'I' : spec.qualities().charAt(i));
                }
                record.setBaseQualityString(quals.toString());
                record.setReadNegativeStrandFlag(spec.negative());
                record.setAttribute("RG", "rg1");
                writer.addAlignment(record);
            }
        }
    }

    /** A metrics table without its comment lines, which carry the command line and the clock. */
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

    static void run(final String name, final List<Read> reads, final String... extra)
            throws Exception {
        run(name, reads, SAMFileHeader.SortOrder.coordinate, "m", extra);
    }

    static void run(final String name, final List<Read> reads, final SAMFileHeader.SortOrder order,
                    final String prefix, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("rrbsmetrics");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path in = dir.resolve("in.bam");
        writeBam(in, reads, order);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final Path out = Files.createDirectory(dir.resolve("out"));
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "M=" + out.resolve(prefix), "R=" + reference));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.analysis.CollectRrbsMetrics()
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
        for (final File file : out.toFile().listFiles()) {
            written.add(file.getName());
        }
        Collections.sort(written);
        emit("files", name, String.join(" ", written));
        for (final String file : written) {
            if (file.endsWith("rrbs_summary_metrics")) {
                emit("summary", name, table(out.resolve(file)));
            } else if (file.endsWith("rrbs_detail_metrics")) {
                emit("detail", name, table(out.resolve(file)));
            }
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // A read that reads the reference back: every CpG is seen and none is converted.
        run("unconverted", List.of(read("a", 1, reference(1, 24))));
        // The same read with every CpG converted, which is what bisulfite does to an unmethylated
        // site: the C of each `CG` is read as a T.
        run("cpg-converted", List.of(read("a", 1, reference(1, 24).replace("CG", "TG"))));
        // And the same read with only the ISOLATED cytosine converted, which is the other branch:
        // the `CA` at offset five becomes `TA` and no CpG moves.
        run("non-cpg-converted", List.of(read("a", 1, reference(1, 24).replace("CA", "TA"))));
        // Both at once, which is what a real conversion looks like.
        run("both-converted", List.of(read("a", 1,
                reference(1, 24).replace("CG", "TG").replace("CA", "TA"))));

        // The last base of the block is never a CpG: the read that ends ON the C of the second
        // `CG` pair does not report it, and the one base longer does.
        run("cpg-at-the-last-base", List.of(read("a", 1, reference(1, 9))));
        run("cpg-one-base-longer", List.of(read("a", 1, reference(1, 10))));

        // The two quality thresholds, which are different numbers: the C's own and its neighbour's.
        final String twenty = "5";  // Phred 20
        final String nineteen = "4";
        final String ten = "+";     // Phred 10
        final String nine = "*";
        run("c-at-the-threshold", List.of(new Read("a", 1, reference(1, 12),
                "I" + twenty + "IIIIIIIIII", false, null)));
        run("c-under-the-threshold", List.of(new Read("a", 1, reference(1, 12),
                "I" + nineteen + "IIIIIIIIII", false, null)));
        run("neighbour-at-the-threshold", List.of(new Read("a", 1, reference(1, 12),
                "II" + ten + "IIIIIIIII", false, null)));
        run("neighbour-under-the-threshold", List.of(new Read("a", 1, reference(1, 12),
                "II" + nine + "IIIIIIIII", false, null)));
        // The isolated cytosine sits at offset five, so the same two thresholds again on the
        // branch that reads its quality off the whole read rather than off the block.
        run("isolated-c-under-the-threshold", List.of(new Read("a", 1, reference(1, 12),
                "IIIII" + nineteen + "IIIIII", false, null)));
        run("isolated-c-neighbour-under-the-threshold", List.of(new Read("a", 1, reference(1, 12),
                "IIIIII" + nine + "IIIII", false, null)));

        // A negative-strand read, whose site is reported one further left than the mirror.
        run("negative-strand", List.of(new Read("a", 1, reference(1, 24), null, true, null)));

        // A read that starts partway into itself, which is where the two quality readings differ:
        // the CpG branch reads the block's own qualities and the non-CpG branch the whole read's.
        run("soft-clipped", List.of(new Read("a", 5, "NNNN" + reference(5, 16), null, false,
                "4S16M")));
        // The same read with the low quality at the WHOLE READ's index of the isolated cytosine
        // and not at the block's, which is the one place the two readings disagree.
        run("soft-clipped-low-quality-off-by-the-clip", List.of(new Read("a", 5,
                "NNNN" + reference(5, 16), "IIIII" + nineteen + "IIIIIIIIIIIIII", false,
                "4S16M")));

        // The two filters, and the boundary of the second: the bound is Math.round(length * rate)
        // and the test is strictly greater.
        run("short-read", List.of(read("a", 1, reference(1, 4))));
        run("short-read-allowed", List.of(read("a", 1, reference(1, 4))),
                "MINIMUM_READ_LENGTH=4");
        // The mismatches are A-to-G, which bisulfite does not explain, so they count where a
        // C-to-T would not. A read of twenty bases has a bound of `Math.round(20 * 0.1)` = 2 and
        // the test is strictly greater, so two are kept and three are dropped.
        run("one-mismatch", List.of(read("a", 1, mutated(20, 1))));
        run("two-mismatches", List.of(read("a", 1, mutated(20, 1, 7))));
        run("three-mismatches", List.of(read("a", 1, mutated(20, 1, 7, 8))));
        // And a read of the same three mismatches with the rate raised, which keeps it.
        run("three-mismatches-allowed", List.of(read("a", 1, mutated(20, 1, 7, 8))),
                "MAX_MISMATCH_RATE=0.2");
        // A C-to-T on a CpG is not a mismatch at all: it is a conversion, which is the whole
        // point of the tool.
        run("four-conversions-are-not-mismatches",
                List.of(read("a", 1, reference(1, 20).replace("CG", "TG"))));

        // A read with no CpG at all, which still contributes non-CpG cytosines: the comment in
        // the reference says these are held back until a CpG is seen, and the code does not.
        run("no-cpg", List.of(read("a", 4, reference(4, 6))));

        // The sequence filter, the accumulation level, and a prefix without its dot.
        run("sequence-filter", List.of(read("a", 1, reference(1, 20))),
                "SEQUENCE_NAMES=chr2");
        run("sequence-filter-matching", List.of(read("a", 1, reference(1, 20))),
                "SEQUENCE_NAMES=chr1");
        run("accumulation-level", List.of(read("a", 1, reference(1, 20))),
                "LEVEL=ALL_READS", "LEVEL=SAMPLE");
        run("prefix-with-a-dot", List.of(read("a", 1, reference(1, 20))),
                SAMFileHeader.SortOrder.coordinate, "m.");

        // A file the header calls queryname sorted, on both paths.
        run("queryname-sorted", List.of(read("a", 1, reference(1, 20))),
                SAMFileHeader.SortOrder.queryname, "m");
        run("queryname-sorted-assumed", List.of(read("a", 1, reference(1, 20))),
                SAMFileHeader.SortOrder.queryname, "m", "ASSUME_SORTED=true");

        // And a file with no reads at all.
        run("empty", List.of());

        System.out.print(buf);
    }
}
