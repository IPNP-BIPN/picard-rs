/*
 * BamToBfq's binary output, taken from the reference.
 *
 * `bfq` is Maq's binary fastq: one record is the read's name, its length, and one byte per base
 * carrying the base in the top two bits and the quality in the bottom six. Everything about it is
 * fixed-width and little-endian, so the file is byte-comparable and the interesting questions are
 * which reads reach it and what each byte says.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE FILE IS GZIP: `IOUtil.openFileForWriting` compresses a `.bfq` the way it compresses a
 *     `.gz`, so the bytes on disk are a DEFLATE stream and reproducing them is the programme's
 *     own hard problem. The payload is recorded beside the file, and it is what a port claims;
 *   - ONE BYTE CARRIES BOTH: `(base << 6) | quality`, so a quality over 63 could not be written
 *     and an `A` at quality 30 is the same byte as a `C` at quality 30 shifted by sixty-four;
 *   - `A`, `C`, `G` AND `T` ARE 0, 1, 2 AND 3, and lower case is the same as upper;
 *   - THE TOOL'S "UNKNOWN BASE" BRANCH IS UNREACHABLE THROUGH A BAM: htsjdk refuses the record
 *     at write time, so the fixture cannot be built and the branch cannot be entered from this
 *     input format at all. The refusal the WRITER makes is recorded in its place;
 *   - AN `N` IS WRITTEN AS AN `A`, and its quality is not the read's: it is ONE inside the seed
 *     region for the first few of them and ZERO after that, and one outside it;
 *   - THE FILES ARE NAMED `<prefix><index>.<read>.bfq`, the index counting chunks from one and
 *     the read being 1 or 2, so a single-end run writes half as many files;
 *   - `--READ_CHUNK_SIZE` SPLITS THE OUTPUT rather than truncating it, and the last chunk is
 *     short;
 *   - `--READS_TO_ALIGN` TRUNCATES IT, and it counts records rather than files;
 *   - `--PAIRED_RUN` DECIDES WHICH READS ARE WRITTEN AT ALL, an unpaired file through a paired
 *     run being an error rather than an empty output;
 *   - `--BASES_TO_WRITE` TRIMS EVERY READ to that length, which changes the record's length field
 *     and not only its bytes;
 *   - `--INCLUDE_NON_PF_READS` DECIDES WHETHER A VENDOR-FAILED READ IS WRITTEN;
 *   - `--CLIP_ADAPTERS` REWRITES THE CLIPPED TAIL AS `A` AT QUALITY ONE rather than shortening
 *     the record, so the length field does not move;
 *   - THE NAME PREFIX IS STRIPPED FROM THE FRONT of every read name, and a name that does not
 *     start with it is written whole;
 *   - AND `--OUTPUT_FILE_PREFIX` IS MUTUALLY EXCLUSIVE WITH `--FLOWCELL_BARCODE` AND `--LANE`,
 *     which together build the same prefix.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     bfq\t<case>/<file>\t<the file's bytes, hex>
 *     plain\t<case>/<file>\t<the same file gunzipped, hex>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: BamToBfqDump
 */

import htsjdk.samtools.*;
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

public class BamToBfqDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static String hex(final byte[] bytes) {
        final StringBuilder text = new StringBuilder();
        for (final byte value : bytes) {
            text.append(String.format("%02x", value));
        }
        return text.toString();
    }

    /** One read of the fixture. */
    record Read(String name, String bases, String qualities, int flags) {}

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
                record.setBaseQualityString(spec.qualities());
                record.setFlags(spec.flags() | 0x4);
                record.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                writer.addAlignment(record);
            }
        }
    }

    /** A pair of unmapped reads, which is what a paired run wants. */
    static List<Read> pair(final String name, final String first, final String second,
                           final int extra) {
        return List.of(
                new Read(name, first, "I".repeat(first.length()), 0x1 | 0x40 | 0x8 | extra),
                new Read(name, second, "I".repeat(second.length()), 0x1 | 0x80 | 0x8 | extra));
    }

    static void run(final String name, final List<Read> reads, final boolean paired,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("bamtobfq");
        final Path in = dir.resolve("in.bam");
        try {
            writeBam(in, reads);
        } catch (final Exception e) {
            // The tool has a branch for a base no code knows, and a BAM cannot carry one: the
            // writer refuses the record before the tool ever sees it. The refusal is recorded
            // here, which is what says the branch is unreachable through this input format.
            emit("error", name, "fixture " + e.getClass().getName() + ":"
                    + String.valueOf(e.getMessage()));
            return;
        }
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final Path out = Files.createDirectory(dir.resolve("out"));
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "ANALYSIS_DIR=" + out, "PAIRED_RUN=" + paired));
        if (Arrays.stream(extra).noneMatch(a -> a.startsWith("OUTPUT_FILE_PREFIX")
                || a.startsWith("FLOWCELL_BARCODE") || a.startsWith("F="))) {
            argv.add("OUTPUT_FILE_PREFIX=run");
        }
        argv.addAll(Arrays.asList(extra));
        final java.io.ByteArrayOutputStream errBytes = new java.io.ByteArrayOutputStream();
        final java.io.PrintStream realErr = System.err;
        try {
            final int code;
            try {
                System.setErr(new java.io.PrintStream(errBytes, true, StandardCharsets.UTF_8));
                code = new picard.fastq.BamToBfq().instanceMain(argv.toArray(new String[0]));
            } finally {
                System.err.flush();
                System.setErr(realErr);
            }
            if (code != 0) {
                emit("error", name, "exit " + code + " " + reason(errBytes.toString(StandardCharsets.UTF_8)));
                return;
            }
        } catch (final Exception e) {
            System.setErr(realErr);
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
            final byte[] raw = Files.readAllBytes(out.resolve(file));
            emit("bfq", name + "/" + file, hex(raw));
            // The file is GZIP, which `IOUtil.openFileForWriting` applies to a `.bfq` the same way
            // it applies it to a `.gz`. Its bytes are deterministic here because the header's
            // timestamp is zero, but reproducing them means reproducing a DEFLATE stream, which is
            // the programme's own hard problem. The payload is what a port can claim, so it is
            // recorded beside the file: the records themselves, uncompressed.
            emit("plain", name + "/" + file, hex(gunzip(raw)));
        }
    }

    /** The file's payload, which is what a port compares without owning a deflater. */
    static byte[] gunzip(final byte[] raw) throws Exception {
        try (final java.util.zip.GZIPInputStream in = new java.util.zip.GZIPInputStream(
                new java.io.ByteArrayInputStream(raw))) {
            return in.readAllBytes();
        }
    }

    /** The first line of a refusal, which the reference prints above its usage. */
    static String reason(final String stderr) {
        for (final String line : stderr.split("\n", -1)) {
            final String trimmed = line.trim();
            if (!trimmed.isEmpty()) {
                return trimmed;
            }
        }
        return "";
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // One pair, which is the shape everything else is measured against.
        run("a-pair", pair("read1", "ACGT", "TGCA", 0), true);
        // The four bases and their lower case, which is the encoding itself.
        run("every-base", pair("read1", "ACGTacgt", "TTTTAAAA", 0), true);
        // An N, whose quality is not the read's.
        run("an-n-inside-the-seed", pair("read1", "ANGT", "ACGT", 0), true);
        run("an-n-outside-the-seed", List.of(
                new Read("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGN",
                        "I".repeat(40), 0x1 | 0x40 | 0x8),
                new Read("read1", "ACGT", "IIII", 0x1 | 0x80 | 0x8)), true);
        // A base no code knows.
        run("an-unknown-base", pair("read1", "ACXT", "ACGT", 0), true);

        // The qualities, which are the low six bits.
        run("a-low-quality", List.of(
                new Read("read1", "ACGT", "!!!!", 0x1 | 0x40 | 0x8),
                new Read("read1", "ACGT", "IIII", 0x1 | 0x80 | 0x8)), true);

        // Two pairs, chunked and truncated.
        final List<Read> two = new ArrayList<>();
        two.addAll(pair("read1", "ACGT", "TGCA", 0));
        two.addAll(pair("read2", "GGGG", "CCCC", 0));
        run("two-pairs", two, true);
        run("two-pairs-chunked", two, true, "READ_CHUNK_SIZE=1");
        run("two-pairs-truncated", two, true, "READS_TO_ALIGN=1");

        // A single-end run, which writes half as many files.
        run("single-end", List.of(
                new Read("read1", "ACGT", "IIII", 0x4),
                new Read("read2", "TTTT", "IIII", 0x4)), false);

        // The trimming and the clipping, which do different things to the length.
        run("bases-to-write", pair("read1", "ACGTACGT", "TGCATGCA", 0), true, "BASES_TO_WRITE=4");
        run("clip-adapters", pair("read1", "ACGTACGT", "TGCATGCA", 0), true, "CLIP_ADAPTERS=true");

        // A vendor-failed read, on both sides of the flag.
        run("a-non-pf-read", pair("read1", "ACGT", "TGCA", 0x200), true);
        run("a-non-pf-read-included", pair("read1", "ACGT", "TGCA", 0x200), true,
                "INCLUDE_NON_PF_READS=true");

        // The name prefix, stripped from the front and left alone when it does not match.
        run("a-stripped-name", pair("RUN:read1", "ACGT", "TGCA", 0), true,
                "READ_NAME_PREFIX=RUN:");
        run("a-name-the-prefix-does-not-match", pair("read1", "ACGT", "TGCA", 0), true,
                "READ_NAME_PREFIX=OTHER:");

        // The two ways of naming the output, which may not be given together.
        run("a-flowcell-and-a-lane", pair("read1", "ACGT", "TGCA", 0), true,
                "F=30PYMAAXX", "L=3");
        run("both-ways-of-naming", pair("read1", "ACGT", "TGCA", 0), true,
                "F=30PYMAAXX", "OUTPUT_FILE_PREFIX=run");

        System.out.print(buf);
    }
}
