/*
 * SortGff's output, taken from the reference.
 *
 * A GFF3 file sorted by contig and then by start. What is measured is the order the comparator
 * produces, what a sequence dictionary changes about it, and what survives the round trip through
 * the codec.
 *
 * Nine behaviours this is built to catch.
 *
 *   - WITHOUT A DICTIONARY THE CONTIGS SORT LEXICOGRAPHICALLY, so `chr10` comes before `chr2`;
 *   - WITH ONE THEY SORT BY ITS OWN ORDER, which is what a dictionary is for here;
 *   - A CONTIG THE DICTIONARY DOES NOT NAME GETS INDEX -1 AND SORTS FIRST, before every contig
 *     the dictionary does name, rather than being refused or left where it was;
 *   - WITHIN A CONTIG THE ORDER IS THE START ALONE, so two features that start together keep the
 *     order they were read in and their ends are never compared;
 *   - THE VERSION DIRECTIVE IS THE CODEC'S OWN AND NOT THE INPUT'S: a file that opens with
 *     `##gff-version 3.1.26` comes back opening with 3.1.25;
 *   - THE COMMENT LINES ARE CARRIED OVER;
 *   - A PARENT IS ALLOWED TO FOLLOW ITS CHILD in the input and the sort puts each wherever its
 *     own coordinates say;
 *   - --nRecordsInMemory CHANGES NOTHING ABOUT THE OUTPUT, only where the sort holds its records;
 *   - AND A FILE WITH NO FEATURES IS REFUSED exactly as a file that is not GFF at all is: the
 *     codec's `canDecode` wants a feature and not only a directive, so an empty file is an
 *     IllegalArgumentException naming the input.
 *
 * Output:
 *
 *     gff\t<name>\t<that input file, escaped>
 *     out\t<case>\t<the sorted file, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SortGffDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static String feature(final String contig, final int start, final int end, final String id) {
        return String.join("\t", contig, "fixture", "gene", Integer.toString(start),
                Integer.toString(end), ".", "+", ".", "ID=" + id);
    }

    static String gff(final List<String> features) {
        final List<String> lines = new ArrayList<>();
        lines.add("##gff-version 3.1.26");
        lines.add("#a comment the sorter carries over");
        lines.addAll(features);
        return String.join("\n", lines) + "\n";
    }

    static Path dictionary(final Path dir, final List<String> contigs) throws Exception {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for (final String contig : contigs) {
            dict.addSequence(new SAMSequenceRecord(contig, 1000000));
        }
        header.setSequenceDictionary(dict);
        final Path path = dir.resolve("ref.dict");
        try (final Writer writer = Files.newBufferedWriter(path)) {
            new SAMTextHeaderCodec().encode(writer, header);
        }
        return path;
    }

    static void run(final String name, final String input, final List<String> contigs,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("sortgff");
        final Path in = dir.resolve("in.gff3");
        Files.writeString(in, input, StandardCharsets.UTF_8);
        final Path out = dir.resolve("out.gff3");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + out));
        if (contigs != null) {
            argv.add("SD=" + dictionary(dir, contigs));
        }
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.annotation.SortGff()
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
        emit("out", name, Files.readString(out).replace(dir.toString(), "<dir>"));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Contigs whose lexicographic order is not their numeric one, and two features that
        // start together on one of them.
        final String mixed = gff(List.of(
                feature("chr2", 500, 600, "b1"),
                feature("chr10", 100, 200, "a1"),
                feature("chr1", 900, 1000, "c1"),
                feature("chr1", 100, 300, "c2"),
                feature("chr1", 100, 200, "c3")));
        emit("gff", "mixed", mixed);
        run("lexicographic", mixed, null);
        run("dictionary-order", mixed, List.of("chr1", "chr2", "chr10"));
        // A dictionary that names only some of the contigs, so the rest get index -1.
        run("dictionary-partial", mixed, List.of("chr2", "chr10"));

        // A parent written after its child.
        final String outOfOrder = gff(List.of(
                String.join("\t", "chr1", "fixture", "exon", "200", "300", ".", "+", ".",
                        "ID=e1;Parent=g1"),
                String.join("\t", "chr1", "fixture", "gene", "100", "400", ".", "+", ".",
                        "ID=g1")));
        emit("gff", "child-first", outOfOrder);
        run("child-before-parent", outOfOrder, null);

        // Attributes the codec has to escape on the way out.
        final String escaped = gff(List.of(
                String.join("\t", "chr1", "fixture", "gene", "100", "200", ".", "+", ".",
                        "ID=g1;Note=a%2Cb;Name=with space")));
        emit("gff", "escaped", escaped);
        run("escaped-attributes", escaped, null);

        // A file with no features at all, which the codec refuses rather than sorting.
        final String empty = gff(List.of());
        emit("gff", "empty", empty);
        run("no-features", empty, null);

        // A file the codec cannot read.
        run("not-a-gff", "this is not a gff\n", null);

        // The in-memory record count, which spills to disk when it is small.
        run("spill-to-disk", mixed, null, "nRecordsInMemory=1");

        System.out.print(buf);
    }
}
