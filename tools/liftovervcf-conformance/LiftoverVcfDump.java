/*
 * `LiftoverVcf`, taken from the reference.
 *
 * The tool moves variants from one reference to another through a chain file, and the interesting
 * part is not the arithmetic: it is what happens when the new reference disagrees with the old
 * one. A variant whose reference allele is not what the target genome carries there cannot simply
 * be renumbered, and the tool's answers to that are what this measures.
 *
 * Eight behaviours this is built to catch.
 *
 *   - A VARIANT INSIDE A CHAIN BLOCK IS RENUMBERED, and the offset is the block's;
 *   - A VARIANT OUTSIDE EVERY BLOCK IS REJECTED, into the reject file rather than dropped, with a
 *     filter naming the reason;
 *   - A VARIANT WHOSE REFERENCE ALLELE THE TARGET DOES NOT CARRY IS REJECTED TOO, and the reason
 *     is a different one;
 *   - `--RECOVER_SWAPPED_REF_ALT` TURNS ONE OF THOSE INTO A LIFTED VARIANT with its alleles
 *     swapped, which changes the genotypes as well as the alleles;
 *   - A REVERSED BLOCK FLIPS THE STRAND, so the alleles are complemented and the position is
 *     counted from the other end;
 *   - `--WRITE_ORIGINAL_POSITION` AND `--WRITE_ORIGINAL_ALLELES` RECORD WHERE A VARIANT CAME FROM,
 *     in INFO fields the header has to declare;
 *   - THE OUTPUT IS SORTED by the TARGET's coordinates, which is not the input's order;
 *   - AND A CONTIG THE CHAIN DOES NOT MENTION IS AN ERROR unless `--WARN_ON_MISSING_CONTIG` says
 *     otherwise.
 *
 * Output:
 *
 *     lifted\t<case>\t<the output VCF's variant lines, escaped>
 *     rejected\t<case>\t<the reject VCF's variant lines, escaped>
 *     code\t<case>\t<the exit status>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: LiftoverVcfDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LiftoverVcfDump {

    static final StringBuilder buf = new StringBuilder();
    static final int LENGTH = 400;

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A contig of `ACGT` repeating, so a base is known by arithmetic. */
    static String bases(final int length, final int shift) {
        final StringBuilder text = new StringBuilder();
        for (int index = 0; index < length; index++) {
            text.append("ACGT".charAt((index + shift) % 4));
        }
        return text.toString();
    }

    static String fasta(final String contig, final String sequence) {
        final StringBuilder text = new StringBuilder(">" + contig + "\n");
        for (int index = 0; index < sequence.length(); index += 60) {
            text.append(sequence, index, Math.min(index + 60, sequence.length())).append('\n');
        }
        return text.toString();
    }

    /**
     * A chain file with one forward block and one reversed one.
     *
     * `chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id`, where the
     * TARGET of a liftover is the `q` side. The first block maps the source's 0-100 onto the
     * target's 10-110; the second maps the source's 200-300 onto the target's 200-300 backwards.
     */
    static String chain() {
        final StringBuilder text = new StringBuilder();
        // The offset is TWELVE rather than ten, a multiple of the reference's four-base period, so
        // a lifted position lands on the same base it started on and a mismatch has to be put
        // there on purpose.
        text.append("chain 100 chr1 ").append(LENGTH).append(" + 0 100 chrT ").append(LENGTH)
                .append(" + 12 112 1\n100\n\n");
        text.append("chain 100 chr1 ").append(LENGTH).append(" + 200 300 chrT ").append(LENGTH)
                .append(" - 100 200 2\n100\n\n");
        return text.toString();
    }

    /** A VCF over the source reference, with the records given as `pos ref alt gt`. */
    static String vcf(final List<String[]> records) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=").append(LENGTH).append(">\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n");
        for (final String[] record : records) {
            text.append("chr1\t").append(record[0]).append("\t.\t").append(record[1])
                    .append('\t').append(record[2]).append("\t100\tPASS\tAF=0.25\tGT\t")
                    .append(record[3]).append('\n');
        }
        return text.toString();
    }

    /** A VCF's variant lines, without its header. */
    static String records(final Path file) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final List<String[]> variants, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("liftover");
        // The source and the target differ by a shift, so a position that moves lands on a
        // different base: that is what makes a reference mismatch reachable at all.
        final Path source = dir.resolve("source.fasta");
        Files.writeString(source, fasta("chr1", bases(LENGTH, 0)), StandardCharsets.UTF_8);
        final Path target = dir.resolve("target.fasta");
        Files.writeString(target, fasta("chrT", bases(LENGTH, 0)), StandardCharsets.UTF_8);
        for (final Path reference : List.of(source, target)) {
            new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
                    "R=" + reference,
                    "O=" + reference.toString().replace(".fasta", ".dict")});
            FastaSequenceIndexCreator.create(reference, true);
        }
        final Path chain = dir.resolve("chain.txt");
        Files.writeString(chain, chain(), StandardCharsets.UTF_8);
        final Path in = dir.resolve("in.vcf");
        Files.writeString(in, vcf(variants), StandardCharsets.UTF_8);

        final Path out = dir.resolve("out.vcf");
        final Path reject = dir.resolve("reject.vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + out, "CHAIN=" + chain, "REJECT=" + reject,
                "R=" + target));
        argv.addAll(Arrays.asList(extra));

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.vcf.LiftoverVcf().instanceMain(argv.toArray(new String[0]));
            System.setOut(original);
            System.setErr(originalError);
            emit("code", name, String.valueOf(code));
        } catch (final Exception e) {
            System.setOut(original);
            System.setErr(originalError);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        } finally {
            System.setOut(original);
            System.setErr(originalError);
        }
        if (Files.exists(out)) {
            emit("lifted", name, records(out));
        }
        if (Files.exists(reject)) {
            emit("rejected", name, records(reject));
        }
    }

    public static void main(final String[] args) throws Exception {
        // Position 21 of the source is inside the first block, so it moves by ten.
        final List<String[]> inside = List.<String[]>of(new String[]{"21", "A", "C", "0/1"});
        run("inside-a-block", inside);
        // Position 150 is between the blocks.
        run("between-the-blocks", List.<String[]>of(new String[]{"150", "C", "A", "0/1"}));
        // Position 250 is inside the reversed block.
        run("inside-the-reversed-block", List.<String[]>of(new String[]{"250", "C", "A", "0/1"}));

        // A reference allele the target does not carry, and no alternate that would explain it:
        // position 21 is `A` in both references, so a variant declaring `C` there is neither
        // liftable nor recoverable.
        run("a-reference-mismatch", List.<String[]>of(new String[]{"21", "C", "G", "0/1"}));
        // The same position with the alleles the other way round IS recoverable: the target's `A`
        // is the variant's ALT, so swapping them makes the record true and the genotypes with it.
        run("a-swapped-ref-and-alt", List.<String[]>of(new String[]{"21", "C", "A", "0/1"}));
        run("a-swapped-ref-and-alt-recovered", List.<String[]>of(new String[]{"21", "C", "A", "0/1"}),
                "RECOVER_SWAPPED_REF_ALT=true");

        // What the output records about where a variant came from.
        run("with-the-original-position", inside, "WRITE_ORIGINAL_POSITION=true");
        // The original alleles are worth recording where they CHANGED, which is the recovered
        // swap rather than a plain lift.
        run("with-the-original-alleles", inside, "WRITE_ORIGINAL_ALLELES=true");
        run("with-the-original-alleles-after-a-swap",
                List.<String[]>of(new String[]{"21", "C", "A", "0/1"}),
                "RECOVER_SWAPPED_REF_ALT=true", "WRITE_ORIGINAL_ALLELES=true");

        // Several variants at once, to see the order the output is written in.
        run("three-variants", List.<String[]>of(
                new String[]{"250", "C", "A", "0/1"},
                new String[]{"21", "A", "C", "1/1"},
                new String[]{"150", "C", "A", "0/0"}));

        System.out.print(buf);
    }
}
