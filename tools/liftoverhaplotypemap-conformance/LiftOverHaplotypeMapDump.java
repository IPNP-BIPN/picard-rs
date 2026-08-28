/*
 * LiftOverHaplotypeMap's output, taken from the reference.
 *
 * A haplotype database whose SNPs are moved onto another reference through a UCSC chain file.
 * What is measured is what survives the move, what is rewritten on the way out, and what the exit
 * code says when a SNP does not make it.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE EXIT CODE FOR A FAILED LIFTOVER IS 101 AND NOT 1, whatever the count of failures;
 *   - A SNP THAT DOES NOT LIFT IS DROPPED AND THE RUN STILL WRITES ITS FILE, so the output is a
 *     shorter map rather than none at all, and a database whose every SNP fails leaves a file
 *     holding its header and its column line alone;
 *   - THE ALLELES ARE CARRIED OVER UNCHANGED ACROSS A NEGATIVE-STRAND CHAIN: the tool moves the
 *     coordinates and never complements the bases, so `chr2:50 A/C` becomes `chrB:851 A/C`;
 *   - AND SO IS THE FREQUENCY, which is the minor allele's whatever the strand;
 *   - THE FREQUENCY IS REFORMATTED ON THE WAY OUT, `0.10` coming back as `0.1`;
 *   - THE ANCHOR COLUMN IS REWRITTEN, not carried: the first row of a block by position gets an
 *     EMPTY anchor and every later row gets that first row's name, whatever the input named;
 *   - THE PANELS ARE CARRIED, comma-separated, and an absent one leaves the column empty;
 *   - THE HEADER OF THE OUTPUT IS THE --SEQUENCE_DICTIONARY and not the input's, M5 and UR
 *     included;
 *   - A DICTIONARY THAT DOES NOT NAME A CONTIG THE CHAIN LIFTS TO IS REFUSED before any SNP is
 *     looked at, while a contig the chain does not cover is not: its SNPs simply fail;
 *   - AND A BLOCK MAY NOT SPAN TWO CONTIGS AT ALL. It is refused when the input is READ, by
 *     HaplotypeBlock rather than by the tool, so the liftover never sees it.
 *
 * Output:
 *
 *     db\t<name>\t<that haplotype database, escaped>
 *     chain\t<name>\t<that chain file, escaped>
 *     out\t<case>\t<the lifted database, escaped>
 *     code\t<case>\t<exit code>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class LiftOverHaplotypeMapDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A chain that lifts chr1:1-500 onto chrA at +100, and chr2:1-300 onto chrB reversed. */
    static final String CHAIN = String.join("\n",
            "chain 1000 chr1 1000 + 0 500 chrA 1000 + 100 600 1",
            "500",
            "",
            "chain 1000 chr2 1000 + 0 300 chrB 1000 - 100 400 2",
            "300",
            "") + "\n";

    /** The same chain with its second block removed, so nothing lifts onto chrB. */
    static final String CHAIN_ONE_CONTIG = String.join("\n",
            "chain 1000 chr1 1000 + 0 500 chrA 1000 + 100 600 1",
            "500",
            "") + "\n";

    static String row(final String chrom, final int pos, final String name, final char major,
                      final char minor, final String maf, final String anchor, final String panels) {
        return String.join("\t", chrom, Integer.toString(pos), name, String.valueOf(major),
                String.valueOf(minor), maf, anchor == null ? "" : anchor,
                panels == null ? "" : panels);
    }

    static String database(final List<String> rows) {
        final List<String> lines = new ArrayList<>();
        lines.add("@HD\tVN:1.6\tSO:coordinate");
        lines.add("@SQ\tSN:chr1\tLN:1000");
        lines.add("@SQ\tSN:chr2\tLN:1000");
        lines.add("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS");
        lines.addAll(rows);
        return String.join("\n", lines) + "\n";
    }

    /** A FASTA of the TO contigs, from which the sequence dictionary is taken. */
    static String fasta(final List<String> contigs) {
        final StringBuilder fasta = new StringBuilder();
        for (final String contig : contigs) {
            fasta.append('>').append(contig).append('\n');
            for (int i = 0; i < 20; i++) {
                fasta.append("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC").append('\n');
            }
        }
        return fasta.toString();
    }

    static void run(final String name, final String database, final String chain,
                    final List<String> contigs) throws Exception {
        final Path dir = Files.createTempDirectory("liftoverhaplotypemap");
        final Path in = dir.resolve("in.txt");
        Files.writeString(in, database, StandardCharsets.UTF_8);
        final Path chainFile = dir.resolve("in.chain");
        Files.writeString(chainFile, chain, StandardCharsets.UTF_8);
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(contigs), StandardCharsets.UTF_8);
        final Path dict = dir.resolve("ref.dict");
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dict});
        final Path out = dir.resolve("out.txt");
        try {
            final int code = new picard.fingerprint.LiftOverHaplotypeMap().instanceMain(new String[]{
                    "I=" + in, "O=" + out, "SD=" + dict, "CHAIN=" + chainFile});
            emit("code", name, Integer.toString(code));
        } catch (final Exception e) {
            Throwable cause = e;
            while (cause.getCause() != null) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":"
                    + String.valueOf(cause.getMessage()).replace(dir.toString(), "<dir>"));
            return;
        }
        if (Files.exists(out)) {
            emit("out", name, Files.readString(out).replace(dir.toString(), "<dir>"));
        } else {
            emit("error", name, "no output file");
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        emit("chain", "two-blocks", CHAIN);
        emit("chain", "one-block", CHAIN_ONE_CONTIG);

        // A block of two SNPs that both lift, whose input anchor is named on the LAST of them.
        final String liftable = database(List.of(
                row("chr1", 100, "rs1", 'A', 'C', "0.10", "rs2", "panelA"),
                row("chr1", 200, "rs2", 'G', 'T', "0.20", "rs2", "panelA,panelB")));
        emit("db", "liftable", liftable);
        run("all-lift", liftable, CHAIN, List.of("chrA", "chrB"));

        // One SNP past the end of the chain's block, which does not lift.
        final String partly = database(List.of(
                row("chr1", 100, "rs1", 'A', 'C', "0.10", "rs1", null),
                row("chr1", 600, "rs2", 'G', 'T', "0.20", "rs1", null)));
        emit("db", "partly", partly);
        run("one-fails", partly, CHAIN, List.of("chrA", "chrB"));

        // A whole block that fails, beside one that does not.
        final String wholeBlockFails = database(List.of(
                row("chr1", 100, "rs1", 'A', 'C', "0.10", "rs1", null),
                row("chr1", 600, "rs9", 'G', 'T', "0.20", "rs9", null),
                row("chr1", 700, "rs8", 'C', 'A', "0.30", "rs9", null)));
        emit("db", "whole-block-fails", wholeBlockFails);
        run("whole-block-fails", wholeBlockFails, CHAIN, List.of("chrA", "chrB"));

        // A SNP that lifts across the negative-strand chain, whose alleles are not complemented.
        final String reversed = database(List.of(
                row("chr2", 50, "rs1", 'A', 'C', "0.10", "rs1", null)));
        emit("db", "reversed", reversed);
        run("negative-strand", reversed, CHAIN, List.of("chrA", "chrB"));

        // One block whose two SNPs sit on two different contigs, which the reader refuses.
        final String twoContigs = database(List.of(
                row("chr2", 50, "rs1", 'A', 'C', "0.10", "rs1", null),
                row("chr1", 100, "rs2", 'G', 'T', "0.20", "rs1", null)));
        emit("db", "two-contigs", twoContigs);
        run("block-across-contigs", twoContigs, CHAIN, List.of("chrA", "chrB"));

        // A dictionary that does not name chrB, which the chain lifts onto.
        run("dictionary-missing-contig", liftable, CHAIN, List.of("chrA"));

        // A chain with no block for chr2 at all: the SNP fails rather than the run being refused.
        run("no-chain-for-contig", reversed, CHAIN_ONE_CONTIG, List.of("chrA"));

        // A database whose every SNP fails.
        run("everything-fails", database(List.of(
                row("chr1", 900, "rs1", 'A', 'C', "0.10", "rs1", null))),
                CHAIN, List.of("chrA", "chrB"));

        System.out.print(buf);
    }
}
