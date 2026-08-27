/*
 * ConvertHaplotypeDatabaseToVcf's VCF, taken from the reference.
 *
 * A haplotype database is a table of major and minor alleles grouped into blocks by an anchor
 * SNP. The VCF is that table read against the reference, which is what is supposed to decide
 * which of the two alleles is REF. What is measured is which one actually is, what the block
 * structure adds, and which tables are refused.
 *
 * Ten behaviours this is built to catch.
 *
 *   - REF IS THE ALLELE THAT DISAGREES WITH THE REFERENCE, in both directions: the tool asks
 *     which allele matches the reference base and then writes THAT one as ALT. A row whose major
 *     allele is the reference base and a row whose minor allele is come out with the same REF and
 *     ALT, which is the wrong way round in each;
 *   - AF IS CONSISTENT WITH THE ALT IT WROTE and not with the table: a row whose major allele
 *     matches carries `AF=1-MAF`, a row whose minor allele matches carries `AF=MAF`, and in both
 *     that number is the frequency of the allele written as ALT;
 *   - A BLOCK OF ONE SNP IS NOT PHASED and carries no phase set;
 *   - A BLOCK OF MORE THAN ONE IS PHASED and every one of its rows carries the same PS;
 *   - THAT PS IS THE POSITION OF THE BLOCK'S FIRST SNP BY POSITION, not of the row named as the
 *     anchor: a block whose anchor sits second is phased on the first all the same;
 *   - A ROW WHOSE MAJOR ALLELE MATCHES, INSIDE A PHASED BLOCK, HAS ITS GENOTYPE REVERSED to
 *     `1|0`, while the same row in a block of one keeps `0/1`;
 *   - THE RECORDS COME OUT IN DICTIONARY ORDER whatever order the table listed them in;
 *   - THE `##source` LINE IS THE LITERAL `HaplotypeMap::writeAsVcf`, while the `##reference` line
 *     the code sets to that same literal is overwritten by the file URI of the FASTA;
 *   - A ROW WHOSE TWO ALLELES BOTH DISAGREE WITH THE REFERENCE IS REFUSED, by a message that
 *     prints the SNP as `<contig>:<position>` and neither allele;
 *   - AND A TABLE WITH NO `@` HEADER, A SHORT ROW, AND AN ANCHOR THAT NAMES NO ROW ARE EACH
 *     REFUSED BY A DIFFERENT MESSAGE, while a table with a header and no rows at all writes a
 *     VCF holding nothing.
 *
 * Output:
 *
 *     db\t<name>\t<that haplotype database, escaped>
 *     reference\t<name>\t<that FASTA, escaped>
 *     out\t<case>\t<the VCF, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class ConvertHaplotypeDatabaseToVcfDump {

    static final StringBuilder buf = new StringBuilder();

    /** chr1 is 60 bases of a repeating pattern, so every position's base is known by index. */
    static final String CHR1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    static final String CHR2 = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC";

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** The reference FASTA, whose two contigs are in the order the dictionary will hold them. */
    static String fasta() {
        return ">chr1\n" + CHR1 + "\n>chr2\n" + CHR2 + "\n";
    }

    /** One row of the table, whose anchor may be its own name or another's. */
    static String row(final String chrom, final int pos, final String name, final char major,
                      final char minor, final String maf, final String anchor) {
        return String.join("\t", chrom, Integer.toString(pos), name, String.valueOf(major),
                String.valueOf(minor), maf, anchor == null ? "" : anchor);
    }

    static String database(final List<String> rows) {
        final List<String> lines = new ArrayList<>();
        lines.add("@HD\tVN:1.6\tSO:coordinate");
        lines.add("@SQ\tSN:chr1\tLN:60");
        lines.add("@SQ\tSN:chr2\tLN:60");
        lines.add("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS");
        lines.addAll(rows);
        return String.join("\n", lines) + "\n";
    }

    static void run(final String name, final String database) throws Exception {
        final Path dir = Files.createTempDirectory("haplotypedbtovcf");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
                "R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path in = dir.resolve("in.txt");
        Files.writeString(in, database, StandardCharsets.UTF_8);
        final Path out = dir.resolve("out.vcf");
        try {
            final int code = new picard.fingerprint.ConvertHaplotypeDatabaseToVcf()
                    .instanceMain(new String[]{"I=" + in, "O=" + out, "R=" + reference});
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

        // chr1 is ACGT repeating from position 1, so position 1 is A, 2 is C, 3 is G, 4 is T.
        emit("reference", "ref", fasta());

        // One SNP whose MAJOR allele is the reference base, which the tool writes as ALT.
        final String plain = database(List.of(row("chr1", 1, "rs1", 'A', 'C', "0.25", "rs1")));
        emit("db", "plain", plain);
        run("major-is-reference", plain);

        // One SNP whose MINOR allele is the reference base. The two alleles come out in the
        // same order as the case above and only the frequency tells the two rows apart.
        final String swapped = database(List.of(row("chr1", 1, "rs1", 'C', 'A', "0.25", "rs1")));
        emit("db", "swapped", swapped);
        run("minor-is-reference", swapped);

        // A block of three, whose middle row's minor allele matches and whose anchor is named
        // SECOND, so the phase set is the FIRST row's position and not the anchor's.
        final String block = database(List.of(
                row("chr1", 1, "rs1", 'A', 'C', "0.10", "rs2"),
                row("chr1", 2, "rs2", 'G', 'C', "0.20", "rs2"),
                row("chr1", 3, "rs3", 'G', 'T', "0.30", "rs2")));
        emit("db", "block", block);
        run("phased-block", block);

        // The same three listed in the opposite order, and on two contigs: the output sorts.
        final String unsorted = database(List.of(
                row("chr2", 5, "rs9", 'G', 'C', "0.40", "rs9"),
                row("chr1", 3, "rs3", 'G', 'A', "0.30", "rs3"),
                row("chr1", 1, "rs1", 'A', 'C', "0.10", "rs1")));
        emit("db", "unsorted", unsorted);
        run("sorted-output", unsorted);

        // A row neither of whose alleles is the reference base, which is C at position 2.
        final String disagreeing = database(List.of(row("chr1", 2, "rs1", 'A', 'T', "0.25", "rs1")));
        emit("db", "disagreeing", disagreeing);
        run("neither-allele-matches", disagreeing);

        // A table with no @ header at all.
        run("no-header", "#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\n");

        // A row with too few fields.
        run("short-row", database(List.of("chr1\t1\trs1\tA")));

        // A block whose anchor is never itself declared as a row.
        run("dangling-anchor", database(List.of(row("chr1", 1, "rs1", 'A', 'C', "0.10", "rsX"))));

        // A table with a header and no rows at all.
        run("no-rows", database(List.of()));

        System.out.print(buf);
    }
}
