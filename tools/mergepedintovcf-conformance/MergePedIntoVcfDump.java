/*
 * MergePedIntoVcf's VCF, taken from the reference.
 *
 * The tool takes a genotyping-array VCF and a zCall PED beside it, and writes the VCF back with
 * the zCall genotype added to each record. Every input is text, so what is measured is the merge
 * and nothing else.
 *
 * Fourteen behaviours this is built to catch.
 *
 *   - THE ORIGINAL GENOTYPE IS KEPT AS `GTA` AND THE ZCALL ONE ADDED AS `GTZ`, so a record whose
 *     two callers disagree carries both;
 *   - THE PED'S ALLELES ARE `A`, `B` OR `0` AND NOT BASES: each is looked up in the record's own
 *     ALLELE_A and ALLELE_B, and a `0` is a no-call;
 *   - THE LOOKED-UP ALLELE IS ALWAYS BUILT AS NON-REFERENCE, so a PED calling `A` where the
 *     record's ALLELE_A is its REF makes an allele the context does not hold and the run fails;
 *   - THE TOOL RECALCULATES AC, AF AND AN, and refuses to write them if the header does not
 *     declare them;
 *   - THE PED'S ALLELES ARE MATCHED BY THE MAP'S SNP NAME AND NOT BY POSITION, the two files
 *     being read in step;
 *   - THE FIRST SIX PED FIELDS ARE IGNORED whatever they hold;
 *   - THE THRESHOLDS FILE ADDS `ZTHRESH_X` AND `ZTHRESH_Y` to the records whose ID it names, and
 *     leaves the others without them;
 *   - A THRESHOLD PAIR OF `NA` IS WRITTEN AS THE MISSING VALUE rather than dropped;
 *   - ONE `NA` OF A PAIR IS REFUSED, by a message saying they must both exist or both not;
 *   - THE ZCALL VERSION AND THE THRESHOLDS FILE'S NAME GO INTO THE HEADER;
 *   - A PED NAMING MORE THAN ONE SAMPLE IS REFUSED, and so is an allele of more than one
 *     character;
 *   - A VCF NAMING MORE THAN ONE SAMPLE IS REFUSED;
 *   - A VCF WHOSE GENOTYPES CARRY NOTHING BUT `GT` CANNOT BE PROCESSED AT ALL: the tool puts its
 *     two new fields into the map `getExtendedAttributes` answers, and that map is IMMUTABLE when
 *     the genotype has no extended attributes, so it throws where it stands;
 *   - EVERY ONE OF THOSE FAILURES IS SWALLOWED: doWork catches Exception, prints the stack trace
 *     and returns ZERO, so the tool reports success and writes no file;
 *   - AND THE THRESHOLDS MAP IS STATIC, so every run after the first still holds every earlier
 *     run's thresholds and writes them onto records its own thresholds file never named. The
 *     effect shows from the second case onward and is plainest in the last, whose file names rs3
 *     alone and whose output carries thresholds on all three.
 *
 * The cases run in the order written, and the last one depends on that order.
 *
 * Output:
 *
 *     in\t<name>\t<that input file, escaped>
 *     out\t<case>\t<the VCF's records and its added header lines, escaped>
 *     code\t<case>\t<exit code>
 *     none\t<case>\t<no output was written>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class MergePedIntoVcfDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /**
     * A single-sample VCF of three SNPs, each with an ID the PED and MAP name too.
     *
     * The genotypes carry a non-standard FORMAT field on purpose. The tool puts GTA and GTZ into
     * the map `Genotype.getExtendedAttributes()` answers, and that map is IMMUTABLE when the
     * genotype has none: a VCF of nothing but GT makes the tool throw where it stands.
     */
    static String vcf(final List<String> samples, final boolean extendedAttribute) {
        final List<String> lines = new ArrayList<>();
        lines.add("##fileformat=VCFv4.2");
        lines.add("##contig=<ID=chr1,length=100000>");
        lines.add("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        lines.add("##INFO=<ID=ALLELE_A,Number=1,Type=String,Description=\"A allele\">");
        lines.add("##INFO=<ID=ALLELE_B,Number=1,Type=String,Description=\"B allele\">");
        lines.add("##FORMAT=<ID=IGC,Number=1,Type=Float,Description=\"GenCall confidence\">");
        // The tool recalculates these three and refuses to write them undeclared.
        lines.add("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">");
        lines.add("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">");
        lines.add("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">");
        final List<String> columns = new ArrayList<>(List.of(
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"));
        columns.addAll(samples);
        lines.add(String.join("\t", columns));
        final String[][] records = {
                {"1000", "rs1", "A", "C"},
                {"2000", "rs2", "G", "T"},
                {"3000", "rs3", "C", "A"},
        };
        final String[] genotypes = {"0/0", "0/1", "1/1"};
        for (int i = 0; i < records.length; i++) {
            // The zCall translation reads ALLELE_A and ALLELE_B off the record: a PED allele of
            // `A` becomes whichever base ALLELE_A names, and `B` whichever ALLELE_B names.
            final String info = "ALLELE_A=" + records[i][2] + ";ALLELE_B=" + records[i][3];
            final List<String> row = new ArrayList<>(List.of(
                    "chr1", records[i][0], records[i][1], records[i][2], records[i][3],
                    ".", ".", info, extendedAttribute ? "GT:IGC" : "GT"));
            for (int s = 0; s < samples.size(); s++) {
                row.add(extendedAttribute ? genotypes[i] + ":0.7" : genotypes[i]);
            }
            lines.add(String.join("\t", row));
        }
        return String.join("\n", lines) + "\n";
    }

    /** The MAP: one line per SNP, the name in the second field. */
    static String map() {
        return String.join("\n",
                "chr1\trs1\t0\t1000",
                "chr1\trs2\t0\t2000",
                "chr1\trs3\t0\t3000") + "\n";
    }

    /** The PED: six ignored fields, then two allele characters per SNP. */
    static String ped(final String alleles) {
        return "FAM\tIND\t0\t0\t0\t-9\t" + alleles + "\n";
    }

    static String thresholds(final List<String> lines) {
        return String.join("\n", lines) + "\n";
    }

    /** The VCF's records and the header lines the tool added. */
    static String interesting(final String vcf) {
        final List<String> kept = new ArrayList<>();
        for (final String line : vcf.split("\n", -1)) {
            if (line.isEmpty()) {
                continue;
            }
            if (line.startsWith("##")) {
                if (line.contains("zCall") || line.contains("ZTHRESH") || line.contains("GTA")
                        || line.contains("GTZ")) {
                    kept.add(line);
                }
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final List<String> samples, final String pedAlleles,
                    final List<String> thresholdLines, final String... extra) throws Exception {
        run(name, samples, pedAlleles, thresholdLines, true, extra);
    }

    static void run(final String name, final List<String> samples, final String pedAlleles,
                    final List<String> thresholdLines, final boolean extendedAttribute,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("mergepedintovcf");
        final Path vcf = dir.resolve("in.vcf");
        Files.writeString(vcf, vcf(samples, extendedAttribute), StandardCharsets.UTF_8);
        final Path pedFile = dir.resolve("in.ped");
        Files.writeString(pedFile, ped(pedAlleles), StandardCharsets.UTF_8);
        final Path mapFile = dir.resolve("in.map");
        Files.writeString(mapFile, map(), StandardCharsets.UTF_8);
        final Path thresholdsFile = dir.resolve("thresholds-" + name + ".txt");
        Files.writeString(thresholdsFile, thresholds(thresholdLines), StandardCharsets.UTF_8);
        final Path out = dir.resolve("out.vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "VCF=" + vcf, "PED=" + pedFile, "MAP=" + mapFile,
                "ZCALL_T_FILE=" + thresholdsFile, "ZCALL_VERSION=1.0", "O=" + out));
        argv.addAll(List.of(extra));
        final int code = new picard.arrays.MergePedIntoVcf()
                .instanceMain(argv.toArray(new String[0]));
        emit("code", name, Integer.toString(code));
        if (!Files.exists(out)) {
            emit("none", name, "no output was written");
            return;
        }
        emit("out", name, interesting(Files.readString(out, StandardCharsets.UTF_8))
                .replace(dir.toString(), "<dir>"));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final List<String> one = List.of("sample1");
        emit("in", "vcf", vcf(one, true));
        emit("in", "map", map());
        emit("in", "ped", ped("B B B B B B"));

        // The PED agreeing with the VCF at every SNP, and one thresholds line.
        run("agreeing", one, "B B B B B B", List.of("rs1\t0.5\t0.6"));

        // The PED disagreeing at rs1, so GTA and GTZ differ there.
        run("disagreeing", one, "B B 0 0 B B", List.of("rs1\t0.5\t0.6"));

        // Thresholds for every SNP, and for none.
        run("all-thresholds", one, "B B B B B B",
                List.of("rs1\t0.5\t0.6", "rs2\t0.7\t0.8", "rs3\t0.9\t1.0"));

        // A pair of NA, which is written as the missing value.
        run("na-thresholds", one, "B B B B B B", List.of("rs2\tNA\tNA"));

        // The first six PED fields, which are ignored: a different family and sex.
        run("other-ped-header", one, "B B B B B B", List.of("rs1\t0.5\t0.6"));

        // One NA of a pair, which is refused.
        run("half-na", one, "B B B B B B", List.of("rs1\tNA\t0.6"));

        // A PED calling the `A` allele where the record's ALLELE_A is its REF, which the tool
        // builds as a NON-reference allele and then cannot find in the context.
        run("a-allele-is-the-reference", one, "A A B B B B", List.of("rs1\t0.5\t0.6"));

        // A PED allele of more than one character.
        run("long-allele", one, "BB B B B B B", List.of("rs1\t0.5\t0.6"));

        // A VCF naming two samples.
        run("two-vcf-samples", List.of("sample1", "sample2"), "B B B B B B",
                List.of("rs1\t0.5\t0.6"));

        // A VCF whose genotypes carry nothing but GT, which the tool cannot process at all.
        run("no-extended-attributes", one, "B B B B B B", List.of("rs1\t0.5\t0.6"), false);

        // The static thresholds map: this run's file names rs3 alone, and the records for rs1 and
        // rs2 still carry the thresholds every run before it put in the map.
        run("static-map-leak", one, "B B B B B B", List.of("rs3\t0.1\t0.2"));

        System.out.print(buf);
    }
}
