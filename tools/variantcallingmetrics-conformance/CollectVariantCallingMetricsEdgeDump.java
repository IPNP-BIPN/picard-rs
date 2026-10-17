/*
 * CollectVariantCallingMetrics on the edge cases of a VCF, taken from the reference.
 *
 * Which counter a site reaches is `VariantContext.getType()`, not the allele lengths alone, and
 * the per-sample rows come out in the order of the accumulator's `HashMap`. Five behaviours this is
 * built to catch:
 *
 *   - A MULTIALLELIC INDEL IS A COMPLEX INDEL, and only a biallelic one sharing its first base is
 *     a simple insertion or deletion;
 *   - A SYMBOLIC, MNP OR MIXED SITE IS NEITHER A SNP NOR AN INDEL, so it reaches no counter, not
 *     even a filtered one, and a site whose ALT is `.` is excluded;
 *   - A HALF-CALLED GENOTYPE (`.|1`) IS MIXED, not a no-call: it is counted for its sample,
 *     though it is neither het nor hom-var;
 *   - THE DBSNP INDEL BITSET COVERS A KNOWN DELETION'S WHOLE SPAN, so an insertion two bases into
 *     it is known although the dbSNP file has nothing at its position;
 *   - AND THE DETAIL ROWS ARE IN HASHMAP ORDER, NA3 before NA1, which the table is NOT sorted to
 *     hide.
 *
 * Output:
 *
 *     in\t<name>\t<that vcf, escaped>
 *     detail\t<case>\t<the detail table, rows in written order, escaped>
 *     summary\t<case>\t<the summary table, escaped>
 */

import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.vcf.VCFCodec;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class CollectVariantCallingMetricsEdgeDump {
    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final List<String> EDGE = List.of(
            "##fileformat=VCFv4.2",
            "##FILTER=<ID=LowQual,Description=\"Low quality\">",
            "##FILTER=<ID=q10,Description=\"Quality below 10\">",
            "##FILTER=<ID=s50,Description=\"Less than 50% of samples\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
            "##FORMAT=<ID=FT,Number=.,Type=String,Description=\"Genotype-level filter\">",
            "##FORMAT=<ID=GF,Number=1,Type=Float,Description=\"A float\">",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
            "##FORMAT=<ID=GS,Number=1,Type=String,Description=\"A string\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype quality\">",
            "##FORMAT=<ID=MI,Number=.,Type=Integer,Description=\"Unbounded ints\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">",
            "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">",
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total alleles\">",
            "##INFO=<ID=BIG,Number=1,Type=Integer,Description=\"Large integer\">",
            "##INFO=<ID=CH,Number=1,Type=Character,Description=\"Char\">",
            "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
            "##INFO=<ID=F2,Number=2,Type=Float,Description=\"Two floats\">",
            "##INFO=<ID=GI,Number=G,Type=Integer,Description=\"Per genotype ints\">",
            "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2\">",
            "##INFO=<ID=MID,Number=.,Type=Integer,Description=\"Unbounded ints\">",
            "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mapping quality\">",
            "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Samples with data\">",
            "##INFO=<ID=RF,Number=R,Type=Float,Description=\"Per allele floats\">",
            "##INFO=<ID=STRS,Number=.,Type=String,Description=\"Strings\">",
            "##contig=<ID=20,length=62435964>",
            "##contig=<ID=21,length=48129895>",
            "##contig=<ID=X,length=155270560>",
            "##reference=file:///ref.fa",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\tNA3",
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tAA=G;AC=3;AF=0.500;AN=6;BIG=70000;CH=Z;DB;DP=14;F2=1.50,2.50;GI=1,2,3;H2;MID=1,300,70000;MQ=60.25;NS=3;RF=0.100,0.200;STRS=a,bb,ccc\tGT:AD:DP:GF:GQ:GS:HQ:MI:PL\t0|0:1,0:1:0.500:48:x:51,51:1:0,10,100\t1|0:3,4:8:1.25:48:yy:51,51:1,2:10,0,100\t1/1:.:5:.:43",
            "20\t17330\t.\tT\tA\t3\tq10\tAF=0.017;DP=11;NS=3\tGT:DP:FT:GQ:HQ\t0|0:3:PASS:49:58,50\t0|1:5:q10:3:65,3\t0/0:3:q10;s50:41",
            "20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tAA=T;AC=1,2;AF=0.333,0.667;DB;DP=10;GI=1,2,3,4,5,6;NS=2\tGT:AD:DP:GQ:HQ:PL\t1|2:1,2,3:6:21:23,27:1,2,3,4,5,6\t2|1:.:0:2:18,2\t2/2:4,5,6:4:35:.:7,8,9,10,11,12",
            "20\t1230237\t.\tT\t.\t47\t.\tAA=T;DP=13;NS=3\tGT:DP:GQ:HQ\t0|0:7:54:56,60\t0|0:4:48:51,51\t0/0:2:61",
            "20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tAA=G;DP=9;MID=5;NS=3\tGT:DP:GQ\t0/1:4:35\t0/2:2:17\t1/1:3:40",
            "20\t1234568\t.\tG\tGA,GAC,GACT\t.\tq10;s50\tBIG=-2147483647;DP=200\tGT:GQ:PL\t0/3:99:0,1,2,3,4,5,6,7,8,9\t.\t1|3",
            "20\t1240000\t.\tAC\tGT\t50\tPASS\t.\tGT\t0/1\t0/0\t1/1",
            "20\t1250000\t.\tA\tC,AT\t50\tPASS\t.\tGT\t0/1\t0/2\t1/2",
            "21\t100\t.\tA\tC\t0\tLowQual\t.\tGT:DP\t./.\t0:3\t.|1:2",
            "21\t200\tid1;id2\tACGT\tA\t0\tPASS\tAF=0.00;DP=0;MQ=1.000e-03\tGT\t0/1\t1/1\t0/0",
            "21\t202\t.\tG\tGT\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0",
            "X\t10\t.\tC\tCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t123.46\tPASS\tCH=.;DP=.;STRS=onlyone\tGT:GF:GS\t1:.:hello\t0:3.14\t.:0.00:z",
            "X\t20\t.\tA\t<DEL>\t99\tPASS\t.\tGT:AD\t0/1:5,5\t0/0:10,0\t1/1:0,10");

    /** The known sites: the input without the insertion inside the deletion's span. */
    static final String KNOWN_ONLY_BY_SPAN = "21\t202\t";

    /** A table without its comments, its rows in the order the tool wrote them. */
    static String table(final String text) {
        final List<String> rows = new ArrayList<>();
        for (final String line : text.split("\n", -1)) {
            if (!line.startsWith("#") && !line.isEmpty()) {
                rows.add(line);
            }
        }
        return String.join("\n", rows);
    }

    static void index(final Path vcf) throws Exception {
        IndexFactory.createLinearIndex(vcf.toFile(), new VCFCodec())
                .writeBasedOnFeatureFile(vcf.toFile());
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Files.createTempDirectory("variantcallingmetricsedge");
        final String vcf = String.join("\n", EDGE) + "\n";
        final List<String> knownLines = new ArrayList<>();
        for (final String line : EDGE) {
            if (!line.startsWith(KNOWN_ONLY_BY_SPAN)) {
                knownLines.add(line);
            }
        }
        final String known = String.join("\n", knownLines) + "\n";
        emit("in", "edge", vcf);
        emit("in", "dbsnp", known);

        final Path in = dir.resolve("in.vcf");
        Files.writeString(in, vcf, StandardCharsets.UTF_8);
        index(in);
        final Path db = dir.resolve("dbsnp.vcf");
        Files.writeString(db, known, StandardCharsets.UTF_8);
        index(db);
        final Path out = dir.resolve("out");
        final int code = new picard.vcf.CollectVariantCallingMetrics()
                .instanceMain(new String[] {"I=" + in, "DBSNP=" + db, "O=" + out});
        if (code != 0) {
            throw new IllegalStateException("exit " + code);
        }
        emit("detail", "edge", table(Files.readString(
                Path.of(out + ".variant_calling_detail_metrics"), StandardCharsets.UTF_8)));
        emit("summary", "edge", table(Files.readString(
                Path.of(out + ".variant_calling_summary_metrics"), StandardCharsets.UTF_8)));
        System.out.print(buf);
    }
}
