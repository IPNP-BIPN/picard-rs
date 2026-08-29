/*
 * `FindMendelianViolations`, taken from the reference.
 *
 * A trio's genotypes have to be possible: a child's alleles come one from each parent, and a call
 * that cannot be explained that way is a violation. The tool counts them, and what makes it more
 * than an arithmetic check is everything it declines to count.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A HOMOZYGOUS CHILD OF TWO HOMOZYGOUS PARENTS OF THE OTHER ALLELE IS A VIOLATION, and the
 *     kind of violation is named: the counts are per parent as well as in total;
 *   - A CALL BELOW `--MIN_GQ` IS NOT COUNTED AT ALL, in either direction: a low-quality trio is
 *     skipped rather than called clean;
 *   - `--MIN_DP` DOES THE SAME with the depth, which is a different field;
 *   - `--MIN_HET_FRACTION` REJECTS A HETEROZYGOUS CALL whose allele depths do not look
 *     heterozygous, so a call the caller made is dropped by the checker;
 *   - THE SEX CHROMOSOMES ARE COUNTED DIFFERENTLY: a male child's X is haploid, so a het call
 *     there is its own kind of violation;
 *   - `--SKIP_CHROMS` LEAVES A CONTIG OUT, the mitochondrion by default;
 *   - AND `--VCF_DIR` WRITES THE OFFENDING RECORDS OUT, one file per violation kind.
 *
 * Output:
 *
 *     metrics\t<case>\t<the table without its comments, escaped>
 *     files\t<case>\t<the violation VCFs written, sorted, space separated>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MendelianViolationsDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;

public class MendelianViolationsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** One record: the contig, the position, and a genotype per member of the trio. */
    record Site(String contig, int position, String father, String mother, String child) {}

    /**
     * A VCF of the trio, with a genotype quality and a depth on every call.
     *
     * The genotype string carries its own GQ, DP and AD where a case wants them: `0/1:20:10:5,5`
     * is the whole FORMAT, and a bare `0/1` gets the defaults.
     */
    static String vcf(final List<Site> sites) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=1000>\n");
        text.append("##contig=<ID=chrX,length=1000>\n");
        text.append("##contig=<ID=chrM,length=1000>\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n");
        text.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
        text.append("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">\n");
        text.append("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfather\tmother\tchild\n");
        for (final Site site : sites) {
            text.append(site.contig()).append('\t').append(site.position())
                    .append("\t.\tA\tC\t100\tPASS\t.\tGT:GQ:DP:AD:PL\t")
                    .append(full(site.father())).append('\t')
                    .append(full(site.mother())).append('\t')
                    .append(full(site.child())).append('\n');
        }
        return text.toString();
    }

    /**
     * A genotype with the defaults filled in where the case did not name them.
     *
     * The likelihoods are not decoration: the tool reads `getPL()` when it decides whether a
     * violation is a de-novo call, and a genotype without them fails inside a worker thread with
     * `Cannot load from int array because the return value of Genotype.getPL() is null`. So every
     * genotype carries a PL consistent with its call.
     */
    static String full(final String genotype) {
        // A case that spelled out its own likelihoods keeps them: the quality the tool reads is
        // derived from the PL, so a case about a low-quality call has to say so THERE and not only
        // in the GQ column.
        if (genotype.split(":").length >= 5) {
            return genotype;
        }
        final String withDepths = genotype.contains(":") ? genotype : genotype + ":60:30:15,15";
        final String call = withDepths.split(":")[0];
        final String likelihoods = call.equals("0/0") ? "0,50,100"
                : call.equals("1/1") ? "100,50,0" : "50,0,50";
        return withDepths + ":" + likelihoods;
    }

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

    static void run(final String name, final List<Site> sites, final String childSex,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("mendelian");
        final Path in = dir.resolve("in.vcf");
        Files.writeString(in, vcf(sites), StandardCharsets.UTF_8);
        // The VCF is QUERIED per contig rather than streamed, so it needs its Tribble index.
        htsjdk.tribble.index.IndexFactory.writeIndex(
                htsjdk.tribble.index.IndexFactory.createLinearIndex(
                        in.toFile(), new htsjdk.variant.vcf.VCFCodec()),
                new java.io.File(in + ".idx"));
        // A pedigree: family, individual, father, mother, sex, phenotype.
        final Path trios = dir.resolve("trios.ped");
        Files.writeString(trios,
                "fam\tfather\t0\t0\t1\t0\nfam\tmother\t0\t0\t2\t0\n"
                        + "fam\tchild\tfather\tmother\t" + childSex + "\t0\n",
                StandardCharsets.UTF_8);
        final Path out = dir.resolve("metrics.txt");
        final Path vcfs = Files.createDirectories(dir.resolve("vcfs"));

        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "TRIOS=" + trios, "O=" + out));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        final boolean wantsVcfs = tail.remove("WITH_VCFS");
        if (wantsVcfs) {
            tail.add("VCF_DIR=" + vcfs);
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.vcf.MendelianViolations.FindMendelianViolations()
                    .instanceMain(argv.toArray(new String[0]));
            System.setOut(original);
            System.setErr(originalError);
            if (code != 0) {
                emit("error", name, "exit " + code);
                return;
            }
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
        emit("metrics", name, table(out));
        if (wantsVcfs) {
            final List<String> written = new ArrayList<>();
            try (final Stream<Path> walk = Files.walk(vcfs)) {
                walk.filter(Files::isRegularFile)
                        .forEach(path -> written.add(path.getFileName().toString()));
            }
            Collections.sort(written);
            emit("files", name, String.join(" ", written));
        }
    }

    public static void main(final String[] args) throws Exception {
        // A child that could have come from its parents, and one that could not.
        run("a-possible-child", List.of(new Site("chr1", 100, "0/0", "0/1", "0/1")), "1");
        run("a-violation", List.of(new Site("chr1", 100, "0/0", "0/0", "1/1")), "1");
        // Which parent the impossible allele came from is counted separately.
        run("a-violation-from-the-father", List.of(new Site("chr1", 100, "0/0", "1/1", "0/0")), "1");

        // A call the checker declines to look at.
        run("a-low-quality-child",
                List.of(new Site("chr1", 100, "0/0", "0/0", "1/1:10:30:0,30:10,5,0")), "1");
        run("a-low-quality-child-with-a-lower-floor",
                List.of(new Site("chr1", 100, "0/0", "0/0", "1/1:10:30:0,30:10,5,0")), "1",
                "MIN_GQ=5");
        run("a-shallow-child",
                List.of(new Site("chr1", 100, "0/0", "0/0", "1/1:60:3:0,3")), "1", "MIN_DP=10");
        // A heterozygous call whose allele depths do not look heterozygous.
        run("a-lopsided-het",
                List.of(new Site("chr1", 100, "0/1:60:30:29,1", "0/0", "0/0")), "1");

        // The sex chromosomes, and the contig that is skipped by default.
        run("a-male-child-on-the-x", List.of(new Site("chrX", 100, "0/0", "0/1", "0/1")), "1");
        run("a-female-child-on-the-x", List.of(new Site("chrX", 100, "0/0", "0/1", "0/1")), "2");
        run("the-mitochondrion", List.of(new Site("chrM", 100, "0/0", "0/0", "1/1")), "1");

        // And the records themselves, written out.
        run("with-the-offending-records",
                List.of(new Site("chr1", 100, "0/0", "0/0", "1/1")), "1", "WITH_VCFS");
        run("with-the-offending-records-and-a-tab-report",
                List.of(new Site("chr1", 100, "0/0", "0/0", "1/1")), "1", "TAB_MODE=true");

        System.out.print(buf);
    }
}
