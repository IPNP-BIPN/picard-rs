/*
 * CollectVariantCallingMetrics' metrics, taken from the reference.
 *
 * The tool walks a VCF against a dbSNP one and counts what it finds, per sample and over the file.
 * Every input is text. What is measured is which variants are counted where, and what the derived
 * columns make of the counts.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE OUTPUT ARGUMENT IS A PREFIX, two fixed extensions being appended;
 *   - A VARIANT IN THE DBSNP FILE IS `IN_DB_SNP` AND ONE THAT IS NOT IS `NOVEL`, and the two
 *     always sum to the total;
 *   - THE TI/TV RATIOS ARE COUNTED SEPARATELY FOR THE TWO, so a file whose known variants are
 *     transitions and whose novel ones are transversions reports very different numbers;
 *   - A FILTERED VARIANT IS COUNTED AS FILTERED AND NOWHERE ELSE, so it reaches neither the
 *     known nor the novel tally;
 *   - INDELS ARE COUNTED APART FROM SNPS, with their own dbSNP and insertion-to-deletion columns;
 *   - A MULTIALLELIC SNP HAS ITS OWN COLUMN and is not counted among the plain ones;
 *   - THE DETAIL FILE HAS ONE ROW PER SAMPLE AND THE SUMMARY ONE ROW FOR THE FILE, and the
 *     summary's totals are not the detail rows' sums when the samples do not all carry every
 *     variant;
 *   - A SAMPLE'S HET/HOMVAR RATIO AND ITS TOTAL HET DEPTH COME FROM ITS OWN GENOTYPES, the depth
 *     being read off the AD field of passing biallelic SNP hets;
 *   - --TARGET_INTERVALS RESTRICTS THE WALK, so a variant outside them is counted nowhere;
 *   - A SINGLETON IS COUNTED WHEN EXACTLY ONE SAMPLE CARRIES THE ALTERNATE;
 *   - THE SUMMARY ROW CARRIES NO SAMPLE_ALIAS, being the file's and not a sample's;
 *   - AND A VCF WITH NO VARIANTS AT ALL STILL WRITES BOTH FILES, their counts nought. Not every
 *     ratio is then NaN: PCT_DBSNP is a division of nought by nought and comes out `?`, while the
 *     TI/TV columns come out `0`, so the two kinds of empty read differently.
 *
 * Output:
 *
 *     in\t<name>\t<that vcf, escaped>
 *     detail\t<case>\t<the detail table, its rows sorted, escaped>
 *     summary\t<case>\t<the summary table, escaped>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CollectVariantCallingMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 100000;

    /** One VCF record: position, ref, alt, filter, and one genotype per sample. */
    record Record(int position, String reference, String alternate, String filter,
                  List<String> genotypes) {}

    static Record snp(final int position, final String reference, final String alternate,
                      final String... genotypes) {
        return new Record(position, reference, alternate, "PASS", Arrays.asList(genotypes));
    }

    static String vcf(final List<String> samples, final List<Record> records) {
        final List<String> lines = new ArrayList<>();
        lines.add("##fileformat=VCFv4.2");
        lines.add("##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">");
        lines.add("##FILTER=<ID=LOW,Description=\"low quality\">");
        lines.add("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        lines.add("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">");
        lines.add("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">");
        final List<String> columns = new ArrayList<>(List.of(
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"));
        columns.addAll(samples);
        lines.add(String.join("\t", columns));
        for (final Record record : records) {
            final List<String> row = new ArrayList<>(List.of(
                    "chr1", Integer.toString(record.position()), ".", record.reference(),
                    record.alternate(), "50", record.filter(), ".", "GT:AD:GQ"));
            for (final String genotype : record.genotypes()) {
                row.add(genotype + ":10,10:99");
            }
            lines.add(String.join("\t", row));
        }
        return String.join("\n", lines) + "\n";
    }

    /** A dbSNP VCF naming the positions given. */
    static String dbsnp(final List<Record> records) {
        final List<String> lines = new ArrayList<>();
        lines.add("##fileformat=VCFv4.2");
        lines.add("##contig=<ID=chr1,length=" + CONTIG_LENGTH + ">");
        lines.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        for (final Record record : records) {
            lines.add(String.join("\t", "chr1", Integer.toString(record.position()),
                    "rs" + record.position(), record.reference(), record.alternate(), ".", ".",
                    "."));
        }
        return String.join("\n", lines) + "\n";
    }

    /** A tribble linear index beside the VCF, which the processor requires. */
    static void index(final Path vcf) throws Exception {
        IndexFactory.createLinearIndex(vcf.toFile(), new VCFCodec())
                .writeBasedOnFeatureFile(vcf.toFile());
    }

    /** A table without its comments, its data rows SORTED. */
    static String table(final String text) {
        final List<String> head = new ArrayList<>();
        final List<String> rows = new ArrayList<>();
        boolean seenHeader = false;
        for (final String line : text.split("\n", -1)) {
            if (line.startsWith("#") || line.isEmpty()) {
                continue;
            }
            if (!seenHeader) {
                head.add(line);
                seenHeader = true;
            } else {
                rows.add(line);
            }
        }
        Collections.sort(rows);
        head.addAll(rows);
        return String.join("\n", head);
    }

    static void run(final String name, final List<String> samples, final List<Record> records,
                    final List<Record> known, final String intervals, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("variantcallingmetrics");
        final Path in = dir.resolve("in.vcf");
        Files.writeString(in, vcf(samples, records), StandardCharsets.UTF_8);
        // The processor wants an index whatever the arguments say, so both inputs get one.
        index(in);
        final Path db = dir.resolve("dbsnp.vcf");
        Files.writeString(db, dbsnp(known), StandardCharsets.UTF_8);
        index(db);
        final Path out = dir.resolve("out");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "DBSNP=" + db, "O=" + out));
        if (intervals != null) {
            final Path list = dir.resolve("targets.interval_list");
            Files.writeString(list, intervals, StandardCharsets.UTF_8);
            argv.add("TI=" + list);
        }
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.vcf.CollectVariantCallingMetrics()
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
        emit("detail", name, table(Files.readString(
                Path.of(out + ".variant_calling_detail_metrics"), StandardCharsets.UTF_8)));
        emit("summary", name, table(Files.readString(
                Path.of(out + ".variant_calling_summary_metrics"), StandardCharsets.UTF_8)));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final List<String> one = List.of("s1");
        final List<String> two = List.of("s1", "s2");

        // Six SNPs: the three that are known are two transitions and one transversion, and the
        // three that are not are one transition and two transversions. The two TI/TV columns are
        // then 2 and 0.5, which is what tells them apart.
        final List<Record> mixed = List.of(
                snp(1000, "A", "G", "0/1"),
                snp(2000, "C", "T", "0/1"),
                snp(3000, "A", "C", "0/1"),
                snp(4000, "G", "T", "0/1"),
                snp(5000, "T", "C", "0/1"),
                snp(6000, "A", "T", "0/1"));
        final List<Record> known = List.of(
                snp(1000, "A", "G"), snp(2000, "C", "T"), snp(3000, "A", "C"));
        emit("in", "mixed", vcf(one, mixed));
        emit("in", "dbsnp", dbsnp(known));
        run("known-and-novel", one, mixed, known, null);

        // Nothing in dbSNP at all, so every variant is novel.
        run("nothing-known", one, mixed, List.of(), null);
        // And everything in it.
        run("everything-known", one, mixed, mixed, null);

        // A filtered variant, which reaches neither tally.
        run("filtered", one, List.of(
                snp(1000, "A", "G", "0/1"),
                new Record(2000, "C", "T", "LOW", List.of("0/1"))), known, null);

        // Indels, counted apart from the SNPs.
        run("indels", one, List.of(
                snp(1000, "A", "G", "0/1"),
                snp(2000, "C", "CTT", "0/1"),
                snp(3000, "GTT", "G", "0/1")), known, null);

        // A multiallelic SNP, which has its own column.
        run("multiallelic", one, List.of(
                snp(1000, "A", "G", "0/1"),
                snp(2000, "C", "T,A", "0/1")), known, null);

        // Two samples, one of which is homozygous reference at half the sites.
        run("two-samples", two, List.of(
                snp(1000, "A", "G", "0/1", "0/0"),
                snp(2000, "C", "T", "0/1", "0/1"),
                snp(3000, "A", "C", "1/1", "0/0"),
                snp(4000, "G", "T", "0/1", "0/1")), known, null);
        // The summary's totals are not the detail rows' sums when a sample is hom-ref somewhere.

        // A target interval list covering the first two variants only.
        run("targeted", one, mixed, known,
                "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH
                        + "\nchr1\t1\t2500\t+\ttargets\n");

        // A VCF with no variants at all.
        run("empty", one, List.of(), known, null);

        System.out.print(buf);
    }
}
