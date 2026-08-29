/*
 * CollectArraysVariantCallingMetrics' three files, taken from the reference.
 *
 * The tool counts an Illumina array's calls, and almost everything it reports comes out of the
 * VCF's HEADER rather than out of its variants: the chip, the sample, the two genders, the call
 * rate the genotyping software reported, the dates, the cluster file and the twenty-three control
 * codes. What the variants decide is a handful of counters and one threshold.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE HEADER IS REQUIRED AND THE MISSING PIECE IS NAMED: a VCF without `autocallVersion` is
 *     refused by an IllegalArgumentException naming the line, so the header is a contract;
 *   - SOME OF IT IS OPTIONAL AND SOME IS NOT, and the difference is not guessable: `zcallVersion`
 *     and `pipelineVersion` may be absent where `clusterFile` and `p95Red` may not;
 *   - AN ABSENT GENDER IS `NotReported` AND NOT AN ERROR, and it renders as its symbol;
 *   - THE TWENTY-THREE CONTROL CODES ARE A FILE OF THEIR OWN, parsed out of header lines whose
 *     value is four fields separated by a pipe;
 *   - A FILTERED ASSAY IS NOT COUNTED, unless its filter is `DUPE`, which is counted as if it had
 *     passed;
 *   - `ZEROED_OUT_ASSAY` IS COUNTED APART, being a filter that means the chip could not read the
 *     assay rather than that the call was poor;
 *   - A CALL IS A GENOTYPE THAT IS CALLED, and an autocall is one whose `GTA` is not
 *     `VCFConstants.EMPTY_GENOTYPE`, which is `./.` and not a single dot. The two spellings do
 *     opposite things: a `GTA` of `./.` is NOT an autocall, and a `GTA` of `.` is one, because a
 *     single dot is a missing attribute and the default the reference asks for is the genotype's
 *     own string. A variant with no GTA at all is an autocall for the same reason;
 *   - `NUM_IN_DB_SNP` COMES FROM THE dbSNP FILE, so the same VCF against two dbSNPs gives two
 *     answers, and `NOVEL_SNPS` is the difference;
 *   - THE CALL RATE THRESHOLD DECIDES A PASS/FAIL COLUMN and nothing else;
 *   - THE THREE FILES ARE NAMED FROM ONE PREFIX with a dot and their own extensions;
 *   - AND THE NUMBER OF PROCESSORS DOES NOT CHANGE A NUMBER: the accumulator's results are merged,
 *     and one processor and two give the same file, which is what makes the multithreading safe
 *     to leave in a golden.
 *
 * Output:
 *
 *     vcf\t<case>\t<the variant lines, escaped>
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     detail\t<case>\t<the detail table without its comments, escaped>
 *     summary\t<case>\t<the summary table without its comments, escaped>
 *     controls\t<case>\t<the control-code table without its comments, escaped>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CollectArraysVariantCallingMetricsDump
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class CollectArraysVariantCallingMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 1000;

    /** The header lines the accumulator reads, in the reference's own spelling. */
    static Map<String, String> headerFields() {
        final Map<String, String> fields = new LinkedHashMap<>();
        fields.put("arrayType", "TestArray-24v1-0_A1");
        fields.put("extendedIlluminaManifestVersion", "1.4");
        fields.put("chipWellBarcode", "1234567890_R01C01");
        fields.put("analysisVersionNumber", "1");
        fields.put("sampleAlias", "sample1");
        fields.put("expectedGender", "Female");
        fields.put("fingerprintGender", "Female");
        fields.put("gtcCallRate", "0.995");
        fields.put("autocallGender", "Female");
        fields.put("autocallDate", "09/21/2016 20:40");
        fields.put("imagingDate", "8/15/2015 7:28:52 AM");
        fields.put("clusterFile", "TestArray-24v1-0_A1.egt");
        fields.put("autocallVersion", "3.0.0");
        fields.put("zcallVersion", "1.0.0");
        fields.put("zcallThresholds", "thresholds.7.txt");
        fields.put("p95Red", "1000");
        fields.put("p95Green", "2000");
        fields.put("scannerName", "N370");
        fields.put("pipelineVersion", "IlluminaGenotypingArray_v1.0");
        return fields;
    }

    /** The twenty-three control codes, whose value is four fields separated by a pipe. */
    static final String[][] CONTROLS = {
            {"DNP(High)", "Staining", "10", "20"},
            {"DNP(Bgnd)", "Staining", "11", "21"},
            {"Biotin(High)", "Staining", "12", "22"},
            {"Biotin(Bgnd)", "Staining", "13", "23"},
            {"Extension(A)", "Extension", "14", "24"},
            {"Extension(T)", "Extension", "15", "25"},
            {"Extension(C)", "Extension", "16", "26"},
            {"Extension(G)", "Extension", "17", "27"},
            {"TargetRemoval", "TargetRemoval", "18", "28"},
            {"Hyb(High)", "Hybridization", "19", "29"},
            {"Hyb(Medium)", "Hybridization", "20", "30"},
            {"Hyb(Low)", "Hybridization", "21", "31"},
            {"String(PM)", "Stringency", "22", "32"},
            {"String(MM)", "Stringency", "23", "33"},
            {"NSB(Bgnd)Red", "Non-SpecificBinding", "24", "34"},
            {"NSB(Bgnd)Purple", "Non-SpecificBinding", "25", "35"},
            {"NSB(Bgnd)Blue", "Non-SpecificBinding", "26", "36"},
            {"NSB(Bgnd)Green", "Non-SpecificBinding", "27", "37"},
            {"NP(A)", "Non-Polymorphic", "28", "38"},
            {"NP(T)", "Non-Polymorphic", "29", "39"},
            {"NP(C)", "Non-Polymorphic", "30", "40"},
            {"NP(G)", "Non-Polymorphic", "31", "41"},
            {"Restore", "Restoration", "32", "42"}};

    /** One variant line of the fixture. */
    record Variant(int position, String reference, String alternate, String genotype,
                   String filter, String gta) {}

    static Variant call(final int position, final String genotype) {
        return new Variant(position, "A", "C", genotype, "PASS", null);
    }

    static String vcf(final List<Variant> variants, final List<String> omitted,
                      final Map<String, String> overrides) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=").append(CONTIG_LENGTH).append(">\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
        text.append("##FORMAT=<ID=GTA,Number=1,Type=String,Description=\"Autocall genotype\">\n");
        text.append("##FILTER=<ID=DUPE,Description=\"Duplicate assay\">\n");
        text.append("##FILTER=<ID=ZEROED_OUT_ASSAY,Description=\"Assay zeroed out\">\n");
        text.append("##FILTER=<ID=TRIALLELIC,Description=\"Triallelic\">\n");
        final Map<String, String> fields = headerFields();
        fields.putAll(overrides);
        for (final Map.Entry<String, String> field : fields.entrySet()) {
            if (omitted.contains(field.getKey())) {
                continue;
            }
            text.append("##").append(field.getKey()).append('=').append(field.getValue())
                    .append('\n');
        }
        for (final String[] control : CONTROLS) {
            text.append("##").append(control[0]).append('=')
                    .append(String.join("|", control)).append('\n');
        }
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n");
        for (final Variant variant : variants) {
            text.append("chr1\t").append(variant.position()).append("\trs")
                    .append(variant.position()).append('\t').append(variant.reference())
                    .append('\t').append(variant.alternate()).append("\t100\t")
                    .append(variant.filter()).append("\t.\tGT:GQ");
            if (variant.gta() != null) {
                text.append(":GTA");
            }
            text.append('\t').append(variant.genotype()).append(":99");
            if (variant.gta() != null) {
                text.append(':').append(variant.gta());
            }
            text.append('\n');
        }
        return text.toString();
    }

    /** A VCF, block-compressed and indexed, which the tool opens with an index requirement. */
    static Path writeVcf(final Path dir, final String name, final String text) throws Exception {
        final Path plain = dir.resolve(name + ".plain.vcf");
        Files.writeString(plain, text, StandardCharsets.UTF_8);
        final Path out = dir.resolve(name + ".vcf.gz");
        try (final htsjdk.variant.vcf.VCFFileReader reader =
                     new htsjdk.variant.vcf.VCFFileReader(plain, false);
             final htsjdk.variant.variantcontext.writer.VariantContextWriter writer =
                     new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder()
                             .setOutputPath(out)
                             .setOption(htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY)
                             .setReferenceDictionary(reader.getFileHeader().getSequenceDictionary())
                             .build()) {
            writer.writeHeader(reader.getFileHeader());
            for (final htsjdk.variant.variantcontext.VariantContext variant : reader) {
                writer.add(variant);
            }
        }
        return out;
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

    static void run(final String name, final List<Variant> variants, final List<String> omitted,
                    final Map<String, String> overrides, final List<Variant> dbsnp,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("arraysmetrics");
        final String text = vcf(variants, omitted, overrides);
        final Path input = writeVcf(dir, "in", text);
        emit("vcf", name, text.lines().filter(line -> !line.startsWith("##"))
                .reduce((a, b) -> a + "\n" + b).orElse(""));
        final Path known = dir.resolve("dbsnp.vcf");
        Files.writeString(known, vcf(dbsnp, List.of(), Map.of()), StandardCharsets.UTF_8);

        final Path out = Files.createDirectory(dir.resolve("out"));
        final List<String> argv = new ArrayList<>(List.of(
                "INPUT=" + input, "OUTPUT=" + out.resolve("m"), "DBSNP=" + known));
        argv.addAll(Arrays.asList(extra));
        // One processor unless the case named another: the count is a scalar argument, so naming
        // it twice is a refusal rather than an override.
        if (extra.length == 0 || Arrays.stream(extra).noneMatch(a -> a.startsWith("NUM_PROCESSORS"))) {
            argv.add("NUM_PROCESSORS=1");
        }
        try {
            final int code = new picard.arrays.CollectArraysVariantCallingMetrics()
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
            final String kind = file.contains("control") ? "controls"
                    : file.contains("detail") ? "detail" : "summary";
            emit(kind, name, table(out.resolve(file)));
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final List<Variant> plain = List.of(
                call(100, "0/0"), call(200, "0/1"), call(300, "1/1"), call(400, "./."));
        final List<Variant> known = List.of(call(100, "0/0"), call(200, "0/1"));

        run("plain", plain, List.of(), Map.of(), known);
        // The same VCF against a dbSNP that knows nothing, which is what NOVEL_SNPS is.
        run("nothing-in-dbsnp", plain, List.of(), Map.of(), List.of());
        // A filtered assay, a duplicate one and a zeroed-out one, which are three different fates.
        run("a-filtered-assay", List.of(
                call(100, "0/0"),
                new Variant(200, "A", "C", "0/1", "TRIALLELIC", null)), List.of(), Map.of(), known);
        run("a-duplicate-assay", List.of(
                call(100, "0/0"),
                new Variant(200, "A", "C", "0/1", "DUPE", null)), List.of(), Map.of(), known);
        run("a-zeroed-out-assay", List.of(
                call(100, "0/0"),
                new Variant(200, "A", "C", "0/1", "ZEROED_OUT_ASSAY", null)),
                List.of(), Map.of(), known);
        // The GTA attribute, which is what separates a call from an autocall.
        run("an-autocall", List.of(
                new Variant(100, "A", "C", "0/1", "PASS", "AC")), List.of(), Map.of(), known);
        // The autocall test compares the GTA against `VCFConstants.EMPTY_GENOTYPE`, which is a
        // single dot, so a GTA of "./." or "NC" is an autocall and only "." is not.
        run("a-call-whose-autocall-is-empty", List.of(
                new Variant(100, "A", "C", "0/1", "PASS", ".")), List.of(), Map.of(), known);
        run("a-call-whose-autocall-is-a-no-call-genotype", List.of(
                new Variant(100, "A", "C", "0/1", "PASS", "./.")), List.of(), Map.of(), known);
        // An indel, which is counted apart from a SNP.
        run("an-indel", List.of(
                new Variant(100, "A", "AC", "0/1", "PASS", null)), List.of(), Map.of(), known);

        // The header: what may be left out and what may not.
        run("without-zcall", plain, List.of("zcallVersion", "zcallThresholds"), Map.of(), known);
        run("without-the-pipeline-version", plain, List.of("pipelineVersion"), Map.of(), known);
        run("without-the-genders", plain, List.of("expectedGender", "fingerprintGender"),
                Map.of(), known);
        run("without-the-autocall-version", plain, List.of("autocallVersion"), Map.of(), known);
        run("without-the-cluster-file", plain, List.of("clusterFile"), Map.of(), known);
        run("without-p95-red", plain, List.of("p95Red"), Map.of(), known);
        run("without-the-call-rate", plain, List.of("gtcCallRate"), Map.of(), known);

        // The threshold, on each side of the call rate the header reports.
        run("under-the-call-rate-threshold", plain, List.of(), Map.of(), known,
                "CALL_RATE_PF_THRESHOLD=0.999");
        run("over-the-call-rate-threshold", plain, List.of(), Map.of(), known,
                "CALL_RATE_PF_THRESHOLD=0.5");

        // And two processors, which must not change a number.
        run("two-processors", plain, List.of(), Map.of(), known, "NUM_PROCESSORS=2");

        System.out.print(buf);
    }
}
