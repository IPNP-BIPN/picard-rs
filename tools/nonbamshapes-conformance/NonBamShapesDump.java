/*
 * Oracle dump for what a screen of the non-BAM tools on large inputs found the port answering
 * differently from Picard 3.5.0, none of which the tools' own suites carried:
 *
 *  - CollectVariantCallingMetrics' dbSNP sets mark a record's every position from start to end,
 *    and only for a SNP or an INDEL (an MNP, a mixed or a symbolic record is in neither); a site is
 *    classified by VariantContext.getType(), so an MNP or a mixed site counts nowhere, a
 *    multiallelic indel is complex, and a filtered MNP is no filtered indel; the detail rows come
 *    out in the HashMap order of the sample names.
 *  - the .idx every plain .vcf output gets from VariantContextWriterBuilder's default
 *    INDEX_ON_THE_FLY, whatever CREATE_INDEX says (GtcToVcf, CombineGenotypingArrayVcfs,
 *    ConvertHaplotypeDatabaseToVcf, LiftoverVcf, ExtractFingerprint, IdentifyContaminant), and its
 *    double properties written by Double.toString (a variance under a thousandth is 3.1E-4).
 *  - ##fileDate and "Generated on:" written as new Date().toString(), not a placeholder.
 *  - BaitDesigner's bait numbers formatted by NumberFormat.getIntegerInstance(): 01,200.
 *  - CheckFingerprint's indexed genotype lookup, which a deletion covering a SNP's position
 *    from the base before does not shadow: the answer is the plain file's.
 *  - CreateSequenceDictionary's default OUTPUT, and its refusal of one that exists.
 *  - CollectIndependentReplicateMetrics' refusal of a paired record without an MC tag.
 *  - ExtractIlluminaBarcodes' compressed barcode files: one gzip member, not BGZF.
 *
 * An index is reported by what it says, not its bytes: the bytes carry the output's path and
 * modification time. A date is reported by its shape.
 *
 * Fixtures come from tools/coverage/MakeFixtures.java, mounted at /work/fixtures.
 *
 *   java -cp picard.jar:. NonBamShapesDump
 */
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class NonBamShapesDump {
    /** The corpus files the cases read, copied into the work directory and carried in the dump. */
    static final String[] CORPUS = {
        "arr.bpm", "arr.egt", "arr_dbsnp.vcf", "arr_dbsnp.vcf.idx", "arr_extended.csv", "arr_flags.egt",
        "arr_manifest.csv", "arr_ref.dict", "arr_ref.fasta", "arr_ref.fasta.fai", "arr_s1.gtc", "arr_s1.vcf",
        "arr_s1.vcf.idx", "arr_s2.vcf", "arr_s2.vcf.idx", "fp_a.bai", "fp_a.bam", "fp_genotypes.vcf",
        "fp_genotypes.vcf.idx", "fp_genotypes_a.vcf", "fp_genotypes_a.vcf.idx", "fp_haplotypes.txt",
        "lift_vcf.chain", "liftover.vcf", "ref.dict", "ref.fasta", "ref.fasta.fai", "replicates.bai",
        "replicates.bam", "replicates.vcf", "ill_barcodes.txt"};
    /** Every input and output lives under this directory, which messages name as `<dir>`. */
    static File W;
    static String F;
    static final String DATE = "[A-Z][a-z]{2} [A-Z][a-z]{2} \\d{2} \\d{2}:\\d{2}:\\d{2} [A-Za-z+\\-0-9:]+ \\d{4}";

    public static void main(final String[] args) throws Exception {
        W = Files.createTempDirectory("shapes").toFile();
        F = W.getAbsolutePath() + "/";
        for (final String name : CORPUS) {
            Files.copy(new File("/work/fixtures/" + name).toPath(), new File(W, name).toPath());
            input(new File(W, name));
        }
        variantCallingMetrics();
        indexes();
        baitNames();
        fingerprintOverlap();
        sequenceDictionary();
        replicatesWithoutMateCigars();
        compressedBarcodes();
    }

    /**
     * ExtractIlluminaBarcodes with COMPRESS_OUTPUTS: each barcode file is one gzip member
     * (`openFileForBufferedWriting`), not BGZF, so its bytes are the reference's to the byte.
     */
    static void compressedBarcodes() throws Exception {
        final java.nio.file.Path run = new File("/work/fixtures/ill_run").toPath();
        try (java.util.stream.Stream<java.nio.file.Path> files = Files.walk(run)) {
            for (final java.nio.file.Path from : (Iterable<java.nio.file.Path>) files.sorted()::iterator) {
                final File to = new File(W, "ill_run/" + run.relativize(from));
                if (Files.isDirectory(from)) {
                    to.mkdirs();
                } else {
                    Files.copy(from, to.toPath());
                    input(to);
                }
            }
        }
        final File out = dir("barcodes");
        tool("barcodes_gzip", new picard.illumina.ExtractIlluminaBarcodes(),
                "BASECALLS_DIR=" + F + "ill_run/Data/Intensities/BaseCalls", "READ_STRUCTURE=2T2B",
                "BARCODE_FILE=" + F + "ill_barcodes.txt", "LANE=1", "OUTPUT_DIR=" + out,
                "METRICS_FILE=" + new File(out, "metrics.txt"), "COMPRESS_OUTPUTS=true");
        final StringBuilder hex = new StringBuilder();
        for (final byte b : Files.readAllBytes(new File(out, "s_1_1101_barcode.txt.gz").toPath())) {
            hex.append(String.format("%02x", b));
        }
        emit("gzip_bytes", "barcodes_gzip", hex.toString());
    }

    /**
     * CreateSequenceDictionary with no OUTPUT (the reference's name with `.dict`), with an OUTPUT
     * that exists (refused before anything is read), and on a reference whose extension is not a
     * FASTA one. Paths are reported relative to the run's directory.
     */
    static void sequenceDictionary() throws Exception {
        final File dir = dir("dict");
        final File fasta = new File(dir, "r.fa");
        Files.copy(new File(F + "ref.fasta").toPath(), fasta.toPath());
        input(fasta);
        tool("dict_default_output", new picard.sam.CreateSequenceDictionary(), "REFERENCE=" + fasta);
        final File dict = new File(dir, "r.dict");
        if (dict.exists()) {
            final List<String> kept = new ArrayList<>();
            for (final String line : Files.readAllLines(dict.toPath())) {
                kept.add(line.replaceAll("\tUR:[^\t]*", ""));
            }
            emit("dict", "dict_default_output", String.join("\n", kept));
        } else {
            emit("dict_absent", "dict_default_output", "");
        }
        tool("dict_existing_output", new picard.sam.CreateSequenceDictionary(), "REFERENCE=" + fasta);
        final File other = new File(dir, "r.seq");
        Files.copy(fasta.toPath(), other.toPath());
        input(other);
        tool("dict_not_fasta_no_output", new picard.sam.CreateSequenceDictionary(), "REFERENCE=" + other);
        tool("dict_not_fasta_output", new picard.sam.CreateSequenceDictionary(), "REFERENCE=" + other,
                "OUTPUT=" + new File(dir, "other.dict"));
        emit("dict_written", "dict_not_fasta_output", Boolean.toString(new File(dir, "other.dict").exists()));
    }

    /**
     * CollectIndependentReplicateMetrics on the corpus's replicates without their MC tags: the
     * duplicate-set sort asks every paired record with a mapped mate for its mate's cigar.
     */
    static void replicatesWithoutMateCigars() throws Exception {
        final File dir = dir("cirm");
        final File bam = new File(dir, "no_mc.bam");
        try (htsjdk.samtools.SamReader reader = htsjdk.samtools.SamReaderFactory.makeDefault().open(new File(F + "replicates.bam"));
             htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory().setCreateIndex(true)
                     .makeBAMWriter(reader.getFileHeader(), true, bam)) {
            for (final htsjdk.samtools.SAMRecord record : reader) {
                record.setAttribute("MC", null);
                writer.addAlignment(record);
            }
        }
        input(bam);
        input(new File(dir, "no_mc.bai"));
        tool("cirm_no_mate_cigar", new picard.analysis.replicates.CollectIndependentReplicateMetrics(),
                "INPUT=" + bam, "VCF=" + F + "replicates.vcf", "OUTPUT=" + new File(dir, "out.txt"),
                "MATRIX_OUTPUT=" + new File(dir, "matrix.txt"));
    }

    /** A callset and a dbSNP over every type a VariantContext can have, three samples. */
    static void variantCallingMetrics() throws Exception {
        final File dir = dir("cvcm");
        final String contigs = "##contig=<ID=chr1,length=10000>\n";
        final String dict = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n";
        final File sd = new File(dir, "ref.dict");
        write(sd, dict);
        final StringBuilder db = new StringBuilder("##fileformat=VCFv4.2\n" + contigs
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        // A SNP; a deletion whose span covers 202; an MNP and a mixed record, which mark nothing;
        // an insertion; a symbolic deletion, which marks nothing either.
        for (final String row : new String[]{
                "100\trsA\tA\tG", "200\trsB\tACGT\tA", "300\trsC\tAC\tGT", "400\trsD\tA\tG,AT",
                "500\trsE\tA\tAT", "800\trsF\tAC\tGT", "1100\trsG\tACG\tT", "1400\trsH\tA\t<DEL>"}) {
            final String[] f = row.split("\t");
            db.append("chr1\t").append(f[0]).append('\t').append(f[1]).append('\t').append(f[2]).append('\t')
                    .append(f[3]).append("\t.\tPASS\t").append(f[3].startsWith("<") ? "END=1500" : ".").append('\n');
        }
        final File dbsnp = new File(dir, "dbsnp.vcf");
        write(dbsnp, db.toString());
        final StringBuilder vcf = new StringBuilder("##fileformat=VCFv4.2\n" + contigs
                + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tfather\tmother\n");
        final String[][] sites = {
            {"100", "A", "G", "PASS"}, {"150", "C", "T,G", "PASS"}, {"202", "T", "TA", "PASS"},
            {"300", "A", "AT", "PASS"}, {"400", "A", "ATT", "PASS"}, {"500", "A", "AT", "PASS"},
            {"600", "AT", "A", "PASS"}, {"700", "A", "ATT,ATTT", "PASS"}, {"800", "AC", "GT", "PASS"},
            {"900", "A", "G,AT", "PASS"}, {"1100", "ACG", "T", "PASS"}, {"1200", "AC", "GT", "LowQual"},
            {"1300", "A", "AT", "LowQual"}, {"1450", "C", "A", "PASS"}, {"1600", "G", "A", "PASS"}};
        final String[] calls = {"0/1:7,5:40", "1/1:0,9:30", "0/0:8,0:50", "0/1:3,4:0", "./.:.:."};
        int n = 0;
        for (final String[] s : sites) {
            vcf.append("chr1\t").append(s[0]).append("\t.\t").append(s[1]).append('\t').append(s[2])
                    .append("\t50\t").append(s[3]).append("\t.\tGT:AD:GQ");
            for (int k = 0; k < 3; k++) {
                String call = calls[(n + k) % calls.length];
                if (s[2].contains(",") && call.startsWith("1/1")) call = "1/2:0,4,5:30";
                if (s[2].contains(",") && call.contains(":7,5:")) call = "0/1:7,5,0:40";
                if (s[2].contains(",") && call.contains(":3,4:")) call = "0/1:3,4,0:0";
                if (s[2].contains(",") && call.contains(":8,0:")) call = "0/0:8,0,0:50";
                vcf.append('\t').append(call);
            }
            vcf.append('\n');
            n++;
        }
        final File input = new File(dir, "calls.vcf");
        write(input, vcf.toString());
        index(input);
        index(dbsnp);
        final String out = new File(dir, "cvcm").getPath();
        tool("cvcm_types", new picard.vcf.CollectVariantCallingMetrics(), "INPUT=" + input, "DBSNP=" + dbsnp,
                "SEQUENCE_DICTIONARY=" + sd, "OUTPUT=" + out);
        metrics("cvcm_types", "detail", new File(out + ".variant_calling_detail_metrics"));
        metrics("cvcm_types", "summary", new File(out + ".variant_calling_summary_metrics"));
    }

    /** The index beside each tool's plain .vcf, and the dates the tools stamp. */
    static void indexes() throws Exception {
        final File dir = dir("idx");

        // Sixty thousand SNPs and one deletion: the feature lengths barely vary, so the variance
        // property is under a thousandth.
        final StringBuilder snps = new StringBuilder("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=70000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (int pos = 60001; pos >= 1; pos--) {
            if (pos == 50) snps.append("chr1\t").append(pos).append("\t.\tAC\tA\t.\tPASS\t.\n");
            else snps.append("chr1\t").append(pos).append("\t.\tA\tG\t.\tPASS\t.\n");
        }
        final File unsorted = new File(dir, "unsorted.vcf");
        write(unsorted, snps.toString());
        final File sorted = new File(dir, "sorted.vcf");
        tool("idx_sortvcf", new picard.vcf.SortVcf(), "INPUT=" + unsorted, "OUTPUT=" + sorted);
        idx("idx_sortvcf", sorted);

        final File gtc = new File(dir, "gtc.vcf");
        tool("idx_gtctovcf", new picard.arrays.GtcToVcf(), "INPUT=" + F + "arr_s1.gtc", "OUTPUT=" + gtc,
                "EXTENDED_ILLUMINA_MANIFEST=" + F + "arr_extended.csv", "CLUSTER_FILE=" + F + "arr.egt",
                "ILLUMINA_BEAD_POOL_MANIFEST_FILE=" + F + "arr.bpm", "R=" + F + "arr_ref.fasta",
                "EXPECTED_GENDER=Female", "SAMPLE_ALIAS=sample1", "PIPELINE_VERSION=1.0", "ANALYSIS_VERSION_NUMBER=1");
        idx("idx_gtctovcf", gtc);
        dates("idx_gtctovcf", gtc, "##fileDate=");

        for (final String create : new String[]{"true", "false"}) {
            final String kase = "idx_combine_create_index_" + create;
            final File combined = new File(dir, "combined_" + create + ".vcf");
            tool(kase, new picard.arrays.CombineGenotypingArrayVcfs(), "INPUT=" + F + "arr_s1.vcf",
                    "INPUT=" + F + "arr_s2.vcf", "OUTPUT=" + combined, "CREATE_INDEX=" + create);
            idx(kase, combined);
        }
        final File combinedDefault = new File(dir, "combined_default.vcf");
        tool("idx_combine_default", new picard.arrays.CombineGenotypingArrayVcfs(), "INPUT=" + F + "arr_s1.vcf",
                "INPUT=" + F + "arr_s2.vcf", "OUTPUT=" + combinedDefault);
        idx("idx_combine_default", combinedDefault);

        final File haplotypes = new File(dir, "haplotypes.vcf");
        tool("idx_haplotypes", new picard.fingerprint.ConvertHaplotypeDatabaseToVcf(), "INPUT=" + F + "fp_haplotypes.txt",
                "OUTPUT=" + haplotypes, "REFERENCE_SEQUENCE=" + F + "ref.fasta");
        idx("idx_haplotypes", haplotypes);

        final File lifted = new File(dir, "lifted.vcf");
        final File rejected = new File(dir, "rejected.vcf");
        tool("idx_liftover", new picard.vcf.LiftoverVcf(), "INPUT=" + F + "liftover.vcf", "OUTPUT=" + lifted,
                "REJECT=" + rejected, "CHAIN=" + F + "lift_vcf.chain", "REFERENCE_SEQUENCE=" + F + "ref.fasta");
        idx("idx_liftover", lifted);
        idx("idx_liftover_reject", rejected);
        final File unsortedLift = new File(dir, "lifted_unsorted.vcf");
        tool("idx_liftover_disable_sort", new picard.vcf.LiftoverVcf(), "INPUT=" + F + "liftover.vcf",
                "OUTPUT=" + unsortedLift, "REJECT=" + new File(dir, "rejected2.vcf"), "CHAIN=" + F + "lift_vcf.chain",
                "REFERENCE_SEQUENCE=" + F + "ref.fasta", "DISABLE_SORT=true");
        idx("idx_liftover_disable_sort", unsortedLift);

        final File extracted = new File(dir, "extracted.vcf");
        tool("idx_extract", new picard.fingerprint.ExtractFingerprint(), "INPUT=" + F + "fp_a.bam", "OUTPUT=" + extracted,
                "HAPLOTYPE_MAP=" + F + "fp_haplotypes.txt", "REFERENCE_SEQUENCE=" + F + "ref.fasta", "SAMPLE_ALIAS=sampleA");
        idx("idx_extract", extracted);
        dates("idx_extract", extracted, "##fileDate=");

        final File contaminant = new File(dir, "contaminant.vcf");
        tool("idx_contaminant", new picard.fingerprint.IdentifyContaminant(), "INPUT=" + F + "fp_a.bam",
                "OUTPUT=" + contaminant, "HAPLOTYPE_MAP=" + F + "fp_haplotypes.txt",
                "REFERENCE_SEQUENCE=" + F + "ref.fasta", "SAMPLE_ALIAS=sampleA");
        idx("idx_contaminant", contaminant);
        dates("idx_contaminant", contaminant, "##fileDate=");

        final File report = new File(dir, "report.txt");
        tool("date_extended_manifest", new picard.arrays.illumina.CreateExtendedIlluminaManifest(), "INPUT=" + F + "arr_manifest.csv",
                "OUTPUT=" + new File(dir, "extended.csv"), "REPORT_FILE=" + report, "BAD_ASSAYS_FILE=" + new File(dir, "bad.txt"),
                "CLUSTER_FILE=" + F + "arr_flags.egt", "DBSNP_FILE=" + F + "arr_dbsnp.vcf", "TARGET_BUILD=37",
                "REFERENCE_SEQUENCE=" + F + "arr_ref.fasta");
        dates("date_extended_manifest", report, "Generated on: ");
    }

    /** One target long enough for 1,225 baits, on a reference of its own. */
    static void baitNames() throws Exception {
        final File dir = dir("baits");
        long seed = 11;
        final StringBuilder fa = new StringBuilder(">chr1\n");
        for (int i = 0; i < 100000; i++) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            fa.append("ACGT".charAt((int) ((seed >>> 33) & 3)));
            if (i % 60 == 59) fa.append('\n');
        }
        fa.append('\n');
        final File fasta = new File(dir, "ref.fasta");
        write(fasta, fa.toString());
        write(new File(dir, "ref.fasta.fai"), "chr1\t100000\t6\t60\t61\n");
        write(new File(dir, "ref.dict"), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n");
        final File targets = new File(dir, "targets.interval_list");
        write(targets, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t99000\t+\tbig\n");
        final File out = new File(dir, "out");
        out.mkdirs();
        tool("baits_thousands", new picard.util.BaitDesigner(), "TARGETS=" + targets, "REFERENCE_SEQUENCE=" + fasta,
                "DESIGN_NAME=design", "OUTPUT_DIRECTORY=" + out);
        final List<String> names = new ArrayList<>();
        for (final String line : Files.readAllLines(new File(out, "design.baits.interval_list").toPath())) {
            if (!line.startsWith("@")) names.add(line);
        }
        emit("baits", "baits_thousands", String.join("\n", names));
        final List<String> headers = new ArrayList<>();
        for (final String line : Files.readAllLines(new File(out, "design.design.fasta").toPath())) {
            if (line.startsWith(">")) headers.add(line);
        }
        emit("fasta_names", "baits_thousands", String.join("\n", headers));
    }

    /** CheckFingerprint against indexed genotypes where a deletion overlaps every other SNP. */
    static void fingerprintOverlap() throws Exception {
        final File dir = dir("fp");
        final StringBuilder out = new StringBuilder();
        int record = 0;
        for (final String line : Files.readAllLines(new File(F + "fp_genotypes.vcf").toPath())) {
            if (line.startsWith("#") || record++ % 2 == 1) {
                out.append(line).append('\n');
                continue;
            }
            final String[] f = line.split("\t", -1);
            final int pos = Integer.parseInt(f[1]);
            final StringBuilder deletion = new StringBuilder();
            deletion.append(f[0]).append('\t').append(pos - 1).append("\tdel").append(pos).append("\tNN\tN\t.\tPASS\t.");
            for (int k = 8; k < f.length; k++) deletion.append('\t').append(f[k]);
            out.append(deletion).append('\n').append(line).append('\n');
        }
        final File genotypes = new File(dir, "overlapping.vcf");
        write(genotypes, out.toString());
        index(genotypes);
        final String prefix = new File(dir, "check").getPath();
        // A VCF as the INPUT: the corpus BAM leaves no haplotype confidently checked either way.
        tool("fp_overlap", new picard.fingerprint.CheckFingerprint(), "INPUT=" + F + "fp_genotypes_a.vcf",
                "GENOTYPES=" + genotypes, "OUTPUT=" + prefix, "HAPLOTYPE_MAP=" + F + "fp_haplotypes.txt",
                "EXPECTED_SAMPLE_ALIAS=sampleA");
        metrics("fp_overlap", "summary", new File(prefix + ".fingerprinting_summary_metrics"));
        metrics("fp_overlap", "detail", new File(prefix + ".fingerprinting_detail_metrics"));
    }

    /** One run; an exception's text with the work directory written as `<dir>`. */
    static void tool(final String kase, final picard.cmdline.CommandLineProgram program, final String... argv) {
        emit("tool", kase, program.getClass().getSimpleName());
        emit("args", kase, String.join(" ", argv).replace(W.getAbsolutePath(), "<dir>"));
        try {
            emit("exit", kase, Integer.toString(program.instanceMain(argv)));
        } catch (final Throwable t) {
            emit("thrown", kase, t.toString().replace(W.getAbsolutePath(), "<dir>"));
        }
    }

    /** A directory of the work directory's, for one group of cases. */
    static File dir(final String name) {
        final File dir = new File(W, name);
        dir.mkdirs();
        return dir;
    }

    /** An input file, carried whole (base64) under its path in the work directory. */
    static void input(final File file) throws Exception {
        final String name = W.toPath().relativize(file.toPath()).toString();
        emit("input", name, java.util.Base64.getEncoder().encodeToString(Files.readAllBytes(file.toPath())));
    }

    /** A metrics file without its `# ` header lines, which name the command line and the date. */
    static void metrics(final String kase, final String which, final File file) throws Exception {
        if (!file.exists()) {
            emit("metrics_absent_" + which, kase, "");
            return;
        }
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readAllLines(file.toPath())) {
            if (!line.startsWith("# ")) kept.add(line);
        }
        emit("metrics_" + which, kase, String.join("\n", kept));
    }

    /** What an index says: its kind, its properties in order, and the sequences it covers. */
    static void idx(final String kase, final File vcf) {
        final File file = new File(vcf.getPath() + ".idx");
        if (!file.exists()) {
            emit("idx_absent", kase, "");
            return;
        }
        final Index index = IndexFactory.loadIndex(file.getPath());
        emit("idx_type", kase, index.getClass().getSimpleName());
        final List<String> properties = new ArrayList<>();
        for (final Map.Entry<String, String> e : index.getProperties().entrySet()) {
            properties.add(e.getKey() + "=" + e.getValue());
        }
        emit("idx_properties", kase, String.join("\n", properties));
        emit("idx_sequences", kase, String.join(",", index.getSequenceNames()));
    }

    /** Whether the line that carries a date carries one in Date.toString()'s shape. */
    static void dates(final String kase, final File file, final String prefix) throws Exception {
        String shape = "no line starting " + prefix;
        for (final String line : Files.readAllLines(file.toPath())) {
            if (line.startsWith(prefix)) {
                shape = line.substring(prefix.length()).matches(DATE) ? "a date" : "not a date";
                break;
            }
        }
        emit("date", kase, shape);
    }

    /** An input VCF's index, which is an input too. */
    static void index(final File vcf) throws Exception {
        final File idx = new File(vcf.getPath() + ".idx");
        IndexFactory.writeIndex(IndexFactory.createLinearIndex(vcf, new htsjdk.variant.vcf.VCFCodec()), idx);
        input(idx);
    }

    /** An input file, written and carried. */
    static void write(final File file, final String text) throws Exception {
        try (PrintStream ps = new PrintStream(file)) {
            ps.print(text);
        }
        input(file);
    }

    static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
