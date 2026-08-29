/*
 * CheckFingerprint's two metrics files, taken from the reference.
 *
 * The tool asks whether the sample a file was taken from is the sample a set of genotypes says it
 * is. The answer is a LOD score per haplotype block and a summary over them, and what decides it
 * is which alleles the reads carry at the map's sites, so the fixture is the same one the
 * `extractfingerprint` suite uses: three sites, two of them in one block.
 *
 * Ten behaviours this is built to catch.
 *
 *   - THE LOD IS A LOG10 LIKELIHOOD RATIO between the sample being the same and being a random
 *     other, so agreement is positive and disagreement is negative;
 *   - THE EXPECTED SAMPLE DEFAULTS TO THE OBSERVED ONE, taken from the input's read groups, so a
 *     run that names neither compares a sample against its own genotypes;
 *   - A SAMPLE THE GENOTYPES DO NOT CARRY IS NOT AN ERROR: the run returns
 *     `EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND`, which is one, and writes nothing;
 *   - AND A RUN WITH NOTHING TO CHECK RETURNS TWO AND WRITES ITS FILES ANYWAY:
 *     `EXIT_CODE_WHEN_NO_VALID_CHECKS` is a code for a run that produced metrics whose every
 *     comparison was inconclusive, so a caller that reads the files without reading the code sees
 *     numbers that answer nothing;
 *   - `--OUTPUT` AND THE TWO EXPLICIT FILES ARE MUTUALLY EXCLUSIVE, and the first builds the other
 *     two by appending its own suffixes;
 *   - THE SUFFIXES ARE `.fingerprinting_summary_metrics` AND `.fingerprinting_detail_metrics`;
 *   - THE DETAIL FILE IS ONE ROW PER HAPLOTYPE BLOCK, so the fixture's three sites are two rows;
 *   - A SITE WITH NO READS CONTRIBUTES NOTHING rather than a zero, which is what makes a run over
 *     a file that covers one block shorter than one that covers both;
 *   - THE GENOTYPES MAY BE A VCF WITH SEVERAL SAMPLES, and `--EXPECTED_SAMPLE_ALIAS` names which;
 *   - THE INPUT MAY ITSELF BE A VCF, in which case `--OBSERVED_SAMPLE_ALIAS` names the sample;
 *   - AND THE DICTIONARIES MUST AGREE, a mismatch between the input's and the genotypes' being
 *     refused by name.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     genotypes\t<case>\t<the genotype VCF's variant lines, escaped>
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     summary\t<case>\t<the summary table without its comments, escaped>
 *     detail\t<case>\t<the detail table without its comments, escaped>
 *     error\t<case>\t<the exit code, and the reason where there is one>
 *
 * Usage: CheckFingerprintDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class CheckFingerprintDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 600;

    /** chr1 is `ACGT` repeating, so the base at a position is known by arithmetic. */
    static String fasta() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < CONTIG_LENGTH; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /**
     * The haplotype map: three sites, the first two in one block and the third alone.
     *
     * Position 101 is A and position 105 is A and position 201 is A in the reference, so each
     * site's major allele is its reference base and its minor allele is C.
     */
    static String database() {
        final List<String> lines = new ArrayList<>();
        lines.add("@HD\tVN:1.6\tSO:coordinate");
        lines.add("@SQ\tSN:chr1\tLN:" + CONTIG_LENGTH);
        lines.add("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS");
        lines.add("chr1\t101\trs1\tA\tC\t0.4\trs1\t");
        lines.add("chr1\t105\trs2\tA\tC\t0.4\trs1\t");
        lines.add("chr1\t201\trs3\tA\tC\t0.3\trs3\t");
        return String.join("\n", lines) + "\n";
    }

    record Read(String name, int start, String bases, String qualities) {}

    static SAMFileHeader header(final List<String> samples) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        int i = 0;
        for (final String sample : samples) {
            final SAMReadGroupRecord group = new SAMReadGroupRecord("rg" + (++i));
            group.setSample(sample);
            group.setLibrary("lib" + i);
            group.setPlatformUnit("unit" + i);
            header.addReadGroup(group);
        }
        return header;
    }

    static void writeBam(final Path bam, final List<String> samples, final List<Read> reads) {
        final SAMFileHeader header = header(samples);
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, false, bam.toFile())) {
            int index = 0;
            for (final Read spec : reads) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName(spec.name());
                record.setFlags(0);
                record.setReferenceName("chr1");
                record.setAlignmentStart(spec.start());
                record.setMappingQuality(60);
                record.setCigarString(spec.bases().length() + "M");
                record.setReadString(spec.bases());
                record.setBaseQualityString(spec.qualities());
                // Round-robin the read groups, so a two-sample file has reads of each.
                record.setAttribute("RG", "rg" + ((index++ % samples.size()) + 1));
                writer.addAlignment(record);
            }
        }
    }

    /** The VCF's records, without its header lines. */
    static String records(final String vcf) {
        final List<String> kept = new ArrayList<>();
        for (final String line : vcf.split("\n", -1)) {
            if (line.startsWith("##") || line.isEmpty()) {
                continue;
            }
            kept.add(line);
        }
        return String.join("\n", kept);
    }

    static String sampleColumn(final String vcf) {
        for (final String line : vcf.split("\n", -1)) {
            if (line.startsWith("#CHROM")) {
                final String[] columns = line.split("\t");
                return columns[columns.length - 1];
            }
        }
        return "";
    }

    /** A genotype VCF over the map's three sites, for the samples named. */
    static String genotypes(final List<String> samples, final List<String> calls) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=").append(CONTIG_LENGTH).append(">\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (final String sample : samples) {
            text.append('\t').append(sample);
        }
        text.append('\n');
        final int[] positions = {101, 105, 201};
        for (int site = 0; site < positions.length; site++) {
            text.append("chr1\t").append(positions[site]).append("\trs").append(site + 1)
                    .append("\tA\tC\t100\tPASS\t.\tGT:PL");
            for (int i = 0; i < samples.size(); i++) {
                final String call = calls.get(i);
                final String likelihoods = call.equals("0/0") ? "0,50,100"
                        : call.equals("0/1") ? "50,0,50" : "100,50,0";
                text.append('\t').append(call).append(':').append(likelihoods);
            }
            text.append('\n');
        }
        return text.toString();
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

    /** One run: a BAM of the reads given, a genotype VCF of the calls given, and the two files. */
    static void run(final String name, final List<String> samples, final List<Read> reads,
                    final List<String> genotypeSamples, final List<String> calls,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("checkfingerprint");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path map = dir.resolve("map.txt");
        Files.writeString(map, database(), StandardCharsets.UTF_8);
        final Path in = dir.resolve("in.bam");
        writeBam(in, samples, reads);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());
        final Path known = dir.resolve("genotypes.vcf");
        final String vcf = genotypes(genotypeSamples, calls);
        Files.writeString(known, vcf, StandardCharsets.UTF_8);
        emit("genotypes", name, vcf.lines().filter(line -> !line.startsWith("##"))
                .reduce((a, b) -> a + "\n" + b).orElse(""));

        final Path out = Files.createDirectory(dir.resolve("out"));
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "G=" + known, "H=" + map, "R=" + reference,
                "OUTPUT=" + out.resolve("check")));
        argv.addAll(Arrays.asList(extra));
        final int code;
        try {
            code = new picard.fingerprint.CheckFingerprint()
                    .instanceMain(argv.toArray(new String[0]));
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
        if (code != 0) {
            // The tool RETURNS a code for a sample it could not find rather than throwing, and it
            // writes nothing, so the code is what a caller sees.
            emit("error", name, "exit " + code);
            return;
        }
        for (final String file : written) {
            emit(file.endsWith("summary_metrics") ? "summary" : "detail", name,
                    table(out.resolve(file)));
        }
    }

    /**
     * Reads covering the sites, every base the one given.
     *
     * Each read covers ONE haplotype block and no more: twenty bases at 95 reach rs1 and rs2, and
     * twenty at 195 reach rs3. A first version used long reads spanning both blocks, and the
     * depths that came back were not the ones it had written: a read is evidence for one block,
     * not for every site it happens to overlap.
     */
    static List<Read> pileup(final int depth, final char base, final char quality) {
        final List<Read> reads = new ArrayList<>();
        for (int i = 0; i < depth; i++) {
            for (final int start : new int[]{95, 195}) {
                final StringBuilder bases = new StringBuilder();
                final StringBuilder quals = new StringBuilder();
                for (int j = 0; j < 20; j++) {
                    bases.append(base);
                    quals.append(quality);
                }
                reads.add(new Read("r" + start + "_" + i, start, bases.toString(),
                        quals.toString()));
            }
        }
        return reads;
    }


    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
        final List<String> one = List.of("sample1");

        // Reads that agree with the genotypes, and reads that do not.
        run("agreeing", one, pileup(10, 'A', 'I'), one, List.of("0/0"));
        run("disagreeing", one, pileup(10, 'C', 'I'), one, List.of("0/0"));
        run("heterozygous-genotypes", one, pileup(10, 'A', 'I'), one, List.of("0/1"));
        run("homozygous-alternate-genotypes", one, pileup(10, 'C', 'I'), one, List.of("1/1"));

        // A shallower pileup, which moves the score without changing its sign.
        run("one-read-each", one, pileup(1, 'A', 'I'), one, List.of("0/0"));

        // A file that covers one block only, which is fewer rows and not a zero.
        final List<Read> firstBlockOnly = new ArrayList<>();
        for (final Read read : pileup(10, 'A', 'I')) {
            if (read.start() == 95) {
                firstBlockOnly.add(read);
            }
        }
        run("one-block-covered", one, firstBlockOnly, one, List.of("0/0"));

        // No reads at all.
        run("no-reads", one, List.of(), one, List.of("0/0"));

        // The sample the genotypes do not carry, which is a code and not an exception.
        run("a-sample-the-genotypes-do-not-have", one, pileup(10, 'A', 'I'),
                List.of("other"), List.of("0/0"));
        // And the same file with the expected sample named explicitly.
        run("an-expected-sample-named", one, pileup(10, 'A', 'I'),
                List.of("other"), List.of("0/0"), "EXPECTED_SAMPLE_ALIAS=other");

        // Two samples in the genotypes, one of which is the one to compare against.
        run("two-samples-in-the-genotypes", one, pileup(10, 'A', 'I'),
                List.of("sample1", "sample2"), List.of("0/0", "1/1"));
        run("the-other-sample-named", one, pileup(10, 'A', 'I'),
                List.of("sample1", "sample2"), List.of("0/0", "1/1"),
                "EXPECTED_SAMPLE_ALIAS=sample2");

        System.out.print(buf);
    }
}
