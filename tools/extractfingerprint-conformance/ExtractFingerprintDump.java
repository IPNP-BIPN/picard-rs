/*
 * ExtractFingerprint's VCF, taken from the reference.
 *
 * The tool pileups a BAM at the sites of a haplotype map and writes, for each site, the phred
 * likelihoods of the three genotypes of the CONTAMINATING sample under an assumed contamination.
 * What is measured is which bases reach the model, what the contamination argument does to it,
 * and how the output is named.
 *
 * Thirteen behaviours this is built to catch.
 *
 *   - THE PLS ARE THE CONTAMINATOR'S AND NOT THE SAMPLE'S by default, which is what the tool is
 *     for: a run of reads that are all one allele reports the OTHER allele's genotype as the
 *     likeliest contaminant;
 *   - --EXTRACT_CONTAMINATION FLIPS THE CONTAMINATION ARGUMENT rather than the output: the tool
 *     computes `1 - CONTAMINATION` when it is not set, so the same number means opposite things
 *     under the two settings;
 *   - A CONTAMINATION OF NOUGHT AND ONE ARE THE TWO EXTREMES, and their PLs are each other's;
 *   - A BASE MATCHING NEITHER ALLELE IS COUNTED AND IGNORED, reaching neither the likelihoods nor
 *     the depth the model uses;
 *   - THE BASE QUALITY WEIGHTS THE EVIDENCE, so the same pileup at quality twenty and at quality
 *     forty gives different PLs, and one under the fingerprinter's floor gives none at all;
 *   - --LOCUS_MAX_READS CAPS THE EVIDENCE, and the cap is not simply per site: forty reads under
 *     a cap of ten report a depth of SIXTEEN at the block of two SNPs and of ten at the block of
 *     one, so what the cap bounds is the block and not the record;
 *   - THE SAMPLE NAME GETS `-contaminant` APPENDED by default, and --SAMPLE_ALIAS replaces it
 *     outright rather than adding to it;
 *   - --EXTRACT_CONTAMINATED DROPS THAT SUFFIX, naming the sample as the header does;
 *   - THE VCF CARRIES ONE RECORD PER REPRESENTATIVE SNP, so a block of three SNPs is one record;
 *   - --EXTRACT_NON_REPRESENTATIVES_TOO WRITES EVERY SNP OF EVERY BLOCK instead;
 *   - A FILE NAMING MORE THAN ONE SAMPLE IS REFUSED, by a message counting the fingerprints;
 *   - A SITE WITH NO READS AT ALL STILL GETS A RECORD, its PLs being nought across the board
 *     rather than the priors;
 *   - AND THE BLOCK'S MINOR-ALLELE FREQUENCY IS THE PRIOR, so two blocks of the same pileup at a
 *     middling quality report different PLs: 0,13,30 where the frequency is 0.4 and 0,12,30 where
 *     it is 0.3.
 *
 * Output:
 *
 *     db\t<name>\t<that haplotype database, escaped>
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     out\t<case>\t<the VCF without its header lines, escaped>
 *     sample\t<case>\t<the sample column's name>
 *     error\t<case>\t<exception class>:<message>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ExtractFingerprintDump {

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

    static void run(final String name, final List<String> samples, final List<Read> reads,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("extractfingerprint");
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

        final Path out = dir.resolve("out.vcf");
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "O=" + out, "H=" + map, "R=" + reference));
        argv.addAll(Arrays.asList(extra));
        try {
            final int code = new picard.fingerprint.ExtractFingerprint()
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
        if (!Files.exists(out)) {
            emit("error", name, "no output");
            return;
        }
        final String vcf = Files.readString(out, StandardCharsets.UTF_8);
        emit("sample", name, sampleColumn(vcf));
        emit("out", name, records(vcf).replace(dir.toString(), "<dir>"));
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
        emit("db", "map", database());

        final List<String> one = List.of("sample1");
        // Ten reads of the major allele at every site, at quality forty.
        final List<Read> major = pileup(10, 'A', 'I');
        // And of the minor allele.
        final List<Read> minor = pileup(10, 'C', 'I');

        run("major-half-contaminated", one, major, "C=0.5");
        run("minor-half-contaminated", one, minor, "C=0.5");
        run("contamination-nought", one, major, "C=0.0");
        run("contamination-one", one, major, "C=1.0");
        // The same number under the other setting, which flips it.
        run("extract-contaminated", one, major, "C=0.5", "EXTRACT_CONTAMINATED=true");
        run("extract-contaminated-nought", one, major, "C=0.0", "EXTRACT_CONTAMINATED=true");

        // A base matching neither allele.
        run("neither-allele", one, pileup(10, 'G', 'I'), "C=0.5");
        // The same pileup at a middling base quality, and at one under the floor.
        run("middling-base-quality", one, pileup(10, 'A', '5'), "C=0.5");
        run("base-quality-under-the-floor", one, pileup(10, 'A', '#'), "C=0.5");
        // A deeper pileup, capped and not.
        run("deep-uncapped", one, pileup(40, 'A', 'I'), "C=0.5");
        run("deep-capped", one, pileup(40, 'A', 'I'), "C=0.5", "LOCUS_MAX_READS=10");

        // The sample name.
        run("sample-alias", one, major, "C=0.5", "SAMPLE_ALIAS=named");
        // Every SNP rather than the representatives.
        run("all-snps", one, major, "C=0.5", "EXTRACT_NON_REPRESENTATIVES_TOO=true");

        // No reads at all, so the PLs are the priors.
        run("no-reads", one, List.of(), "C=0.5");
        // Two samples in one file.
        run("two-samples", List.of("sample1", "sample2"), major, "C=0.5");

        System.out.print(buf);
    }
}
