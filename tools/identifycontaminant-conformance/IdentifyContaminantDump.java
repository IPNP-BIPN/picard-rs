/*
 * IdentifyContaminant's VCF, taken from the reference.
 *
 * The tool is ExtractFingerprint with one argument negated and one default changed, and nothing
 * else: it copies its own arguments across, sets `EXTRACT_CONTAMINATION = !EXTRACT_CONTAMINATED`,
 * and calls the other tool's doWork. The measurement puts the same fixture through it so the
 * difference between the two files is exactly that negation.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE DEFAULT IS THE OPPOSITE ONE: `EXTRACT_CONTAMINATED` is false by default, so
 *     `EXTRACT_CONTAMINATION` is TRUE, and a default run reports the CONTAMINANT where
 *     ExtractFingerprint's default run reports the contaminated sample;
 *   - THE CONTAMINATION ARGUMENT IS THEREFORE NOT FLIPPED by default. At 0.5 the flip is
 *     numerically a no-op, so only the sample's NAME tells the two tools apart there; at nought
 *     and at one the PLs are each other's, which is where the negation shows;
 *   - --EXTRACT_CONTAMINATED=true RESTORES THE OTHER TOOL'S DEFAULT, the suffix going with it;
 *   - THE SAMPLE NAME GAINS `-contaminant` BY DEFAULT here, for the same reason;
 *   - --SAMPLE_ALIAS STILL REPLACES THE NAME outright;
 *   - --LOCUS_MAX_READS DEFAULTS TO TWO HUNDRED rather than fifty, so a hundred reads pass
 *     uncapped here and the same hundred under an explicit fifty report seventy-nine at the block
 *     of two SNPs and fifty at the block of one;
 *   - THE MODEL IS THE SAME ONE, so every other case answers as ExtractFingerprint does under the
 *     matching setting;
 *   - AND A FILE NAMING MORE THAN ONE SAMPLE IS REFUSED by the same message.
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

public class IdentifyContaminantDump {

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
        final Path dir = Files.createTempDirectory("identifycontaminant");
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
            final int code = new picard.fingerprint.IdentifyContaminant()
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
        final List<Read> major = pileup(10, 'A', 'I');
        final List<Read> minor = pileup(10, 'C', 'I');

        // The default, which is the opposite of ExtractFingerprint's.
        run("major-default", one, major, "C=0.5");
        run("minor-default", one, minor, "C=0.5");
        // And the setting that restores the other tool's default.
        run("extract-contaminated", one, major, "C=0.5", "EXTRACT_CONTAMINATED=true");

        run("contamination-nought", one, major, "C=0.0");
        run("contamination-one", one, major, "C=1.0");

        // The two tools' LOCUS_MAX_READS defaults are fifty and two hundred, so a pileup of a
        // hundred is capped by one and not by the other.
        run("hundred-reads", one, pileup(100, 'A', 'I'), "C=0.5");
        run("hundred-reads-capped", one, pileup(100, 'A', 'I'), "C=0.5", "LOCUS_MAX_READS=50");

        run("sample-alias", one, major, "C=0.5", "SAMPLE_ALIAS=named");
        run("neither-allele", one, pileup(10, 'G', 'I'), "C=0.5");
        run("no-reads", one, List.of(), "C=0.5");
        run("two-samples", List.of("sample1", "sample2"), major, "C=0.5");

        System.out.print(buf);
    }
}
