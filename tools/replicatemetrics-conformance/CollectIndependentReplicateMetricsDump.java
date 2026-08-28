/*
 * CollectIndependentReplicateMetrics' counters, taken from the reference.
 *
 * The tool asks how many of a library's duplicates are not duplicates at all: it looks at duplicate
 * sets of size two and three that cover a HETEROZYGOUS site of the sample, and splits them by
 * whether the reads in a set agree on the allele. A set that disagrees cannot have come from one
 * molecule, so it counts the independent replication the whole estimate is built on.
 *
 * Twelve behaviours this is built to catch.
 *
 *   - THE SITES COME FROM THE VCF AND THE READS FROM THE BAM, and only heterozygous sites of the
 *     named sample are looked at, so a homozygous site contributes nothing at all;
 *   - THE SAMPLE MAY BE OMITTED ONLY WHEN THE VCF HAS EXACTLY ONE, and a two-sample VCF without
 *     `--SAMPLE` is refused by an IllegalArgumentException naming the count;
 *   - `--MINIMUM_GQ` IS A FILTER ON THE SITE and not on a read: a site under it is not a site;
 *   - `--MINIMUM_MQ` IS A FILTER ON THE READ, applied before the sets are built, so a read under
 *     it is not in the set that would otherwise have held it;
 *   - `--MINIMUM_BQ` IS A FILTER ON THE BASE, so a read whose base at the site is poor is in the
 *     set and contributes no allele;
 *   - UNPAIRED READS ARE FILTERED BY DEFAULT, `--FILTER_UNPAIRED_READS` being true, so a fixture
 *     of single reads counts nothing until it is turned off;
 *   - A DUPLICATE SET OF TWO AND ONE OF THREE ARE COUNTED SEPARATELY, and a set of four is
 *     neither, its reads counted only in `nReadsInBigSets`;
 *   - A SET WHOSE READS DISAGREE ON THE ALLELE is what the estimate is built on, and a set whose
 *     reads agree is filed under the allele they agree on;
 *   - A THIRD ALLELE AT A SITE TAKES THE WHOLE SITE OUT, counted in `nThreeAllelesSites`;
 *   - THE UMI IS READ FROM `--BARCODE_TAG` AND ITS QUALITY FROM `--BARCODE_BQ`, and a barcode
 *     with any base under `--MINIMUM_BARCODE_BQ` is not used;
 *   - `--MATRIX_OUTPUT` IS OPTIONAL, and the confusion matrix it writes is a second file rather
 *     than a section of the first;
 *   - `--STOP_AFTER` COUNTS DUPLICATE SETS, not reads and not sites;
 *   - AND A RUN THAT EXAMINED NOTHING REPORTS ONE THREE-ALLELE SITE. The counter is incremented
 *     in the tail of `doWork` whenever the loop ended with no locus pending a merge, which is the
 *     state a run over an empty site list is in from the start: a homozygous site, a site under
 *     the quality floor and a file of unpaired reads all report `nThreeAllelesSites=1` having
 *     looked at nothing at all.
 *
 * Output:
 *
 *     sam\t<case>\t<the input as sam, without its header, escaped>
 *     vcf\t<case>\t<the variant lines, escaped>
 *     metrics\t<case>\t<the metrics table without its comments, escaped>
 *     matrix\t<case>\t<the confusion matrix file without its comments, escaped>
 *     files\t<case>\t<the basenames written, sorted, space separated>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CollectIndependentReplicateMetricsDump
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

public class CollectIndependentReplicateMetricsDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int CONTIG_LENGTH = 300;
    /** The site every read is built around, one-based, whose reference base is an `A`. */
    static final int SITE = 101;

    static String bases() {
        final StringBuilder bases = new StringBuilder();
        for (int i = 0; i < CONTIG_LENGTH; i++) {
            bases.append("ACGT".charAt(i % 4));
        }
        return bases.toString();
    }

    static String fasta() {
        final String bases = bases();
        final StringBuilder fasta = new StringBuilder(">chr1\n");
        for (int i = 0; i < bases.length(); i += 60) {
            fasta.append(bases, i, Math.min(i + 60, bases.length())).append('\n');
        }
        return fasta.toString();
    }

    /** One read pair of a duplicate set: where it starts, the base it reads at the site, its UMI. */
    record Pair(String name, int start, char allele, String umi, int mappingQuality,
                char baseQuality, String umiQuality) {}

    static Pair pair(final String name, final char allele, final String umi) {
        return new Pair(name, SITE - 10, allele, umi, 60, 'I', null);
    }

    static SAMFileHeader header() {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", CONTIG_LENGTH));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMReadGroupRecord group = new SAMReadGroupRecord("rg1");
        group.setSample("sample1");
        group.setLibrary("lib1");
        header.addReadGroup(group);
        return header;
    }

    /** The bases a read of twenty covers from `start`, with the site's base replaced. */
    static String readBases(final int start, final char allele) {
        final StringBuilder bases = new StringBuilder(bases().substring(start - 1, start - 1 + 20));
        bases.setCharAt(SITE - start, allele);
        return bases.toString();
    }

    static void writeBam(final Path bam, final List<Pair> pairs, final boolean paired) {
        final SAMFileHeader header = header();
        final List<SAMRecord> records = new ArrayList<>();
        for (final Pair spec : pairs) {
            final SAMRecord first = new SAMRecord(header);
            first.setReadName(spec.name());
            first.setReferenceName("chr1");
            first.setAlignmentStart(spec.start());
            first.setMappingQuality(spec.mappingQuality());
            first.setCigarString("20M");
            first.setReadString(readBases(spec.start(), spec.allele()));
            final StringBuilder quals = new StringBuilder();
            for (int i = 0; i < 20; i++) {
                quals.append(i == SITE - spec.start() ? spec.baseQuality() : 'I');
            }
            first.setBaseQualityString(quals.toString());
            first.setAttribute("RG", "rg1");
            if (spec.umi() != null) {
                first.setAttribute("RX", spec.umi());
                first.setAttribute("QX", spec.umiQuality() == null
                        ? "I".repeat(spec.umi().length()) : spec.umiQuality());
            }
            if (paired) {
                first.setFlags(0x1 | 0x2 | 0x40 | 0x20);
                first.setMateReferenceName("chr1");
                first.setMateAlignmentStart(spec.start() + 60);
                first.setInferredInsertSize(80);
                // The duplicate-set iterator asks each read for its mate's cigar, and refuses a
                // file that does not carry one, so the fixture writes the tag the way a marked
                // file has it.
                first.setAttribute("MC", "20M");
                final SAMRecord second = new SAMRecord(header);
                second.setReadName(spec.name());
                second.setReferenceName("chr1");
                second.setAlignmentStart(spec.start() + 60);
                second.setMappingQuality(spec.mappingQuality());
                second.setCigarString("20M");
                second.setReadString(bases().substring(spec.start() + 59, spec.start() + 79));
                second.setBaseQualityString("I".repeat(20));
                second.setFlags(0x1 | 0x2 | 0x80 | 0x10);
                second.setMateReferenceName("chr1");
                second.setMateAlignmentStart(spec.start());
                second.setInferredInsertSize(-80);
                second.setAttribute("MC", "20M");
                second.setAttribute("RG", "rg1");
                if (spec.umi() != null) {
                    second.setAttribute("RX", spec.umi());
                    second.setAttribute("QX", spec.umiQuality() == null
                            ? "I".repeat(spec.umi().length()) : spec.umiQuality());
                }
                records.add(second);
            } else {
                first.setFlags(0);
            }
            records.add(first);
        }
        records.sort((a, b) -> Integer.compare(a.getAlignmentStart(), b.getAlignmentStart()));
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(header, true, bam.toFile())) {
            for (final SAMRecord record : records) {
                writer.addAlignment(record);
            }
        }
    }

    /** A VCF of one site, with the genotype and quality the case asks for. */
    static String vcf(final String genotype, final int gq, final String alternate,
                      final List<String> samples) {
        final StringBuilder text = new StringBuilder();
        text.append("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=chr1,length=").append(CONTIG_LENGTH).append(">\n");
        text.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        text.append("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (final String sample : samples) {
            text.append('\t').append(sample);
        }
        text.append('\n');
        text.append("chr1\t").append(SITE).append("\t.\t")
                .append(bases().charAt(SITE - 1)).append('\t').append(alternate)
                .append("\t100\tPASS\t.\tGT:GQ");
        for (int i = 0; i < samples.size(); i++) {
            text.append('\t').append(genotype).append(':').append(gq);
        }
        text.append('\n');
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

    static void run(final String name, final List<Pair> pairs, final String genotype, final int gq,
                    final boolean paired, final List<String> samples, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("replicatemetrics");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary()
                .instanceMain(new String[]{"R=" + reference, "O=" + dir.resolve("ref.dict")});
        FastaSequenceIndexCreator.create(reference, true);
        final Path in = dir.resolve("in.bam");
        writeBam(in, pairs, paired);
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(in.toFile())) {
            for (final SAMRecord record : reader) {
                sam.append(record.getSAMString());
            }
        }
        emit("sam", name, sam.toString());
        final Path sites = dir.resolve("sites.vcf");
        final String text = vcf(genotype, gq, "C", samples);
        Files.writeString(sites, text, StandardCharsets.UTF_8);
        emit("vcf", name, text.lines().filter(line -> !line.startsWith("##"))
                .reduce((a, b) -> a + "\n" + b).orElse(""));

        final Path out = Files.createDirectory(dir.resolve("out"));
        final List<String> argv = new ArrayList<>(List.of(
                "I=" + in, "V=" + sites, "O=" + out.resolve("m.independent_replicate_metrics"),
                "R=" + reference));
        argv.addAll(Arrays.asList(extra));
        for (int i = 0; i < argv.size(); i++) {
            argv.set(i, argv.get(i).replace("<out>", out.toString()));
        }
        final java.io.ByteArrayOutputStream errBytes = new java.io.ByteArrayOutputStream();
        final java.io.PrintStream realErr = System.err;
        try {
            final int code;
            try {
                System.setErr(new java.io.PrintStream(errBytes, true, StandardCharsets.UTF_8));
                code = new picard.analysis.replicates.CollectIndependentReplicateMetrics()
                        .instanceMain(argv.toArray(new String[0]));
            } finally {
                System.err.flush();
                System.setErr(realErr);
            }
            if (code != 0) {
                // A validation failure is a return code and not an exception, and its reason is
                // the last line under a usage the golden has no room for.
                emit("error", name, "exit " + code + " " + reason(errBytes.toString(StandardCharsets.UTF_8)));
                return;
            }
        } catch (final Exception e) {
            System.setErr(realErr);
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
            emit(file.endsWith("matrix") ? "matrix" : "metrics", name, table(out.resolve(file)));
        }
    }

    /** The last line of a refusal, which the reference prints under its whole usage. */
    static String reason(final String stderr) {
        final List<String> kept = new ArrayList<>();
        for (final String line : stderr.split("\n", -1)) {
            final String trimmed = line.trim();
            if (!trimmed.isEmpty()) {
                kept.add(trimmed);
            }
        }
        // The reason is the FIRST line: the reference prints it and then its whole usage, of
        // which the last line is the argument list's rather than anything about the failure.
        return kept.isEmpty() ? "" : kept.get(0);
    }

    static List<String> one() {
        return List.of("sample1");
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        final char reference = bases().charAt(SITE - 1);
        final char alternate = 'C';

        // A duplicate set of two whose reads disagree on the allele, and one whose reads agree.
        run("doubleton-disagreeing", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, one());
        run("doubleton-agreeing-on-the-reference", List.of(
                pair("a", reference, "AAAA"), pair("b", reference, "CCCC")),
                "0/1", 99, true, one());
        run("doubleton-agreeing-on-the-alternate", List.of(
                pair("a", alternate, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, one());
        // The same UMI on both reads of a disagreeing set, which is the interesting confusion.
        run("doubleton-disagreeing-same-umi", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "AAAA")),
                "0/1", 99, true, one());

        // A set of three, and a set of four which is neither a double nor a triple.
        run("tripleton", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC"),
                pair("c", alternate, "GGGG")), "0/1", 99, true, one());
        run("a-set-of-four", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC"),
                pair("c", alternate, "GGGG"), pair("d", reference, "TTTT")),
                "0/1", 99, true, one());

        // A third allele at the site takes the whole site out.
        run("a-third-allele", List.of(
                pair("a", reference, "AAAA"), pair("b", 'G', "CCCC")),
                "0/1", 99, true, one());

        // The site's own filters: a homozygous genotype and a low quality one.
        run("a-homozygous-site", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "1/1", 99, true, one());
        run("a-low-quality-site", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 10, true, one());
        run("a-low-quality-site-allowed", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 10, true, one(), "GQ=5");

        // The read's own filters, and the base's.
        run("a-low-mapping-quality-read", List.of(
                new Pair("a", SITE - 10, reference, "AAAA", 10, 'I', null),
                pair("b", alternate, "CCCC")), "0/1", 99, true, one());
        run("a-low-base-quality", List.of(
                new Pair("a", SITE - 10, reference, "AAAA", 60, '$', null),
                pair("b", alternate, "CCCC")), "0/1", 99, true, one());
        run("a-low-base-quality-allowed", List.of(
                new Pair("a", SITE - 10, reference, "AAAA", 60, '$', null),
                pair("b", alternate, "CCCC")), "0/1", 99, true, one(), "BQ=2");

        // Unpaired reads, which the default filters out.
        run("unpaired-reads", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, false, one());
        run("unpaired-reads-kept", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, false, one(), "FUR=false");

        // The barcode: its tag, its quality tag, and the floor under which it is not used.
        run("a-low-quality-barcode", List.of(
                new Pair("a", SITE - 10, reference, "AAAA", 60, 'I', "!!!!"),
                pair("b", alternate, "CCCC")), "0/1", 99, true, one());
        run("another-barcode-tag", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, one(), "BARCODE_TAG=BX");

        // The second file, and the sample rules.
        run("a-matrix", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, one(), "MO=<out>/m.matrix");
        run("two-samples-without-a-name", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, List.of("sample1", "sample2"));
        run("two-samples-with-a-name", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, List.of("sample1", "sample2"), "ALIAS=sample2");
        // The short name is ALIAS and not S, which is worth a case of its own: `S=` is not an
        // option this tool has.
        run("a-sample-by-the-wrong-short-name", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, one(), "S=sample1");
        run("a-sample-that-is-not-there", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, one(), "ALIAS=nobody");

        // And the stop, which counts duplicate sets.
        run("stop-after-one-set", List.of(
                pair("a", reference, "AAAA"), pair("b", alternate, "CCCC")),
                "0/1", 99, true, one(), "STOP_AFTER=1");

        System.out.print(buf);
    }
}
