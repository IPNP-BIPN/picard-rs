/*
 * `GtcToVcf`, taken from the reference.
 *
 * The tool turns a chip's genotype calls into a VCF. It reads four files to do it: the calls
 * themselves, the bead pool manifest that says what each locus IS, the cluster file that says how
 * the call was made, and the extended manifest that says where the locus sits on the target build
 * and which alleles it has there.
 *
 * Seven behaviours this is built to catch.
 *
 *   - A CALL BECOMES A GENOTYPE against the build-37 alleles, not against the chip's A and B, and
 *     where the reference base is NEITHER of the chip's two alleles the record is TRIALLELIC: the
 *     reference base becomes the REF, both chip alleles become ALTs, and the filter says so;
 *   - A NO-CALL IS WRITTEN, as `./.` rather than as a missing record;
 *   - THE RECORD CARRIES THE CHIP'S OWN NUMBERS in its FORMAT fields: the normalized intensities,
 *     the score, and the cluster's mean and deviation;
 *   - THE HEADER CARRIES THE RUN'S IDENTITY, the sample alias, the pipeline version and the
 *     analysis version the command line gave;
 *   - THE OUTPUT IS SORTED by the target build's coordinates rather than by the chip's order;
 *   - A LOCUS THE EXTENDED MANIFEST FLAGS IS DROPPED, not written with a filter: `ILLUMINA_FLAGGED`
 *     and `DUPE` each remove their record from the output entirely;
 *   - AND A SAMPLE'S GENDER IS READ OR GIVEN, and disagreeing with the call is not an error.
 *
 * Output:
 *
 *     vcf\t<case>\t<the variant lines, escaped>
 *     header\t<case>\t<the header lines that name the sample and the versions, escaped>
 *     code\t<case>\t<the exit status>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: GtcToVcfDump
 */

import htsjdk.samtools.reference.FastaSequenceIndexCreator;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GtcToVcfDump {

    static final StringBuilder buf = new StringBuilder();
    static final int LENGTH = 5000;

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** A two-contig reference of `ACGT` repeating, so a locus's reference base is arithmetic. */
    static String fasta() {
        final StringBuilder out = new StringBuilder();
        for (final String contig : new String[]{"1", "2"}) {
            out.append('>').append(contig).append('\n');
            final StringBuilder bases = new StringBuilder();
            for (int index = 0; index < LENGTH; index++) {
                bases.append("ACGT".charAt(index % 4));
            }
            for (int index = 0; index < bases.length(); index += 60) {
                out.append(bases, index, Math.min(index + 60, bases.length())).append('\n');
            }
        }
        return out.toString();
    }

    /** The reference base at a position, which the manifest has to agree with. */
    static String base(final int position) {
        return String.valueOf("ACGT".charAt((position - 1) % 4));
    }

    static String lines(final Path file, final boolean variants) throws Exception {
        final List<String> kept = new ArrayList<>();
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.isEmpty()) {
                continue;
            }
            final boolean isHeader = line.startsWith("#");
            if (variants == !isHeader) {
                // The header lines worth comparing are the ones the command line decided; the rest
                // are htsjdk's own and are the same for every run.
                if (!variants && !(line.startsWith("##fileformat") || line.contains("sampleAlias")
                        || line.contains("pipelineVersion") || line.contains("analysisVersion")
                        || line.startsWith("#CHROM"))) {
                    continue;
                }
                kept.add(line);
            }
        }
        return String.join("\n", kept);
    }

    static void run(final String name, final List<MakeExtendedManifest.Row> rows,
                    final MakeGtc.Sample sample, final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("gtctovcf");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        // The tool refuses any reference whose dictionary does not declare GRCh37: it is written
        // for build 37 and says so rather than guessing.
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
                "R=" + reference, "O=" + dir.resolve("ref.dict"),
                "GENOME_ASSEMBLY=GRCh37", "SPECIES=Homo sapiens"});
        FastaSequenceIndexCreator.create(reference, true);

        final List<MakeBpm.Locus> loci = new ArrayList<>();
        final List<String> names = new ArrayList<>();
        for (final MakeExtendedManifest.Row row : rows) {
            loci.add(new MakeBpm.Locus(row.name(), row.snp(), row.chrom(), row.position(),
                    row.addressA(), row.addressB(), row.addressB() == 0 ? 0 : 1, 1,
                    "TOP", "TOP", row.refStrand()));
            names.add(row.name());
        }
        final Path bpm = MakeBpm.write(dir.resolve("fixture.bpm"), "fixture.bpm", loci);
        final Path egt = MakeEgt.write(dir.resolve("fixture.egt"), "fixture.bpm", names);
        final Path manifest = MakeExtendedManifest.write(dir.resolve("extended.csv"), rows);
        final java.util.Set<Integer> unique = new java.util.TreeSet<>();
        for (final MakeBpm.Locus locus : loci) {
            unique.add(locus.normalizationId() + 100 * locus.assayType());
        }
        final Path gtc = MakeGtc.write(dir.resolve("sample.gtc"), sample, "fixture.egt",
                "fixture.bpm", unique.size());
        final Path out = dir.resolve("out.vcf");

        final List<String> argv = new ArrayList<>(List.of(
                "INPUT=" + gtc, "OUTPUT=" + out, "R=" + reference,
                "EXTENDED_ILLUMINA_MANIFEST=" + manifest, "CLUSTER_FILE=" + egt,
                "ILLUMINA_BEAD_POOL_MANIFEST_FILE=" + bpm,
                "SAMPLE_ALIAS=sample1", "ANALYSIS_VERSION_NUMBER=1", "PIPELINE_VERSION=1.0"));
        final List<String> tail = new ArrayList<>(Arrays.asList(extra));
        if (tail.stream().noneMatch(argument -> argument.startsWith("EXPECTED_GENDER="))) {
            tail.add("EXPECTED_GENDER=Female");
        }
        argv.addAll(tail);

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.arrays.GtcToVcf().instanceMain(argv.toArray(new String[0]));
            System.setOut(original);
            System.setErr(originalError);
            emit("code", name, String.valueOf(code));
            if (code != 0) {
                // A refused command line prints its reason rather than throwing it.
                final List<String> reasons = new ArrayList<>();
                for (final String line : said.toString(StandardCharsets.UTF_8).split("\n", -1)) {
                    if (line.startsWith("ERROR:") || line.contains("Exception")) {
                        reasons.add(line.replace(dir.toString(), "<dir>"));
                    }
                }
                emit("error", name, String.join(" | ", reasons));
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
        emit("header", name, lines(out, false));
        emit("vcf", name, lines(out, true));
    }

    public static void main(final String[] args) throws Exception {
        // Four loci on two contigs, called AA, AB, BB and no-call.
        final List<MakeExtendedManifest.Row> four = List.of(
                new MakeExtendedManifest.Row("rs1", "[A/G]", "1", 1001, 11, 0, "+",
                        base(1001), "A", "G", "PASS"),
                new MakeExtendedManifest.Row("rs2", "[T/C]", "1", 2001, 12, 13, "+",
                        base(2001), "T", "C", "PASS"),
                new MakeExtendedManifest.Row("rs3", "[A/C]", "2", 3001, 14, 0, "+",
                        base(3001), "A", "C", "PASS"),
                new MakeExtendedManifest.Row("rs4", "[A/T]", "2", 4001, 15, 16, "+",
                        base(4001), "A", "T", "PASS"));
        run("four-loci", four, MakeGtc.fixture("sample1"));

        // Every call the same, so the genotypes are the manifest's answer rather than the chip's.
        run("all-homozygous-reference", four, new MakeGtc.Sample("sample1",
                List.of(1, 1, 1, 1), List.of(1000, 2000, 3000, 4000),
                List.of(1100, 2100, 3100, 4100), List.of(0.9f, 0.9f, 0.9f, 0.9f), 1.0f));
        run("all-no-calls", four, new MakeGtc.Sample("sample1",
                List.of(0, 0, 0, 0), List.of(1000, 2000, 3000, 4000),
                List.of(1100, 2100, 3100, 4100), List.of(0f, 0f, 0f, 0f), 0.0f));

        // A locus the extended manifest flags as failed.
        final List<MakeExtendedManifest.Row> flagged = new ArrayList<>(four);
        // The flag is an enum of the reference's own, so a locus is flagged with one of its
        // names: `ILLUMINA_FLAGGED` is the one that says the chip's maker called the assay bad.
        flagged.set(1, new MakeExtendedManifest.Row("rs2", "[T/C]", "1", 2001, 12, 13, "+",
                base(2001), "T", "C", "ILLUMINA_FLAGGED"));
        run("a-flagged-locus", flagged, MakeGtc.fixture("sample1"));
        final List<MakeExtendedManifest.Row> duped = new ArrayList<>(four);
        duped.set(2, new MakeExtendedManifest.Row("rs3", "[A/C]", "2", 3001, 14, 0, "+",
                base(3001), "A", "C", "DUPE"));
        run("a-duplicate-locus", duped, MakeGtc.fixture("sample1"));

        // The gender the command line declares.
        run("a-male-sample", four, MakeGtc.fixture("sample1"), "EXPECTED_GENDER=Male");

        System.out.print(buf);
    }
}
