/*
 * `CreateExtendedIlluminaManifest`, taken from the reference.
 *
 * The tool takes Illumina's own manifest and works out, for every locus, where it sits on the
 * target build and which alleles it has there. That is the file `GtcToVcf` reads, and the answer
 * is not a copy of the input: a locus whose probe does not match the reference, or whose position
 * is not on the build, is FLAGGED rather than dropped, and the flag is what the downstream tool
 * acts on.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE OUTPUT IS THE INPUT PLUS SEVEN COLUMNS, and the header carries the versions and the
 *     files the run was given;
 *   - A LOCUS WHOSE ALLELES MATCH THE REFERENCE IS `PASS`, with the reference base as
 *     `build37RefAllele`;
 *   - AN rsID COMES FROM THE dbSNP FILE, matched by position, and a locus dbSNP does not carry
 *     keeps the manifest's own name;
 *   - A LOCUS ON A CONTIG THE REFERENCE DOES NOT HAVE IS FLAGGED, not refused;
 *   - THE REPORT COUNTS THE FLAGS, one row per kind;
 *   - AND THE BAD ASSAYS FILE LISTS what was flagged, which is the same information the other way
 *     round.
 *
 * Output:
 *
 *     manifest\t<case>\t<the assay rows of the written manifest, escaped>
 *     report\t<case>\t<the report file, escaped>
 *     bad\t<case>\t<the bad assays file, escaped>
 *     code\t<case>\t<the exit status>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: CreateExtendedManifestDump
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

public class CreateExtendedManifestDump {

    static final StringBuilder buf = new StringBuilder();
    static final int LENGTH = 5000;

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

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

    static String base(final int position) {
        return String.valueOf("ACGT".charAt((position - 1) % 4));
    }

    /** A dbSNP VCF over the positions given, each with an rs name of its own. */
    static String dbsnp(final List<int[]> sites) {
        final StringBuilder text = new StringBuilder("##fileformat=VCFv4.2\n");
        text.append("##contig=<ID=1,length=").append(LENGTH).append(">\n");
        text.append("##contig=<ID=2,length=").append(LENGTH).append(">\n");
        text.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (final int[] site : sites) {
            text.append(site[0]).append('\t').append(site[1]).append("\trs")
                    .append(site[1]).append('\t').append(base(site[1]))
                    .append("\tG\t100\tPASS\t.\n");
        }
        return text.toString();
    }

    /** The assay rows of a manifest, without its header sections. */
    static String assays(final Path file) throws Exception {
        final List<String> kept = new ArrayList<>();
        boolean inAssays = false;
        for (final String line : Files.readString(file, StandardCharsets.UTF_8).split("\n", -1)) {
            if (line.startsWith("[Assay]")) {
                inAssays = true;
                continue;
            }
            if (line.startsWith("[Controls]")) {
                break;
            }
            if (inAssays && !line.isEmpty()) {
                kept.add(line);
            }
        }
        return String.join("\n", kept);
    }

    /** The plain manifest rows a case's loci make. */
    static List<MakeExtendedManifest.Row> rows(final List<MakeBpm.Locus> loci) {
        final List<MakeExtendedManifest.Row> rows = new ArrayList<>();
        for (final MakeBpm.Locus locus : loci) {
            final String[] alleles = locus.snp().replace("[", "").replace("]", "").split("/");
            rows.add(new MakeExtendedManifest.Row(locus.name(), locus.snp(), locus.chrom(),
                    locus.position(), locus.addressA(), locus.addressB(), locus.refStrand(),
                    "", alleles[0], alleles[1], ""));
        }
        return rows;
    }

    static void run(final String name, final List<MakeBpm.Locus> loci, final List<int[]> known,
                    final String... extra) throws Exception {
        final Path dir = Files.createTempDirectory("extendedmanifest");
        final Path reference = dir.resolve("ref.fasta");
        Files.writeString(reference, fasta(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
                "R=" + reference, "O=" + dir.resolve("ref.dict"),
                "GENOME_ASSEMBLY=GRCh37", "SPECIES=Homo sapiens"});
        FastaSequenceIndexCreator.create(reference, true);

        // The tool's INPUT is Illumina's own manifest CSV rather than the binary bead pool file:
        // the extension is what it WRITES, and the columns it adds are the ones it computes.
        final Path manifest = MakeExtendedManifest.writePlain(dir.resolve("fixture.csv"),
                rows(loci));
        final Path bpm = MakeBpm.write(dir.resolve("fixture.bpm"), "fixture.bpm", loci);
        final List<String> names = new ArrayList<>();
        for (final MakeBpm.Locus locus : loci) {
            names.add(locus.name());
        }
        final Path egt = MakeEgt.write(dir.resolve("fixture.egt"), "fixture.bpm", names);
        final Path sites = dir.resolve("dbsnp.vcf");
        Files.writeString(sites, dbsnp(known), StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.writeIndex(
                htsjdk.tribble.index.IndexFactory.createLinearIndex(
                        sites.toFile(), new htsjdk.variant.vcf.VCFCodec()),
                new java.io.File(sites + ".idx"));

        final Path out = dir.resolve("extended.csv");
        final Path report = dir.resolve("report.txt");
        final Path bad = dir.resolve("bad.txt");
        final List<String> argv = new ArrayList<>(List.of(
                "INPUT=" + manifest, "OUTPUT=" + out, "REPORT_FILE=" + report,
                "BAD_ASSAYS_FILE=" + bad, "CLUSTER_FILE=" + egt, "DBSNP_FILE=" + sites,
                "R=" + reference, "TARGET_BUILD=37"));
        argv.addAll(Arrays.asList(extra));

        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final int code = new picard.arrays.illumina.CreateExtendedIlluminaManifest()
                    .instanceMain(argv.toArray(new String[0]));
            System.setOut(original);
            System.setErr(originalError);
            emit("code", name, String.valueOf(code));
            if (code != 0) {
                final List<String> reasons = new ArrayList<>();
                for (final String line : said.toString(StandardCharsets.UTF_8).split("\n", -1)) {
                    if (line.startsWith("ERROR:")) {
                        reasons.add(line);
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
        emit("manifest", name, assays(out));
        if (Files.exists(report)) {
            // The report opens with a clock and names its input by path, neither of which is the
            // measurement: what is kept is the counts under them.
            final List<String> kept = new ArrayList<>();
            for (final String line : Files.readString(report, StandardCharsets.UTF_8)
                    .split("\n", -1)) {
                if (line.isEmpty() || line.startsWith("Generated on")) {
                    continue;
                }
                kept.add(line.replace(dir.toString(), "<dir>"));
            }
            emit("report", name, String.join("\n", kept));
        }
        if (Files.exists(bad)) {
            emit("bad", name, assays(bad));
        }
    }

    public static void main(final String[] args) throws Exception {
        // Four loci whose reference bases the fixture reference actually carries.
        final List<MakeBpm.Locus> four = List.of(
                new MakeBpm.Locus("rs1", "[A/G]", "1", 1001, 11, 0, 0, 1, "TOP", "TOP", "+"),
                new MakeBpm.Locus("rs2", "[T/C]", "1", 2001, 12, 13, 1, 1, "BOT", "BOT", "-"),
                new MakeBpm.Locus("rs3", "[A/C]", "2", 3001, 14, 0, 0, 2, "TOP", "PLUS", "+"),
                new MakeBpm.Locus("rs4", "[A/T]", "2", 4001, 15, 16, 2, 2, "PLUS", "TOP", "+"));
        run("four-loci", four, List.of(new int[]{1, 1001}, new int[]{2, 3001}));
        // Nothing in dbSNP at all.
        run("no-known-sites", four, List.of());
        // A locus on a contig the reference does not have.
        final List<MakeBpm.Locus> offContig = new ArrayList<>(four);
        offContig.set(0, new MakeBpm.Locus("rs1", "[A/G]", "9", 1001, 11, 0, 0, 1,
                "TOP", "TOP", "+"));
        run("a-contig-that-is-not-there", offContig, List.of());
        // A locus past the end of its contig.
        final List<MakeBpm.Locus> offEnd = new ArrayList<>(four);
        offEnd.set(1, new MakeBpm.Locus("rs2", "[T/C]", "1", 99000, 12, 13, 1, 1,
                "BOT", "BOT", "-"));
        run("a-position-past-the-end", offEnd, List.of());

        System.out.print(buf);
    }
}
