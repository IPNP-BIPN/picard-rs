import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * The genotyping-array corpus of tools/coverage/MakeArraysFixtures.java at chip scale, for
 * tools/benchmark/screen_nonbam.sh: the same writers (ArraysBpm, ArraysEgt, ArraysGtc,
 * ArraysManifest) over LOCI loci on two contigs, so a run measures the work and not the startup.
 * Compiled with tools/coverage/*.java inside the oracle container: `MakeBigArrays <dir> <loci>`.
 */
public class MakeBigArrays {

    public static void main(final String[] args) throws Exception {
        final Path dir = Paths.get(args[0]);
        final int loci = Integer.parseInt(args[1]);
        Files.createDirectories(dir);
        // Twenty bases a locus on each of two contigs, every locus at a base of the ACGT pattern.
        final int length = loci * 10 + 2000;
        final Path reference = dir.resolve("arr_ref.fasta");
        final StringBuilder fasta = new StringBuilder();
        for (final String contig : new String[]{"1", "2"}) {
            fasta.append('>').append(contig).append('\n');
            for (int index = 0; index < length; index += 60) {
                for (int b = index; b < Math.min(index + 60, length); b++) {
                    fasta.append("ACGT".charAt(b % 4));
                }
                fasta.append('\n');
            }
        }
        Files.writeString(reference, fasta.toString(), StandardCharsets.UTF_8);
        new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
                "R=" + reference, "O=" + dir.resolve("arr_ref.dict"),
                "GENOME_ASSEMBLY=GRCh37", "SPECIES=Homo sapiens"});
        htsjdk.samtools.reference.FastaSequenceIndexCreator.create(reference, true);

        final String[] snps = {"[A/G]", "[T/C]", "[A/C]", "[A/T]"};
        final List<ArraysManifest.Row> rows = new ArrayList<>();
        final List<ArraysBpm.Locus> bpmLoci = new ArrayList<>();
        final List<String> names = new ArrayList<>();
        final List<ArraysManifest.Row> plain = new ArrayList<>();
        final java.util.Set<Integer> unique = new java.util.TreeSet<>();
        for (int i = 0; i < loci; i++) {
            final String chrom = i < loci / 2 ? "1" : "2";
            final int position = 1001 + 20 * (i < loci / 2 ? i : i - loci / 2);
            final String snp = snps[i % 4];
            final int assayType = i % 3 == 1 ? 1 : 0;
            final int addressA = 11 + 2 * i;
            final int addressB = assayType == 1 ? 12 + 2 * i : 0;
            final String name = "rs" + (i + 1);
            final String[] alleles = snp.replace("[", "").replace("]", "").split("/");
            final String flag = i % 97 == 5 ? "DUPE" : i % 101 == 7 ? "ILLUMINA_FLAGGED" : "PASS";
            rows.add(new ArraysManifest.Row(name, snp, chrom, position, addressA, addressB, "+",
                    MakeArraysFixtures.base(position), alleles[0], alleles[1], flag));
            final int normalization = 1 + i % 40;
            bpmLoci.add(new ArraysBpm.Locus(name, snp, chrom, position, addressA, addressB,
                    assayType, normalization, "TOP", "TOP", "+"));
            unique.add(normalization + 100 * assayType);
            names.add(name);
            plain.add(new ArraysManifest.Row(name, snp, chrom, position, addressA, addressB, "+",
                    "", alleles[0], alleles[1], ""));
        }
        final Path bpm = ArraysBpm.write(dir.resolve("arr.bpm"), "fixture.bpm", bpmLoci);
        final Path egt = ArraysEgt.write(dir.resolve("arr.egt"), "fixture.bpm", names);
        final Path manifest = ArraysManifest.write(dir.resolve("arr_extended.csv"), rows);
        ArraysManifest.writePlain(dir.resolve("arr_manifest.csv"), plain);

        final StringBuilder sites = new StringBuilder("##fileformat=VCFv4.2\n");
        sites.append("##contig=<ID=1,length=").append(length).append(">\n");
        sites.append("##contig=<ID=2,length=").append(length).append(">\n");
        sites.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (int i = 0; i < loci; i += 2) {
            final ArraysManifest.Row row = rows.get(i);
            sites.append(row.chrom()).append('\t').append(row.position()).append("\trs").append(i)
                    .append('\t').append(MakeArraysFixtures.base(row.position())).append("\tG\t100\tPASS\t.\n");
        }
        final Path dbsnp = dir.resolve("arr_dbsnp.vcf");
        Files.writeString(dbsnp, sites.toString(), StandardCharsets.UTF_8);
        htsjdk.tribble.index.IndexFactory.writeIndex(
                htsjdk.tribble.index.IndexFactory.createLinearIndex(
                        dbsnp.toFile(), new htsjdk.variant.vcf.VCFCodec()),
                new java.io.File(dbsnp + ".idx"));

        for (int s = 1; s <= 2; s++) {
            final List<Integer> genotypes = new ArrayList<>();
            final List<Integer> x = new ArrayList<>();
            final List<Integer> y = new ArrayList<>();
            final List<Float> scores = new ArrayList<>();
            for (int i = 0; i < loci; i++) {
                genotypes.add((i * s + i / 7) % 4);
                x.add(500 + (i * 37 * s) % 5000);
                y.add(500 + (i * 53 + s) % 5000);
                scores.add(0.3f + ((i * s) % 7) / 10f);
            }
            final String sample = "sample" + s;
            final Path gtc = ArraysGtc.write(dir.resolve("arr_s" + s + ".gtc"),
                    new ArraysGtc.Sample(sample, genotypes, x, y, scores, 0.97f),
                    "fixture.egt", "fixture.bpm", unique.size());
            MakeArraysFixtures.gtcToVcf(gtc, dir.resolve("arr_s" + s + ".vcf"), reference, manifest,
                    egt, bpm, sample, s == 1 ? "Female" : "Male");
        }
    }
}
