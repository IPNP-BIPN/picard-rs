/*
 * An extended Illumina manifest, which is a CSV rather than a binary file.
 *
 * `CreateExtendedIlluminaManifest` writes one and `GtcToVcf` reads it: it is the ordinary Illumina
 * manifest plus the seven columns that say where each locus sits on the target build and which
 * alleles it has there. Its shape is a header of `key,value` rows ending at `Loci Count`, then an
 * `[Assay]` line, then a column header and a row per locus.
 *
 * The parser refuses a column name it does not know, so every name here is one of its own.
 *
 * Usage: MakeExtendedManifest <file>
 */

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class MakeExtendedManifest {

    /** The columns, in the order the file writes them. */
    static final List<String> COLUMNS = List.of(
            "IlmnID", "Name", "IlmnStrand", "SNP", "AddressA_ID", "AlleleA_ProbeSeq",
            "AddressB_ID", "AlleleB_ProbeSeq", "GenomeBuild", "Chr", "MapInfo", "Ploidy",
            "Species", "Source", "SourceVersion", "SourceStrand", "SourceSeq", "TopGenomicSeq",
            "BeadSetID", "Exp_Clusters", "RefStrand", "Intensity_Only",
            "build37Chr", "build37Pos", "build37RefAllele", "build37AlleleA", "build37AlleleB",
            "build37Rsid", "build37Flag");

    /** One locus's row, with the build-37 columns the extension adds. */
    record Row(String name, String snp, String chrom, int position, int addressA, int addressB,
               String refStrand, String refAllele, String alleleA, String alleleB, String flag) {}

    static String row(final Row row) {
        final List<String> values = new ArrayList<>(List.of(
                row.name() + "_ilmn", row.name(), "TOP", row.snp(),
                String.valueOf(row.addressA()), "ACGTACGT",
                row.addressB() == 0 ? "" : String.valueOf(row.addressB()),
                row.addressB() == 0 ? "" : "ACGTACGA",
                "37", row.chrom(), String.valueOf(row.position()), "diploid", "Homo sapiens",
                "source", "1", "TOP", "ACGTACGTACGT", "ACGT", "1", "3", row.refStrand(), "0",
                row.chrom(), String.valueOf(row.position()), row.refAllele(), row.alleleA(),
                row.alleleB(), row.name(), row.flag()));
        return String.join(",", values);
    }

    static String text(final List<Row> rows) {
        final StringBuilder out = new StringBuilder();
        out.append("Illumina, Inc.\n");
        out.append("[Heading]\n");
        out.append("Descriptor File Name,fixture.bpm\n");
        out.append("Assay Format,Infinium HTS\n");
        out.append("Date Manufactured,1/1/2020\n");
        out.append("Loci Count ,").append(rows.size()).append('\n');
        out.append("[Assay]\n");
        out.append(String.join(",", COLUMNS)).append('\n');
        for (final Row row : rows) {
            out.append(row(row)).append('\n');
        }
        out.append("[Controls]\n");
        return out.toString();
    }

    static Path write(final Path file, final List<Row> rows) throws IOException {
        Files.createDirectories(file.getParent());
        Files.writeString(file, text(rows), StandardCharsets.UTF_8);
        return file;
    }

    public static void main(final String[] args) throws Exception {
        write(Paths.get(args[0]), List.of(
                new Row("rs1", "[A/G]", "1", 1000, 11, 0, "+", "A", "A", "G", "PASS")));
        System.out.println("wrote " + args[0]);
    }
}
