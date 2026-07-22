/*
 * Oracle dump harness for CollectQualityYieldMetricsSNVQ conformance in picard-rs.
 *
 * Emits the corpus as escaped TSV to stdout: a `sam` row (the exact input) and a `metrics` row
 * (the .txt file the oracle produced) per case. The committed corpus is
 * `java ... SnvqDump | gzip > tests/data/snvq.txt.gz`, regenerated and compared in CI.
 *
 * The tool reads per-base SNV qualities from the tags qa/qc/qg/qt (q + lowercased ACGT), FASTQ
 * encoded, and counts base and SNVQ yield at the 20/30/40 thresholds. It is entirely integer
 * counting plus ratio derivation, so no transcendental math is involved. The reads here are
 * unmapped (the tool sets usesNoRefReads), and exercise: a plain PF read, a vendor-fail read, an
 * N base (which is unequal to all four alt bases so all four SNVQs count), secondary and
 * supplementary reads (excluded by default), and a mix of quality thresholds.
 *
 * Run under the pinned en_US locale (htsjdk-rs decision 0011):
 *   java -Duser.language=en -Duser.country=US -cp picard-fat.jar:. SnvqDump
 */
import htsjdk.samtools.*;
import picard.analysis.CollectQualityYieldMetricsSNVQ;

import java.io.File;
import java.nio.file.Files;

public class SnvqDump {
    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("s1");
        header.addReadGroup(rg);

        final File sam = File.createTempFile("snvq-", ".sam");
        sam.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(header, false, sam);

        // A PF read, all four bases present, qualities and SNV qualities spanning the thresholds.
        w.addAlignment(read(header, "r1", "ACGTACGT",
                new int[] {19, 20, 29, 30, 39, 40, 41, 45},
                "IIIIIIII", "5555IIII", "##!!5555", "IIII####", 0));
        // A vendor-fail read: counted in TOTAL_* but not PF_*.
        w.addAlignment(read(header, "r2", "ACGT",
                new int[] {40, 40, 40, 40},
                "IIII", "IIII", "IIII", "IIII", SAMFlag.READ_FAILS_VENDOR_QUALITY_CHECK.intValue()));
        // An N base is unequal to all four alt bases, so all four SNVQs are counted at that position.
        w.addAlignment(read(header, "r3", "ANGT",
                new int[] {40, 40, 40, 40},
                "IIII", "IIII", "IIII", "IIII", 0));
        // Secondary and supplementary reads are excluded by default.
        w.addAlignment(read(header, "r4", "ACGT",
                new int[] {40, 40, 40, 40},
                "IIII", "IIII", "IIII", "IIII", SAMFlag.SECONDARY_ALIGNMENT.intValue()));
        w.addAlignment(read(header, "r5", "ACGT",
                new int[] {40, 40, 40, 40},
                "IIII", "IIII", "IIII", "IIII", SAMFlag.SUPPLEMENTARY_ALIGNMENT.intValue()));
        w.close();

        final File out = File.createTempFile("snvq-", ".txt");
        out.deleteOnExit();
        final int rc = new CollectQualityYieldMetricsSNVQ().instanceMain(new String[] {
                "INPUT=" + sam.getAbsolutePath(),
                "OUTPUT=" + out.getAbsolutePath(),
        });
        if (rc != 0) {
            System.err.println("tool exited " + rc);
            System.exit(rc);
        }

        emit("sam", "basic", sam);
        emit("metrics", "basic", out);
    }

    private static SAMRecord read(final SAMFileHeader header, final String name, final String bases,
                                  final int[] quals, final String qa, final String qc,
                                  final String qg, final String qt, final int extraFlags) {
        final SAMRecord r = new SAMRecord(header);
        r.setReadName(name);
        // Mapped with a trivial cigar: the SAM validator forbids secondary/supplementary flags on
        // unmapped reads, and the collector counts bases regardless of alignment.
        r.setReferenceName("chr1");
        r.setAlignmentStart(100);
        r.setCigarString(bases.length() + "M");
        r.setReadBases(bases.getBytes());
        final byte[] q = new byte[quals.length];
        for (int i = 0; i < quals.length; i++) q[i] = (byte) quals[i];
        r.setBaseQualities(q);
        r.setAttribute("qa", qa);
        r.setAttribute("qc", qc);
        r.setAttribute("qg", qg);
        r.setAttribute("qt", qt);
        r.setAttribute("RG", "rg1");
        r.setFlags(r.getFlags() | extraFlags);
        return r;
    }

    private static void emit(final String kind, final String kase, final File file) throws Exception {
        final String payload = new String(Files.readAllBytes(file.toPath()))
                .replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + payload);
    }
}
