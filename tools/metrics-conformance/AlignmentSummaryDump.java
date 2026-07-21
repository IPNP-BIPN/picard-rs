/*
 * Builds BAMs and reference FASTAs, runs Picard's CollectAlignmentSummaryMetrics, and emits all
 * three so the port can be compared against the real output rather than against an expectation.
 *
 * Output, one line per case:
 *   bam     <TAB> <case> <TAB> <hex of the input BAM>
 *   fasta   <TAB> <case> <TAB> <the reference bases, one contig>
 *   metrics <TAB> <case> <TAB> <metrics file, \n and \t escaped>
 *
 * The case list is chosen around the places the metric definitions and the code disagree:
 * multi-block alignments (where BAD_CYCLES bins by block offset), supplementary records (counted
 * for bases and not for reads), clipping at both ends, and the paired/unpaired row rules.
 */

import htsjdk.samtools.*;
import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AlignmentSummaryDump {

    static final int REF_LENGTH = 200;

    static String hex(byte[] b) {
        StringBuilder sb = new StringBuilder();
        for (byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    static String esc(String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static SAMFileHeader header() {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", REF_LENGTH));
        h.setSequenceDictionary(d);
        h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("sample1");
        rg.setLibrary("lib1");
        rg.setPlatform("ILLUMINA");
        h.addReadGroup(rg);
        return h;
    }

    static SAMRecord read(SAMFileHeader h, String name, int start, String cigar, String bases,
                          int flags, int mapq) {
        SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setFlags(flags);
        r.setReferenceIndex(0);
        r.setAlignmentStart(start);
        r.setMappingQuality(mapq);
        r.setCigarString(cigar);
        r.setReadString(bases);
        byte[] q = new byte[bases.length()];
        Arrays.fill(q, (byte) 30);
        r.setBaseQualities(q);
        r.setAttribute("RG", "rg1");
        if ((flags & 0x4) != 0) {
            r.setReadUnmappedFlag(true);
            r.setReferenceIndex(0);
            r.setAlignmentStart(start);
        }
        return r;
    }

    static String repeat(char c, int n) {
        char[] a = new char[n];
        Arrays.fill(a, c);
        return new String(a);
    }

    /** A reference of A's with C's at the given 0-based positions. */
    static String reference(int... mismatchAt) {
        char[] bases = repeat('A', REF_LENGTH).toCharArray();
        for (int p : mismatchAt) bases[p] = 'C';
        return new String(bases);
    }

    static void emit(String name, List<SAMRecord> records, String refBases) throws Exception {
        SAMFileHeader h = header();

        File bam = File.createTempFile("asm", ".bam");
        try (SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(false)
                .makeBAMWriter(h, true, bam)) {
            for (SAMRecord r : records) w.addAlignment(r);
        }

        File fasta = File.createTempFile("asm", ".fasta");
        try (PrintWriter p = new PrintWriter(fasta)) {
            p.println(">chr1");
            for (int i = 0; i < refBases.length(); i += 60) {
                p.println(refBases.substring(i, Math.min(i + 60, refBases.length())));
            }
        }
        File fai = new File(fasta.getPath() + ".fai");
        try (PrintWriter p = new PrintWriter(fai)) {
            p.printf("chr1\t%d\t6\t60\t61%n", refBases.length());
        }
        File dict = new File(fasta.getPath().replaceAll("\\.fasta$", "") + ".dict");
        try (PrintWriter p = new PrintWriter(dict)) {
            p.println("@HD\tVN:1.6\tSO:unsorted");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", refBases.length());
        }

        File out = File.createTempFile("asm", ".txt");
        int rc = new picard.analysis.CollectAlignmentSummaryMetrics().instanceMain(new String[] {
                "INPUT=" + bam.getPath(),
                "OUTPUT=" + out.getPath(),
                "REFERENCE_SEQUENCE=" + fasta.getPath(),
        });
        if (rc != 0) {
            System.err.println("case " + name + " exited " + rc);
            return;
        }

        System.out.println("bam\t" + name + "\t" + hex(Files.readAllBytes(bam.toPath())));
        System.out.println("fasta\t" + name + "\t" + refBases);
        System.out.println("metrics\t" + name + "\t"
                + esc(new String(Files.readAllBytes(out.toPath()))));

        bam.delete(); fasta.delete(); fai.delete(); dict.delete(); out.delete();
    }

    public static void main(String[] args) throws Exception {
        SAMFileHeader h = header();

        // One perfectly matching unpaired read: the floor.
        emit("perfect", List.of(read(h, "r1", 1, "20M", repeat('A', 20), 0, 60)),
             reference());

        // The BAD_CYCLES probe, as a BAM this time. Two mismatches in different alignment
        // blocks that share a block offset, so they collapse into one bad cycle.
        emit("collide", List.of(read(h, "r1", 1, "10M5D10M", repeat('A', 20), 0, 60)),
             reference(3, 18));

        // The same read with the second mismatch at a different block offset.
        emit("distinct", List.of(read(h, "r1", 1, "10M5D10M", repeat('A', 20), 0, 60)),
             reference(3, 20));

        // Three alignment blocks, so the collision can happen twice over.
        emit("three_blocks",
             List.of(read(h, "r1", 1, "5M3D5M3D5M", repeat('A', 15), 0, 60)),
             reference(1, 9, 17));

        // An insertion: the read advances and the reference does not, which shifts the
        // read-to-block offset the other way.
        emit("insertion", List.of(read(h, "r1", 1, "8M4I8M", repeat('A', 20), 0, 60)),
             reference(2, 12));

        // Adjacent match operators, which getAlignmentBlocks does not merge.
        emit("split_match", List.of(read(h, "r1", 1, "10M10M", repeat('A', 20), 0, 60)),
             reference(3, 13));

        // Clipping, at each end and on each strand.
        emit("soft_clip_3prime",
             List.of(read(h, "r1", 1, "15M5S", repeat('A', 20), 0, 60)), reference());
        emit("soft_clip_5prime",
             List.of(read(h, "r1", 1, "5S15M", repeat('A', 20), 0, 60)), reference());
        emit("soft_clip_reverse",
             List.of(read(h, "r1", 1, "5S15M", repeat('A', 20), 0x10, 60)), reference());
        emit("hard_clip", List.of(read(h, "r1", 1, "15M5H", repeat('A', 15), 0, 60)),
             reference());

        // No-calls in the read, which reach the bad-cycle histogram by the other branch.
        emit("no_calls", List.of(read(h, "r1", 1, "20M", "AAAANAAAAAAAAANAAAAA".substring(0, 20),
             0, 60)), reference());

        // An unmapped read takes the whole-read no-call path instead of the block path.
        emit("unmapped", List.of(read(h, "r1", 1, "*", "AAAANAAAAAAAAANAAAAA", 0x4, 0)),
             reference());

        // A low mapping quality read is aligned but not high quality.
        emit("low_mapq", List.of(read(h, "r1", 1, "20M", repeat('A', 20), 0, 5)),
             reference(3));

        // Vendor-failed: counted in TOTAL_READS and not in PF_READS. Paired with a passing read
        // so that finish() does not throw.
        emit("vendor_failed", List.of(
                read(h, "r1", 1, "20M", repeat('A', 20), 0, 60),
                read(h, "r2", 1, "20M", repeat('A', 20), 0x200, 60)),
             reference());

        // Supplementary: excluded from read counts, included in base counts.
        emit("supplementary", List.of(
                read(h, "r1", 1, "20M", repeat('A', 20), 0, 60),
                read(h, "r2", 30, "20M", repeat('A', 20), 0x800, 60)),
             reference());

        // Secondary: excluded from both.
        emit("secondary", List.of(
                read(h, "r1", 1, "20M", repeat('A', 20), 0, 60),
                read(h, "r2", 30, "20M", repeat('A', 20), 0x100, 60)),
             reference());

        // A proper FR pair, which produces the FIRST_OF_PAIR / SECOND_OF_PAIR / PAIR rows and
        // no UNPAIRED row.
        SAMRecord p1 = read(h, "p", 1, "20M", repeat('A', 20), 0x1 | 0x2 | 0x40 | 0x20, 60);
        p1.setMateReferenceIndex(0);
        p1.setMateAlignmentStart(100);
        p1.setInferredInsertSize(119);
        SAMRecord p2 = read(h, "p", 100, "20M", repeat('A', 20), 0x1 | 0x2 | 0x80 | 0x10, 60);
        p2.setMateReferenceIndex(0);
        p2.setMateAlignmentStart(1);
        p2.setInferredInsertSize(-119);
        emit("proper_pair", List.of(p1, p2), reference(3));

        // An improper pair in RF orientation: chimeric under the default expectation of FR.
        SAMRecord c1 = read(h, "c", 1, "20M", repeat('A', 20), 0x1 | 0x40 | 0x10, 60);
        c1.setMateReferenceIndex(0);
        c1.setMateAlignmentStart(100);
        c1.setInferredInsertSize(119);
        SAMRecord c2 = read(h, "c", 100, "20M", repeat('A', 20), 0x1 | 0x80 | 0x20, 60);
        c2.setMateReferenceIndex(0);
        c2.setMateAlignmentStart(1);
        c2.setInferredInsertSize(-119);
        emit("chimeric_orientation", List.of(c1, c2), reference());

        // A pair carrying MQ, which switches the chimera denominator to the two-ended test.
        SAMRecord m1 = read(h, "m", 1, "20M", repeat('A', 20), 0x1 | 0x2 | 0x40 | 0x20, 60);
        m1.setMateReferenceIndex(0);
        m1.setMateAlignmentStart(100);
        m1.setInferredInsertSize(119);
        m1.setAttribute("MQ", 5);
        SAMRecord m2 = read(h, "m", 100, "20M", repeat('A', 20), 0x1 | 0x2 | 0x80 | 0x10, 60);
        m2.setMateReferenceIndex(0);
        m2.setMateAlignmentStart(1);
        m2.setInferredInsertSize(-119);
        m2.setAttribute("MQ", 60);
        emit("mate_mapping_quality", List.of(m1, m2), reference());

        // An adapter read: unmapped, so the mapping-quality exemption does not apply.
        emit("adapter", List.of(read(h, "r1", 1, "*",
                "AATGATACGGCGACCACCGAGATCT", 0x4, 0)), reference());

        // The XN noise tag.
        SAMRecord noise = read(h, "r1", 1, "20M", repeat('A', 20), 0, 60);
        noise.setAttribute("XN", 1);
        emit("noise_read", List.of(noise), reference());

        // The two divergences fulcrumgenomics/riker documents in its ERRATA against Picard.
        // Riker is an independent Rust reimplementation of these same tools whose stated goal is
        // functional rather than byte equivalence, so its errata is a list of the places a
        // careful reimplementer chooses to differ. Each is turned into a corpus case here, so
        // that what this port reproduces is pinned rather than argued about.

        // 1. "Picard computes mean_aligned_read_length over all PF reads, including unmapped
        //    reads which contribute zero to the sum." One mapped 20-base read and one unmapped
        //    20-base read: Picard's mean is 10, riker's would be 20.
        emit("riker_mean_aligned_dilution", List.of(
                read(h, "mapped", 1, "20M", repeat('A', 20), 0, 60),
                read(h, "unmapped", 1, "*", repeat('A', 20), 0x4, 0)),
             reference());

        // 2. "Picard counts all mapped, paired, non-proper reads as improperly paired, including
        //    reads whose mate is unmapped." Riker requires both mates mapped.
        SAMRecord halfMapped = read(h, "half", 1, "20M", repeat('A', 20), 0x1 | 0x40 | 0x8, 60);
        halfMapped.setMateReferenceIndex(0);
        halfMapped.setMateAlignmentStart(1);
        SAMRecord unmappedMate = read(h, "half", 1, "*", repeat('A', 20), 0x1 | 0x80 | 0x4, 0);
        unmappedMate.setMateReferenceIndex(0);
        unmappedMate.setMateAlignmentStart(1);
        emit("riker_improper_pair_unmapped_mate", List.of(halfMapped, unmappedMate),
             reference());

        // Reads of different lengths, so the read-length histogram has something to summarise.
        List<SAMRecord> mixed = new ArrayList<>();
        int[] lengths = {10, 15, 20, 20, 25, 30, 30, 30, 40, 50};
        for (int i = 0; i < lengths.length; i++) {
            mixed.add(read(h, "r" + i, 1 + i, lengths[i] + "M", repeat('A', lengths[i]), 0, 60));
        }
        emit("mixed_lengths", mixed, reference(5, 25, 45));

        // No case runs off the end of the reference. The `refIndex + i < refLength` guard in
        // collectQualityData looks like a bounds check on real data, but htsjdk's writer
        // rejects such a record at STRICT stringency with CIGAR_MAPS_OFF_REFERENCE, so it
        // cannot be reached through a BAM whose header agrees with the FASTA. It is reproduced
        // in the port anyway, unexercised, because a BAM read at a laxer stringency could
        // reach it.
    }
}
