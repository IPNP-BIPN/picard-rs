/*
 * Generates a benchmark BAM and its reference FASTA.
 *
 * Usage: MakeBenchmarkBam <read count> <output prefix>
 *
 * The reads are synthetic but not uniform: a fixed-seed mixture of read lengths, mapping
 * qualities, strands, pairing states, indels, soft clips and mismatches, so the collector takes
 * every branch rather than the cheapest one. A benchmark over a million identical perfect reads
 * measures branch prediction, not the tool.
 */

import htsjdk.samtools.*;
import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Random;

public class MakeBenchmarkBam {

    static final int REF_LENGTH = 1_000_000;

    public static void main(String[] args) throws Exception {
        final int n = Integer.parseInt(args[0]);
        final String prefix = args[1];
        final Random rng = new Random(20260721L);

        // A reference with a mismatch every 97 bases, so reads hit them at varying offsets.
        char[] refBases = new char[REF_LENGTH];
        Arrays.fill(refBases, 'A');
        for (int i = 0; i < REF_LENGTH; i += 97) refBases[i] = 'C';
        String reference = new String(refBases);

        try (PrintWriter p = new PrintWriter(prefix + ".fasta")) {
            p.println(">chr1");
            for (int i = 0; i < REF_LENGTH; i += 60) {
                p.println(reference.substring(i, Math.min(i + 60, REF_LENGTH)));
            }
        }
        try (PrintWriter p = new PrintWriter(prefix + ".fasta.fai")) {
            p.printf("chr1\t%d\t6\t60\t61%n", REF_LENGTH);
        }
        try (PrintWriter p = new PrintWriter(prefix + ".dict")) {
            p.println("@HD\tVN:1.6\tSO:unsorted");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", REF_LENGTH);
        }

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

        // Coordinate-sorted output, so the writer is not the thing being measured on read-back.
        int[] starts = new int[n];
        for (int i = 0; i < n; i++) starts[i] = 1 + rng.nextInt(REF_LENGTH - 1000);
        Arrays.sort(starts);

        String[] cigars = {"100M", "50M5D50M", "40M10I50M", "10S90M", "90M10S", "5H95M",
                           "30M100N70M", "100M", "100M", "60M2D40M"};

        try (SAMFileWriter w = new SAMFileWriterFactory().setCreateIndex(false)
                .makeBAMWriter(h, true, new File(prefix + ".bam"))) {
            for (int i = 0; i < n; i++) {
                SAMRecord r = new SAMRecord(h);
                r.setReadName("read" + i);
                r.setReferenceIndex(0);
                r.setAlignmentStart(starts[i]);

                String cigar = cigars[i % cigars.length];
                r.setCigarString(cigar);
                int readLen = 0;
                for (CigarElement e : r.getCigar().getCigarElements()) {
                    if (e.getOperator().consumesReadBases()) readLen += e.getLength();
                }

                byte[] bases = new byte[readLen];
                for (int j = 0; j < readLen; j++) {
                    int roll = rng.nextInt(100);
                    bases[j] = roll < 2 ? (byte) 'N' : roll < 8 ? (byte) 'T' : (byte) 'A';
                }
                r.setReadBases(bases);
                byte[] q = new byte[readLen];
                for (int j = 0; j < readLen; j++) q[j] = (byte) (10 + rng.nextInt(30));
                r.setBaseQualities(q);
                r.setMappingQuality(rng.nextInt(100) < 10 ? rng.nextInt(20) : 60);

                int flags = 0;
                if (rng.nextInt(100) < 50) {
                    // Paired, mostly proper, first or second of pair.
                    flags |= 0x1;
                    if (rng.nextInt(100) < 90) flags |= 0x2;
                    flags |= (i % 2 == 0) ? 0x40 : 0x80;
                    r.setMateReferenceIndex(0);
                    r.setMateAlignmentStart(Math.max(1, starts[i] + 200));
                    r.setInferredInsertSize(rng.nextInt(100) < 2 ? 200000 : 300);
                    if (rng.nextInt(100) < 5) flags |= 0x8;   // mate unmapped
                    if (rng.nextInt(100) < 50) flags |= 0x20; // mate reverse
                }
                if (rng.nextInt(100) < 50) flags |= 0x10;     // reverse
                if (rng.nextInt(1000) < 5) flags |= 0x200;    // vendor failed
                if (rng.nextInt(1000) < 10) flags |= 0x800;   // supplementary
                if (rng.nextInt(1000) < 5) flags |= 0x100;    // secondary
                r.setFlags(flags);
                r.setAttribute("RG", "rg1");
                if (rng.nextInt(100) < 20) r.setAttribute("MQ", rng.nextInt(60));
                if (rng.nextInt(1000) < 5) r.setAttribute("SA", "chr1,1,+,50M,60,0;");
                w.addAlignment(r);
            }
        }
        System.err.println("wrote " + n + " reads to " + prefix + ".bam");
    }
}
