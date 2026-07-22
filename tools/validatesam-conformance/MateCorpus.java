/*
 * Oracle corpus generator for ValidateSamFile mate-pair validation (VERBOSE mode, SAM input, no
 * reference) conformance. Ports the checks in SamFileValidator.validateMateFields /
 * PairEndInfo.validateMates: a clean matched pair, a paired read whose mate is absent, and the
 * mate-field mismatches (alignment start, CIGAR via the MC tag, mate negative-strand flag) and the
 * both-marked-first-of-pair case. The verbose output is raw and has no timestamp or banner.
 *
 *   java -cp picard-fat.jar:. MateCorpus > corpus.txt
 */
import java.io.*;
import java.nio.file.*;

public class MateCorpus {
    static StringBuilder buf = new StringBuilder();

    static void emit(String k, String c, String p) {
        buf.append(k).append('\t').append(c).append('\t')
           .append(p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")).append('\n');
    }

    static void run(String c, String sam) throws Exception {
        File d = Files.createTempDirectory("m").toFile();
        d.deleteOnExit();
        File in = new File(d, "in.sam");
        try (PrintStream p = new PrintStream(in)) { p.print(sam); }
        File o = new File(d, "o.txt");
        new picard.sam.ValidateSamFile().instanceMain(new String[]{
            "INPUT=" + in, "OUTPUT=" + o, "MODE=VERBOSE"});
        emit("input", c, sam);
        emit("output", c, new String(Files.readAllBytes(o.toPath())));
    }

    public static void main(String[] a) throws Exception {
        String h = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:100\n"
                 + "@RG\tID:rg1\tSM:s\tPL:illumina\n";

        // Clean matched pair.
        run("pair_clean", h
            + "p\t99\tchr1\t10\t60\t4M\t=\t20\t14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\tMC:Z:4M\n"
            + "p\t147\tchr1\t20\t60\t4M\t=\t10\t-14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\tMC:Z:4M\n");

        // Paired read whose mate is absent.
        run("mate_not_found", h
            + "p\t99\tchr1\t10\t60\t4M\t=\t20\t14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        // Mate alignment-start mismatch: mate is actually at 30, read1 records it at 20.
        run("mismatch_mate_start", h
            + "p\t99\tchr1\t10\t60\t4M\t=\t20\t14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n"
            + "p\t147\tchr1\t30\t60\t4M\t=\t10\t-14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        // Mate CIGAR mismatch via the MC tag.
        run("mismatch_mate_cigar", h
            + "p\t99\tchr1\t10\t60\t4M\t=\t20\t14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\tMC:Z:5M\n"
            + "p\t147\tchr1\t20\t60\t4M\t=\t10\t-14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\tMC:Z:4M\n");

        // Mate negative-strand flag mismatch: read1 marks its mate reverse, but read2 is forward.
        run("mismatch_mate_neg_strand", h
            + "p\t99\tchr1\t10\t60\t4M\t=\t20\t14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n"
            + "p\t131\tchr1\t20\t60\t4M\t=\t10\t-14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        // Both mates marked as first of pair.
        run("mates_same_end", h
            + "p\t99\tchr1\t10\t60\t4M\t=\t20\t14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n"
            + "p\t83\tchr1\t20\t60\t4M\t=\t10\t-14\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");

        System.out.print(buf);
    }
}
