/*
 * Full-file oracle for MergeBamAlignment (unpaired and paired, coordinate output) against Picard
 * 3.4.0. Emits the reference .dict, the unmapped and aligned SAM, and the complete output SAM
 * (header + records). The port reads the same committed .dict, so the @SQ (including M5 and the
 * absolute UR path) matches byte-for-byte without canonicalization.
 *
 *   java -cp picard-fat.jar:. MergeBamAlignmentFullDump | gzip > merge_bam_alignment_full.txt.gz
 */
import java.io.*; import java.nio.file.*; import java.util.*;
public class MergeBamAlignmentFullDump {
  static final String REF = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static void run(String c, String unmapped, String aligned) throws Exception {
    File d = Files.createTempDirectory("mbaf").toFile();
    File ref = new File(d, "ref.fasta");
    try (PrintStream p = new PrintStream(ref)) { p.println(">chr1"); p.println(REF); }
    File dict = new File(d, "ref.dict");
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
      "R=" + ref.getAbsolutePath(), "O=" + dict.getAbsolutePath()});
    File u = new File(d, "u.sam"); Files.write(u.toPath(), unmapped.getBytes());
    File a = new File(d, "a.sam"); Files.write(a.toPath(), aligned.getBytes());
    File o = new File(d, "o.sam");
    int rc = new picard.sam.MergeBamAlignment().instanceMain(new String[]{
      "UNMAPPED_BAM=" + u.getAbsolutePath(), "ALIGNED_BAM=" + a.getAbsolutePath(),
      "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(), "OUTPUT=" + o.getAbsolutePath()});
    emit("dict", c, new String(Files.readAllBytes(dict.toPath())));
    emit("unmapped", c, unmapped);
    emit("aligned", c, aligned);
    emit("rc", c, String.valueOf(rc));
    emit("full", c, new String(Files.readAllBytes(o.toPath())));
  }
  public static void main(String[] x) throws Exception {
    String uHdr = "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n";
    String aHdr = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";
    run("unpaired",
        uHdr + "a\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tab:Z:keepme\n" +
               "b\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n",
        aHdr + "a\t0\tchr1\t20\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n" +
               "b\t0\tchr1\t5\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n");
    run("paired",
        uHdr + "p1\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n" +
               "p1\t141\t*\t0\t0\t*\t*\t0\t0\tTTTTGGGG\tIIIIIIII\tRG:Z:rg1\n",
        aHdr + "p1\t99\tchr1\t1\t60\t8M\t=\t20\t27\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n" +
               "p1\t147\tchr1\t20\t60\t8M\t=\t1\t-27\tCCCCAAAA\tIIIIIIII\tRG:Z:rg1\n");
    // Unmapped-read passthrough: read 'a' aligns, read 'z' has no alignment (unmapped in the aligned
    // BAM) and stays unmapped in the output, carrying only PG/RG, sorted after the mapped read.
    run("unmapped_passthrough",
        uHdr + "a\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n" +
               "z\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGCCCC\tJJJJJJJJ\tRG:Z:rg1\n",
        aHdr + "a\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n" +
               "z\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGCCCC\tJJJJJJJJ\tRG:Z:rg1\n");
    System.out.print(buf);
  }
}
