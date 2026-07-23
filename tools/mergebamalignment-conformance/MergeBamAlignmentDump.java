/*
 * Oracle dump for MergeBamAlignment's merged record (unpaired, single primary hit, coordinate output)
 * against Picard 3.4.0. Emits an escaped TSV of the merged record data line(s), one per case:
 *   unmapped\t<case>\t<unmapped sam>
 *   aligned\t<case>\t<aligned sam>
 *   record\t<case>\t<merged record data line>
 * The header (@SQ/@RG/@PG) is a later slice and is not compared. The reference is a fixed 40 bp chr1.
 *
 *   java -cp picard-fat.jar:. MergeBamAlignmentDump | gzip > merge_bam_alignment.txt.gz
 */
import java.io.*; import java.nio.file.*; import java.util.*;
public class MergeBamAlignmentDump {
  static final String REF = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // chr1, 40 bp
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static void run(String c, String unmapped, String aligned) throws Exception {
    File d = Files.createTempDirectory("mba").toFile();
    File ref = new File(d, "ref.fasta");
    try (PrintStream p = new PrintStream(ref)) { p.println(">chr1"); p.println(REF); }
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
      "R=" + ref.getAbsolutePath(), "O=" + new File(d, "ref.dict").getAbsolutePath()});
    File u = new File(d, "u.sam"); Files.write(u.toPath(), unmapped.getBytes());
    File a = new File(d, "a.sam"); Files.write(a.toPath(), aligned.getBytes());
    File o = new File(d, "o.sam");
    int rc = new picard.sam.MergeBamAlignment().instanceMain(new String[]{
      "UNMAPPED_BAM=" + u.getAbsolutePath(), "ALIGNED_BAM=" + a.getAbsolutePath(),
      "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(), "OUTPUT=" + o.getAbsolutePath()});
    emit("unmapped", c, unmapped);
    emit("aligned", c, aligned);
    emit("rc", c, String.valueOf(rc));
    for (String line : new String(Files.readAllBytes(o.toPath())).split("\n")) {
      if (!line.startsWith("@") && !line.isEmpty()) emit("record", c, line);
    }
  }
  public static void main(String[] x) throws Exception {
    String uHdr = "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n";
    String aHdr = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";
    // Perfect forward match.
    run("forward_match",
        uHdr + "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tab:Z:keepme\n",
        aHdr + "r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\tX0:i:1\n");
    // One mismatch at read position 2 (C->A), low quality there (#=2).
    run("forward_mismatch",
        uHdr + "r2\t4\t*\t0\t0\t*\t*\t0\t0\tAAGTACGT\t#IIIIIII\tRG:Z:rg1\n",
        aHdr + "r2\t0\tchr1\t1\t60\t8M\t*\t0\t0\tAAGTACGT\t#IIIIIII\tRG:Z:rg1\n");
    // Negative-strand hit over ref[0:6]=ACGTAC; original read is its reverse-complement GTACGT.
    run("reverse",
        uHdr + "r3\t4\t*\t0\t0\t*\t*\t0\t0\tGTACGT\tABCDEF\tRG:Z:rg1\n",
        aHdr + "r3\t16\tchr1\t1\t60\t6M\t*\t0\t0\tACGTAC\tFEDCBA\tRG:Z:rg1\n");
    System.out.print(buf);
  }
}
