/*
 * Oracle corpus generator for ValidateSamFile's reference-based NM value check (INVALID_TAG_NM),
 * VERBOSE mode. Each case carries the reference FASTA, the SAM input, and the exact verbose output.
 */
import java.io.*; import java.nio.file.*;
public class NmValCorpus {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static void run(String c,String fasta,String sam) throws Exception {
    File dir=Files.createTempDirectory("nmv").toFile(); dir.deleteOnExit();
    File fa=new File(dir,"ref.fasta"); try(PrintStream p=new PrintStream(fa)){p.print(fasta);}
    try(PrintStream p=new PrintStream(new File(dir,"ref.fasta.fai"))){p.print("chr1\t40\t6\t40\t41\n");}
    try(PrintStream p=new PrintStream(new File(dir,"ref.dict"))){p.print("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:40\n");}
    File in=new File(dir,"in.sam"); try(PrintStream p=new PrintStream(in)){p.print(sam);}
    File o=new File(dir,"o.txt");
    new picard.sam.ValidateSamFile().instanceMain(new String[]{
      "INPUT="+in,"OUTPUT="+o,"MODE=VERBOSE","REFERENCE_SEQUENCE="+fa.getAbsolutePath()});
    emit("fasta",c,fasta); emit("input",c,sam); emit("output",c,new String(Files.readAllBytes(o.toPath())));
  }
  public static void main(String[] a) throws Exception {
    String fasta=">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
    String h="@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:40\n@RG\tID:rg1\tSM:s\tPL:illumina\n";
    run("nm_correct", fasta, h+"a\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\n");
    run("nm_wrong_low", fasta, h+"a\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACCTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:0\n");
    run("nm_wrong_high", fasta, h+"a\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACCTACGT\tIIIIIIII\tRG:Z:rg1\tNM:i:5\n");
    // deletion: 4M2D4M, 4M2D4M read: reality NM = 2 deleted + 4 mismatched = 6
    run("nm_del_tag2", fasta, h+"a\t0\tchr1\t1\t60\t4M2D4M\t*\t0\t0\tACGTCGTA\tIIIIIIII\tRG:Z:rg1\tNM:i:2\n");
    // 
    run("nm_del_tag0", fasta, h+"a\t0\tchr1\t1\t60\t4M2D4M\t*\t0\t0\tACGTCGTA\tIIIIIIII\tRG:Z:rg1\tNM:i:0\n");
    System.out.print(buf);
  }
}
