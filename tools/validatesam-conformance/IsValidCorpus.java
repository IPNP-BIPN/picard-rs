import java.io.*; import java.nio.file.*;
public class IsValidCorpus {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static void run(String c,String sam) throws Exception {
    File d=Files.createTempDirectory("v").toFile(); d.deleteOnExit();
    File in=new File(d,"in.sam"); try(PrintStream p=new PrintStream(in)){p.print(sam);}
    File o=new File(d,"o.txt");
    new picard.sam.ValidateSamFile().instanceMain(new String[]{"INPUT="+in,"OUTPUT="+o,"MODE=VERBOSE"});
    emit("input",c,sam); emit("output",c,new String(Files.readAllBytes(o.toPath())));
  }
  public static void main(String[] a) throws Exception {
    String h="@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\tSM:s\tPL:illumina\n";
    run("proper_pair_unpaired", h+"a\t2\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    run("mate_unmapped_unpaired", h+"a\t8\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    run("mate_neg_unpaired", h+"a\t32\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    run("first_of_pair_unpaired", h+"a\t64\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    run("second_of_pair_unpaired", h+"a\t128\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    run("mapq_on_unmapped", h+"a\t4\t*\t0\t60\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    run("unmapped_secondary", h+"a\t260\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    run("unmapped_supplementary", h+"a\t2052\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    run("mapped_no_cigar", h+"a\t0\tchr1\t10\t60\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    run("rg_not_found", h+"a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:other\tNM:i:0\n");
    System.out.print(buf);
  }
}
