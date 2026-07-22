/* Oracle corpus generator for ValidateSamFile SUMMARY mode. */
import java.io.*; import java.nio.file.*;
public class SummaryCorpus {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static void run(String c,String sam) throws Exception {
    File d=Files.createTempDirectory("s").toFile(); d.deleteOnExit();
    File in=new File(d,"in.sam"); try(PrintStream p=new PrintStream(in)){p.print(sam);}
    File o=new File(d,"o.txt");
    new picard.sam.ValidateSamFile().instanceMain(new String[]{"INPUT="+in,"OUTPUT="+o,"MODE=SUMMARY"});
    emit("input",c,sam); emit("output",c,new String(Files.readAllBytes(o.toPath())));
  }
  public static void main(String[] a) throws Exception {
    String sq="@SQ\tSN:chr1\tLN:100\n"; String rg="@RG\tID:rg1\tSM:s\tPL:illumina\n";
    run("summary_clean", "@HD\tVN:1.6\tSO:coordinate\n"+sq+rg+"a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    run("summary_errors", "@HD\tVN:1.6\tSO:coordinate\n"+sq
      +"a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
      +"b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\t*\n");
    run("summary_flags", "@HD\tVN:1.6\tSO:coordinate\n"+sq+rg
      +"a\t2\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n"
      +"b\t64\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tNM:i:0\n");
    System.out.print(buf);
  }
}
