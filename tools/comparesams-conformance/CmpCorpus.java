import java.io.*; import java.nio.file.*;
public class CmpCorpus {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static void run(String c,String sam1,String sam2) throws Exception {
    File d=Files.createTempDirectory("c").toFile(); d.deleteOnExit();
    File a=new File(d,"a.sam"); try(PrintStream p=new PrintStream(a)){p.print(sam1);}
    File b=new File(d,"b.sam"); try(PrintStream p=new PrintStream(b)){p.print(sam2);}
    File o=new File(d,"o.tsv");
    ByteArrayOutputStream bout=new ByteArrayOutputStream(); PrintStream old=System.out; System.setOut(new PrintStream(bout));
    int rc=-99;
    try { rc=new picard.sam.CompareSAMs().instanceMain(new String[]{a.getAbsolutePath(),b.getAbsolutePath(),"O="+o.getAbsolutePath()}); }
    catch(Throwable t){ System.setOut(old); emit("input1",c,sam1); emit("input2",c,sam2); emit("verdict",c,"THREW:"+t.getClass().getSimpleName()); return; }
    System.setOut(old);
    String verdict=bout.toString().trim();
    // Canonicalize the report: strip banner (## StringHeader / # lines), replace paths.
    String rep=new String(Files.readAllBytes(o.toPath()))
      .replace(a.getAbsolutePath(),"LEFT").replace(b.getAbsolutePath(),"RIGHT");
    StringBuilder cr=new StringBuilder();
    for(String line: rep.split("\n",-1)){ if(line.startsWith("## htsjdk")||line.startsWith("# ")) continue; cr.append(line).append("\n"); }
    emit("input1",c,sam1); emit("input2",c,sam2); emit("verdict",c,verdict+" rc="+rc); emit("report",c,cr.toString());
  }
  public static void main(String[] x) throws Exception {
    String h="@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n";
    String a="a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    String bb="b\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    run("equal", h+a+bb, h+a+bb);
    run("differ_pos", h+a, h+"a\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    run("unmapped_left", h+"a\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n", h+a);
    run("unmapped_both", h+"a\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n", h+"a\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
    run("missing_left_trailing", h+a, h+a+bb);
    run("dup_differ", h+a, h+"a\t1024\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    run("strand_differs_still_match", h+a, h+"a\t16\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    System.out.print(buf);
  }
}
