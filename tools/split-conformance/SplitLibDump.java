import java.io.*; import java.nio.file.*; import java.util.*;
public class SplitLibDump {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static void run(String c,String sam) throws Exception {
    File d=Files.createTempDirectory("splitlib").toFile(); d.deleteOnExit();
    File in=new File(d,"in.sam"); try(PrintStream p=new PrintStream(in)){p.print(sam);}
    File out=new File(d,"out"); out.mkdirs();
    int rc=new picard.sam.SplitSamByLibrary().instanceMain(new String[]{
      "I="+in.getAbsolutePath(), "O="+out.getAbsolutePath()});
    emit("input",c,sam); emit("rc",c,String.valueOf(rc));
    File[] files=out.listFiles((dir,n)->n.endsWith(".sam"));
    if(files!=null){
      Arrays.sort(files, Comparator.comparing(File::getName));
      for(File f: files){
        String base=f.getName().substring(0, f.getName().length()-4); // strip .sam
        emit("file:"+base,c,new String(Files.readAllBytes(f.toPath())));
      }
    }
  }
  public static void main(String[] x) throws Exception {
    String hd="@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n";
    // two libraries, one record each
    run("two_libs", hd+"@RG\tID:rg1\tLB:libA\tSM:s\n@RG\tID:rg2\tLB:libB\tSM:s\n"
      +"a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
      +"b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n");
    // a read group with no LB, and a record with no RG -> unknown
    run("with_unknown", hd+"@RG\tID:rg1\tLB:libA\tSM:s\n@RG\tID:rg2\tSM:s\n"
      +"a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
      +"b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n"
      +"c\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    // two read groups share a library -> both @RG in that file
    run("shared_lib", hd+"@RG\tID:rg1\tLB:libA\tSM:s\n@RG\tID:rg2\tLB:libA\tSM:s\n"
      +"a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
      +"b\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n");
    // a library with no records still gets an (empty-body) file
    run("empty_lib", hd+"@RG\tID:rg1\tLB:libA\tSM:s\n@RG\tID:rg2\tLB:libB\tSM:s\n"
      +"a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    // library name with a reserved char -> makeFileNameSafe replaces it
    run("safe_name", hd+"@RG\tID:rg1\tLB:lib/A:1\tSM:s\n"
      +"a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    System.out.print(buf);
  }
}
