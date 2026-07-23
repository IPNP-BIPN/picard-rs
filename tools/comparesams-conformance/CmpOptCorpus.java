import java.io.*; import java.nio.file.*;
public class CmpOptCorpus {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static void run(String c,String opts,String sam1,String sam2) throws Exception {
    File d=Files.createTempDirectory("c").toFile(); d.deleteOnExit();
    File a=new File(d,"a.sam"); try(PrintStream p=new PrintStream(a)){p.print(sam1);}
    File b=new File(d,"b.sam"); try(PrintStream p=new PrintStream(b)){p.print(sam2);}
    File o=new File(d,"o.tsv");
    java.util.List<String> args=new java.util.ArrayList<>();
    args.add(a.getAbsolutePath()); args.add(b.getAbsolutePath()); args.add("O="+o.getAbsolutePath());
    if(!opts.isEmpty()) for(String e:opts.split(" ")) args.add(e);
    ByteArrayOutputStream bout=new ByteArrayOutputStream(); PrintStream old=System.out; System.setOut(new PrintStream(bout));
    int rc=new picard.sam.CompareSAMs().instanceMain(args.toArray(new String[0]));
    System.setOut(old);
    String rep=new String(Files.readAllBytes(o.toPath())).replace(a.getAbsolutePath(),"LEFT").replace(b.getAbsolutePath(),"RIGHT");
    StringBuilder cr=new StringBuilder();
    for(String line: rep.split("\n",-1)){ if(line.startsWith("## htsjdk")||line.startsWith("# ")) continue; cr.append(line).append("\n"); }
    emit("opts",c,opts); emit("input1",c,sam1); emit("input2",c,sam2);
    emit("verdict",c,bout.toString().trim()+" rc="+rc); emit("report",c,cr.toString());
  }
  public static void main(String[] x) throws Exception {
    String h="@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n";
    String lowa="a\t0\tchr1\t10\t2\t4M\t*\t0\t0\tACGT\tIIII\n";
    String lowb="a\t0\tchr1\t50\t2\t4M\t*\t0\t0\tACGT\tIIII\n";
    run("lenient_lowmq_match","LENIENT_LOW_MQ_ALIGNMENT=true",h+lowa,h+lowb);
    run("lenient_lowmq_strict","",h+lowa,h+lowb);
    // both reads unknown mapping quality (255) at different positions: differ strictly, match leniently.
    String unka="a\t0\tchr1\t10\t255\t4M\t*\t0\t0\tACGT\tIIII\n";
    String unkb="a\t0\tchr1\t50\t255\t4M\t*\t0\t0\tACGT\tIIII\n";
    run("lenient_unknown_match","LENIENT_UNKNOWN_MQ_ALIGNMENT=true",h+unka,h+unkb);
    run("lenient_unknown_strict","",h+unka,h+unkb);
    String m1="a\t0\tchr1\t10\t40\t4M\t*\t0\t0\tACGT\tIIII\n";
    String m2="a\t0\tchr1\t10\t30\t4M\t*\t0\t0\tACGT\tIIII\n";
    run("compare_mq_single","COMPARE_MQ=true",h+m1,h+m2);
    // two reads with distinct mq pairs -> histogram with two sorted bins
    String pa="p\t0\tchr1\t10\t40\t4M\t*\t0\t0\tACGT\tIIII\nq\t0\tchr1\t20\t10\t4M\t*\t0\t0\tACGT\tIIII\n";
    String pb="p\t0\tchr1\t10\t35\t4M\t*\t0\t0\tACGT\tIIII\nq\t0\tchr1\t20\t10\t4M\t*\t0\t0\tACGT\tIIII\n";
    run("compare_mq_multi","COMPARE_MQ=true",h+pa,h+pb);
    System.out.print(buf);
  }
}
