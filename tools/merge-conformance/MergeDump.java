import java.io.*; import java.nio.file.*; import java.util.*;
public class MergeDump {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static void run(String c,String so,String... sams) throws Exception { run(c,so,false,sams); }
  static void run(String c,String so,boolean msd,String... sams) throws Exception {
    File d=Files.createTempDirectory("mg").toFile();
    List<String> args=new ArrayList<>();
    for(int i=0;i<sams.length;i++){File f=new File(d,"in"+i+".sam");try(PrintStream p=new PrintStream(f)){p.print(sams[i]);}args.add("I="+f.getAbsolutePath());}
    File o=new File(d,"o.sam"); args.add("O="+o.getAbsolutePath()); args.add("SORT_ORDER="+so);
    args.add("MERGE_SEQUENCE_DICTIONARIES="+msd);
    int rc=new picard.sam.MergeSamFiles().instanceMain(args.toArray(new String[0]));
    emit("so",c,so); emit("msd",c,String.valueOf(msd)); for(int i=0;i<sams.length;i++) emit("input",c,sams[i]);
    emit("rc",c,String.valueOf(rc));
    emit("output",c,new String(Files.readAllBytes(o.toPath())));
  }
  public static void main(String[] x) throws Exception {
    String h="@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s\tLB:lib1\n";
    String a=h+"a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\nc\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    String b=h+"b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\nd\t0\tchr1\t40\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    run("coord_two",  "coordinate", a, b);
    // equal-coordinate reads across files: both at chr1:10, different names -> comparator orders them.
    String e1=h+"x\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    String e2=h+"y\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    run("equal_coord", "coordinate", e1, e2);
    // three files
    run("coord_three","coordinate", h+"a\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n",
                                    h+"b\t0\tchr1\t15\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n",
                                    h+"c\t0\tchr1\t25\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    // queryname sort order (inputs queryname-sorted)
    String hq="@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s\tLB:lib1\n";
    run("queryname_two","queryname", hq+"a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\nc\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n",
                                     hq+"b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\nd\t0\tchr1\t40\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    // shared header with a @PG and a @CO
    String hp="@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s\tLB:lib1\n@PG\tID:p1\tPN:tool\tVN:1\n@CO\ta comment\n";
    run("with_pg_co","coordinate", hp+"a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n",
                                   hp+"b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n");
    // distinct read groups across files (no collision): the merged header unions both @RG.
    String sq="@SQ\tSN:chr1\tLN:1000\n";
    String ua="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s1\tLB:lib1\n"+
      "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    String ub="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg2\tSM:s2\tLB:lib2\n"+
      "b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg2\n";
    run("distinct_read_groups","coordinate", ua, ub);
    // colliding read groups: same ID:rg1, different content -> rg1, rg1.1; the 2nd file's records
    // have their RG:Z rewritten to rg1.1.
    String ca="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s1\tLB:lib1\n"+
      "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    String cb="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s2\tLB:lib2\n"+
      "b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    run("collide_read_groups","coordinate", ca, cb);
    // three colliding read groups -> rg1, rg1.1, rg1.2.
    String cc="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s3\tLB:lib3\n"+
      "c\t0\tchr1\t30\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n";
    run("collide_read_groups_three","coordinate", ca, cb, cc);
    // colliding program groups (no PP): same ID:p1, different content -> p1, p1.1; records' PG:Z
    // rewritten.
    String pa="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s\tLB:lib1\n@PG\tID:p1\tPN:tool\tVN:1\n"+
      "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tPG:Z:p1\n";
    String pb="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s\tLB:lib1\n@PG\tID:p1\tPN:tool\tVN:2\n"+
      "b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tPG:Z:p1\n";
    run("collide_program_groups","coordinate", pa, pb);
    // MERGE_SEQUENCE_DICTIONARIES: file 0 knows chr1,chr2; file 1 knows chr2,chr3 -> superset
    // chr1,chr2,chr3, and file 1's records are remapped by name.
    String da="@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n"+
      "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"+
      "b\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    String db="@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chr3\tLN:3000\n"+
      "c\t0\tchr2\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\n"+
      "d\t0\tchr3\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    run("merge_dictionaries","coordinate", true, da, db);
    // @PG PP-chain, identical in both files -> deduped to one p1 (root) and one p2 (PP:p1).
    String ga="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s\tLB:lib1\n@PG\tID:p1\tPN:a\n@PG\tID:p2\tPN:b\tPP:p1\n"+
      "a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tPG:Z:p2\n";
    String gb="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@RG\tID:rg1\tSM:s\tLB:lib1\n@PG\tID:p1\tPN:a\n@PG\tID:p2\tPN:b\tPP:p1\n"+
      "b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\tPG:Z:p2\n";
    run("chain_dedup","coordinate", ga, gb);
    // @PG PP-chain where the root p1 differs -> p1, p1.1; each p2 rechains to its own parent, so the
    // two p2 collide and the child is renamed p2.3 (shared counter is at 2 after the roots).
    String ha="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@PG\tID:p1\tPN:a\tVN:1\n@PG\tID:p2\tPN:b\tPP:p1\n"+
      "r0\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2\n";
    String hb="@HD\tVN:1.6\tSO:coordinate\n"+sq+"@PG\tID:p1\tPN:a\tVN:2\n@PG\tID:p2\tPN:b\tPP:p1\n"+
      "r1\t0\tchr1\t20\t60\t4M\t*\t0\t0\tACGT\tIIII\tPG:Z:p2\n";
    run("chain_collide","coordinate", ha, hb);
    System.out.print(buf);
  }
}
