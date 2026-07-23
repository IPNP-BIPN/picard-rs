import java.io.*; import java.nio.file.*; import java.util.*;
public class SplitDump {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  // opts: space-separated KEY=VALUE for SplitSamByNumberOfReads.
  static void run(String c,String opts,String sam) throws Exception {
    File d=Files.createTempDirectory("split").toFile(); d.deleteOnExit();
    File in=new File(d,"in.sam"); try(PrintStream p=new PrintStream(in)){p.print(sam);}
    File out=new File(d,"out"); out.mkdirs();
    List<String> args=new ArrayList<>();
    args.add("I="+in.getAbsolutePath()); args.add("O="+out.getAbsolutePath());
    if(!opts.isEmpty()) for(String e:opts.split(" ")) args.add(e);
    int rc=new picard.sam.SplitSamByNumberOfReads().instanceMain(args.toArray(new String[0]));
    emit("opts",c,opts); emit("input",c,sam); emit("rc",c,String.valueOf(rc));
    // shards are named shard_0001.sam, shard_0002.sam, ... in order.
    File[] shards=out.listFiles((dir,n)->n.endsWith(".sam"));
    Arrays.sort(shards, Comparator.comparing(File::getName));
    for(File s: shards){ emit("shard",c,new String(Files.readAllBytes(s.toPath()))); }
  }
  public static void main(String[] x) throws Exception {
    String h="@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n";
    StringBuilder recs=new StringBuilder();
    String[] names={"a","b","c","d","e","f"};
    for(int i=0;i<names.length;i++) recs.append(names[i]).append("\t0\tchr1\t").append(1+i).append("\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    run("n_files_2","SPLIT_TO_N_FILES=2",h+recs);
    run("n_files_3","SPLIT_TO_N_FILES=3",h+recs);
    run("n_reads_2","SPLIT_TO_N_READS=2",h+recs);
    // queryname group straddles the boundary and stays together: a, b, b, c with N_READS=2.
    // reads_per_file=2, but the second "b" is not a name change so shard 1 keeps all three.
    String g=h+"a\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\nb\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
      +"b\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\tIIII\nc\t0\tchr1\t4\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    run("group_kept","SPLIT_TO_N_READS=2",g);
    // total-reads override changes readsPerFile.
    run("total_override","SPLIT_TO_N_FILES=2 TOTAL_READS=6",h+recs);
    System.out.print(buf);
  }
}
