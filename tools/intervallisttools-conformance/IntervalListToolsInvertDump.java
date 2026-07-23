/*
 * Oracle dump for IntervalListTools invert-dependent paths (SUBTRACT / SYMDIFF / the INVERT option)
 * against Picard 3.4.0. These all route through IntervalList.invert, which needs contig lengths.
 * Emits the inputs, action, invert flag, and output (with @PG, stripped by the conformance).
 *   java -cp picard-fat.jar:. IntervalListToolsInvertDump | gzip -n > interval_list_tools_invert.txt.gz
 */
import java.io.*;
import java.nio.file.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class IntervalListToolsInvertDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")).append('\n');
  }

  static File dir, dict;

  static File makeDict(String contig, int length) throws Exception {
    StringBuilder seq = new StringBuilder();
    for (int i = 0; i < length; i++) seq.append("ACGT".charAt(i % 4));
    File fasta = new File(dir, contig + ".fasta");
    Files.write(fasta.toPath(), (">" + contig + "\n" + seq + "\n").getBytes());
    File d = new File(dir, contig + ".dict");
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
        "R=" + fasta.getAbsolutePath(), "O=" + d.getAbsolutePath()});
    return d;
  }

  static File makeInput(String name, Interval... intervals) throws Exception {
    SAMSequenceDictionary sd = SAMSequenceDictionaryExtractor.extractDictionary(dict.toPath());
    IntervalList list = new IntervalList(new SAMFileHeader(sd));
    for (Interval iv : intervals) list.add(iv);
    File f = new File(dir, name + ".interval_list");
    list.write(f);
    return f;
  }

  static void run(String c, String action, boolean invert, File in1, File second) throws Exception {
    File out = new File(dir, c + ".out.interval_list");
    java.util.List<String> args = new java.util.ArrayList<>();
    args.add("I=" + in1.getAbsolutePath());
    if (second != null) args.add("SECOND_INPUT=" + second.getAbsolutePath());
    args.add("O=" + out.getAbsolutePath());
    args.add("ACTION=" + action);
    args.add("INVERT=" + invert);
    new picard.util.IntervalListTools().instanceMain(args.toArray(new String[0]));
    emit("input1", c, new String(Files.readAllBytes(in1.toPath())));
    emit("input2", c, second == null ? "" : new String(Files.readAllBytes(second.toPath())));
    emit("action", c, action);
    emit("invert", c, String.valueOf(invert));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }

  public static void main(String[] x) throws Exception {
    dir = Files.createTempDirectory("ilti").toFile();
    dict = makeDict("chr1", 100);

    // INVERT option with default ACTION=CONCAT: complement of {10-20, 50-60} on chr1[1..100].
    File a = makeInput("a",
        new Interval("chr1", 10, 20, false, "A"),
        new Interval("chr1", 50, 60, false, "C"));
    run("invert_concat", "CONCAT", true, a, null);

    // SUBTRACT: INPUT {1-50} minus SECOND_INPUT {20-30} = {1-19, 31-50}.
    File big = makeInput("big", new Interval("chr1", 1, 50, false, "A"));
    File mid = makeInput("mid", new Interval("chr1", 20, 30, false, "B"));
    run("subtract", "SUBTRACT", false, big, mid);

    // SYMDIFF: {1-30} xor {20-50} = {1-19, 31-50}.
    File left = makeInput("left", new Interval("chr1", 1, 30, false, "A"));
    File right = makeInput("right", new Interval("chr1", 20, 50, false, "B"));
    run("symdiff", "SYMDIFF", false, left, right);

    System.out.print(buf);
  }
}
