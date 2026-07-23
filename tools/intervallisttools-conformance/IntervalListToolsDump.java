/*
 * Oracle dump for IntervalListTools (CONCAT / UNION slice) against Picard 3.4.0. For each case it
 * builds a shared "from" .dict and two input interval_lists on it, runs IntervalListTools with the
 * given ACTION/SORT/UNIQUE/DONT_MERGE_ABUTTING, and emits the two inputs, the option values, and the
 * output interval_list. The output carries a @PG whose CL is the command line (non-reproducible); the
 * conformance and this dump strip @PG lines before comparing.
 *   java -cp picard-fat.jar:. IntervalListToolsDump | gzip -n > interval_list_tools.txt.gz
 */
import java.io.*;
import java.nio.file.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class IntervalListToolsDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")).append('\n');
  }

  static File dir;
  static File dict;

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

  static void run(String c, File in1, File in2, String action, boolean sort, boolean unique,
                  boolean dontMergeAbutting) throws Exception {
    File out = new File(dir, c + ".out.interval_list");
    new picard.util.IntervalListTools().instanceMain(new String[]{
        "I=" + in1.getAbsolutePath(),
        "I=" + in2.getAbsolutePath(),
        "O=" + out.getAbsolutePath(),
        "ACTION=" + action,
        "SORT=" + sort,
        "UNIQUE=" + unique,
        "DONT_MERGE_ABUTTING=" + dontMergeAbutting});
    emit("input1", c, new String(Files.readAllBytes(in1.toPath())));
    emit("input2", c, new String(Files.readAllBytes(in2.toPath())));
    emit("action", c, action);
    emit("sort", c, String.valueOf(sort));
    emit("unique", c, String.valueOf(unique));
    emit("dont_merge_abutting", c, String.valueOf(dontMergeAbutting));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }

  public static void main(String[] x) throws Exception {
    dir = Files.createTempDirectory("ilt").toFile();
    dict = makeDict("chr1", 100);

    // input1: 1-10 A, 50-60 C ; input2: 5-15 B (overlaps A), 61-70 D (abuts C at 60|61).
    File in1 = makeInput("in1",
        new Interval("chr1", 1, 10, false, "A"),
        new Interval("chr1", 50, 60, false, "C"));
    File in2 = makeInput("in2",
        new Interval("chr1", 5, 15, false, "B"),
        new Interval("chr1", 61, 70, false, "D"));

    // CONCAT, default SORT=true, UNIQUE=false: sorted concatenation, no merging.
    run("concat_sorted", in1, in2, "CONCAT", true, false, false);
    // CONCAT, SORT=false: input file order preserved (in1 then in2).
    run("concat_unsorted", in1, in2, "CONCAT", false, false, false);
    // UNION: sort + uniqued(merge abutting, concat names). 1-15 "A|B", 50-70 "C|D".
    run("union", in1, in2, "UNION", true, true, false);
    // CONCAT + UNIQUE + DONT_MERGE_ABUTTING: merge overlaps but keep abutting 50-60 / 61-70 apart.
    run("unique_no_abut", in1, in2, "CONCAT", true, true, true);

    System.out.print(buf);
  }
}
