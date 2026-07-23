/*
 * Oracle dump for IntervalListTools set-ops (INTERSECT / OVERLAPS) against Picard 3.4.0. INTERSECT
 * intersects two INPUTs; OVERLAPS keeps whole INPUT intervals that overlap any SECOND_INPUT interval.
 * Emits the inputs, action, and output (with @PG, stripped by the conformance).
 *   java -cp picard-fat.jar:. IntervalListToolsSetOpsDump | gzip -n > interval_list_tools_setops.txt.gz
 */
import java.io.*;
import java.nio.file.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class IntervalListToolsSetOpsDump {
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

  static void runIntersect(String c, File in1, File in2) throws Exception {
    File out = new File(dir, c + ".out.interval_list");
    new picard.util.IntervalListTools().instanceMain(new String[]{
        "I=" + in1.getAbsolutePath(), "I=" + in2.getAbsolutePath(),
        "O=" + out.getAbsolutePath(), "ACTION=INTERSECT"});
    emit("input1", c, new String(Files.readAllBytes(in1.toPath())));
    emit("input2", c, new String(Files.readAllBytes(in2.toPath())));
    emit("action", c, "INTERSECT");
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }

  static void runOverlaps(String c, File in1, File second) throws Exception {
    File out = new File(dir, c + ".out.interval_list");
    new picard.util.IntervalListTools().instanceMain(new String[]{
        "I=" + in1.getAbsolutePath(), "SECOND_INPUT=" + second.getAbsolutePath(),
        "O=" + out.getAbsolutePath(), "ACTION=OVERLAPS"});
    emit("input1", c, new String(Files.readAllBytes(in1.toPath())));
    emit("input2", c, new String(Files.readAllBytes(second.toPath())));
    emit("action", c, "OVERLAPS");
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }

  public static void main(String[] x) throws Exception {
    dir = Files.createTempDirectory("ilts").toFile();
    dict = makeDict("chr1", 100);

    // INTERSECT: in1 = A(1-20), C(50-60); in2 = B(10-30), D(55-70).
    // B.intersect(A)=10-20 "B intersection A"; D.intersect(C)=55-60 "D intersection C".
    File a = makeInput("a",
        new Interval("chr1", 1, 20, false, "A"),
        new Interval("chr1", 50, 60, false, "C"));
    File b = makeInput("b",
        new Interval("chr1", 10, 30, false, "B"),
        new Interval("chr1", 55, 70, false, "D"));
    runIntersect("intersect", a, b);

    // OVERLAPS: keep whole INPUT intervals overlapping any SECOND_INPUT interval.
    // in1 = A(1-20), C(50-60), E(80-90); second = B(10-30). Only A overlaps.
    File a2 = makeInput("a2",
        new Interval("chr1", 1, 20, false, "A"),
        new Interval("chr1", 50, 60, false, "C"),
        new Interval("chr1", 80, 90, false, "E"));
    File second = makeInput("second",
        new Interval("chr1", 10, 30, false, "B"));
    runOverlaps("overlaps", a2, second);

    System.out.print(buf);
  }
}
