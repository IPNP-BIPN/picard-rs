/*
 * Oracle dump for IntervalListTools PADDING and BREAK_BANDS_AT_MULTIPLES_OF against Picard 3.4.0.
 * PADDING pads each input interval (clamped to [1, contigLength]) before reducing; BREAK_BANDS splits
 * the final intervals at band multiples. Emits the input, options, and output (@PG stripped by the
 * conformance).
 *   java -cp picard-fat.jar:. IntervalListToolsPadBreakDump | gzip -n > interval_list_tools_padbreak.txt.gz
 */
import java.io.*;
import java.nio.file.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class IntervalListToolsPadBreakDump {
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

  static void run(String c, File in1, int padding, int breakBands) throws Exception {
    File out = new File(dir, c + ".out.interval_list");
    new picard.util.IntervalListTools().instanceMain(new String[]{
        "I=" + in1.getAbsolutePath(),
        "O=" + out.getAbsolutePath(),
        "ACTION=CONCAT",
        "PADDING=" + padding,
        "BREAK_BANDS_AT_MULTIPLES_OF=" + breakBands});
    emit("input1", c, new String(Files.readAllBytes(in1.toPath())));
    emit("padding", c, String.valueOf(padding));
    emit("break_bands", c, String.valueOf(breakBands));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }

  public static void main(String[] x) throws Exception {
    dir = Files.createTempDirectory("iltp").toFile();
    dict = makeDict("chr1", 100);

    // PADDING=10: A(5-20) -> 1-30 (start clamped to 1); C(95-100) -> 85-100 (end clamped to 100).
    File pad = makeInput("pad",
        new Interval("chr1", 5, 20, false, "A"),
        new Interval("chr1", 95, 100, false, "C"));
    run("padding", pad, 10, 0);

    // BREAK_BANDS=10: A(5-25) -> A.1(5-9), A.2(10-19), A.3(20-25).
    File brk = makeInput("brk", new Interval("chr1", 5, 25, false, "A"));
    run("break_bands", brk, 0, 10);

    System.out.print(buf);
  }
}
