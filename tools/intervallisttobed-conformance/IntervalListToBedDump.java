/*
 * Oracle dump for IntervalListToBed against Picard 3.4.0. Builds a reference .dict, synthesizes an
 * interval list whose header carries that dict, runs IntervalListToBed with the given SCORE/SORT, and
 * emits the interval list, SCORE, SORT, and the BED output.
 *   java -cp picard-fat.jar:. IntervalListToBedDump | gzip -n > interval_list_to_bed.txt.gz
 */
import java.io.*; import java.nio.file.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
public class IntervalListToBedDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static void run(String c, int score, boolean sort, Interval[] ivs) throws Exception {
    File d = Files.createTempDirectory("ib").toFile();
    // A minimal reference so CreateSequenceDictionary gives us a real dict for the list header.
    String fasta = ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nAAAACCCCGGGGTTTTACGT\n>chr10\nGGGGCCCCAAAATTTTACGT\n";
    File in = new File(d, "in.fasta"); Files.write(in.toPath(), fasta.getBytes());
    File dict = new File(d, "in.dict");
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
      "R=" + in.getAbsolutePath(), "O=" + dict.getAbsolutePath()});

    SAMSequenceDictionary sd = SAMSequenceDictionaryExtractor.extractDictionary(dict.toPath());
    SAMFileHeader header = new SAMFileHeader(sd);
    IntervalList il = new IntervalList(header);
    for (Interval iv : ivs) il.add(iv);
    File ilFile = new File(d, "in.interval_list");
    il.write(ilFile);

    File out = new File(d, "out.bed");
    int rc = new picard.util.IntervalListToBed().instanceMain(new String[]{
      "INPUT=" + ilFile.getAbsolutePath(),
      "OUTPUT=" + out.getAbsolutePath(),
      "SCORE=" + score,
      "SORT=" + sort});

    emit("interval_list", c, new String(Files.readAllBytes(ilFile.toPath())));
    emit("score", c, String.valueOf(score));
    emit("sort", c, String.valueOf(sort));
    emit("rc", c, String.valueOf(rc));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }
  public static void main(String[] x) throws Exception {
    // Default SCORE/SORT: dictionary order (chr1, chr2, chr10) beats string order, with a strand tie.
    run("sorted_default", 500, true, new Interval[]{
      new Interval("chr10", 3, 9,  false, "j"),
      new Interval("chr2",  5, 8,  true,  "b"),
      new Interval("chr1",  1, 4,  false, "a"),
      new Interval("chr1",  1, 4,  true,  "a_rev"),
    });
    // SORT off keeps file order; a non-default SCORE too.
    run("unsorted", 1000, false, new Interval[]{
      new Interval("chr10", 3, 9, false, "j"),
      new Interval("chr1",  1, 4, false, "a"),
    });
    System.out.print(buf);
  }
}
