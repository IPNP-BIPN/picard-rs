/*
 * Oracle dump for LiftOverIntervalList against Picard 3.4.0. For each case it builds a "from" .dict
 * and a "to" .dict (via CreateSequenceDictionary), a UCSC chain file, and an input interval_list on
 * the "from" build, then runs LiftOverIntervalList and emits the input interval_list, the "to" .dict
 * used as SEQUENCE_DICTIONARY, the chain text, MIN_LIFTOVER_PCT, the return code, and the output
 * interval_list.
 *   java -cp picard-fat.jar:. LiftOverIntervalListDump | gzip -n > liftover_interval_list.txt.gz
 */
import java.io.*;
import java.nio.file.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class LiftOverIntervalListDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")).append('\n');
  }

  static File dir;

  /** Build a single-line FASTA of the given length and return its CreateSequenceDictionary .dict file. */
  static File makeDict(String contig, int length) throws Exception {
    StringBuilder seq = new StringBuilder();
    for (int i = 0; i < length; i++) seq.append("ACGT".charAt(i % 4));
    File fasta = new File(dir, contig + ".fasta");
    Files.write(fasta.toPath(), (">" + contig + "\n" + seq + "\n").getBytes());
    File dict = new File(dir, contig + ".dict");
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
        "R=" + fasta.getAbsolutePath(), "O=" + dict.getAbsolutePath()});
    return dict;
  }

  /** Build an input interval_list on the "from" dictionary with the given intervals. */
  static File makeInput(String name, File fromDict, Interval... intervals) throws Exception {
    SAMSequenceDictionary sd = SAMSequenceDictionaryExtractor.extractDictionary(fromDict.toPath());
    IntervalList list = new IntervalList(new SAMFileHeader(sd));
    for (Interval iv : intervals) list.add(iv);
    File f = new File(dir, name + ".interval_list");
    list.write(f);
    return f;
  }

  static void run(String c, File input, File toDict, String chainText, double minPct) throws Exception {
    File chain = new File(dir, c + ".chain");
    Files.write(chain.toPath(), chainText.getBytes());
    File out = new File(dir, c + ".out.interval_list");
    int rc = new picard.util.LiftOverIntervalList().instanceMain(new String[]{
        "I=" + input.getAbsolutePath(),
        "O=" + out.getAbsolutePath(),
        "SD=" + toDict.getAbsolutePath(),
        "CHAIN=" + chain.getAbsolutePath(),
        "MIN_LIFTOVER_PCT=" + minPct});
    emit("input", c, new String(Files.readAllBytes(input.toPath())));
    emit("sd", c, new String(Files.readAllBytes(toDict.toPath())));
    emit("chain", c, chainText);
    emit("min_pct", c, String.valueOf(minPct));
    emit("rc", c, String.valueOf(rc));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }

  public static void main(String[] x) throws Exception {
    dir = Files.createTempDirectory("lo").toFile();
    File chr1 = makeDict("chr1", 100);   // "from" build
    File chr2 = makeDict("chr2", 200);   // "to" build (used for every case)

    // same_strand: chr1[0,100) -> chr2[50,150) '+', single block. Two intervals fully inside.
    File in1 = makeInput("same",
        chr1,
        new Interval("chr1", 11, 20, false, "A"),
        new Interval("chr1", 31, 40, false, "B"));
    run("same_strand", in1, chr2,
        "chain 1000 chr1 100 + 0 100 chr2 200 + 50 150 1\n100\n\n", 0.95);

    // neg_strand: chr1[0,100) -> chr2[0,100) '-' (toSequenceSize 200). One interval, flipped.
    File in2 = makeInput("neg",
        chr1,
        new Interval("chr1", 11, 20, false, "A"));
    run("neg_strand", in2, chr2,
        "chain 1000 chr1 100 + 0 100 chr2 200 - 0 100 2\n100\n\n", 0.95);

    // with_reject: chr1[0,50) -> chr2[0,50) '+'. One interval lifts, one is outside and rejected.
    File in3 = makeInput("reject",
        chr1,
        new Interval("chr1", 11, 20, false, "A"),
        new Interval("chr1", 71, 80, false, "B"));
    run("with_reject", in3, chr2,
        "chain 1000 chr1 100 + 0 50 chr2 200 + 0 50 3\n50\n\n", 0.95);

    System.out.print(buf);
  }
}
