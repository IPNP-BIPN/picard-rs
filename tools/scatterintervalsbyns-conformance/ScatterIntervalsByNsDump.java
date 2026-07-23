/*
 * Oracle dump for ScatterIntervalsByNs against Picard 3.4.0. Builds the reference .dict and .fai (via
 * htsjdk), runs ScatterIntervalsByNs with the given OUTPUT_TYPE/MAX_TO_MERGE, and emits the reference
 * FASTA, the .dict, the options, and the interval list.
 *   java -cp picard-fat.jar:. ScatterIntervalsByNsDump | gzip -n > scatter_intervals_by_ns.txt.gz
 */
import java.io.*; import java.nio.file.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
public class ScatterIntervalsByNsDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static String FASTA, DICT;
  static File dir, ref, dict;
  static void setup(String fasta) throws Exception {
    FASTA = fasta;
    dir = Files.createTempDirectory("sn").toFile();
    ref = new File(dir, "ref.fasta"); Files.write(ref.toPath(), fasta.getBytes());
    dict = new File(dir, "ref.dict");
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
      "R=" + ref.getAbsolutePath(), "O=" + dict.getAbsolutePath()});
    FastaSequenceIndexCreator.create(ref.toPath(), true);
    DICT = new String(Files.readAllBytes(dict.toPath()));
  }
  static void run(String c, String outputType, int maxToMerge) throws Exception {
    File out = new File(dir, c + ".interval_list");
    int rc = new picard.util.ScatterIntervalsByNs().instanceMain(new String[]{
      "R=" + ref.getAbsolutePath(),
      "OUTPUT=" + out.getAbsolutePath(),
      "OUTPUT_TYPE=" + outputType,
      "MAX_TO_MERGE=" + maxToMerge});
    emit("reference", c, FASTA);
    emit("dict", c, DICT);
    emit("output_type", c, outputType);
    emit("max_to_merge", c, String.valueOf(maxToMerge));
    emit("rc", c, String.valueOf(rc));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }
  public static void main(String[] x) throws Exception {
    // chr1 has an ACGT / N(4) / ACGT / N(1) / ACGT pattern (uniform 20-base lines for faidx); chr2 is
    // all called. With MAX_TO_MERGE=1 the 4-N run splits and the 1-N run folds its flanks.
    setup(">chr1\nACGTACGTACGTACGTACGT\n" +
          "NNNNACGTACGTNACGTACG\n" +
          ">chr2\nACGTACGTACGTACGTACGT\n");
    run("both_merge1", "BOTH", 1);
    run("n_only", "N", 1);
    run("acgt_merge5", "ACGT", 5);
    System.out.print(buf);
  }
}
