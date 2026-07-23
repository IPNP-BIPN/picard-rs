/*
 * Oracle dump for BedToIntervalList against Picard 3.4.0. Builds a reference .dict, runs
 * BedToIntervalList on a BED file with the given SORT/UNIQUE/KEEP_LENGTH_ZERO_INTERVALS flags, and
 * emits the .dict, the BED, the flags, and the resulting interval list.
 *   java -cp picard-fat.jar:. BedToIntervalListDump | gzip -n > bed_to_interval_list.txt.gz
 */
import java.io.*; import java.nio.file.*;
public class BedToIntervalListDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static String FASTA =
    ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nAAAACCCCGGGGTTTTACGT\n>chr10\nGGGGCCCCAAAATTTTACGT\n";
  static File dir, dict;
  static void setup() throws Exception {
    dir = Files.createTempDirectory("b2i").toFile();
    File in = new File(dir, "in.fasta"); Files.write(in.toPath(), FASTA.getBytes());
    dict = new File(dir, "in.dict");
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
      "R=" + in.getAbsolutePath(), "O=" + dict.getAbsolutePath()});
  }
  static void run(String c, String bed, boolean sort, boolean unique, boolean keepZero) throws Exception {
    File bedFile = new File(dir, c + ".bed"); Files.write(bedFile.toPath(), bed.getBytes());
    File out = new File(dir, c + ".interval_list");
    int rc = new picard.util.BedToIntervalList().instanceMain(new String[]{
      "INPUT=" + bedFile.getAbsolutePath(),
      "SEQUENCE_DICTIONARY=" + dict.getAbsolutePath(),
      "OUTPUT=" + out.getAbsolutePath(),
      "SORT=" + sort, "UNIQUE=" + unique, "KEEP_LENGTH_ZERO_INTERVALS=" + keepZero});
    emit("dict", c, new String(Files.readAllBytes(dict.toPath())));
    emit("bed", c, bed);
    emit("sort", c, String.valueOf(sort));
    emit("unique", c, String.valueOf(unique));
    emit("keep_zero", c, String.valueOf(keepZero));
    emit("rc", c, String.valueOf(rc));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }
  public static void main(String[] x) throws Exception {
    setup();
    // Default flags: sort by dictionary order (chr1, chr2, chr10), empty name, a - strand, a length
    // zero feature that must be dropped.
    run("default",
        "chr10\t2\t9\tj\t0\t+\n" +
        "chr2\t0\t5\t\t0\t-\n" +
        "chr1\t0\t4\ta\t0\t+\n" +
        "chr1\t5\t5\tzero\t0\t+\n", true, false, false);
    // UNIQUE merges the two overlapping chr1 features and concatenates their names.
    run("unique",
        "chr1\t0\t10\ta\t0\t+\n" +
        "chr1\t5\t15\tb\t0\t+\n", true, true, false);
    // KEEP_LENGTH_ZERO_INTERVALS keeps the zero-length feature at chr1:1.
    run("keep_zero",
        "chr1\t0\t0\tz\t0\t+\n" +
        "chr2\t3\t7\ty\t0\t+\n", true, false, true);
    System.out.print(buf);
  }
}
