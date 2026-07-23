/*
 * Oracle dump harness for SamToFastqWithTags (the per-tag-group FASTQ output) conformance in
 * picard-rs. Runs the tool for a few cases and emits an escaped TSV to stdout:
 *   input\t<case>\t<sam text>
 *   out\t<case>:<filename>\t<file contents>
 * The base read FASTQ (out.fastq) is the already-ported SamToFastq path and is not dumped here.
 *
 *   java -cp picard-fat.jar:. SamToFastqWithTagsDump | gzip > sam_to_fastq_with_tags.txt.gz
 */
import java.io.*; import java.nio.file.*; import java.util.*;
public class SamToFastqWithTagsDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")).append('\n');
  }
  static void run(String c, String sam, String... extraArgs) throws Exception { run(c, sam, false, extraArgs); }
  static void run(String c, String sam, boolean paired, String... extraArgs) throws Exception {
    File d = Files.createTempDirectory("stfwt").toFile();
    File s = new File(d, "in.sam"); Files.write(s.toPath(), sam.getBytes());
    File o = new File(d, "out.fastq");
    List<String> args = new ArrayList<>();
    args.add("I=" + s.getAbsolutePath());
    args.add("FASTQ=" + o.getAbsolutePath());
    if (paired) args.add("SECOND_END_FASTQ=" + new File(d, "out2.fastq").getAbsolutePath());
    args.addAll(Arrays.asList(extraArgs));
    int rc = new picard.sam.SamToFastqWithTags().instanceMain(args.toArray(new String[0]));
    emit("input", c, sam);
    emit("rc", c, String.valueOf(rc));
    File[] files = d.listFiles();
    Arrays.sort(files, Comparator.comparing(File::getName));
    for (File f : files) {
      String n = f.getName();
      if (n.endsWith(".fastq") && !n.equals("out.fastq") && !n.equals("out2.fastq")) {
        emit("out", c + ":" + n, new String(Files.readAllBytes(f.toPath())));
      }
    }
  }
  public static void main(String[] x) throws Exception {
    String sam = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n" +
      "r1\t4\t*\t0\t0\t*\t*\t0\t0\tAACCGGTT\tIIIIIIII\tCR:Z:ACGT\tCY:Z:FFFF\tCB:Z:TTGG\tUR:Z:CC\tUY:Z:!!\n" +
      "r2\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGCCCC\tJJJJJJJJ\tCR:Z:TGCA\tCY:Z:####\tCB:Z:AAAA\tUR:Z:GG\tUY:Z:@@\n";
    // A single tag group with quality, and a two-tag group with concatenated quality.
    run("two_groups", sam,
        "SEQUENCE_TAG_GROUP=CR", "QUALITY_TAG_GROUP=CY",
        "SEQUENCE_TAG_GROUP=CB,UR", "QUALITY_TAG_GROUP=CY,UY");
    // A group with no quality group at all: quality is filled with '~'.
    run("no_quality", sam, "SEQUENCE_TAG_GROUP=CB,UR");
    // Paired input (SECOND_END_FASTQ): the tag FASTQ carries both ends (name/1, name/2) in one file.
    String pairedSam = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n" +
      "p1\t77\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\tCR:Z:ACGT\tCY:Z:FFFF\n" +
      "p1\t141\t*\t0\t0\t*\t*\t0\t0\tCCCC\tJJJJ\tCR:Z:TTTT\tCY:Z:####\n";
    run("paired", pairedSam, true, "SEQUENCE_TAG_GROUP=CR", "QUALITY_TAG_GROUP=CY");
    System.out.print(buf);
  }
}
