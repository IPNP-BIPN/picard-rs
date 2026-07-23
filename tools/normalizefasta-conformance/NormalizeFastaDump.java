/*
 * Oracle dump for NormalizeFasta against Picard 3.4.0. Emits input FASTA, the options, and the
 * normalized output, per case. Output is plain text so every row compares raw.
 *   java -cp picard-fat.jar:. NormalizeFastaDump | gzip > normalize_fasta.txt.gz
 */
import java.io.*; import java.nio.file.*;
public class NormalizeFastaDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static void run(String c, String fasta, int lineLength, boolean truncate) throws Exception {
    File d = Files.createTempDirectory("nf").toFile();
    File in = new File(d, "in.fasta"); Files.write(in.toPath(), fasta.getBytes());
    File out = new File(d, "out.fasta");
    int rc = new picard.reference.NormalizeFasta().instanceMain(new String[]{
      "I=" + in.getAbsolutePath(), "O=" + out.getAbsolutePath(),
      "LINE_LENGTH=" + lineLength, "TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE=" + truncate});
    emit("input", c, fasta);
    emit("line_length", c, String.valueOf(lineLength));
    emit("truncate", c, String.valueOf(truncate));
    emit("rc", c, String.valueOf(rc));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }
  public static void main(String[] x) throws Exception {
    String fasta = ">chr1 a description here\nACGTacgtACGTNNNNacgt\nACGTACGT\n>chr2\nTTTTTTTTTTGGGGGGGGGGCCCCC\n>empty\n";
    run("default", fasta, 100, false);
    run("wrap10", fasta, 10, false);
    run("truncate", fasta, 100, true);
    System.out.print(buf);
  }
}
