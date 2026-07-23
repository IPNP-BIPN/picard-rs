/*
 * Oracle dump for NonNFastaSize against Picard 3.4.0. Builds the required .dict and .fai (via
 * htsjdk), runs NonNFastaSize (whole-genome, no INTERVALS), and emits the input FASTA and the count.
 *   java -cp picard-fat.jar:. NonNFastaSizeDump | gzip > non_n_fasta_size.txt.gz
 */
import java.io.*; import java.nio.file.*;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
public class NonNFastaSizeDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static void run(String c, String fasta) throws Exception {
    File d = Files.createTempDirectory("nn").toFile();
    File in = new File(d, "in.fasta"); Files.write(in.toPath(), fasta.getBytes());
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
      "R=" + in.getAbsolutePath(), "O=" + new File(d, "in.dict").getAbsolutePath()});
    FastaSequenceIndexCreator.create(in.toPath(), true);
    File out = new File(d, "count.txt");
    int rc = new picard.reference.NonNFastaSize().instanceMain(new String[]{
      "I=" + in.getAbsolutePath(), "O=" + out.getAbsolutePath()});
    emit("input", c, fasta);
    emit("rc", c, String.valueOf(rc));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }
  public static void main(String[] x) throws Exception {
    // Uniform line width per sequence (faidx requirement); mixes case and Ns.
    run("mixed", ">chr1\nACGTacgtACGTNNNNacgt\nACGTACGT\n>chr2\nTTTTTTTTTTGGGGGGGGGGCCCCC\n");
    run("all_n", ">c\nNNNNNNNNNN\nnnnnn\n");
    System.out.print(buf);
  }
}
