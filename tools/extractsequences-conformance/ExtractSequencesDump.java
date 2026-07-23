/*
 * Oracle dump for ExtractSequences against Picard 3.4.0. Builds the reference .dict and .fai (via
 * htsjdk), builds an interval list whose header matches that dict, runs ExtractSequences, and emits
 * the interval list, the reference FASTA, the LINE_LENGTH, and the extracted FASTA.
 *   java -cp picard-fat.jar:. ExtractSequencesDump | gzip -n > extract_sequences.txt.gz
 */
import java.io.*; import java.nio.file.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
public class ExtractSequencesDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static void run(String c, String fasta, int lineLength, Interval[] ivs) throws Exception {
    File d = Files.createTempDirectory("es").toFile();
    File in = new File(d, "in.fasta"); Files.write(in.toPath(), fasta.getBytes());
    File dict = new File(d, "in.dict");
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
      "R=" + in.getAbsolutePath(), "O=" + dict.getAbsolutePath()});
    FastaSequenceIndexCreator.create(in.toPath(), true);

    SAMSequenceDictionary sd = SAMSequenceDictionaryExtractor.extractDictionary(dict.toPath());
    SAMFileHeader header = new SAMFileHeader(sd);
    IntervalList il = new IntervalList(header);
    for (Interval iv : ivs) il.add(iv);
    File ilFile = new File(d, "in.interval_list");
    il.write(ilFile);

    File out = new File(d, "out.fasta");
    int rc = new picard.reference.ExtractSequences().instanceMain(new String[]{
      "INTERVAL_LIST=" + ilFile.getAbsolutePath(),
      "R=" + in.getAbsolutePath(),
      "O=" + out.getAbsolutePath(),
      "LINE_LENGTH=" + lineLength});

    emit("interval_list", c, new String(Files.readAllBytes(ilFile.toPath())));
    emit("reference", c, fasta);
    emit("line_length", c, String.valueOf(lineLength));
    emit("rc", c, String.valueOf(rc));
    emit("output", c, new String(Files.readAllBytes(out.toPath())));
  }
  public static void main(String[] x) throws Exception {
    // chr1 40 bases (two 20-base lines), chr2 12 bases (one line): uniform width per contig for faidx.
    String fasta = ">chr1\nACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGT\n>chr2\nAAAACCCCGGGG\n";
    // Positive strand spanning a wrap, a negative strand (reverse complement), a whole small contig.
    run("mixed", fasta, 8, new Interval[]{
      new Interval("chr1", 1, 16, false, "fwd16"),
      new Interval("chr1", 5, 12, true,  "rc8"),
      new Interval("chr2", 1, 12, false, "chr2full"),
    });
    // A length that is an exact multiple of LINE_LENGTH must not gain a trailing blank line.
    run("exact_multiple", fasta, 4, new Interval[]{
      new Interval("chr1", 1, 8, false, "eight"),
    });
    System.out.print(buf);
  }
}
