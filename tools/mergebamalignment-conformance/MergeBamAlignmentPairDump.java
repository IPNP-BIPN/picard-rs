/*
 * `MergeBamAlignment`'s treatment of a PAIR, taken from the reference.
 *
 * A pair is not two records merged twice. The mate fields of each end are written from the other
 * end's alignment, the insert size is computed from both, the proper-pair flag is decided by the
 * merger rather than copied from the aligner, and where the two ends overlap the tool clips one of
 * them so the same base is not counted twice.
 *
 * Eight behaviours this is built to catch.
 *
 *   - THE MATE FIELDS COME FROM THE OTHER END, so a pair the aligner left with empty mate
 *     references comes out with them filled in;
 *   - THE INSERT SIZE IS COMPUTED, and it is signed by which end starts first;
 *   - THE PROPER-PAIR FLAG IS THE MERGER'S, not the aligner's, unless
 *     `--ALIGNER_PROPER_PAIR_FLAGS` says otherwise;
 *   - `--ADD_MATE_CIGAR` PUTS THE OTHER END'S CIGAR IN `MC`;
 *   - OVERLAPPING ENDS ARE CLIPPED, softly by default and hard on request, so the pair does not
 *     count the same base twice;
 *   - AN END WITH NO ALIGNMENT IS CARRIED as unmapped and still given its mate's coordinates;
 *   - A PAIR WHOSE ENDS LAND ON DIFFERENT CONTIGS is kept, and is not a proper pair;
 *   - AND `--PAIRED_RUN` IS INERT, being deprecated and read nowhere.
 *
 * Output:
 *
 *     record\t<case>\t<each merged record's data line>
 *     rc\t<case>\t<the exit status>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MergeBamAlignmentPairDump
 */
import java.io.*;
import java.nio.file.*;
import java.util.*;

public class MergeBamAlignmentPairDump {
  /** chr1 and chr2, forty bases each, `ACGT` repeating. */
  static final String REF = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
  static final StringBuilder buf = new StringBuilder();

  static void emit(final String kind, final String name, final String payload) {
    buf.append(kind).append('\t').append(name).append('\t')
       .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
       .append('\n');
  }

  static void run(final String name, final String unmapped, final String aligned,
                  final String... extra) throws Exception {
    final File dir = Files.createTempDirectory("mbap").toFile();
    final File reference = new File(dir, "ref.fasta");
    try (PrintStream out = new PrintStream(reference)) {
      out.println(">chr1");
      out.println(REF);
      out.println(">chr2");
      out.println(REF);
    }
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
        "R=" + reference.getAbsolutePath(), "O=" + new File(dir, "ref.dict").getAbsolutePath()});
    final File unmappedFile = new File(dir, "u.sam");
    Files.write(unmappedFile.toPath(), unmapped.getBytes());
    final File alignedFile = new File(dir, "a.sam");
    Files.write(alignedFile.toPath(), aligned.getBytes());
    final File output = new File(dir, "o.sam");

    final List<String> argv = new ArrayList<>(Arrays.asList(
        "UNMAPPED_BAM=" + unmappedFile.getAbsolutePath(),
        "ALIGNED_BAM=" + alignedFile.getAbsolutePath(),
        "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
        "OUTPUT=" + output.getAbsolutePath(),
        "SORT_ORDER=queryname"));
    argv.addAll(Arrays.asList(extra));

    final PrintStream originalOut = System.out;
    final PrintStream originalError = System.err;
    final ByteArrayOutputStream said = new ByteArrayOutputStream();
    final int code;
    try {
      System.setOut(new PrintStream(said, true, "UTF-8"));
      System.setErr(new PrintStream(said, true, "UTF-8"));
      code = new picard.sam.MergeBamAlignment().instanceMain(argv.toArray(new String[0]));
      System.setOut(originalOut);
      System.setErr(originalError);
    } catch (final Exception e) {
      System.setOut(originalOut);
      System.setErr(originalError);
      Throwable cause = e;
      while (cause.getCause() != null && cause.getCause() != cause) {
        cause = cause.getCause();
      }
      emit("error", name, cause.getClass().getName() + ":"
          + String.valueOf(cause.getMessage()).replace(dir.getAbsolutePath(), "<dir>"));
      return;
    } finally {
      System.setOut(originalOut);
      System.setErr(originalError);
    }
    emit("rc", name, String.valueOf(code));
    if (!output.exists()) {
      return;
    }
    for (final String line : new String(Files.readAllBytes(output.toPath())).split("\n", -1)) {
      if (line.startsWith("@") || line.isEmpty()) {
        continue;
      }
      emit("record", name, line);
    }
  }

  public static void main(final String[] args) throws Exception {
    final String unmappedHeader =
        "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n";
    final String alignedHeader =
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n@SQ\tSN:chr2\tLN:40\n"
        + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";
    // The unmapped pair: two ends of one template, neither placed.
    final String unmappedPair =
        "p\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        + "p\t141\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";

    // One aligned end, with the mate fields the aligner happened to write. They are WRONG on
    // purpose: the mate is on the other contig at position 39 with an insert size of 999, none of
    // which is true, and the merger is what decides what the record says in the end.
    final String first =
        "p\t65\tchr1\t%d\t60\t%s\tchr2\t39\t999\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";
    final String second =
        "p\t129\tchr1\t%d\t60\t%s\tchr2\t39\t-999\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";

    // The mate fields the aligner wrote are replaced by the ones the other end implies.
    run("mate-fields-from-the-other-end", unmappedHeader + unmappedPair,
        alignedHeader + String.format(first, 1, "8M") + String.format(second, 25, "8M"));

    // The same pair with the ends swapped, to see which end the insert size is signed for.
    run("the-second-end-first", unmappedHeader + unmappedPair,
        alignedHeader + String.format(first, 25, "8M") + String.format(second, 1, "8M"));

    // A pair the aligner called proper across two contigs, and the flag that keeps its word.
    final String acrossContigs =
        "p\t67\tchr1\t1\t60\t8M\tchr2\t1\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        + "p\t131\tchr2\t1\t60\t8M\tchr1\t1\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";
    run("across-two-contigs", unmappedHeader + unmappedPair, alignedHeader + acrossContigs);
    run("across-two-contigs-with-the-aligners-flags", unmappedHeader + unmappedPair,
        alignedHeader + acrossContigs, "ALIGNER_PROPER_PAIR_FLAGS=true");

    // The mate's cigar, and the run that does not add it.
    run("without-the-mate-cigar", unmappedHeader + unmappedPair,
        alignedHeader + String.format(first, 1, "8M") + String.format(second, 25, "8M"),
        "ADD_MATE_CIGAR=false");

    // Two ends that overlap, clipped softly, not at all, and hard.
    final String overlapping = alignedHeader
        + "p\t65\tchr1\t5\t60\t8M\tchr1\t3\t-10\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        + "p\t145\tchr1\t3\t60\t8M\tchr1\t5\t10\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";
    run("overlapping-ends", unmappedHeader + unmappedPair, overlapping);
    run("overlapping-ends-unclipped", unmappedHeader + unmappedPair, overlapping,
        "CLIP_OVERLAPPING_READS=false");
    run("overlapping-ends-hard-clipped", unmappedHeader + unmappedPair, overlapping,
        "HARD_CLIP_OVERLAPPING_READS=true");

    // One end with no alignment at all: the other end says so in its flags.
    run("one-end-unaligned", unmappedHeader + unmappedPair,
        alignedHeader
        + "p\t73\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"
        + "p\t141\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n");

    // A soft-clipped end, to see what the insert size is measured between.
    run("a-soft-clipped-end", unmappedHeader + unmappedPair,
        alignedHeader + String.format(first, 1, "2S6M") + String.format(second, 25, "8M"));

    // And the deprecated argument, which the tool reads nowhere.
    run("the-deprecated-paired-run", unmappedHeader + unmappedPair,
        alignedHeader + String.format(first, 1, "8M") + String.format(second, 25, "8M"),
        "PAIRED_RUN=false");

    System.out.print(buf);
  }
}
