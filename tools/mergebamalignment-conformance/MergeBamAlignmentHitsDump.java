/*
 * `MergeBamAlignment`'s choice of a PRIMARY alignment, taken from the reference.
 *
 * An aligner may report a read several times. The merger has to choose which of those alignments
 * is the primary one, mark the rest as secondary, and decide whether to carry them at all. Four
 * strategies make that choice, and they do not agree.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE ALIGNER'S OWN PRIMARY IS KEPT where it named one, whatever the strategy;
 *   - `BestMapq` PICKS THE HIGHEST MAPPING QUALITY where the aligner named none;
 *   - `EarliestFragment` PICKS THE ALIGNMENT THAT MAPS THE EARLIEST BASE OF THE READ, which is not
 *     the earliest position on the reference;
 *   - `MostDistant` PICKS THE PAIR WITH THE LARGEST INSERT, and falls back to the mapping quality
 *     where every pairing would be chimeric;
 *   - `--INCLUDE_SECONDARY_ALIGNMENTS` DECIDES WHETHER THE REST ARE WRITTEN, and the primary is
 *     written either way;
 *   - AND `EarliestFragment` REFUSES A PAIRED RUN, which is a refusal rather than a choice.
 *
 * Output:
 *
 *     record\t<case>\t<each merged record's data line>
 *     rc\t<case>\t<the exit status>
 *     refusal\t<case>\t<the reason a refused command line printed>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MergeBamAlignmentHitsDump
 */
import java.io.*;
import java.nio.file.*;
import java.util.*;

public class MergeBamAlignmentHitsDump {
  static final int LENGTH = 200;
  static final StringBuilder buf = new StringBuilder();

  static String reference() {
    final StringBuilder bases = new StringBuilder();
    for (int index = 0; index < LENGTH; index++) {
      bases.append("ACGT".charAt(index % 4));
    }
    return bases.toString();
  }

  static String window(final int start, final int length) {
    return reference().substring(start - 1, start - 1 + length);
  }

  static void emit(final String kind, final String name, final String payload) {
    buf.append(kind).append('\t').append(name).append('\t')
       .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
       .append('\n');
  }

  static void run(final String name, final String unmapped, final String aligned,
                  final String... extra) throws Exception {
    final File dir = Files.createTempDirectory("mbah2").toFile();
    final File reference = new File(dir, "ref.fasta");
    try (PrintStream out = new PrintStream(reference)) {
      out.println(">chr1");
      out.println(reference());
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
    if (code != 0) {
      final List<String> lines = new ArrayList<>();
      for (final String line : said.toString("UTF-8").split("\n", -1)) {
        if (!line.trim().isEmpty()) {
          lines.add(line.replace(dir.getAbsolutePath(), "<dir>"));
        }
      }
      final List<String> reasons = new ArrayList<>();
      for (final String line : lines) {
        if (line.startsWith("ERROR")) {
          reasons.add(line);
        }
      }
      if (reasons.isEmpty() && !lines.isEmpty()) {
        reasons.add(lines.get(lines.size() - 1));
      }
      emit("refusal", name, String.join("\n", reasons));
      return;
    }
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
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:" + LENGTH + "\n"
        + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";
    final String bases = window(41, 20);
    final String qualities = "I".repeat(20);
    final String unmappedRead =
        unmappedHeader + "r\t4\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n";

    // Three alignments of one read, none of them named primary by the aligner: the first has the
    // best mapping quality, the second maps the earliest base of the read, the third neither.
    // The best mapping quality is on the alignment that does NOT map the read's first base, so the
    // two strategies part company: one wants the quality, the other the earliest base.
    final String threeHits = alignedHeader
        + "r\t256\tchr1\t121\t60\t5S15M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n"
        + "r\t256\tchr1\t61\t30\t20M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n"
        + "r\t256\tchr1\t41\t10\t10S10M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n";
    run("three-hits-best-mapq", unmappedRead, threeHits);
    run("three-hits-earliest-fragment", unmappedRead, threeHits,
        "PRIMARY_ALIGNMENT_STRATEGY=EarliestFragment");
    run("three-hits-best-end-mapq", unmappedRead, threeHits,
        "PRIMARY_ALIGNMENT_STRATEGY=BestEndMapq");
    run("three-hits-most-distant", unmappedRead, threeHits,
        "PRIMARY_ALIGNMENT_STRATEGY=MostDistant");
    run("three-hits-without-the-secondaries", unmappedRead, threeHits,
        "INCLUDE_SECONDARY_ALIGNMENTS=false");

    // The same three with the aligner naming one primary, which every strategy keeps.
    final String namedPrimary = alignedHeader
        + "r\t256\tchr1\t121\t60\t5S15M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n"
        + "r\t0\tchr1\t61\t30\t20M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n"
        + "r\t256\tchr1\t41\t10\t10S10M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n";
    run("the-aligners-own-primary", unmappedRead, namedPrimary);
    run("the-aligners-own-primary-under-another-strategy", unmappedRead, namedPrimary,
        "PRIMARY_ALIGNMENT_STRATEGY=EarliestFragment");

    // A pair with two alignments each, and the strategies that read a pair differently.
    final String unmappedPair =
        unmappedHeader + "p\t77\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n"
        + "p\t141\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t" + qualities + "\tRG:Z:rg1\n";
    final String pairedHits = alignedHeader
        + "p\t321\tchr1\t41\t60\t20M\tchr1\t81\t60\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\tHI:i:1\n"
        + "p\t385\tchr1\t81\t60\t20M\tchr1\t41\t-60\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\tHI:i:1\n"
        + "p\t321\tchr1\t41\t30\t20M\tchr1\t161\t140\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\tHI:i:2\n"
        + "p\t385\tchr1\t161\t30\t20M\tchr1\t41\t-140\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\tHI:i:2\n";
    run("a-pair-with-two-hits", unmappedPair, pairedHits);
    run("a-pair-with-two-hits-most-distant", unmappedPair, pairedHits,
        "PRIMARY_ALIGNMENT_STRATEGY=MostDistant");
    // And the strategy that will not read a pair at all.
    run("a-pair-under-earliest-fragment", unmappedPair, pairedHits,
        "PRIMARY_ALIGNMENT_STRATEGY=EarliestFragment");

    System.out.print(buf);
  }
}
