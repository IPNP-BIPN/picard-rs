/*
 * `MergeBamAlignment`'s merged HEADER, taken from the reference.
 *
 * The records are measured next door. This is the header they arrive under, which is built from
 * three files rather than one: the sequences come from the reference dictionary, the read groups
 * from the UNMAPPED bam, and the programs from the ALIGNED one, merged through
 * `SamFileHeaderMerger`.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE SEQUENCES COME FROM THE DICTIONARY, so the `@SQ` lines carry the dictionary's `M5` and
 *     `UR` and not the aligner's bare `SN`/`LN`;
 *   - THE `UR` IS THE PATH THE COMMAND LINE GAVE, canonicalised or not: the same file named the
 *     long way round is what says which;
 *   - THE READ GROUPS COME FROM THE UNMAPPED BAM, so an aligner that rewrote them is ignored, and
 *     an unmapped bam with none leaves the output with none;
 *   - THE PROGRAMS COME FROM THE ALIGNED BAM, chained by `PP` and followed by the tool's own;
 *   - TWO ALIGNED FILES DECLARING THE SAME PROGRAM ID ARE MERGED rather than refused, and the
 *     second one's id is rewritten;
 *   - `--PROGRAM_RECORD_ID` AND ITS FELLOWS ADD A PROGRAM of the caller's own;
 *   - THE FOUR `--PROGRAM_GROUP` ARGUMENTS COME TOGETHER OR NOT AT ALL, and the inputs' own
 *     comments are kept or dropped;
 *   - THE OUTPUT'S SORT ORDER IS THE ONE ASKED FOR, and it is on the `@HD`;
 *   - AND AN ALIGNED HEADER THAT DISAGREES WITH THE DICTIONARY IS REFUSED, in the reference's own
 *     words.
 *
 * Output:
 *
 *     header\t<case>\t<the output's header lines, escaped>
 *     rc\t<case>\t<the exit status>
 *     refusal\t<case>\t<the `ERROR:` lines a refused command line printed>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MergeBamAlignmentHeaderDump
 */
import java.io.*;
import java.nio.file.*;
import java.util.*;

public class MergeBamAlignmentHeaderDump {
  static final String REF = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // chr1, 40 bp
  static final StringBuilder buf = new StringBuilder();

  static void emit(final String kind, final String name, final String payload) {
    buf.append(kind).append('\t').append(name).append('\t')
       .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
       .append('\n');
  }

  /** One run: the two inputs, and whatever else the case asks for. */
  static void run(final String name, final String unmapped, final String aligned,
                  final String... extra) throws Exception {
    runFrom(name, unmapped, new String[]{aligned}, false, extra);
  }

  /**
   * One run, with as many aligned files as the case wants and the reference named by a relative
   * path where it asks for one.
   */
  static void runFrom(final String name, final String unmapped, final String[] aligned,
                      final boolean relativeReference, final String... extra) throws Exception {
    final File dir = Files.createTempDirectory("mbah").toFile();
    final File reference = new File(dir, "ref.fasta");
    try (PrintStream out = new PrintStream(reference)) {
      out.println(">chr1");
      out.println(REF);
    }
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
        "R=" + reference.getAbsolutePath(), "O=" + new File(dir, "ref.dict").getAbsolutePath()});

    final File unmappedFile = new File(dir, "u.sam");
    Files.write(unmappedFile.toPath(), unmapped.getBytes());
    final List<String> argv = new ArrayList<>();
    argv.add("UNMAPPED_BAM=" + unmappedFile.getAbsolutePath());
    for (int index = 0; index < aligned.length; index++) {
      final File file = new File(dir, "a" + index + ".sam");
      Files.write(file.toPath(), aligned[index].getBytes());
      argv.add("ALIGNED_BAM=" + file.getAbsolutePath());
    }
    // The same file named the long way round is what says whether the `UR` is canonicalised: the
    // path is absolute either way, and only one of the two is the shortest spelling of it.
    argv.add("REFERENCE_SEQUENCE=" + (relativeReference
        ? theLongWayRound(dir, reference) : reference.getAbsolutePath()));
    final File output = new File(dir, "o.sam");
    argv.add("OUTPUT=" + output.getAbsolutePath());
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
          + String.valueOf(cause.getMessage())
              .replace(dir.getAbsolutePath(), "<dir>"));
      return;
    } finally {
      System.setOut(originalOut);
      System.setErr(originalError);
    }
    emit("rc", name, String.valueOf(code));
    if (code != 0) {
      // A refused command line prints its reason rather than throwing it, and the stream carries
      // the whole usage in front of it. The parser's own refusals are `ERROR:` lines; the tool's
      // own are printed bare, after the usage, so the last line is the sentence it refused with.
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
    final List<String> header = new ArrayList<>();
    for (final String line : new String(Files.readAllBytes(output.toPath())).split("\n", -1)) {
      if (!line.startsWith("@")) {
        continue;
      }
      // The temp directory is the run's own, and the version of the tool is the container's: what
      // is measured is which lines are there and what they carry.
      header.add(line.replace(dir.getAbsolutePath(), "<dir>"));
    }
    emit("header", name, String.join("\n", header));
  }

  /** The same file, by a path that goes through a directory and back out of it. */
  static String theLongWayRound(final File dir, final File reference) throws IOException {
    final File detour = new File(dir, "sub");
    Files.createDirectories(detour.toPath());
    return detour.getAbsolutePath() + "/../" + reference.getName();
  }

  public static void main(final String[] args) throws Exception {
    final String unmappedHeader =
        "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n";
    final String alignedHeader =
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n"
        + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";
    final String unmappedRead =
        "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";
    final String alignedRead =
        "r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n";

    run("the-whole-header", unmappedHeader + unmappedRead, alignedHeader + alignedRead);

    // The reference named the long way round, which is what says whether the UR is canonicalised.
    runFrom("a-reference-named-the-long-way-round", unmappedHeader + unmappedRead,
        new String[]{alignedHeader + alignedRead}, true);

    // A read group the aligner rewrote, and an unmapped bam that declares none at all.
    run("a-read-group-the-aligner-rewrote", unmappedHeader + unmappedRead,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n"
        + "@RG\tID:rg1\tSM:other\tLB:otherlib\tPL:ILLUMINA\n"
        + "@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n" + alignedRead);
    run("no-read-group-at-all",
        "@HD\tVN:1.6\tSO:queryname\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n"
        + "@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n"
        + "r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n");

    // Two programs, chained, and the tool's own after them.
    run("a-chain-of-programs", unmappedHeader + unmappedRead,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n"
        + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n"
        + "@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n"
        + "@PG\tID:samtools\tPN:samtools\tVN:1.9\tPP:bwa\tCL:samtools view\n" + alignedRead);

    // Two aligned files whose programs collide.
    runFrom("two-aligned-files-with-the-same-program-id", unmappedHeader + unmappedRead,
        new String[]{
            alignedHeader + alignedRead,
            "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n"
            + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n"
            + "@PG\tID:bwa\tPN:bwa\tVN:2.0\tCL:bwa mem -a\n"
            + "r2\t0\tchr1\t9\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n"},
        false);

    // A program of the caller's own, and a comment.
    run("a-program-from-the-command-line", unmappedHeader + unmappedRead,
        alignedHeader + alignedRead,
        "PROGRAM_RECORD_ID=mine", "PROGRAM_GROUP_NAME=miner",
        "PROGRAM_GROUP_VERSION=3.0", "PROGRAM_GROUP_COMMAND_LINE=mine --do-it");
    // The four arguments that make a program record come together or not at all.
    run("half-a-program-from-the-command-line", unmappedHeader + unmappedRead,
        alignedHeader + alignedRead, "PROGRAM_RECORD_ID=mine");
    run("comments-in-the-inputs",
        "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n"
        + "@CO\tfrom the unmapped\n" + unmappedRead,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:40\n"
        + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n"
        + "@CO\tfrom the aligned\n" + alignedRead);

    // The sort order the output is written in.
    run("a-queryname-output", unmappedHeader + unmappedRead, alignedHeader + alignedRead,
        "SORT_ORDER=queryname");
    run("an-unsorted-output", unmappedHeader + unmappedRead, alignedHeader + alignedRead,
        "SORT_ORDER=unsorted");

    // An aligned header that disagrees with the dictionary.
    run("an-aligned-header-that-disagrees", unmappedHeader + unmappedRead,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:41\n"
        + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n"
        + alignedRead);
    run("an-aligned-header-with-another-contig", unmappedHeader + unmappedRead,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr2\tLN:40\n"
        + "@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n"
        + "r1\t0\tchr2\t1\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rg1\n");

    System.out.print(buf);
  }
}
