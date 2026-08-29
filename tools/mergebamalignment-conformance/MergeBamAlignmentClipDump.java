/*
 * `MergeBamAlignment`'s clipping and unmapping, taken from the reference.
 *
 * Three things can shorten or remove an alignment on the way through: the contig it runs off the
 * end of, the adapter the unmapped bam marked in it, and the contamination filter that decides an
 * alignment is too short to believe.
 *
 * Seven behaviours this is built to catch.
 *
 *   - AN ALIGNMENT THAT RUNS OFF THE END OF ITS CONTIG IS SOFT-CLIPPED to fit, and the record is
 *     kept rather than dropped;
 *   - THE CLIP IS COUNTED FROM THE READ rather than from the reference, so a read that already
 *     ends in a soft clip is not clipped twice;
 *   - `--CLIP_ADAPTERS` CLIPS WHAT `XT` MARKED, from that base to the three-prime end;
 *   - AND IT IS ON BY DEFAULT, so turning it off is what keeps the adapter aligned;
 *   - `--UNMAP_CONTAMINANT_READS` UNMAPS AN ALIGNMENT WITH TOO FEW UNCLIPPED BASES, counted
 *     against `--MIN_UNCLIPPED_BASES`, and only where the alignment is clipped at BOTH ends;
 *   - THE THRESHOLD IS THE WHOLE TEST, so lowering it keeps the same read;
 *   - AND `--UNMAPPED_READ_STRATEGY` DECIDES WHAT AN UNMAPPED CONTAMINANT REMEMBERS: nothing, a
 *     `PA` tag beside its alignment, or a `PA` tag instead of it.
 *
 * Output:
 *
 *     record\t<case>\t<each merged record's data line>
 *     rc\t<case>\t<the exit status>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MergeBamAlignmentClipDump
 */
import java.io.*;
import java.nio.file.*;
import java.util.*;

public class MergeBamAlignmentClipDump {
  /** chr1, two hundred bases of `ACGT` repeating. */
  static final int LENGTH = 200;
  static final StringBuilder buf = new StringBuilder();

  static String reference() {
    final StringBuilder bases = new StringBuilder();
    for (int index = 0; index < LENGTH; index++) {
      bases.append("ACGT".charAt(index % 4));
    }
    return bases.toString();
  }

  /** The reference's own bases over a window, which a read copies before it is edited. */
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
    final File dir = Files.createTempDirectory("mbac").toFile();
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
    final String bases = window(41, 40);
    final String qualities = "I".repeat(40);

    // An alignment that runs ten bases off the end of its contig.
    // The parser refuses a cigar that maps off the reference, so the tool has to be told to let it
    // in before it can be asked what it does with it.
    run("off-the-end-of-the-contig",
        unmappedHeader + "r\t4\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n",
        alignedHeader + "r\t0\tchr1\t171\t60\t40M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n",
        "VALIDATION_STRINGENCY=SILENT");
    // The same alignment already ending in a soft clip, which the clip has to account for.
    run("off-the-end-with-a-soft-clip",
        unmappedHeader + "r\t4\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n",
        alignedHeader + "r\t0\tchr1\t171\t60\t35M5S\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n",
        "VALIDATION_STRINGENCY=SILENT");
    // And one that fits, for the difference.
    run("inside-the-contig",
        unmappedHeader + "r\t4\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n",
        alignedHeader + "r\t0\tchr1\t41\t60\t40M\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n");

    // An adapter the unmapped bam marked, clipped and kept.
    final String marked = unmappedHeader + "r\t4\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t"
        + qualities + "\tRG:Z:rg1\tXT:i:31\n";
    final String alignedRead = alignedHeader + "r\t0\tchr1\t41\t60\t40M\t*\t0\t0\t" + bases
        + "\t" + qualities + "\tRG:Z:rg1\n";
    run("an-adapter-marked-in-the-unmapped-bam", marked, alignedRead);
    run("an-adapter-left-alone", marked, alignedRead, "CLIP_ADAPTERS=false");

    // A contaminant: an alignment with too few unclipped bases to believe.
    // The filter counts soft-clip BLOCKS as well as aligned bases, and wants two of them, so an
    // alignment clipped at one end only is never a contaminant however short it is.
    final String shortAlignment = alignedHeader
        + "r\t0\tchr1\t41\t60\t15S10M15S\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n";
    final String clippedOneEnd = alignedHeader
        + "r\t0\tchr1\t41\t60\t10M30S\t*\t0\t0\t" + bases + "\t" + qualities
        + "\tRG:Z:rg1\n";
    final String plain = unmappedHeader + "r\t4\t*\t0\t0\t*\t*\t0\t0\t" + bases + "\t"
        + qualities + "\tRG:Z:rg1\n";
    run("a-short-alignment-kept", plain, shortAlignment);
    run("clipped-at-one-end-only", plain, clippedOneEnd, "UNMAP_CONTAMINANT_READS=true");
    run("a-short-alignment-unmapped", plain, shortAlignment, "UNMAP_CONTAMINANT_READS=true");
    run("a-short-alignment-under-a-lower-threshold", plain, shortAlignment,
        "UNMAP_CONTAMINANT_READS=true", "MIN_UNCLIPPED_BASES=8");

    // What an unmapped contaminant remembers.
    for (final String strategy : new String[]{"MOVE_TO_TAG", "COPY_TO_TAG", "DO_NOT_CHANGE_INVALID"}) {
      run("a-contaminant-" + strategy.toLowerCase(), plain, shortAlignment,
          "UNMAP_CONTAMINANT_READS=true", "UNMAPPED_READ_STRATEGY=" + strategy);
    }

    System.out.print(buf);
  }
}
