/*
 * Two refusals the port used to answer with "not ported", measured on the reference.
 *
 * `FilterSamReads` with FILTER=includePairedIntervals keeps a paired, primary read when it or its
 * mate overlaps an interval, the mate's extent coming from its MC tag, and WRITE_READS_FILES leaves
 * a `.reads` file beside the output for the input and for the output. The aligned filters walk a
 * queryname-sorted file two records at a time, so a template of three records is an assertion
 * about the record after a first-of-pair rather than a group.
 *
 * `MergeBamAlignment` takes ALIGNED_BAM more than once, and READ1_ALIGNED_BAM and READ2_ALIGNED_BAM
 * more than once each. The aligned files are merged by a MergingSamRecordIterator under one
 * SamFileHeaderMerger, in QUERYNAME order for ALIGNED_BAM and in COORDINATE order for the two ends,
 * so the order the files are named in, their dictionaries, and colliding program ids all show.
 *
 * Behaviours this is built to catch.
 *
 *   - A READ IS KEPT FOR ITS MATE'S OVERLAP, and a secondary, a supplementary or an unpaired read
 *     never is;
 *   - AN INTERVAL ON A CONTIG THE LIST'S HEADER DOES NOT NAME IS IGNORED, and a bad line is the
 *     exception doWork catches, logs and turns into exit status 1;
 *   - THE SORT ORDER IS CHECKED AFTER SORT_ORDER OVERWRITES IT, so a queryname file asked for in
 *     coordinate order is accepted;
 *   - A MAPPED MATE WITHOUT MC is SAMUtils.getMateAlignmentEnd's exception;
 *   - THE .reads FILE HOLDS htsjdk 5's SAMRecord.toString(), the SAM line;
 *   - TWO ALIGNED FILES ARE ONE STREAM, whichever is named first, and their program ids collide;
 *   - ALIGNED FILES WHOSE DICTIONARIES DISAGREE are an AssertionError for ALIGNED_BAM and a
 *     SequenceListsDifferException for the two ends.
 *
 * Output:
 *
 *     argv\t<case>\t<the command line, one argument a line, the case's directory as <dir>>
 *     record\t<case>\t<each output data line>
 *     header\t<case>\t<each output @PG line>
 *     reads\t<case>\t<a .reads file's name and contents>
 *     rc\t<case>\t<the exit status>
 *     log\t<case>\t<what the run logged at ERROR, and the exception it logged>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: FilterAndMergeDump
 */
import java.io.*;
import java.nio.file.*;
import java.util.*;

public class FilterAndMergeDump {
  static final int LENGTH = 1000;
  static final StringBuilder buf = new StringBuilder();
  static File base;

  /** Each case's files, in a directory named after the case, so a run can be replayed. */
  static File caseDirectory(final String name) throws IOException {
    if (base == null) {
      base = Files.createTempDirectory("filterandmerge").toFile();
    }
    final File dir = new File(base, name);
    if (!dir.mkdirs()) {
      throw new IOException("cannot create " + dir);
    }
    return dir;
  }

  static String reference() {
    final StringBuilder bases = new StringBuilder();
    for (int index = 0; index < LENGTH; index++) {
      bases.append("ACGTTGCA".charAt((index * 7 + index / 13) % 8));
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

  static String clean(final String text, final File dir) {
    return text.replace(dir.getAbsolutePath(), "<dir>")
        .replaceAll("^(ERROR|WARNING|INFO)\\t\\d{4}-\\d\\d-\\d\\d \\d\\d:\\d\\d:\\d\\d", "$1");
  }

  /** One run of a tool, with what it said captured: the status, or the throwable that escaped. */
  static Integer invoke(final String name, final File dir, final picard.cmdline.CommandLineProgram tool,
                        final List<String> argv) throws Exception {
    final PrintStream originalOut = System.out;
    final PrintStream originalError = System.err;
    final ByteArrayOutputStream said = new ByteArrayOutputStream();
    Integer code = null;
    try {
      final PrintStream capture = new PrintStream(said, true, "UTF-8");
      System.setOut(capture);
      System.setErr(capture);
      // `Log` keeps the stream it was loaded with, so it is pointed at the capture as well.
      htsjdk.samtools.util.Log.setGlobalPrintStream(capture);
      code = tool.instanceMain(argv.toArray(new String[0]));
    } catch (final Throwable e) {
      Throwable cause = e;
      while (cause.getCause() != null && cause.getCause() != cause) {
        cause = cause.getCause();
      }
      System.setOut(originalOut);
      System.setErr(originalError);
      emit("error", name, cause.getClass().getName() + ":" + clean(String.valueOf(cause.getMessage()), dir));
      return null;
    } finally {
      System.setOut(originalOut);
      System.setErr(originalError);
      htsjdk.samtools.util.Log.setGlobalPrintStream(originalError);
    }
    emit("rc", name, String.valueOf(code));
    if (code != 0) {
      final List<String> logged = new ArrayList<>();
      for (final String line : said.toString("UTF-8").split("\n", -1)) {
        final String kept = clean(line, dir);
        if (kept.startsWith("ERROR") || kept.startsWith("htsjdk.") || kept.startsWith("java.")
            || kept.startsWith("picard.")) {
          logged.add(kept);
        }
      }
      emit("log", name, String.join("\n", logged));
    }
    return code;
  }

  static void emitArguments(final String name, final File dir, final List<String> argv) {
    emit("argv", name, clean(String.join("\n", argv), dir));
  }

  static void records(final String name, final File output, final boolean programs) throws IOException {
    if (!output.exists()) {
      return;
    }
    if (output.getName().endsWith(".bam")) {
      try (htsjdk.samtools.SamReader reader = htsjdk.samtools.SamReaderFactory.makeDefault().open(output)) {
        for (final htsjdk.samtools.SAMRecord record : reader) {
          emit("record", name, record.getSAMString().trim());
        }
      }
      return;
    }
    for (final String line : new String(Files.readAllBytes(output.toPath())).split("\n", -1)) {
      if (line.startsWith("@PG") && programs) {
        emit("header", name, line);
      } else if (!line.startsWith("@") && !line.isEmpty()) {
        emit("record", name, line);
      }
    }
  }

  static File write(final File dir, final String file, final String text) throws IOException {
    final File out = new File(dir, file);
    Files.write(out.toPath(), text.getBytes());
    return out;
  }

  static String line(final String name, final int flag, final int pos, final String cigar,
                     final String mateRef, final int matePos, final int insert, final int start,
                     final String tags) {
    final String bases = window(start, 20);
    return name + "\t" + flag + "\tchr1\t" + pos + "\t" + ((flag & 4) != 0 ? 0 : 60) + "\t" + cigar + "\t" + mateRef + "\t" + matePos
        + "\t" + insert + "\t" + bases + "\t" + "IIIIIIIIIIIIIIIIIIII" + (tags.isEmpty() ? "" : "\t" + tags) + "\n";
  }

  // ---------------------------------------------------------------------------------- FilterSamReads

  static final String COORDINATE_HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:" + LENGTH + "\n";

  /** Pairs on both sides of an interval at 300-400, and the reads it never keeps, in coordinate order. */
  static String coordinateInput(final boolean withMateCigar) {
    return coordinateInput(withMateCigar, true);
  }

  static String coordinateInput(final boolean withMateCigar, final boolean sorted) {
    final String mc = withMateCigar ? "MC:Z:20M" : "";
    final List<String> lines = new ArrayList<>(Arrays.asList(
        line("inside", 99, 310, "20M", "=", 600, 310, 310, "MC:Z:20M"),
        line("mateInside", 99, 100, "20M", "=", 350, 270, 100, mc),
        line("outside", 99, 120, "20M", "=", 500, 400, 120, "MC:Z:20M"),
        line("single", 0, 320, "20M", "*", 0, 0, 320, ""),
        line("secondary", 353, 330, "20M", "=", 700, 390, 330, "MC:Z:20M"),
        line("supplementary", 2145, 340, "20M", "=", 700, 380, 340, "MC:Z:20M"),
        line("mateUnmapped", 73, 360, "20M", "=", 360, 0, 360, ""),
        line("mateUnmapped", 133, 360, "*", "=", 360, 0, 360, "MC:Z:20M"),
        line("endsAtStart", 99, 281, "20M", "=", 401, 140, 281, "MC:Z:20M"),
        line("mateInside", 147, 350, "20M", "=", 100, -270, 350, "MC:Z:20M"),
        line("endsAtStart", 147, 401, "20M", "=", 281, -140, 401, "MC:Z:20M"),
        line("outside", 147, 500, "20M", "=", 120, -400, 500, "MC:Z:20M"),
        line("inside", 147, 600, "20M", "=", 310, -310, 600, "MC:Z:20M")));
    if (sorted) {
      lines.sort(Comparator.comparingInt(l -> Integer.parseInt(l.split("\t")[3])));
    }
    return COORDINATE_HEADER + String.join("", lines);
  }

  static final String INTERVALS = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:" + LENGTH + "\n"
      + "chr1\t300\t400\t+\ttarget\n";

  static void filter(final String name, final String input, final String intervals,
                     final String... extra) throws Exception {
    final File dir = caseDirectory(name);
    final File in = write(dir, "in.sam", input);
    // An OUTPUT among the extra arguments names a BAM instead of the default SAM.
    final boolean bam = Arrays.stream(extra).anyMatch(a -> a.startsWith("OUTPUT="));
    final File output = new File(dir, bam ? "out.bam" : "out.sam");
    final List<String> argv = new ArrayList<>(Arrays.asList("INPUT=" + in.getAbsolutePath()));
    if (!bam) {
      argv.add("OUTPUT=" + output.getAbsolutePath());
    }
    if (intervals != null) {
      argv.add("INTERVAL_LIST=" + write(dir, "targets.interval_list", intervals).getAbsolutePath());
    }
    for (final String arg : extra) {
      argv.add(arg.replace("<dir>", dir.getAbsolutePath()));
    }
    emitArguments(name, dir, argv);
    final Integer code = invoke(name, dir, new picard.sam.FilterSamReads(), argv);
    if (code == null) {
      return;
    }
    records(name, output, false);
    final File[] reads = dir.listFiles((d, f) -> f.endsWith(".reads"));
    Arrays.sort(reads);
    for (final File file : reads) {
      emit("reads", name, file.getName() + "\n" + new String(Files.readAllBytes(file.toPath())));
    }
  }

  static void filterSamReads() throws Exception {
    filter("intervals", coordinateInput(true), INTERVALS, "FILTER=includePairedIntervals");
    filter("intervals_reads_files", coordinateInput(true), INTERVALS, "FILTER=includePairedIntervals",
        "WRITE_READS_FILES=true");
    filter("intervals_out_of_order", coordinateInput(true, false), INTERVALS, "FILTER=includePairedIntervals");
    filter("intervals_out_of_order_bam", coordinateInput(true, false), INTERVALS, "FILTER=includePairedIntervals",
        "OUTPUT=<dir>/out.bam");
    filter("intervals_bam", coordinateInput(true), INTERVALS, "FILTER=includePairedIntervals",
        "OUTPUT=<dir>/out.bam", "WRITE_READS_FILES=true");
    filter("intervals_missing_list", coordinateInput(true), null, "FILTER=includePairedIntervals",
        "INTERVAL_LIST=<dir>/missing.interval_list");
    filter("intervals_unknown_contig_and_empty", coordinateInput(true),
        INTERVALS.replace("target\n", "target\nchr9\t1\t10\t+\tgone\nchr1\t700\t699\t+\tempty\n"),
        "FILTER=includePairedIntervals");
    filter("intervals_no_mate_cigar", coordinateInput(false), INTERVALS, "FILTER=includePairedIntervals");
    final String queryname = coordinateInput(true).replace("SO:coordinate", "SO:queryname");
    filter("intervals_queryname", queryname, INTERVALS, "FILTER=includePairedIntervals");
    filter("intervals_queryname_sorted_to_coordinate", queryname, INTERVALS,
        "FILTER=includePairedIntervals", "SORT_ORDER=coordinate");
    filter("intervals_no_header", coordinateInput(true), "chr1\t300\t400\t+\ttarget\n",
        "FILTER=includePairedIntervals");
    filter("intervals_bad_strand", coordinateInput(true), INTERVALS.replace("+\ttarget", ".\ttarget"),
        "FILTER=includePairedIntervals");
    filter("intervals_beyond_contig", coordinateInput(true), INTERVALS.replace("400", "1400"),
        "FILTER=includePairedIntervals");
    filter("intervals_start_after_end", coordinateInput(true), INTERVALS.replace("300\t400", "300\t200"),
        "FILTER=includePairedIntervals");
    filter("intervals_four_fields", coordinateInput(true), INTERVALS.replace("\ttarget", ""),
        "FILTER=includePairedIntervals");

    // Queryname order. A template with a secondary after its pair, one with the secondary between
    // the two ends, and a lone first-of-pair at the end of the file.
    final String qn = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:" + LENGTH + "\n";
    final String clean = qn
        + line("a", 99, 100, "20M", "=", 200, 120, 100, "")
        + line("a", 147, 200, "20M", "=", 100, -120, 200, "")
        + line("b", 73, 300, "20M", "=", 300, 0, 300, "")
        + line("b", 133, 300, "*", "=", 300, 0, 300, "")
        + line("c", 77, 0, "*", "*", 0, 0, 1, "").replace("\tchr1\t0\t0\t", "\t*\t0\t0\t")
        + line("c", 141, 0, "*", "*", 0, 0, 1, "").replace("\tchr1\t0\t0\t", "\t*\t0\t0\t");
    filter("aligned_include", clean, null, "FILTER=includeAligned", "WRITE_READS_FILES=true");
    filter("aligned_exclude", clean, null, "FILTER=excludeAligned");
    filter("aligned_secondary_after_pair", qn
        + line("a", 99, 100, "20M", "=", 200, 120, 100, "")
        + line("a", 147, 200, "20M", "=", 100, -120, 200, "")
        + line("a", 355, 500, "20M", "=", 200, 0, 500, "")
        + line("b", 99, 300, "20M", "=", 400, 120, 300, "")
        + line("b", 147, 400, "20M", "=", 300, -120, 400, ""), null, "FILTER=includeAligned");
    filter("aligned_secondary_between", qn
        + line("a", 99, 100, "20M", "=", 200, 120, 100, "")
        + line("a", 355, 500, "20M", "=", 200, 0, 500, "")
        + line("a", 147, 200, "20M", "=", 100, -120, 200, ""), null, "FILTER=excludeAligned");
    filter("aligned_lone_first_at_end", qn
        + line("a", 99, 100, "20M", "=", 200, 120, 100, "")
        + line("a", 147, 200, "20M", "=", 100, -120, 200, "")
        + line("z", 73, 300, "20M", "=", 300, 0, 300, ""), null, "FILTER=excludeAligned");
    filter("aligned_names_differ", qn
        + line("a", 99, 100, "20M", "=", 200, 120, 100, "")
        + line("b", 147, 200, "20M", "=", 100, -120, 200, ""), null, "FILTER=includeAligned");
  }

  // ------------------------------------------------------------------------------- MergeBamAlignment

  static final String UNMAPPED_HEADER = "@HD\tVN:1.6\tSO:queryname\n@RG\tID:rg1\tSM:s\tLB:lib1\tPL:ILLUMINA\n";

  static String unmappedPair(final String name, final int start) {
    return name + "\t77\t*\t0\t0\t*\t*\t0\t0\t" + window(start, 20) + "\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n"
        + name + "\t141\t*\t0\t0\t*\t*\t0\t0\t" + reverseComplement(window(start + 150, 20))
        + "\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n";
  }

  static String reverseComplement(final String bases) {
    final StringBuilder out = new StringBuilder();
    for (int index = bases.length() - 1; index >= 0; index--) {
      out.append("TGCA".charAt("ACGT".indexOf(bases.charAt(index))));
    }
    return out.toString();
  }

  static String alignedEnd(final String name, final boolean first, final int start, final String suffix) {
    final int pos = first ? start : start + 150;
    final int matePos = first ? start + 150 : start;
    final int flag = first ? 99 : 147;
    return name + suffix + "\t" + flag + "\tchr1\t" + pos + "\t60\t20M\t=\t" + matePos + "\t"
        + (first ? 170 : -170) + "\t" + window(pos, 20) + "\tIIIIIIIIIIIIIIIIIIII\n";
  }

  static String alignedHeader(final String order, final int length, final String program) {
    return "@HD\tVN:1.6\tSO:" + order + "\n@SQ\tSN:chr1\tLN:" + length + "\n" + program;
  }

  static final String[] NAMES = {"r1", "r2", "r3", "r4", "r5", "r6"};
  static final int[] STARTS = {500, 100, 300, 50, 700, 200};

  static String aligned(final int[] which, final boolean firsts, final boolean seconds, final String suffix1,
                        final String suffix2) {
    final StringBuilder out = new StringBuilder();
    for (final int index : which) {
      if (firsts) out.append(alignedEnd(NAMES[index], true, STARTS[index], suffix1));
      if (seconds) out.append(alignedEnd(NAMES[index], false, STARTS[index], suffix2));
    }
    return out.toString();
  }

  static void merge(final String name, final String[] files, final String[] kinds,
                    final String... extra) throws Exception {
    final File dir = caseDirectory(name);
    final File reference = new File(dir, "ref.fasta");
    try (PrintStream out = new PrintStream(reference)) {
      out.println(">chr1");
      out.println(reference());
    }
    new picard.sam.CreateSequenceDictionary().instanceMain(new String[]{
        "R=" + reference.getAbsolutePath(), "O=" + new File(dir, "ref.dict").getAbsolutePath()});
    final StringBuilder unmapped = new StringBuilder(UNMAPPED_HEADER);
    for (int index = 0; index < NAMES.length; index++) {
      unmapped.append(unmappedPair(NAMES[index], STARTS[index]));
    }
    final File output = new File(dir, "o.sam");
    final List<String> argv = new ArrayList<>(Arrays.asList(
        "UNMAPPED_BAM=" + write(dir, "u.sam", unmapped.toString()).getAbsolutePath(),
        "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
        "OUTPUT=" + output.getAbsolutePath()));
    for (int index = 0; index < files.length; index++) {
      argv.add(kinds[index] + "=" + write(dir, "a" + index + ".sam", files[index]).getAbsolutePath());
    }
    argv.addAll(Arrays.asList(extra));
    emitArguments(name, dir, argv);
    if (invoke(name, dir, new picard.sam.MergeBamAlignment(), argv) != null) {
      records(name, output, true);
    }
  }

  static void mergeBamAlignment() throws Exception {
    final String bwa = "@PG\tID:bwa\tPN:bwa\tVN:1.0\tCL:bwa mem\n";
    final String bwaOther = "@PG\tID:bwa\tPN:bwa\tVN:2.0\tCL:bwa mem -M\n";
    final String[] two = {"ALIGNED_BAM", "ALIGNED_BAM"};
    final String first = alignedHeader("queryname", LENGTH, bwa) + aligned(new int[]{0, 2, 4}, true, true, "", "");
    final String second = alignedHeader("queryname", LENGTH, bwa) + aligned(new int[]{1, 3, 5}, true, true, "", "");
    // r1, r3, r5 against r2, r4, r6: one file's names interleave with the other's.
    merge("aligned_two", new String[]{first, second}, two);
    merge("aligned_two_reversed", new String[]{second, first}, two, "SORT_ORDER=queryname");
    merge("aligned_two_programs_collide", new String[]{first,
        alignedHeader("queryname", LENGTH, bwaOther) + aligned(new int[]{1, 3, 5}, true, true, "", "")
            .replace("IIII\n", "IIII\tPG:Z:bwa\n")}, two, "SORT_ORDER=queryname");
    merge("aligned_two_one_program", new String[]{alignedHeader("queryname", LENGTH, "")
        + aligned(new int[]{0, 2, 4}, true, true, "", ""), second}, two);
    merge("aligned_three_one_empty", new String[]{first, alignedHeader("queryname", LENGTH, bwa), second},
        new String[]{"ALIGNED_BAM", "ALIGNED_BAM", "ALIGNED_BAM"});
    merge("aligned_same_read_in_both", new String[]{first,
        alignedHeader("queryname", LENGTH, bwa) + aligned(new int[]{0}, true, true, "", "").replace("\t60\t", "\t30\t")},
        two);
    merge("aligned_two_unsorted", new String[]{alignedHeader("queryname", LENGTH, bwa)
        + aligned(new int[]{4, 0, 2}, true, true, "", ""), second}, two, "SORT_ORDER=queryname");
    merge("aligned_two_dictionaries_differ", new String[]{first,
        alignedHeader("queryname", LENGTH + 1, bwa) + aligned(new int[]{1, 3, 5}, true, true, "", "")}, two);
    merge("aligned_two_read_group_twice", new String[]{first,
        alignedHeader("queryname", LENGTH, bwa + "@RG\tID:x\tSM:a\n@RG\tID:x\tSM:b\n")
            + aligned(new int[]{1, 3, 5}, true, true, "", "")}, two);

    final String[] ends = {"READ1_ALIGNED_BAM", "READ1_ALIGNED_BAM", "READ2_ALIGNED_BAM", "READ2_ALIGNED_BAM"};
    final String end1a = alignedHeader("queryname", LENGTH, bwa) + aligned(new int[]{0, 2, 4}, true, false, "/1", "");
    final String end1b = alignedHeader("queryname", LENGTH, bwa) + aligned(new int[]{1, 3, 5}, true, false, "/1", "");
    final String end2a = alignedHeader("queryname", LENGTH, bwa) + aligned(new int[]{0, 2, 4}, false, true, "", "/2");
    final String end2b = alignedHeader("queryname", LENGTH, bwa) + aligned(new int[]{1, 3, 5}, false, true, "", "/2");
    // Each end file in name order, merged by coordinate, which is not the order the merge needs.
    merge("ends_two_each", new String[]{end1a, end1b, end2a, end2b}, ends, "SORT_ORDER=queryname");
    merge("ends_two_each_coordinate_output", new String[]{end1a, end1b, end2a, end2b}, ends);
    final String c1a = alignedHeader("coordinate", LENGTH, bwa) + aligned(new int[]{2, 0, 4}, true, false, "/1", "");
    final String c1b = alignedHeader("coordinate", LENGTH, bwa) + aligned(new int[]{3, 1, 5}, true, false, "/1", "");
    final String c2a = alignedHeader("coordinate", LENGTH, bwa) + aligned(new int[]{2, 0, 4}, false, true, "", "/2");
    final String c2b = alignedHeader("coordinate", LENGTH, bwa) + aligned(new int[]{3, 1, 5}, false, true, "", "/2");
    merge("ends_coordinate_files", new String[]{c1a, c1b, c2a, c2b}, ends, "SORT_ORDER=queryname");
    merge("ends_one_and_two", new String[]{end1a, end2a, end2b},
        new String[]{"READ1_ALIGNED_BAM", "READ2_ALIGNED_BAM", "READ2_ALIGNED_BAM"}, "SORT_ORDER=queryname");
    merge("ends_dictionaries_differ", new String[]{end1a,
        alignedHeader("queryname", LENGTH + 1, bwa) + aligned(new int[]{1, 3, 5}, true, false, "/1", ""), end2a, end2b},
        ends);
    merge("ends_programs_collide", new String[]{end1a, end1b, end2a,
        alignedHeader("queryname", LENGTH, bwaOther) + aligned(new int[]{1, 3, 5}, false, true, "", "/2")},
        ends, "SORT_ORDER=queryname");
  }

  public static void main(final String[] args) throws Exception {
    filterSamReads();
    mergeBamAlignment();
    System.out.print(buf);
  }
}
