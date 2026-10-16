/*
 * The header lines a Picard tool writes about its own run, taken from the reference.
 *
 * Every metrics file opens with `CommandLineProgram.getMetricsFile`'s two default headers, the
 * command line and `Started on: <new Date()>`, and the tools that add a `SAMProgramRecord` put the
 * same command line in its `CL`. The conformance suites strip both lines (`strip_line_prefixes`,
 * `strip_pg`), because they carry paths and a clock reading; this dump is what compares them.
 *
 * picard-rs is run the way Picard's users run Picard, `NAME=value`. That form is what
 * `CommandLineSyntaxTranslater.isLegacyPicardStyle` recognises, and for it
 * `CommandLineProgram.getCommandLineParserForArgs` hands out Barclay's
 * `LegacyCommandLineArgumentParser`, whose command line is built in `checkNumArguments`: the tool's
 * simple name, every option that was SET as ` NAME=<field value>`, a three-space separator, and
 * every option left unset whose default string is not `null` as ` NAME=<default>`, both groups in
 * the parser's own declaration order. So the line depends on three things per tool, and the dump
 * emits all three:
 *
 *   - THE OPTION DEFINITIONS, in the parser's order, with the name, the short name, the default
 *     string, whether the option is a collection, its element type and whether it is optional.
 *     These are read off the parser itself (`optionDefinitions`), not off the source, because a
 *     default is an expression (`Defaults.COMPRESSION_LEVEL`, `getClass().getSimpleName()`) and
 *     only the parser knows what it rendered to;
 *   - THE LINE FOR A SET OF COMMAND LINES, parsed and never run: every required option alone, every
 *     option set (short names, lower-case keys, values the type's constructor normalises), every
 *     optional one set to `null`, and every option set to a second value;
 *   - THE LINES OF A REAL RUN for the seven tools whose `doWork` adds a `@PG` carrying
 *     `getCommandLine()`, and for the metrics files those runs write: MarkDuplicates and its three
 *     variants (`AbstractMarkDuplicatesCommandLineProgram.getChainedPgIds`), DownsampleSam
 *     (`CommandLineProgram.getPGRecord`), PositionBasedDownsampleSam and IntervalListTools (their
 *     own `doWork`), with CollectQualityYieldMetrics for a metrics-only tool.
 *
 * Output:
 *
 *     decl\t<Tool>/<index>\t<name>\t<short>\t<collection>\t<element type>\t<optional>\t<default>
 *     positional\t<Tool>\t<element type>\t<min>\t<max>
 *     version\t<Tool>\t<getVersion()>
 *     argv\t<Tool>/<case>\t<the arguments, one per line>
 *     cl\t<Tool>/<case>\t<getCommandLine() after parseArguments, or FAILED>
 *     input\t<file>\t<the input file's text>
 *     run\t<run>\t<the arguments, one per line>
 *     exit\t<run>\t<instanceMain's return>
 *     pg\t<run>\t<the output header's @PG lines>
 *     readpg\t<run>\t<each output record's name and PG tag>
 *     metrics\t<run>/<file>\t<the metrics file's four header lines>
 *
 * The runs read and write under /work/hl, which the Rust side maps to its own directory.
 */

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.PicardCommandLine;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.ParameterizedType;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

public class HeaderLinesDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** Every tool that writes `getMetricsFile()`'s headers or a `@PG` holding `getCommandLine()`. */
    static final String[] TOOLS = {
            "AccumulateVariantCallingMetrics", "CalculateFingerprintMetrics", "CheckFingerprint",
            "ClusterCrosscheckMetrics", "CollectAlignmentSummaryMetrics",
            "CollectArraysVariantCallingMetrics", "CollectBaseDistributionByCycle",
            "CollectDuplicateMetrics", "CollectGcBiasMetrics", "CollectHiSeqXPfFailMetrics",
            "CollectIlluminaBasecallingMetrics", "CollectIlluminaLaneMetrics",
            "CollectIndependentReplicateMetrics", "CollectInsertSizeMetrics",
            "CollectJumpingLibraryMetrics", "CollectMultipleMetrics", "CollectOxoGMetrics",
            "CollectQualityYieldMetrics", "CollectQualityYieldMetricsFlow",
            "CollectQualityYieldMetricsSNVQ", "CollectRnaSeqMetrics", "CollectRrbsMetrics",
            "CollectSamErrorMetrics", "CollectSequencingArtifactMetrics", "CollectHsMetrics",
            "CollectTargetedPcrMetrics", "CollectUmiPrevalenceMetrics",
            "CollectVariantCallingMetrics", "CollectWgsMetrics", "CollectRawWgsMetrics",
            "CollectWgsMetricsWithNonZeroCoverage", "CompareMetrics", "CompareSAMs",
            "ConvertSequencingArtifactToOxoG", "CreateBafRegressMetricsFile",
            "CreateVerifyIDIntensityContaminationMetricsFile", "CrosscheckFingerprints",
            "CrosscheckReadGroupFingerprints", "DownsampleSam", "EstimateLibraryComplexity",
            "ExtractIlluminaBarcodes", "IlluminaBasecallsToFastq", "IlluminaBasecallsToSam",
            "FindMendelianViolations", "GenotypeConcordance", "MarkDuplicates",
            "MarkDuplicatesWithMateCigar", "SimpleMarkDuplicatesWithMateCigar",
            "UmiAwareMarkDuplicatesWithMateCigar", "MarkIlluminaAdapters", "MeanQualityByCycle",
            "QualityScoreDistribution", "PositionBasedDownsampleSam", "IntervalListTools",
    };

    static Map<String, Class<?>> classes;

    static Object read(final Object target, final String name) throws Exception {
        Class<?> c = target.getClass();
        while (c != null) {
            try {
                final Field f = c.getDeclaredField(name);
                f.setAccessible(true);
                return f.get(target);
            } catch (final NoSuchFieldException e) {
                c = c.getSuperclass();
            }
        }
        throw new NoSuchFieldException(name);
    }

    static CommandLineProgram tool(final String name) throws Exception {
        return (CommandLineProgram) classes.get(name).getDeclaredConstructor().newInstance();
    }

    /**
     * The parser `instanceMain` would use for a legacy-style command line. The argument is only
     * there for `isLegacyPicardStyle` to see a `NAME=value`.
     */
    static CommandLineParser legacyParser(final CommandLineProgram program) {
        return program.getCommandLineParserForArgs(new String[]{"PROBE=probe"});
    }

    /** The option's element class: a collection's type argument, a primitive's wrapper. */
    static Class<?> elementType(final Field field, final boolean collection) {
        Class<?> type = field.getType();
        if (collection) {
            final ParameterizedType generic = (ParameterizedType) field.getGenericType();
            type = (Class<?>) generic.getActualTypeArguments()[0];
        }
        if (type == int.class) return Integer.class;
        if (type == long.class) return Long.class;
        if (type == double.class) return Double.class;
        if (type == float.class) return Float.class;
        if (type == boolean.class) return Boolean.class;
        if (type == short.class) return Short.class;
        if (type == byte.class) return Byte.class;
        if (type == char.class) return Character.class;
        return type;
    }

    record Option(String name, String shortName, boolean collection, Class<?> type,
                  boolean optional, int minElements, String defaultValue, Set<String> mutex) { }

    @SuppressWarnings("unchecked")
    static List<Option> options(final CommandLineParser parser) throws Exception {
        final List<Option> out = new ArrayList<>();
        for (final Object def : (List<?>) read(parser, "optionDefinitions")) {
            final boolean collection = (Boolean) read(def, "isCollection");
            out.add(new Option(
                    (String) read(def, "name"),
                    (String) read(def, "shortName"),
                    collection,
                    elementType((Field) read(def, "field"), collection),
                    (Boolean) read(def, "optional"),
                    (Integer) read(def, "minElements"),
                    (String) read(def, "defaultValue"),
                    (Set<String>) read(def, "mutuallyExclusive")));
        }
        return out;
    }

    /** Parse, never run: the line `parseArgs` stores, or FAILED when the parser refuses. */
    static String commandLine(final String name, final List<String> argv) throws Exception {
        final CommandLineProgram program = tool(name);
        final String[] args = argv.toArray(new String[0]);
        final CommandLineParser parser = program.getCommandLineParserForArgs(args);
        final ByteArrayOutputStream messages = new ByteArrayOutputStream();
        final PrintStream sink = new PrintStream(messages);
        boolean ok;
        try {
            ok = parser.parseArguments(sink, args);
        } catch (final RuntimeException e) {
            sink.println(e);
            ok = false;
        }
        if (!ok) {
            // On standard error, which the suite discards: what the parser said, for whoever
            // runs this by hand to see why a candidate was dropped.
            sink.flush();
            final String said = messages.toString();
            System.err.println(name + " " + argv + ": " + said.substring(0, Math.min(400, said.length())));
        }
        return ok ? parser.getCommandLine() : "FAILED";
    }

    /** A value of the option's type, in one of two flavours, as its constructor will read it. */
    static List<String> values(final Option o, final int flavour, final boolean plain) {
        final Class<?> t = o.type();
        final String lower = o.name().toLowerCase(Locale.ROOT);
        if (t.isEnum()) {
            final Object[] constants = t.getEnumConstants();
            final String first = ((Enum<?>) constants[0]).name();
            final String last = ((Enum<?>) constants[constants.length - 1]).name();
            return flavour == 0 ? List.of(first, last) : List.of(last);
        }
        if (t == Boolean.class) return List.of(flavour == 0 ? (plain ? "true" : "TRUE") : "false");
        if (t == Integer.class || t == Long.class || t == Short.class || t == Byte.class) {
            if (plain) return List.of("7");
            return flavour == 0 ? List.of("+0042", "3") : List.of("-3");
        }
        if (t == Double.class || t == Float.class) {
            if (plain) return List.of("0.5");
            return flavour == 0 ? List.of("0.250", "1E7") : List.of("1e-3");
        }
        if (t == File.class) {
            if (plain) return List.of("/work/in/" + lower);
            return flavour == 0 ? List.of("dir//" + lower + "/", "./b") : List.of("/abs//x/");
        }
        if (t == Character.class) return List.of("c");
        if (plain) return List.of("/work/in/" + lower);
        return flavour == 0 ? List.of("dir//" + lower + ".x", "b c") : List.of("s");
    }

    /** `name=value` pairs for one option, one per element. */
    static List<String> assignments(final Option o, final String key, final List<String> values) {
        final List<String> out = new ArrayList<>();
        final List<String> chosen = o.collection() ? values : values.subList(0, 1);
        for (final String v : chosen) out.add(key + "=" + v);
        return out;
    }

    static List<String> positionals(final CommandLineParser parser, final boolean plain) throws Exception {
        final Field field = (Field) read(parser, "positionalArguments");
        if (field == null) return List.of();
        final int min = (Integer) read(parser, "minPositionalArguments");
        final List<String> out = new ArrayList<>();
        for (int i = 0; i < min; i++) out.add(plain ? "/work/in/p" + i : "pos//p" + i + "/");
        return out;
    }

    /** Add each candidate to the command line only if the parser still accepts it. */
    static List<String> grow(final String name, final List<String> base,
                             final List<List<String>> candidates, final List<String> tail) throws Exception {
        final List<String> argv = new ArrayList<>(base);
        for (final List<String> candidate : candidates) {
            final List<String> trial = new ArrayList<>(argv);
            trial.addAll(candidate);
            final List<String> full = new ArrayList<>(trial);
            full.addAll(tail);
            if (!commandLine(name, full).equals("FAILED")) argv.addAll(candidate);
        }
        argv.addAll(tail);
        return argv;
    }

    static boolean required(final Option o) {
        return o.collection() ? o.minElements() > 0 : !o.optional();
    }

    static void declarations(final String name) throws Exception {
        final CommandLineProgram program = tool(name);
        final CommandLineParser parser = legacyParser(program);
        final List<Option> options = options(parser);
        for (int i = 0; i < options.size(); i++) {
            final Option o = options.get(i);
            emit("decl", name + "/" + i, String.join("\t", o.name(), o.shortName(),
                    String.valueOf(o.collection()), o.type().getName(), String.valueOf(o.optional()),
                    o.defaultValue()));
        }
        final Field positional = (Field) read(parser, "positionalArguments");
        if (positional != null) {
            emit("positional", name, String.join("\t", elementType(positional, true).getName(),
                    String.valueOf(read(parser, "minPositionalArguments")),
                    String.valueOf(read(parser, "maxPositionalArguments"))));
        }
        emit("version", name, String.valueOf(parser.getVersion()));

        // A required option is required "unless any of" its mutually exclusive partners is set,
        // and setting both is refused, so the first of a pair stands for both.
        final List<String> requiredArgs = new ArrayList<>();
        final Set<String> set = new HashSet<>();
        for (final Option o : options) {
            if (!required(o) || o.mutex().stream().anyMatch(set::contains)) continue;
            set.add(o.name());
            final List<String> plain = new ArrayList<>();
            for (int i = 0; i < Math.max(1, o.minElements()); i++) {
                plain.add(values(o, 0, true).get(0) + (i == 0 ? "" : String.valueOf(i)));
            }
            requiredArgs.addAll(o.collection() ? assignments(o, o.name(), plain)
                    : assignments(o, o.name(), plain.subList(0, 1)));
        }
        final List<String> plainTail = positionals(parser, true);
        final List<String> trickyTail = positionals(parser, false);

        final Map<String, List<String>> cases = new LinkedHashMap<>();
        final List<String> requiredOnly = new ArrayList<>(requiredArgs);
        requiredOnly.addAll(plainTail);
        cases.put("required", requiredOnly);

        // Every option, last declared first, with the short name where there is one and every
        // third key in lower case: the line comes back in declaration order and in long names.
        final List<List<String>> every = new ArrayList<>();
        for (int i = options.size() - 1; i >= 0; i--) {
            final Option o = options.get(i);
            if (required(o)) continue;
            String key = o.shortName().isEmpty() ? o.name() : o.shortName();
            if (i % 3 == 0) key = key.toLowerCase(Locale.ROOT);
            every.add(assignments(o, key, values(o, 0, false)));
        }
        cases.put("every", grow(name, requiredArgs, every, trickyTail));

        final List<List<String>> nulls = new ArrayList<>();
        for (final Option o : options) {
            if (!required(o)) nulls.add(List.of(o.name() + "=null"));
        }
        cases.put("nulls", grow(name, requiredArgs, nulls, plainTail));

        final List<List<String>> second = new ArrayList<>();
        for (final Option o : options) {
            if (!required(o)) second.add(assignments(o, o.name(), values(o, 1, false)));
        }
        cases.put("second", grow(name, requiredArgs, second, plainTail));

        for (final Map.Entry<String, List<String>> c : cases.entrySet()) {
            emit("argv", name + "/" + c.getKey(), String.join("\n", c.getValue()));
            emit("cl", name + "/" + c.getKey(), commandLine(name, c.getValue()));
        }
    }

    // ---------------------------------------------------------------------------------------
    // The runs.
    // ---------------------------------------------------------------------------------------

    static final String DIR = "/work/hl";

    /**
     * Fragments on one contig, coordinate-sorted, with a header whose `@PG` ids collide with the
     * new records' (`MarkDuplicates`, `DownsampleSam`, `1`) and a chain (`samtools` follows `bwa`),
     * and reads that carry `PG:Z:bwa`, `PG:Z:samtools` or no `PG` at all. The names are Illumina
     * names so that the optical-duplicate finder and PositionBasedDownsampleSam can read a tile.
     */
    static String inputSam() {
        final StringBuilder sam = new StringBuilder();
        sam.append("@HD\tVN:1.6\tSO:coordinate\n");
        sam.append("@SQ\tSN:chr1\tLN:10000\n");
        sam.append("@RG\tID:rg1\tSM:sample\tLB:lib\tPL:ILLUMINA\n");
        sam.append("@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n");
        sam.append("@PG\tID:samtools\tPN:samtools\tPP:bwa\n");
        sam.append("@PG\tID:MarkDuplicates\tPN:MarkDuplicates\n");
        sam.append("@PG\tID:DownsampleSam\tPN:DownsampleSam\n");
        sam.append("@PG\tID:1\tPN:other\n");
        final String[] pgs = {"bwa", "samtools", null};
        for (int i = 0; i < 24; i++) {
            final int start = 100 + (i / 2) * 37;
            final String name = "H0:1:FC:1:" + (1101 + i % 3) + ":" + (1000 + 101 * i) + ":" + (2000 + 53 * i);
            final int flag = (i % 4 == 3) ? 16 : 0;
            sam.append(name).append('\t').append(flag).append("\tchr1\t").append(start)
                    .append("\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII")
                    .append("\tRG:Z:rg1\tRX:Z:").append(i % 5 == 0 ? "AAAA" : "CCGG");
            final String pg = pgs[i % 3];
            if (pg != null) sam.append("\tPG:Z:").append(pg);
            sam.append('\n');
        }
        return sam.toString();
    }

    static String intervalList() {
        return "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n@PG\tID:1\tPN:other\n"
                + "chr1\t100\t200\t+\ta\nchr1\t300\t400\t+\tb\n";
    }

    static void run(final String id, final String name, final String output, final String... args) throws Exception {
        final List<String> argv = List.of(args);
        emit("run", id, name + "\n" + String.join("\n", argv));
        int code;
        try {
            code = tool(name).instanceMain(args);
        } catch (final Throwable t) {
            code = -1;
            emit("thrown", id, t.toString());
        }
        emit("exit", id, String.valueOf(code));
        if (output != null && new File(output).exists()) {
            final StringBuilder pg = new StringBuilder();
            if (output.endsWith(".interval_list")) {
                for (final String line : Files.readAllLines(Path.of(output))) {
                    if (line.startsWith("@PG")) pg.append(line).append('\n');
                }
            } else {
                try (SamReader reader = SamReaderFactory.makeDefault()
                        .validationStringency(ValidationStringency.SILENT).open(new File(output))) {
                    final java.io.StringWriter text = new java.io.StringWriter();
                    new SAMTextHeaderCodec().encode(text, reader.getFileHeader());
                    for (final String line : text.toString().split("\n")) {
                        if (line.startsWith("@PG")) pg.append(line).append('\n');
                    }
                    final StringBuilder tags = new StringBuilder();
                    for (final SAMRecord rec : reader) {
                        tags.append(rec.getReadName()).append(' ')
                                .append(rec.getStringAttribute("PG")).append('\n');
                    }
                    emit("readpg", id, tags.toString());
                }
            }
            emit("pg", id, pg.toString());
        }
        for (final String arg : args) {
            final String value = arg.substring(arg.indexOf('=') + 1);
            if (value.endsWith(".metrics") && new File(value).exists()) {
                final List<String> lines = Files.readAllLines(Path.of(value));
                emit("metrics", id + "/" + new File(value).getName(),
                        String.join("\n", lines.subList(0, Math.min(4, lines.size()))));
            }
        }
    }

    static void runs() throws Exception {
        Files.createDirectories(Path.of(DIR));
        final String sam = inputSam();
        emit("input", "in.sam", sam);
        Files.writeString(Path.of(DIR, "in.sam"), sam);
        final String intervals = intervalList();
        emit("input", "in.interval_list", intervals);
        Files.writeString(Path.of(DIR, "in.interval_list"), intervals);

        // A BAM of the same records, which every tool reads.
        try (SamReader reader = SamReaderFactory.makeDefault().open(new File(DIR, "in.sam"));
             SAMFileWriter writer = new SAMFileWriterFactory()
                     .makeBAMWriter(reader.getFileHeader(), true, new File(DIR, "in.bam"))) {
            for (final SAMRecord rec : reader) writer.addAlignment(rec);
        }
        final String in = "I=" + DIR + "/in.bam";

        run("md", "MarkDuplicates", DIR + "/md.bam", in, "O=" + DIR + "/md.bam", "M=" + DIR + "/md.metrics");
        run("md-named", "MarkDuplicates", DIR + "/md-named.bam", in, "O=" + DIR + "/md-named.bam",
                "M=" + DIR + "/md-named.metrics", "PROGRAM_RECORD_ID=bwa", "PROGRAM_GROUP_NAME=MD",
                "PROGRAM_GROUP_VERSION=9", "PROGRAM_GROUP_COMMAND_LINE=marked by hand");
        run("md-nopg", "MarkDuplicates", DIR + "/md-nopg.bam", in, "O=" + DIR + "/md-nopg.bam",
                "M=" + DIR + "/md-nopg.metrics", "PROGRAM_RECORD_ID=null");
        run("md-notag", "MarkDuplicates", DIR + "/md-notag.bam", in, "O=" + DIR + "/md-notag.bam",
                "M=" + DIR + "/md-notag.metrics", "ADD_PG_TAG_TO_READS=false");
        run("mdwmc", "MarkDuplicatesWithMateCigar", DIR + "/mdwmc.bam", in, "O=" + DIR + "/mdwmc.bam",
                "M=" + DIR + "/mdwmc.metrics");
        run("simple", "SimpleMarkDuplicatesWithMateCigar", DIR + "/simple.bam", in,
                "O=" + DIR + "/simple.bam", "M=" + DIR + "/simple.metrics");
        run("umi", "UmiAwareMarkDuplicatesWithMateCigar", DIR + "/umi.bam", in, "O=" + DIR + "/umi.bam",
                "M=" + DIR + "/umi.metrics", "UMI_METRICS_FILE=" + DIR + "/umi-umi.metrics");
        run("downsample", "DownsampleSam", DIR + "/downsample.bam", in, "O=" + DIR + "/downsample.bam",
                "P=0.5", "RANDOM_SEED=7", "METRICS_FILE=" + DIR + "/downsample.metrics");
        // The output of the run above, read again with the same seed: its header now holds a
        // DownsampleSam record with no `rs` (the input's), which reseeds from that record's
        // hashCode, and the one just written, whose `rs` is the seed that reseeding lands on, so
        // the seed moves a second time.
        run("downsample-again", "DownsampleSam", DIR + "/downsample-again.bam", "I=" + DIR + "/downsample.bam",
                "O=" + DIR + "/downsample-again.bam", "P=0.5", "RANDOM_SEED=7");
        run("position", "PositionBasedDownsampleSam", DIR + "/position.bam", in,
                "O=" + DIR + "/position.bam", "FRACTION=0.5");
        run("intervals", "IntervalListTools", DIR + "/out.interval_list", "I=" + DIR + "/in.interval_list",
                "O=" + DIR + "/out.interval_list");
        run("yield", "CollectQualityYieldMetrics", null, in, "O=" + DIR + "/yield.metrics");
    }

    public static void main(final String[] args) throws Exception {
        // Every command-line program, by simple name, as `PicardCommandLine` lists them.
        classes = new TreeMap<>();
        PicardCommandLine.processAllCommandLinePrograms(
                Collections.singletonList("picard"), (clazz, properties) -> classes.put(clazz.getSimpleName(), clazz));
        for (final String name : TOOLS) {
            if (!classes.containsKey(name)) throw new IllegalStateException("no tool " + name);
            declarations(name);
        }
        runs();
        System.out.print(buf);
    }
}
