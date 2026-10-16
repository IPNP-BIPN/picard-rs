import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

/**
 * What Picard prints, and exits with, for a command line it does not run.
 *
 * Each probe is a real `java -jar picard.jar <Tool> <words>` process, because what a user sees is
 * the whole process: the syntax banner htsjdk's `Log` prints before parsing, the parser's message,
 * the usage after it, which stream each goes to, and the exit status. Running the tools in this JVM
 * would record what the parser returns and miss all four.
 *
 * The probes:
 *
 * * `legacy_unknown`, `legacy_lowercase`: a `NAME=value` the tool does not declare. The legacy
 *   parser upper-cases a name before it looks it up, so the message names `PICARD_RS_PROBE` both
 *   times, and the banner keeps the case the user typed;
 * * `legacy_word`: a word without `=` after a `NAME=value` every tool declares, which only a tool
 *   with positional arguments takes, so it is probed only on the others;
 * * `word`: the same word alone, which is not legacy syntax (no word has `=`), so Barclay's parser
 *   refuses it instead;
 * * `help_short`, `help_long`, `version`: `-h`, `--help` and `--version`, which Barclay's parser
 *   takes, since a `-` word turns the legacy parser off;
 * * `barclay_unknown`, `barclay_embedded`: an undeclared `--NAME value` and `--NAME=value`, which
 *   Barclay refuses with two different messages;
 * * `mixed_help`: a `NAME=value` with `-h`, which Picard reads as Barclay syntax after warning.
 *
 * Two fields are a clock: the second htsjdk's `Log` stamps a line with, and the millisecond log4j
 * stamps one with. Both are replaced, and nothing else is.
 *
 * Rows: `run <tool> <probe> <exit> <stdout> <stderr>`, the streams with `\`, newline, tab and
 * carriage return escaped.
 */
public class CommandLineDump {
    static final String[][] PROBES = {
        {"legacy_unknown", "PICARD_RS_PROBE=1"},
        {"legacy_lowercase", "picard_rs_probe=1"},
        {"legacy_word", "VERBOSITY=INFO", "picard_rs_probe"},
        {"word", "picard_rs_probe"},
        {"help_short", "-h"},
        {"help_long", "--help"},
        {"version", "--version"},
        {"barclay_unknown", "--PICARD_RS_PROBE", "1"},
        {"barclay_embedded", "--PICARD_RS_PROBE=1"},
        {"mixed_help", "PICARD_RS_PROBE=1", "-h"},
    };

    static final Pattern LOG_TIME =
        Pattern.compile("(?m)^(INFO|WARNING|ERROR|DEBUG)\t\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}\t");
    static final Pattern LOG4J_TIME = Pattern.compile("(?m)^\\d{2}:\\d{2}:\\d{2}\\.\\d{3} ");

    public static void main(final String[] args) throws Exception {
        final List<String> jvm = new ArrayList<>();
        jvm.add("java");
        for (final String flag : ManagementFactory.getRuntimeMXBean().getInputArguments()) {
            if (flag.startsWith("-D")) {
                jvm.add(flag);
            }
        }
        jvm.add("-jar");
        jvm.add(System.getenv("PICARD_JAR"));

        final ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        final List<Future<String>> rows = new ArrayList<>();
        for (final Map.Entry<String, Class<?>> entry : Declarations.tools().entrySet()) {
            final String tool = entry.getKey();
            final boolean positional = Declarations.takesPositional(Declarations.legacy(entry.getValue()));
            for (final String[] probe : PROBES) {
                if ((probe[0].equals("legacy_word") || probe[0].equals("word")) && positional) {
                    continue;
                }
                final List<String> command = new ArrayList<>(jvm);
                command.add(tool);
                command.addAll(Arrays.asList(probe).subList(1, probe.length));
                rows.add(pool.submit(() -> run(tool, probe[0], command)));
            }
        }
        System.out.println("# CommandLineDump: what Picard prints and exits with for a command line it does not run");
        for (final Future<String> row : rows) {
            System.out.println(row.get());
        }
        pool.shutdown();
    }

    static String run(final String tool, final String probe, final List<String> command) throws Exception {
        final Process process = new ProcessBuilder(command).start();
        process.getOutputStream().close();
        final Future<byte[]> out = Executors.newSingleThreadExecutor().submit(() -> drain(process.getInputStream()));
        final byte[] err = drain(process.getErrorStream());
        if (!process.waitFor(300, TimeUnit.SECONDS)) {
            process.destroyForcibly();
            throw new IllegalStateException(tool + " " + probe + " did not exit");
        }
        return String.join("\t", "run", tool, probe, Integer.toString(process.exitValue()),
            escape(clock(new String(out.get(), StandardCharsets.UTF_8))),
            escape(clock(new String(err, StandardCharsets.UTF_8))));
    }

    static byte[] drain(final InputStream stream) throws Exception {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        stream.transferTo(bytes);
        return bytes.toByteArray();
    }

    static String clock(final String text) {
        final String logged = LOG_TIME.matcher(text).replaceAll("$1\t1970-01-01 00:00:00\t");
        return LOG4J_TIME.matcher(logged).replaceAll("00:00:00.000 ");
    }

    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\n", "\\n").replace("\t", "\\t").replace("\r", "\\r");
    }
}
