/*
 * Coverage-guided fuzzing of the reference, in one warm JVM.
 *
 * The covering arrays cover PAIRS of argument values; they say nothing about which branches of
 * Picard those pairs reach. This starts from the array's rows, mutates one argument at a time
 * within the array's own value domains, measures Picard's probes with JaCoCo, and keeps a mutant
 * that reaches a probe no earlier case reached or ends in an outcome no earlier case ended in.
 *
 * The tool runs in process through PicardCommandLine.instanceMain, which returns where main would
 * call System.exit, so the JVM and the agent stay warm across iterations.
 *
 * What this driver decides is only what to KEEP. Whether the port agrees on a kept case is decided
 * by run_fuzz.py, which re-runs both sides out of process through run_array.py, with the same
 * canonicalization every covering-array row gets. An in-process outcome can carry state a fresh
 * JVM would not have, so it guides the search and is never compared.
 *
 * Usage: FuzzDriver <tool> <seed file> <iterations> <corpus out> <rng seed> <timeout seconds>
 *
 * The seed file is tab-separated, one record per line:
 *   positional <word>...       appended last, as run_array.run_oracle does
 *   lists      <name>...       a value of these is words, each `--NAME word`
 *   stdin      [<file>]        a corpus file fed to standard input, for a filter
 *   domain     <name> <v>...   every value the array gives the argument; what a mutation draws from
 *   row        --NAME=value... one seed command line
 */

import org.jacoco.agent.rt.IAgent;
import org.jacoco.agent.rt.RT;
import org.jacoco.core.data.ExecutionData;
import org.jacoco.core.tools.ExecFileLoader;
import picard.cmdline.PicardCommandLine;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class FuzzDriver {

    /** Reaches the protected entry point without the System.exit that main would call. */
    static class Runner extends PicardCommandLine {
        int run(String[] argv) {
            return instanceMain(argv);
        }
    }

    /** A mutation that leaves the argument off; run_array.OMITTED is the array's own spelling. */
    static final String UNSET = "<unset>";
    static final String OMITTED = "<absent>";
    static final String POSITIONAL = "[NA - Positional]";
    static final IAgent AGENT = RT.getAgent();
    static final Path OUT = Path.of("/work/out");
    static final Path FIXTURES = Path.of("/work/fixtures");

    static String tool;
    static List<String> positional = new ArrayList<>();
    static Set<String> lists = new HashSet<>();
    static byte[] stdin;
    static Map<String, List<String>> domains = new LinkedHashMap<>();

    public static void main(String[] args) throws Exception {
        tool = args[0];
        Path seedFile = Path.of(args[1]);
        int iterations = Integer.parseInt(args[2]);
        Path corpusOut = Path.of(args[3]);
        Random rng = new Random(Long.parseLong(args[4]));
        long timeout = Long.parseLong(args[5]);

        List<List<String>> corpus = new ArrayList<>();
        for (String line : Files.readAllLines(seedFile)) {
            if (line.isBlank()) continue;
            String[] f = line.split("\t", -1);
            List<String> rest = new ArrayList<>(Arrays.asList(f).subList(1, f.length));
            rest.removeIf(String::isEmpty);
            switch (f[0]) {
                case "positional": positional.addAll(rest); break;
                case "lists": lists.addAll(rest); break;
                case "stdin": if (!rest.isEmpty()) stdin = Files.readAllBytes(FIXTURES.resolve(rest.get(0))); break;
                case "domain": domains.put(rest.get(0), rest.subList(1, rest.size())); break;
                case "row": corpus.add(rest); break;
                default: throw new IllegalArgumentException("unknown seed record: " + f[0]);
            }
        }
        if (corpus.isEmpty()) throw new IllegalArgumentException("empty seed corpus");

        Set<String> covered = new HashSet<>();
        Set<String> outcomes = new HashSet<>();
        List<String> kept = new ArrayList<>();
        ExecutorService runner = Executors.newSingleThreadExecutor(r -> {
            Thread t = new Thread(r, "fuzz-case");
            t.setDaemon(true);
            return t;
        });

        int seeds = corpus.size(), fromSeeds = 0, run = 0;
        boolean timedOut = false;
        // The seeds first, in order: they are the baseline a mutant has to beat.
        for (int i = 0; i < seeds + iterations && !timedOut; i++) {
            List<String> c = i < seeds ? corpus.get(i) : mutate(corpus.get(rng.nextInt(corpus.size())), rng);
            String outcome = evaluate(c, runner, timeout);
            run++;
            if (outcome.startsWith("TIMEOUT")) timedOut = true;
            Set<String> probes = probesCovered();
            int before = covered.size();
            covered.addAll(probes);
            int fresh = covered.size() - before;
            boolean novel = outcomes.add(outcome);
            if (fresh > 0 || novel) {
                // What the replay runs: the pairs a command line carries, without the placeholders.
                String given = c.stream().filter(p -> !isUnset(value(p))).collect(Collectors.joining("\t"));
                kept.add(given + "\t#\t" + outcome + "\t" + fresh);
                // A mutant that reached something new becomes a parent: that is what makes the
                // search guided rather than random.
                if (i >= seeds) corpus.add(c);
            }
            if (i == seeds - 1) fromSeeds = covered.size();
            System.out.printf("%-5d new=%-5d total=%-7d %s%n", i, fresh, covered.size(),
                outcome.substring(0, Math.min(100, outcome.length())));
        }

        try (PrintWriter out = new PrintWriter(corpusOut.toFile())) {
            for (String line : kept) out.println(line);
        }
        System.out.printf(
            "SUMMARY tool=%s seeds=%d iterations=%d run=%d probes_from_seeds=%d probes_total=%d "
                + "kept=%d outcomes=%d timed_out=%b%n",
            tool, seeds, iterations, run, fromSeeds, covered.size(), kept.size(), outcomes.size(), timedOut);
        // A case that timed out may still be running, and a Picard tool that threw can leave a
        // non-daemon executor behind; nothing else should wait for either.
        System.exit(0);
    }

    static String name(String pair) { int eq = pair.indexOf('='); return eq < 0 ? pair : pair.substring(0, eq); }
    static String value(String pair) { int eq = pair.indexOf('='); return eq < 0 ? "" : pair.substring(eq + 1); }
    static boolean isUnset(String v) { return v.equals(UNSET) || v.equals(OMITTED); }

    /** `--NAME=value` pairs to argv, the way run_array.run_oracle writes them. */
    static String[] argv(List<String> pairs) {
        List<String> out = new ArrayList<>();
        out.add(tool);
        List<String> words = positional;
        for (String pair : pairs) {
            String name = name(pair), value = value(pair);
            if (isUnset(value)) continue;
            if (name.equals(POSITIONAL)) {
                words = Arrays.asList(value.split(" "));
            } else if (lists.contains(name)) {
                for (String word : value.split(" ")) { out.add(name); out.add(word); }
            } else { out.add(name); out.add(value); }
        }
        out.addAll(words);
        return out.toArray(new String[0]);
    }

    static String evaluate(List<String> c, ExecutorService runner, long timeout) throws Exception {
        clearOut();
        AGENT.reset();
        String[] command = argv(c);
        PrintStream realOut = System.out, realErr = System.err;
        InputStream realIn = System.in;
        PrintStream quiet = new PrintStream(new ByteArrayOutputStream());
        String outcome;
        System.setOut(quiet);
        System.setErr(quiet);
        if (stdin != null) System.setIn(new ByteArrayInputStream(stdin));
        try {
            Future<Integer> f = runner.submit(() -> new Runner().run(command));
            outcome = "EXIT=" + f.get(timeout, TimeUnit.SECONDS);
        } catch (TimeoutException e) {
            outcome = "TIMEOUT";
        } catch (Throwable t) {
            Throwable root = t;
            while (root.getCause() != null && root.getCause() != root) root = root.getCause();
            // The class, not the message: a message carries per-case values and paths, and would
            // make every rejection novel.
            outcome = "THROWN=" + root.getClass().getName();
        } finally {
            System.setOut(realOut);
            System.setErr(realErr);
            System.setIn(realIn);
        }
        return outcome;
    }

    static void clearOut() throws IOException {
        Files.createDirectories(OUT);
        try (Stream<Path> walk = Files.walk(OUT)) {
            walk.sorted(Comparator.reverseOrder()).filter(p -> !p.equals(OUT)).forEach(p -> p.toFile().delete());
        }
    }

    /**
     * Picard's covered probes, as class#probe. Read straight from the agent's execution data:
     * running JaCoCo's Analyzer over the jar would re-read thousands of classes per case. A probe
     * sits at a branch boundary, so this is a proxy for branch coverage: a guidance signal, stated
     * here so the count is not quoted as a coverage percentage.
     */
    static Set<String> probesCovered() throws Exception {
        ExecFileLoader loader = new ExecFileLoader();
        loader.load(new ByteArrayInputStream(AGENT.getExecutionData(false)));
        Set<String> covered = new HashSet<>();
        for (ExecutionData data : loader.getExecutionDataStore().getContents()) {
            boolean[] probes = data.getProbes();
            for (int i = 0; i < probes.length; i++) if (probes[i]) covered.add(data.getName() + "#" + i);
        }
        return covered;
    }

    static boolean isNumeric(String v) {
        try { Double.parseDouble(v); return !v.isEmpty(); } catch (NumberFormatException e) { return false; }
    }

    /**
     * One edit per mutant: a mutant that changes several things and reaches a new probe says
     * nothing about which edit did it. Only arguments the array varies are touched; a held fixture
     * value is what lets the tool run at all.
     */
    static List<String> mutate(List<String> parent, Random rng) {
        List<String> args = new ArrayList<>(parent);
        List<Integer> mutable = new ArrayList<>();
        for (int i = 0; i < args.size(); i++) if (domains.containsKey(name(args.get(i)))) mutable.add(i);
        if (mutable.isEmpty()) return args;
        int at = mutable.get(rng.nextInt(mutable.size()));
        String name = name(args.get(at)), value = value(args.get(at));
        List<String> domain = domains.get(name);
        int op = rng.nextInt(4);
        if (op == 0 || op == 1) {
            // Another value the array itself uses: always one the inventory considers valid.
            args.set(at, name + "=" + domain.get(rng.nextInt(domain.size())));
        } else if (op == 2 && isNumeric(value)) {
            // Numeric edges, where off-by-one branches live; the array passes declared values
            // only. Kept small: a large count can make one case run for the whole budget.
            String[] edges = {"0", "1", "-1", "2", "100"};
            args.set(at, name + "=" + edges[rng.nextInt(edges.length)]);
        } else {
            // Absent and defaulted are different paths through Barclay.
            args.set(at, name + "=" + UNSET);
        }
        return args;
    }
}
