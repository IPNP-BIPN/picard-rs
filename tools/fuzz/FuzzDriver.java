/*
 * Coverage-guided differential fuzzing of the reference, in one warm JVM.
 *
 * The covering arrays cover interactions between arguments; they say nothing about which *branches
 * of the reference* those interactions reach. A tool can be pairwise-covered and still have most
 * of its logic untouched, because a branch usually needs a particular argument together with a
 * particular shape of input. This closes that gap the way the plan commits to: mutate, measure the
 * reference's branch coverage, and keep whatever reaches something new.
 *
 * Why a Java driver rather than a Python loop of `docker run`: each run of the tool has to happen
 * in the same JVM as the coverage agent, and a fresh container per iteration costs seconds. Here
 * the JVM stays warm, the tool runs in process through PicardCommandLine, and JaCoCo's runtime API
 * hands back the execution data after each iteration.
 *
 * Usage: FuzzDriver <fixtures dir> <seed corpus file> <iterations> <out dir>
 *
 * The seed corpus is one command line per line: `<Tool>\t--ARG=value ...`, which is what
 * run_array.py writes from a covering array. Mutation is seeded and deterministic, so a run that
 * finds something can be replayed exactly.
 */

import org.jacoco.agent.rt.IAgent;
import org.jacoco.agent.rt.RT;
import org.jacoco.core.data.ExecutionData;
import org.jacoco.core.tools.ExecFileLoader;
import picard.cmdline.PicardCommandLine;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class FuzzDriver {

    /** Reaches the protected entry point without the System.exit that main would call. */
    static class Runner extends PicardCommandLine {
        int run(String[] argv) {
            Object code = instanceMain(argv);
            return code instanceof Integer ? (Integer) code : (code == null ? 0 : 1);
        }
    }

    /** One command line and what the reference did with it. */
    static class Case {
        final String tool;
        final List<String> args;
        String outcome = "";
        int newBranches = 0;

        Case(String tool, List<String> args) {
            this.tool = tool;
            this.args = args;
        }

        String commandLine() {
            return tool + "\t" + String.join(" ", args);
        }
    }

    static final IAgent AGENT = RT.getAgent();

    public static void main(String[] args) throws Exception {
        String fixtures = args[0];
        Path seedFile = Path.of(args[1]);
        int iterations = Integer.parseInt(args[2]);
        Path outDir = Path.of(args[3]);
        Files.createDirectories(outDir);

        List<Case> corpus = new ArrayList<>();
        for (String line : Files.readAllLines(seedFile)) {
            if (line.isBlank() || line.startsWith("#")) continue;
            String[] parts = line.split("\t", 2);
            List<String> argv = new ArrayList<>(List.of(parts[1].trim().split("\\s+")));
            corpus.add(new Case(parts[0], argv));
        }
        if (corpus.isEmpty()) throw new IllegalArgumentException("empty seed corpus");

        Random rng = new Random(20260729L);
        Set<String> coveredBranches = new HashSet<>();
        Set<String> seenOutcomes = new HashSet<>();
        List<Case> interesting = new ArrayList<>();

        // The seeds first: they establish the baseline coverage that a mutant has to beat.
        final int seedCount = corpus.size();
        int index = 0;
        for (Case seed : corpus) {
            evaluate(seed, fixtures, outDir, coveredBranches, seenOutcomes, interesting, index++);
        }
        int fromSeeds = coveredBranches.size();

        for (int i = 0; i < iterations; i++) {
            Case parent = corpus.get(rng.nextInt(corpus.size()));
            Case mutant = mutate(parent, rng, fixtures);
            boolean kept = evaluate(
                mutant, fixtures, outDir, coveredBranches, seenOutcomes, interesting, index++);
            // A mutant that reached something new becomes a parent: that is what makes the search
            // guided rather than random.
            if (kept) corpus.add(mutant);
        }

        try (PrintWriter out = new PrintWriter(outDir.resolve("corpus.txt").toFile())) {
            out.println("# tool\targuments\toutcome\tnew_branches");
            for (Case c : interesting) {
                out.printf("%s\t%s\t%d%n", c.commandLine(), c.outcome, c.newBranches);
            }
        }
        System.out.printf(
            "SUMMARY seeds=%d iterations=%d branches_from_seeds=%d branches_total=%d "
                + "interesting=%d distinct_outcomes=%d%n",
            seedCount,
            iterations,
            fromSeeds,
            coveredBranches.size(),
            interesting.size(),
            seenOutcomes.size());
    }

    /**
     * Run one case, measure the reference's branch coverage, and keep it when it reaches a branch
     * nothing has reached before or produces an outcome nothing has produced before.
     */
    static boolean evaluate(
            Case c,
            String fixtures,
            Path outDir,
            Set<String> coveredBranches,
            Set<String> seenOutcomes,
            List<Case> interesting,
            int index)
            throws Exception {
        String output = outDir.resolve("fuzz-" + index + ".txt").toString();
        List<String> argv = new ArrayList<>();
        argv.add(c.tool);
        for (String a : c.args) {
            // Every case writes its own output, so one iteration cannot read another's leftovers.
            if (a.startsWith("--OUTPUT=")) argv.addAll(List.of("--OUTPUT", output));
            else if (a.contains("=")) {
                int eq = a.indexOf('=');
                argv.add(a.substring(0, eq));
                argv.add(a.substring(eq + 1));
            } else argv.add(a);
        }
        if (argv.stream().noneMatch("--OUTPUT"::equals)) argv.addAll(List.of("--OUTPUT", output));

        AGENT.reset();
        PrintStream realOut = System.out, realErr = System.err;
        PrintStream quiet = new PrintStream(new ByteArrayOutputStream());
        String outcome;
        try {
            System.setOut(quiet);
            System.setErr(quiet);
            int code = new Runner().run(argv.toArray(new String[0]));
            outcome = "EXIT=" + code + " sha=" + digest(Path.of(output), c.tool);
        } catch (Throwable t) {
            Throwable root = t;
            while (root.getCause() != null) root = root.getCause();
            outcome = "THROWN=" + root.getClass().getSimpleName() + ": " + root.getMessage();
        } finally {
            System.setOut(realOut);
            System.setErr(realErr);
            new File(output).delete();
        }
        c.outcome = outcome.replace(fixtures, "<FIXTURES>").replace("\n", "\\n");

        Set<String> branches = branchesCovered();
        int before = coveredBranches.size();
        coveredBranches.addAll(branches);
        c.newBranches = coveredBranches.size() - before;

        boolean novelOutcome = seenOutcomes.add(c.outcome);
        boolean kept = c.newBranches > 0 || novelOutcome;
        if (kept) interesting.add(c);
        System.out.printf(
            "%-5d new_branches=%-4d total=%-6d %s%n",
            index, c.newBranches, coveredBranches.size(), c.outcome.substring(0, Math.min(90, c.outcome.length())));
        return kept;
    }

    /**
     * The reference's covered probes, as class#probe identifiers.
     *
     * This reads the agent's execution data directly rather than running JaCoCo's Analyzer over
     * picard.jar. Two reasons, and the second is the one that decided it:
     *
     *   1. The Analyzer re-reads several thousand classes on every iteration, which would dominate
     *      the loop; the execution data is already in memory.
     *   2. picard.jar bundles its own, older asm, and it wins the classpath, so the Analyzer dies
     *      with "Unsupported class file major version 61". Putting JaCoCo's asm first would mean
     *      running the reference against a classpath the oracle contract does not describe, and the
     *      contract is the point of the oracle.
     *
     * A JaCoCo probe sits at a branch boundary, so the set of covered probes is a proxy for branch
     * coverage rather than the thing itself: it is a guidance signal and a stopping condition, not
     * a published coverage figure. Stated here so it is not quoted as one.
     */
    static Set<String> branchesCovered() throws Exception {
        ExecFileLoader loader = new ExecFileLoader();
        loader.load(new ByteArrayInputStream(AGENT.getExecutionData(false)));
        Set<String> covered = new HashSet<>();
        for (ExecutionData data : loader.getExecutionDataStore().getContents()) {
            // Only the reference's own code counts. htsjdk is measured by its own repository, and
            // counting it here would let an htsjdk-only probe masquerade as Picard coverage.
            if (!data.getName().startsWith("picard/")) continue;
            boolean[] probes = data.getProbes();
            for (int i = 0; i < probes.length; i++) {
                if (probes[i]) covered.add(data.getName() + "#" + i);
            }
        }
        return covered;
    }

    static boolean isNumeric(String value) {
        if (value.isEmpty()) return false;
        try {
            Double.parseDouble(value);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    /**
     * The output's fingerprint, over the canonicalized content.
     *
     * Hashing the raw file would make every comparison with the port fail for a reason that is not
     * the tool: a metrics file carries the command line (with per-run paths) and a start time, and
     * the conformance suites already declare both as canonicalized. The same two rules apply here,
     * so a divergence in a fingerprint is a divergence in the metrics.
     */
    static String digest(Path path, String tool) throws Exception {
        if (!Files.exists(path)) return "no-output";
        StringBuilder body = new StringBuilder();
        for (String line : Files.readAllLines(path)) {
            if (line.startsWith("# " + tool) || line.startsWith("# Started on:")) continue;
            body.append(line).append('\n');
        }
        MessageDigest md = MessageDigest.getInstance("SHA-256");
        byte[] hash = md.digest(body.toString().getBytes(java.nio.charset.StandardCharsets.UTF_8));
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 8; i++) sb.append(String.format("%02x", hash[i]));
        return sb.toString();
    }

    /**
     * One edit per mutant, deliberately: a mutant that changes five things and reaches a new
     * branch says nothing about which of the five did it, and the minimizer would have to undo
     * the other four before the case is readable.
     */
    static Case mutate(Case parent, Random rng, String fixtures) {
        List<String> args = new ArrayList<>(parent.args);
        List<Integer> mutable = new ArrayList<>();
        for (int i = 0; i < args.size(); i++) {
            if (args.get(i).contains("=") && !args.get(i).startsWith("--OUTPUT")) mutable.add(i);
        }
        if (mutable.isEmpty()) return new Case(parent.tool, args);

        int at = mutable.get(rng.nextInt(mutable.size()));
        String arg = args.get(at);
        String name = arg.substring(0, arg.indexOf('='));
        String value = arg.substring(arg.indexOf('=') + 1);

        String mutated;
        switch (rng.nextInt(4)) {
            case 0:
                mutated = value.equals("true") ? "false" : value.equals("false") ? "true" : value;
                break;
            case 1:
                // Numeric edge: zero, one, or a large value. Numbers are where off-by-one branches
                // live, and the covering array only ever passes the declared default.
                //
                // Only where the current value is already a number. Replacing a path with "0"
                // produces "Cannot read non-existent file", an outcome the fuzzer reaches on its
                // first mutation and learns nothing further from: it is a rejection by the argument
                // parser, not a branch of the tool.
                if (isNumeric(value)) {
                    mutated = new String[] {"0", "1", "2", "1000000"}[rng.nextInt(4)];
                } else mutated = value;
                break;
            case 2:
                // Swap the input for another fixture: a branch usually needs an argument *and* a
                // shape of input, and the array holds the input fixed per row.
                if (name.equals("--INPUT")) {
                    String[] inputs = {"small.bam", "small.sam", "queryname.bam", "unmapped.bam"};
                    mutated = fixtures + "/" + inputs[rng.nextInt(inputs.length)];
                } else mutated = value;
                break;
            default:
                // Drop the argument: absent and defaulted are different code paths in Barclay.
                args.remove(at);
                return new Case(parent.tool, args);
        }
        args.set(at, name + "=" + mutated);
        return new Case(parent.tool, args);
    }
}
