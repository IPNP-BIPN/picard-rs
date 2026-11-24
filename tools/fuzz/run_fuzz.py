#!/usr/bin/env python3
"""Coverage-guided differential fuzzing: seed from the covering array, fuzz the reference, replay.

1. **Seed** from the tool's covering array: its rows already cover every pair of values, so the
   search starts where combinatorial coverage stops. Each argument's DOMAIN is every value the
   array gives it, which is what a mutation draws from.
2. **Fuzz** inside the pinned container (FuzzDriver.java): one warm JVM, JaCoCo counting Picard's
   probes, a mutant kept when it reaches a probe or an outcome nothing reached before.
3. **Replay** every kept case on both sides out of process, through run_array.py with the flags the
   manifest's `coverage` entry gives the tool, so a case is compared exactly as a covering-array
   row is. A divergence is **minimized**: arguments the array varies are dropped one at a time
   while both sides, re-run, still disagree and the reference still answers what it answered on
   the full case. Only the first MINIMIZE_PER_TOOL divergences of a tool are; the rest keep their
   full command lines.

This is gatk-rs's `tools/fuzz/` (which began as this repository's driver) moved back onto
`PicardCommandLine` and onto this repository's `run_array.py`.

    python3 tools/fuzz/run_fuzz.py --tool CollectQualityYieldMetrics --iterations 30
    python3 tools/fuzz/run_fuzz.py --tool CollectQualityYieldMetrics --iterations 30 \\
        --port target/release/CollectQualityYieldMetrics
    python3 tools/fuzz/run_fuzz.py --shard 0/8 --iterations 30 --port-dir target/release
    python3 tools/fuzz/run_fuzz.py --merge tools/fuzz/findings   # -> tools/fuzz/measured.json

Findings land in tools/fuzz/findings/ and are published by CI rather than committed; what is
committed is `measured.json`, frozen from a CI run's artefact like every other measurement
(decision 0008).
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "coverage"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "conformance"))
import run_array  # noqa: E402
import compare as comparator  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FINDINGS = REPO / "tools" / "fuzz" / "findings"
MEASURED = REPO / "tools" / "fuzz" / "measured.json"
# Only Picard's own code is counted; htsjdk is measured by its own repository.
INCLUDES = "picard.*"
RNG_SEED = 20261105
MINIMIZE_PER_TOOL = 1


def coverage_entry(tool):
    """The manifest's `coverage` entry for a tool: its port binary and how its rows are compared."""
    for entry in comparator.load_manifest()["coverage"]["tools"]:
        if entry["tool"] == tool:
            return entry
    raise SystemExit(f"{tool} has no coverage entry in tools/conformance/manifest.json")


def seed_file(tool, array, path):
    """The array's rows, plus what the driver needs to write them as `run_array.as_cli` would."""
    rows = [run_array.row_arguments(row, array["excluded"]) for row in array["array"]]
    domains = {}
    for row in array["array"]:
        for pair in row["arguments"]:
            name, _, value = pair.partition("=")
            domains.setdefault(name, [])
            if value not in domains[name]:
                domains[name].append(value)
    positional, lists = run_array.tool_spec(tool)
    _, stdin_spec = run_array.tool_extras(tool)
    lines = [
        "\t".join(["positional", *positional]),
        "\t".join(["lists", *lists]),
        # The driver feeds the first file every case; the replay takes them in turn by row, as the
        # array does, and both start from the first.
        "\t".join(["stdin", *([Path(run_array.stdin_file(stdin_spec, 0)).name] if stdin_spec else [])]),
    ]
    lines += ["\t".join(["domain", name, *values]) for name, values in domains.items()]
    lines += ["\t".join(["row", *args]) for args in rows]
    path.write_text("\n".join(lines) + "\n")
    return len(rows), set(domains)


def fuzz(tool, workdir, fixtures, iterations, timeout):
    """Run the driver in the container; returns its SUMMARY fields and the kept cases."""
    harness = workdir / "fuzz"
    command = (
        # A private, writable copy of the corpus: a tool that asserts its input is writable
        # (`$writable` in fixtures.json) would refuse every case on the read-only mount, and a tool
        # that does write beside its input must not touch the corpus the replay reads.
        "mkdir -p /work/tmp /work/out /tmp/fz && cp -r /work/corpus /work/fixtures "
        "&& cp /harness/FuzzDriver.java /tmp/fz/ "
        '&& javac -cp "$ORACLE_CP:$JACOCO_AGENT:$JACOCO_ANALYSIS" -d /tmp/fz /tmp/fz/FuzzDriver.java '
        f"&& java -javaagent:$JACOCO_AGENT=output=none,includes={INCLUDES} "
        '-cp "/tmp/fz:$ORACLE_CP:$JACOCO_AGENT:$JACOCO_ANALYSIS" '
        f"FuzzDriver {tool} /work/fuzz/seed.txt {iterations} /work/fuzz/corpus.txt {RNG_SEED} {timeout}"
    )
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", run_array.PLATFORM,
            "-v", f"{REPO}/tools/fuzz:/harness:ro",
            "-v", f"{fixtures}:/work/corpus:ro",
            "-v", f"{harness}:/work/fuzz",
            "-w", "/work", run_array.IMAGE, command,
        ],
        capture_output=True,
        text=True,
    )
    print(result.stdout[-3000:])
    summary = [line for line in result.stdout.splitlines() if line.startswith("SUMMARY ")]
    if result.returncode != 0 or not summary:
        print(result.stderr[-3000:])
        raise SystemExit(f"{tool}: the fuzz driver failed")
    fields = dict(word.split("=", 1) for word in summary[-1].split()[1:])
    kept = []
    for line in (harness / "corpus.txt").read_text().splitlines():
        args, _, rest = line.partition("\t#\t")
        reached, fresh = rest.rsplit("\t", 1)
        kept.append({
            "arguments": [a for a in args.split("\t") if a],
            "reached": reached,
            "new_probes": int(fresh),
        })
    return fields, kept


def both_sides(tool, args, port, workdir, fixtures, entry):
    """(reference, port) outcomes of one case, each `exit code + canonical answer`, exactly as
    run_array.py's `main` computes and compares a row."""
    on_stdout = "stderr" if entry.get("output") == "stderr" else entry.get("output") == "stdout"
    strip_pg = bool(entry.get("strip_program_records"))
    output_name = entry.get("output_name") or "output.txt"
    exit_is_result = bool(entry.get("exit_code_is_a_result"))
    run_array.EXTRA_PREFIXES[:] = entry.get("strip_line_prefixes", [])

    def side(run):
        try:
            code, text, error = run()
        except run_array.RowTimeout as hung:
            return f"TIMEOUT: the {hung} ran past {run_array.ROW_TIMEOUT} s"
        return f"{code}\t{run_array.outcome(code, text, error, tool, exit_is_result)}"

    reference = side(lambda: run_array.run_oracle(
        tool, args, workdir, on_stdout, strip_pg, output_name, fixtures))
    if port is None:
        return reference, None
    ours = side(lambda: run_array.run_port(
        port, args, workdir, on_stdout, strip_pg, output_name, fixtures, tool))
    return reference, ours


RERUNS = 6


def raced(tool, args, ours, workdir, fixtures, entry):
    """Whether the reference, run again, answers what the port answered. FifoBuffer's two threads
    die on a negative IO_SIZE and log it in whichever order they race, so one command line has two
    reference answers and no port can give both. A miss is 1 in 2**RERUNS for a 50/50 race; an
    answer the reference never gives is still a divergence."""
    for _ in range(RERUNS):
        reference, _ = both_sides(tool, args, None, workdir, fixtures, entry)
        if reference == ours:
            return True
    return False


def minimize(tool, args, varied, port, workdir, fixtures, entry, original):
    """Drop varied arguments while the divergence survives AND the reference still answers what it
    answered on the full case, re-running both sides each step. Without the second condition
    dropping `--INPUT` "minimizes" every divergence into the reference's missing-input refusal.
    """
    current = list(args)
    changed = True
    while changed:
        changed = False
        for arg in list(current):
            if arg.split("=", 1)[0] not in varied:
                continue
            candidate = [a for a in current if a != arg]
            reference, ours = both_sides(tool, candidate, port, workdir, fixtures, entry)
            if reference == original and reference != ours:
                current, changed = candidate, True
                break
    return current, *both_sides(tool, current, port, workdir, fixtures, entry)


def run(options):
    entry = coverage_entry(options.tool)
    t = entry.get("t", options.t)
    array_path = run_array.ARRAYS / f"{options.tool}.t{t}.json"
    if not array_path.exists():
        raise SystemExit(f"no array at {array_path}")
    array = json.loads(array_path.read_text())
    port = options.port or (str(Path(options.port_dir) / entry["port"]) if options.port_dir else None)

    # mkdtemp, removed below, for the reason shard() gives.
    workdir = Path(tempfile.mkdtemp(prefix="picard-fuzz-"))
    (workdir / "fuzz").mkdir()
    try:
        if options.fixtures_dir:
            fixtures = Path(options.fixtures_dir)
        else:
            fixtures = run_array.build_fixtures(workdir)
        seeds, varied = seed_file(options.tool, array, workdir / "fuzz" / "seed.txt")
        print(f"{options.tool}: {seeds} seed rows from the t={t} array")
        fields, kept = fuzz(options.tool, workdir, fixtures, options.iterations, options.timeout)

        divergences = []
        minimized = 0
        if port:
            for case in kept:
                reference, ours = both_sides(options.tool, case["arguments"], port, workdir, fixtures, entry)
                case["matched"] = reference == ours
                if not case["matched"] and raced(options.tool, case["arguments"], ours, workdir, fixtures, entry):
                    # The port's answer is one the reference gives too, on another run.
                    case["matched"] = case["raced"] = True
                if case["matched"]:
                    continue
                if minimized < MINIMIZE_PER_TOOL:
                    minimized += 1
                    minimal, ref_min, ours_min = minimize(
                        options.tool, case["arguments"], varied, port, workdir, fixtures, entry, reference)
                else:
                    minimal, ref_min, ours_min = case["arguments"], reference, ours
                if all(d["arguments"] != minimal for d in divergences):
                    divergences.append({"arguments": minimal, "reference": ref_min, "port": ours_min})
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    summary = {
        "tool": options.tool,
        "t": t,
        "iterations": options.iterations,
        "seeds": seeds,
        "probes_from_seeds": int(fields["probes_from_seeds"]),
        "probes_total": int(fields["probes_total"]),
        "kept": len(kept),
        "timed_out": fields["timed_out"] == "true",
    }
    if port:
        summary["replayed"] = len(kept)
        summary["diverged"] = sum(not case["matched"] for case in kept)
        summary["raced"] = sum(bool(case.get("raced")) for case in kept)
        summary["divergence_classes"] = len(divergences)
        summary["minimized"] = minimized

    FINDINGS.mkdir(parents=True, exist_ok=True)
    (FINDINGS / f"{options.tool}.json").write_text(json.dumps(
        {"summary": summary, "kept": kept, "divergences": divergences}, indent=1) + "\n")
    print(json.dumps(summary))
    for d in divergences[:5]:
        print(f"  diverges: {' '.join(d['arguments'])}")
        print(f"    reference={d['reference'][:200]!r}")
        print(f"    port     ={d['port'][:200]!r}")
    return 0


def merge(directory):
    """Fold per-tool findings into the committed measurement."""
    tools = {}
    for path in sorted(Path(directory).rglob("*.json")):
        summary = json.loads(path.read_text())["summary"]
        if "failed" in summary:
            print(f"not merged, the driver failed: {summary['tool']}: {summary['failed'][:120]}")
            continue
        tools[summary.pop("tool")] = summary
    MEASURED.write_text(json.dumps({
        "$comment": [
            "Produced by tools/fuzz/run_fuzz.py --merge from a CI run's fuzz-findings artefacts.",
            "`probes_*` count Picard's JaCoCo probes reached (a proxy for branches, not a",
            "percentage); `diverged` counts kept cases on which the port, re-run out of process,",
            "did not answer exactly as the reference did; `divergence_classes` counts their distinct",
            "command lines, minimized for the first `minimized` of them only.",
        ],
        "tools": tools,
    }, indent=1, sort_keys=True) + "\n")
    print(f"wrote {MEASURED} with {len(tools)} tools")
    return 0


def shard(options):
    """One CI slice of every tool the manifest's coverage list measures. A tool the driver cannot
    run is recorded as failed and the slice goes on; `--merge` then leaves it out."""
    k, n = (int(x) for x in options.shard.split("/"))
    tools = sorted(e["tool"] for e in comparator.load_manifest()["coverage"]["tools"]
                   if e["status"] == "measured")
    mine = tools[k::n]
    failed = []
    # mkdtemp rather than TemporaryDirectory: the container writes the fixtures as root, and from
    # Python 3.14 TemporaryDirectory's cleanup raises on them even with ignore_cleanup_errors.
    tmp = tempfile.mkdtemp(prefix="picard-fuzz-fixtures-")
    try:
        if options.fixtures_dir is None:
            options.fixtures_dir = str(run_array.build_fixtures(Path(tmp)))
        for tool in mine:
            options.tool = tool
            try:
                run(options)
            except (SystemExit, Exception) as error:  # noqa: BLE001 -- one tool must not end the slice
                failed.append(tool)
                FINDINGS.mkdir(parents=True, exist_ok=True)
                (FINDINGS / f"{tool}.json").write_text(json.dumps(
                    {"summary": {"tool": tool, "failed": str(error)[:500]}}, indent=1) + "\n")
                print(f"{tool}: FAILED {error}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"shard {k}/{n}: {len(mine) - len(failed)} of {len(mine)} fuzzed; failed: {' '.join(failed) or 'none'}")
    return 1 if mine and len(failed) == len(mine) else 0


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tool")
    parser.add_argument("--t", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--timeout", type=int, default=120, help="seconds per in-process case")
    parser.add_argument("--port", help="the port binary (linux/amd64); omit to fuzz the reference alone")
    parser.add_argument("--port-dir", help="where every tool's port binary is, named as the manifest names it")
    parser.add_argument("--fixtures-dir", help="a corpus run_array.py --build-fixtures already built (its fixtures/)")
    parser.add_argument("--merge", metavar="DIR", help="merge per-tool findings into measured.json")
    parser.add_argument("--shard", metavar="K/N", help="fuzz the K-th of N slices of the measured tools")
    options = parser.parse_args(argv)
    if options.merge:
        return merge(options.merge)
    if options.shard:
        return shard(options)
    if not options.tool:
        parser.error("--tool or --shard is required")
    return run(options)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
