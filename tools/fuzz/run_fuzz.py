#!/usr/bin/env python3
"""Drive the coverage-guided fuzzer, then replay what it found against the port.

Three steps, and the third is the one that produces evidence:

1. **Seed** from a covering array. The array's rows are the best starting corpus there is: they
   already cover every pair of arguments, so the fuzzer starts where combinatorial coverage stops
   instead of rediscovering it.
2. **Fuzz** inside the pinned container, one warm JVM, JaCoCo measuring the reference's own probes.
   A mutant is kept when it reaches a probe nothing has reached or produces an outcome nothing has
   produced.
3. **Replay** every kept case against the port and compare the canonicalized output. A divergence
   is then **minimized**: arguments are dropped one at a time for as long as the divergence
   survives, so what gets reported is the smallest command line that still shows it rather than the
   forty-argument line that happened to find it.

    python3 tools/fuzz/run_fuzz.py --tool CollectQualityYieldMetrics --iterations 200
    python3 tools/fuzz/run_fuzz.py --tool CollectAlignmentSummaryMetrics --iterations 200 \\
        --port target/release/collect-alignment-summary-metrics

Findings are written under tools/fuzz/findings/ and are not committed: they are produced in an
emulated container here, and the bit-identity contract accepts results from real x86-64 CI
(decision 0008). CI publishes them as an artefact.
"""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "coverage"))
import run_array  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
FINDINGS = REPO / "tools" / "fuzz" / "findings"
IMAGE = "picard-rs-oracle:3.4.0"
PLATFORM = "linux/amd64"


def seeds_from_array(tool, t, out_path):
    """Write the covering array's rows as the fuzzer's seed corpus."""
    array_path = REPO / "tools" / "coverage" / "arrays" / f"{tool}.t{t}.json"
    if not array_path.exists():
        raise SystemExit(f"no array at {array_path}; generate it in gatk-rs first")
    array = json.load(open(array_path))
    lines = []
    for row in array["array"]:
        args = run_array.row_arguments(row, array["excluded"])
        lines.append(f"{tool}\t" + " ".join(args))
    out_path.write_text("\n".join(lines) + "\n")
    return len(lines)


def fuzz(workdir, seed_file, iterations):
    """Run the driver in the container. Its stdout is the per-iteration coverage log."""
    command = (
        'cp /harness/FuzzDriver.java . '
        '&& javac -cp "$ORACLE_CP:$JACOCO_AGENT:$JACOCO_ANALYSIS" -d . FuzzDriver.java '
        '&& java -javaagent:$JACOCO_AGENT=output=none '
        '-cp ".:$ORACLE_CP:$JACOCO_AGENT:$JACOCO_ANALYSIS" '
        f'FuzzDriver /work/fixtures /work/out/{seed_file.name} {iterations} /work/out'
    )
    result = subprocess.run(
        [
            "docker", "run", "--rm", "--platform", PLATFORM,
            "-v", f"{REPO}/tools/fuzz:/harness:ro",
            "-v", f"{workdir / 'fixtures'}:/work/fixtures:ro",
            "-v", f"{workdir / 'out'}:/work/out",
            "-w", "/work", IMAGE, command,
        ],
        capture_output=True,
        text=True,
    )
    print(result.stdout[-4000:])
    if result.returncode != 0:
        print(result.stderr[-2000:])
        raise SystemExit("the fuzz driver failed")
    return result.stdout


def canonical_digest(path, tool):
    """The port's fingerprint, under the same two rules the driver applies to the reference."""
    if not path.exists():
        return "no-output"
    body = []
    for line in path.read_text().splitlines():
        if line.startswith(f"# {tool}") or line.startswith("# Started on:"):
            continue
        body.append(line)
    text = ("\n".join(body) + "\n") if body else ""
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def run_port(binary, tool, args, workdir):
    """Run the port on one case and return its outcome in the driver's vocabulary."""
    out_dir = workdir / "port"
    out_dir.mkdir(exist_ok=True)
    output = out_dir / "out.txt"
    if output.exists():
        output.unlink()

    argv = [str(binary)]
    for arg in args:
        name, _, value = arg.partition("=")
        value = value.replace("/work/fixtures", str(workdir / "fixtures"))
        if name == "--OUTPUT":
            value = str(output)
        argv.append(f"{name.lstrip('-')}={value}")
    result = subprocess.run(argv, capture_output=True, text=True)
    if result.returncode != 0:
        return f"EXIT={result.returncode} sha=no-output"
    return f"EXIT=0 sha={canonical_digest(output, tool)}"


def oracle_outcome_of(tool, args, workdir):
    """Re-run the reference on a candidate and express the result the way the driver does."""
    code, text, error = run_array.run_oracle(tool, args, workdir)
    if code != 0:
        return f"EXIT={code} sha=no-output"
    produced = workdir / "out" / "output.txt"
    return f"EXIT=0 sha={canonical_digest(produced, tool)}"


def minimize(binary, tool, args, workdir):
    """Drop arguments while the divergence survives, re-running *both* sides at each step.

    The first version of this re-ran only the port and compared against the outcome the reference
    produced for the *original* command line. That is wrong in a way that is easy to miss and
    produces confident nonsense: dropping `--REFERENCE_SEQUENCE` changes what the reference would
    do too, so the reduced case can be reported as diverging when the two sides actually agree on
    it. Correctness costs an oracle run per step, which is a few seconds each and worth it: a
    minimized case is what a reader will trust.
    """
    current = list(args)
    changed = True
    while changed:
        changed = False
        for arg in list(current):
            if arg.startswith("--INPUT") or arg.startswith("--OUTPUT"):
                continue  # the tool cannot run without them
            candidate = [a for a in current if a != arg]
            oracle = oracle_outcome_of(tool, candidate, workdir)
            port = run_port(binary, tool, candidate, workdir)
            if oracle != port:
                current = candidate
                changed = True
                break
    oracle = oracle_outcome_of(tool, current, workdir)
    return current, oracle, run_port(binary, tool, current, workdir)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tool", required=True)
    ap.add_argument("--t", type=int, default=2)
    ap.add_argument("--iterations", type=int, default=100)
    ap.add_argument("--port", help="path to the port binary; omit to fuzz the reference alone")
    args = ap.parse_args(argv)

    workdir = Path(tempfile.mkdtemp(prefix="fuzz-"))
    (workdir / "out").mkdir()
    try:
        print(f"building fixtures in {workdir}")
        run_array.build_fixtures(workdir)

        seed_file = workdir / "out" / "seed.txt"
        seeded = seeds_from_array(args.tool, args.t, seed_file)
        print(f"seeded {seeded} rows from the t={args.t} covering array")

        log = fuzz(workdir, seed_file, args.iterations)
        summary = [l for l in log.splitlines() if l.startswith("SUMMARY")]
        corpus = workdir / "out" / "corpus.txt"
        kept = [
            l.rstrip("\n").split("\t")
            for l in corpus.read_text().splitlines()
            if l.strip() and not l.startswith("#")
        ]

        FINDINGS.mkdir(parents=True, exist_ok=True)
        report = FINDINGS / f"{args.tool}.t{args.t}.txt"

        divergences = []
        if args.port:
            for tool, arg_string, outcome, _new in kept:
                case_args = arg_string.split()
                port_outcome = run_port(args.port, tool, case_args, workdir)
                if port_outcome != outcome:
                    minimal, min_oracle, min_port = minimize(
                        args.port, tool, case_args, workdir
                    )
                    divergences.append((minimal, min_oracle, min_port))

        with open(report, "w") as fh:
            fh.write(f"# {args.tool} t={args.t} iterations={args.iterations}\n")
            for line in summary:
                fh.write(f"# {line}\n")
            for tool, arg_string, outcome, new in kept:
                fh.write(f"keep\t{tool}\t{arg_string}\t{outcome}\t{new}\n")
            for minimal, oracle_outcome, port_outcome in divergences:
                fh.write(
                    f"diverge\t{args.tool}\t{' '.join(minimal)}\t"
                    f"oracle={oracle_outcome}\tport={port_outcome}\n"
                )

        print("\n".join(summary))
        print(f"kept {len(kept)} interesting cases")
        if args.port:
            print(f"port diverges on {len(divergences)} of {len(kept)}")
            for minimal, oracle_outcome, port_outcome in divergences[:5]:
                print(f"  {' '.join(a for a in minimal if not a.startswith('--OUTPUT'))}")
                print(f"    oracle={oracle_outcome}")
                print(f"    port  ={port_outcome}")
        print(f"wrote {report}")
        return 0
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
