#!/usr/bin/env python3
"""Measure what each picard-rs tool costs against Picard 3.5.0, at scale, side by side.

    python3 tools/speed/bench.py --port /target/dist --port-volume picperf-target \\
        --derived ../scale-derived --out speed.json [--tool SortSam ...]
    python3 tools/speed/bench.py --table speed.json            # the ranked Markdown table
    python3 tools/speed/bench.py --compare before.json after.json

Ported from gatk-rs's tools/speed/bench.py, with its rules:

* **Both sides in the same container**, the speed image (tools/speed/Dockerfile), for the host's
  own architecture, so neither side is emulated. `timer.py` times inside it, so the container's
  start is in neither number.
* **Same resources.** One cgroup per row (`--memory`, default 6 GiB) and the same CPU set
  (`--cpus`, default 4) for both sides, and `RAYON_NUM_THREADS` set to that count for the port.
  Picard runs as users get it, `java -jar picard.jar`, with the JVM's default heap (a quarter of
  the cgroup). A side that the cgroup kills is recorded as such, not retried with more.
* **Three numbers**, CPU (user + system, all threads), peak RSS and wall clock, per side, the
  median of the repetitions, and their ratios port / Picard.
* **Bytes.** The md5 of every file the port wrote is kept (`timer.py`), so `--compare` says whether
  an optimisation changed a byte. Byte identity with Picard is the conformance suites' claim, on
  x86-64; on arm64 Picard has no GKL and its compressed bytes are the JDK's, so they are not used.

Rows are `tools/speed/rows.txt`: one tool command line each, over the inputs
`make_scale_inputs.py` derives from the chr20 scale data.
"""

import argparse
import re
import json
import os
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

SPEED = Path(__file__).resolve().parent
IMAGE = "picard-rs-speed:3.5.0"
VARIABLES = {
    "$I": "/derived/bench.bam", "$R": "/derived/bench.fasta", "$D": "/derived/bench.dict",
    "$Q": "/derived/qname.bam", "$T": "/derived/targets.interval_list", "$B": "/derived",
    "$G": "/derived/big", "$N": "/derived/nonbam", "$S": "/data", "$P": "/derived/prep",
    "$O": "/work/out", "$W": "/work/w", "$A": "/derived/nonbam/arrays",
    "$RUN": "/derived/nonbam/ill_big/Data/Intensities/BaseCalls",
    "$R384": "/derived/nonbam/ill_384/Data/Intensities/BaseCalls",
}


def expand(text):
    for name in sorted(VARIABLES, key=len, reverse=True):
        text = text.replace(name, VARIABLES[name])
    return text


def read_rows(path):
    """`(label, tool, arguments, prep)` per row; a `prep` line runs once with Picard, unmeasured."""
    rows, prep = [], []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("prep "):
            prep.append(expand(line[5:]))
            continue
        label, _, arguments = line.partition(" ")
        rows.append((label, label.split("@")[0], expand(arguments)))
    return rows, prep


def median(values):
    return statistics.median(values) if values else None


def ratio(a, b):
    return a / b if a is not None and b else None


def docker(options, work, command, derived_writable=False):
    return [
        "docker", "run", "--rm",
        "--memory", options.memory, "--memory-swap", options.memory,
        "--cpuset-cpus", f"0-{options.cpus - 1}",
        "-e", f"RAYON_NUM_THREADS={options.cpus}",
        "-v", f"{options.data}:/data:ro",
        "-v", f"{options.derived}:/derived{'' if derived_writable else ':ro'}",
        "-v", f"{options.port_volume}:{options.port_mount}:ro",
        "-v", f"{SPEED}:/harness:ro",
        "-v", f"{work}:/work",
        IMAGE, command,
    ]


def run_prep(options, prep):
    """Each prep command once, in order, when its marker is missing: Picard makes those inputs."""
    marks = Path(options.derived) / "prep"
    marks.mkdir(exist_ok=True)
    for index, command in enumerate(prep):
        mark = marks / f".done{index}"
        if mark.exists():
            continue
        # A directory argument Picard writes into must exist before it runs.
        dirs = " ".join(a.split("=", 1)[1] for a in command.split() if a.startswith("OUTPUT_DIR="))
        with tempfile.TemporaryDirectory() as work:
            done = subprocess.run(docker(options, work, f"mkdir -p /derived/prep {dirs} && cd /derived/prep "
                                         f"&& java -jar $PICARD_JAR {command}",
                                         derived_writable=True), capture_output=True)
        if done.returncode == 0:
            mark.touch()
        else:
            # The rows reading it fail and say so; the next run tries again.
            print(f"prep failed, its rows will too: {command[:100]}", file=sys.stderr)


def time_row(options, tool, arguments):
    """`timer.py` once in the container over one row; its raw JSON."""
    # `$W/<name>`: a writable copy of /derived/<name>, made before each run and not timed.
    copies = sorted(set(re.findall(r"/work/w/([\w.]+)", arguments)))
    spec = {
        "setup": " && ".join(f"cp /derived/{name} /work/w/" for name in copies) or None,
        "reps": options.reps,
        "timeout": options.timeout,
        "oracle": f"java -jar $PICARD_JAR {tool} {arguments}",
        "port": f"{options.port}/{tool} {arguments}",
    }
    with tempfile.TemporaryDirectory() as work:
        (Path(work) / "spec.json").write_text(json.dumps(spec))
        result = subprocess.run(docker(options, work, "python3 /harness/timer.py /work/spec.json"),
                                capture_output=True, text=True)
    if result.returncode != 0:
        return {"error": (result.stderr or result.stdout)[-400:]}
    return json.loads(result.stdout.strip().splitlines()[-1])


def summarise(raw):
    """Medians over the summarised repetitions (one onwards) and the port / Picard ratios."""
    if "cold" not in raw:
        return raw

    def side(runs):
        runs = runs[1:] or runs
        return {
            "wall": median([r["wall"] for r in runs]),
            "cpu": median([r["cpu"] for r in runs]),
            "rss_mib": median([r["rss_kib"] for r in runs]) / 1024,
            "exit": max((r["exit"] for r in runs), key=abs),
            "timed_out": any(r["timed_out"] for r in runs),
        }

    oracle, port = side(raw["cold"]["oracle"]), side(raw["cold"]["port"])
    return {
        "picard": oracle,
        "port": port,
        "cpu_ratio": ratio(port["cpu"], oracle["cpu"]),
        "rss_ratio": ratio(port["rss_mib"], oracle["rss_mib"]),
        "wall_ratio": ratio(port["wall"], oracle["wall"]),
        "port_outputs": raw.get("port_outputs", {}),
        "errors": raw.get("errors", {}),
        "raw": raw["cold"],
    }


def verdict(entry):
    if "picard" not in entry:
        return "harness error"
    if entry["port"]["exit"] != 0:
        return "port failed"
    if entry["picard"]["exit"] != 0:
        return "Picard failed"
    cpu_ok = entry["cpu_ratio"] is not None and entry["cpu_ratio"] <= 0.5
    rss_ok = entry["rss_ratio"] is not None and entry["rss_ratio"] <= 0.25
    return "at target" if cpu_ok and rss_ok else " / ".join(
        part for part, ok in (("CPU", cpu_ok), ("RSS", rss_ok)) if not ok) + " over"


def table(document):
    """Markdown, worst RSS ratio first."""
    rows = document["rows"]

    def key(item):
        entry = item[1]
        return -(entry.get("rss_ratio") or 99)

    lines = ["| Row | Picard CPU s | port CPU s | CPU ratio | Picard RSS MiB | port RSS MiB | "
             "RSS ratio | Picard wall s | port wall s | Verdict |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for label, entry in sorted(rows.items(), key=key):
        if "picard" not in entry:
            lines.append(f"| {label} | | | | | | | | | harness error |")
            continue
        p, r = entry["picard"], entry["port"]

        def f(value, digits=2):
            return "" if value is None else f"{value:.{digits}f}"

        lines.append(
            f"| {label} | {f(p['cpu'], 1)} | {f(r['cpu'], 1)} | {f(entry['cpu_ratio'])} | "
            f"{f(p['rss_mib'], 0)} | {f(r['rss_mib'], 0)} | {f(entry['rss_ratio'])} | "
            f"{f(p['wall'], 1)} | {f(r['wall'], 1)} | {verdict(entry)} |")
    return "\n".join(lines)


def compare(first, second):
    """Per row present in both: the ratio moves and whether any port output's md5 changed."""
    a, b = (json.loads(Path(p).read_text())["rows"] for p in (first, second))
    changed = 0
    for label in sorted(set(a) & set(b)):
        x, y = a[label], b[label]
        if "picard" not in x or "picard" not in y:
            continue
        bytes_note = "same bytes" if x["port_outputs"] == y["port_outputs"] else "BYTES CHANGED"
        changed += bytes_note != "same bytes"
        print(f"{label:48s} cpu {x['port']['cpu']:8.1f} -> {y['port']['cpu']:8.1f} s   "
              f"rss {x['port']['rss_mib']:7.0f} -> {y['port']['rss_mib']:7.0f} MiB   {bytes_note}")
    return changed


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--port", default="/target/dist", help="binary directory inside the container")
    parser.add_argument("--port-volume", default="picperf-target", help="volume or host dir holding it")
    parser.add_argument("--port-mount", default="/target")
    parser.add_argument("--data", default=str(Path.home() / "GATK/scale-data"))
    parser.add_argument("--derived", default=str(Path.home() / "GATK/scale-derived"))
    parser.add_argument("--rows", default=str(SPEED / "rows.txt"))
    parser.add_argument("--tool", action="append", help="row label (or tool) to run; repeatable")
    parser.add_argument("--skip", action="append", default=[], help="row label to leave out")
    parser.add_argument("--reps", type=int, default=1, help="summarised repetitions after the warm one")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--memory", default="6g")
    parser.add_argument("--no-prep", action="store_true", help="skip the Picard-made inputs")
    parser.add_argument("--out", help="JSON to write (rows already in it are kept)")
    parser.add_argument("--table", help="print the ranked table of this JSON")
    parser.add_argument("--compare", nargs=2)
    options = parser.parse_args(argv)
    if options.table:
        print(table(json.loads(Path(options.table).read_text())))
        return 0
    if options.compare:
        return 1 if compare(*options.compare) else 0
    options.derived = str(Path(options.derived).resolve())
    options.data = str(Path(options.data).resolve())
    if "/" in options.port_volume:
        options.port_volume = str(Path(options.port_volume).resolve())
    rows, prep = read_rows(options.rows)
    if not options.no_prep:
        run_prep(options, prep)
    out = Path(options.out) if options.out else None
    document = json.loads(out.read_text()) if out and out.exists() else {"rows": {}}
    document["settings"] = {"cpus": options.cpus, "memory": options.memory, "reps": options.reps,
                            "image": IMAGE, "arch": os.uname().machine}
    for label, tool, arguments in rows:
        if options.tool and label not in options.tool and tool not in options.tool:
            continue
        if label in options.skip:
            continue
        entry = summarise(time_row(options, tool, arguments))
        entry["command"] = f"{tool} {arguments}"
        document["rows"][label] = entry
        if "picard" in entry:
            print(f"{label:48s} cpu {entry['picard']['cpu']:8.1f} / {entry['port']['cpu']:8.1f} s  "
                  f"rss {entry['picard']['rss_mib']:7.0f} / {entry['port']['rss_mib']:7.0f} MiB  "
                  f"wall {entry['picard']['wall']:7.1f} / {entry['port']['wall']:7.1f} s  {verdict(entry)}",
                  flush=True)
        else:
            print(f"{label:48s} {entry.get('error', '')[:200]}", flush=True)
        if out:
            out.write_text(json.dumps(document, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
