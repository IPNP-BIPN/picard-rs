#!/usr/bin/env python3
"""The port's side of StrictValidationDump: the same BAMs, the same runs, the same rows.

    python3 tools/strict-conformance/port_side.py target/release DUMP > port.txt
    diff <(grep -v '^#' DUMP) port.txt

DUMP is the dump's output (the CI artefact, or a local container run); its `bam` rows are the
inputs, written back to disk byte for byte, so both sides read the same files. Each run is the
binary named after the tool under the dump's arguments, and its row is the "Ignoring SAM
validation error" lines it printed, then what `Exception in thread "main"` named, or the exit.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

TOOLS = {
    "CollectQualityYieldMetrics": ("collect-quality-yield-metrics", ["OUTPUT={d}/out.cqym.txt"]),
    "SortSam": ("sort-sam", ["OUTPUT={d}/out.sorted.bam", "SORT_ORDER=queryname"]),
    "MarkDuplicates": (
        "mark-duplicates",
        ["OUTPUT={d}/out.md.bam", "METRICS_FILE={d}/out.md.txt"],
    ),
    "CleanSam": ("clean-sam", ["OUTPUT={d}/out.clean.bam"]),
}
PREFIX = 'Exception in thread "main" '


def esc(p):
    return p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")


def main():
    bins, dump = Path(sys.argv[1]), Path(sys.argv[2])
    rows = [line.rstrip("\n").split("\t", 2) for line in dump.read_text().splitlines()]
    with tempfile.TemporaryDirectory() as d:
        for kind, name, payload in rows:
            if kind == "bam":
                Path(d, name + ".bam").write_bytes(bytes.fromhex(payload))
                print(f"bam\t{name}\t{payload}")
                continue
            case, tool, stringency = name.split("/")
            binary, extra = TOOLS[tool]
            args = [str(bins / binary), f"INPUT={d}/{case}.bam"]
            args += [a.format(d=d) for a in extra] + [f"VALIDATION_STRINGENCY={stringency}"]
            run = subprocess.run(args, capture_output=True, text=True)
            lines = [l for l in run.stderr.splitlines() if l.startswith("Ignoring SAM validation error")]
            thrown = [l[len(PREFIX):] for l in run.stderr.splitlines() if l.startswith(PREFIX)]
            if thrown:
                lines.append("exception: " + thrown[0])
            elif run.returncode == 0:
                lines.append("ok")
            else:
                tail = run.stderr.strip().splitlines()[-1:] or [""]
                lines.append(f"exit {run.returncode}: {tail[0]}")
            print(f"run\t{name}\t{esc(chr(10).join(lines))}")


if __name__ == "__main__":
    main()
