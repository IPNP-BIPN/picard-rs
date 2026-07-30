#!/usr/bin/env python3
"""Turn `run_array.py`'s summary lines into a committable measurement.

The runner prints what it found and writes a corpus; nothing until now kept the *number*, so the
dashboard's argument-coverage column has been empty since the arrays were generated. This reads the
runner's own summary back and writes `tools/coverage/measured.json`, which the dashboard reads and
CI re-derives.

What is recorded, per tool:

* `rows`: the array's size at the declared strength;
* `accepted` / `rejected`: how the *reference* answered. A rejected row is not a failed run: the
  reference refuses `FLOW_MODE=true` and refuses a queryname-sorted input without `ASSUME_SORTED`,
  and both are behaviour a port has to reproduce, message included;
* `distinct_outputs`: how many different outputs the accepted rows produced. This is the number
  that says whether the array is testing anything. An array whose accepted rows all produce one
  output covers its argument pairs without observing them, and the corpus cannot tell a port that
  implements those arguments from one that ignores them;
* `matched`: rows where the port answered exactly as the reference did, rejections included;
* `share`: `matched / rows`, which is what the dashboard prints.

The parse is deliberately literal: it reads the runner's printed summary rather than recomputing
anything, so the committed number and the number a reader sees in the log cannot drift apart.

Usage: measure.py --log run.log [--log run2.log ...] --out tools/coverage/measured.json
"""

import argparse
import json
import re
import sys
from pathlib import Path

# `<tool>: rows=16 accepted by the reference=12 rejected=4 distinct outputs=3 distinct rejections=2`
SUMMARY = re.compile(
    r"^(?P<tool>\w+): rows=(?P<rows>\d+) accepted by the reference=(?P<accepted>\d+) "
    r"rejected=(?P<rejected>\d+) distinct outputs=(?P<distinct>\d+) "
    r"distinct rejections=(?P<rejections>\d+)\s*$"
)
# `  port matches the reference on 12/16 rows (75% of the pairwise surface)`
PORT = re.compile(r"^\s*port matches the reference on (?P<matched>\d+)/(?P<rows>\d+) rows")
# `# <tool> t=2 rows=16 from <file>`, the corpus header, for the strength.
STRENGTH = re.compile(r"^# (?P<tool>\w+) t=(?P<t>\d+) rows=\d+")


def parse(text):
    """One log's measurements, keyed by tool, in the order the runner reported them."""
    measured = {}
    current = None
    for line in text.splitlines():
        summary = SUMMARY.match(line)
        if summary:
            current = summary.group("tool")
            measured[current] = {
                "rows": int(summary.group("rows")),
                "accepted": int(summary.group("accepted")),
                "rejected": int(summary.group("rejected")),
                "distinct_outputs": int(summary.group("distinct")),
                "distinct_rejections": int(summary.group("rejections")),
            }
            continue
        port = PORT.match(line)
        if port and current:
            measured[current]["matched"] = int(port.group("matched"))
    return measured


def strength_of(tool, corpus_dir):
    """The `t` the corpus was produced at, read from the corpus the same run wrote."""
    for path in sorted(corpus_dir.glob(f"{tool}.t*.txt")):
        header = path.open().readline()
        match = STRENGTH.match(header)
        if match and match.group("tool") == tool:
            return int(match.group("t"))
    return None


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--log", action="append", required=True, help="a run_array.py log")
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--corpus",
        default=str(Path(__file__).resolve().parent / "corpus"),
        help="where the corpora the same run wrote live",
    )
    args = ap.parse_args(argv)

    corpus_dir = Path(args.corpus)
    tools = {}
    for log in args.log:
        for tool, entry in parse(Path(log).read_text()).items():
            if "matched" not in entry:
                raise SystemExit(
                    f"{log}: {tool} was run without --port, so nothing about the port was "
                    "measured; the column this feeds would be a claim with no measurement behind it"
                )
            entry["t"] = strength_of(tool, corpus_dir)
            entry["share"] = entry["matched"] / entry["rows"] if entry["rows"] else 0.0
            tools[tool] = entry

    if not tools:
        raise SystemExit("no summary line in any log; run_array.py prints one per tool")

    out = {
        "$comment": [
            "Produced by tools/coverage/measure.py from run_array.py's own summary, in the",
            "pinned container on real x86-64. `share` is the fraction of the array's rows on",
            "which the port answered exactly as the reference did, rejections included.",
            "",
            "It is not a byte-identity claim over the argument surface: an array covers pairs of",
            "argument values, and `distinct_outputs` says how many of those pairs the corpus can",
            "actually observe. A tool whose accepted rows all produce one output has an array",
            "that covers its arguments without testing them.",
        ],
        "tools": dict(sorted(tools.items())),
    }
    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
    for tool, entry in sorted(tools.items()):
        print(
            f"{tool}: t={entry['t']} rows={entry['rows']} matched={entry['matched']} "
            f"({entry['share']:.0%}) distinct outputs={entry['distinct_outputs']}"
        )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
