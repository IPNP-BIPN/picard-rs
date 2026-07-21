#!/usr/bin/env python3
"""Classify Picard's metrics collectors by the machinery each one needs.

The calibration gate rests on the claim that members of an archetype cost a small delta after
the first. Decision 0001 measured a 3.1x size spread among three members of the metrics
archetype and concluded that a delta from one arbitrary pair would be an artefact of the pair.

This stratifies instead. It is deliberately mechanical: every signal is a symbol present in the
pinned Picard source, not a judgement about how hard something looks. A stratum is the set of
machinery a tool needs, and the delta claim is only made *within* a stratum.

Usage:  python3 tools/stratify/stratify.py picard/src/main/java/picard
"""

import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Each signal is (name, regex). Presence of the symbol in the file is the whole test.
SIGNALS = [
    # Traversal shape.
    ("single_pass", re.compile(r"extends\s+SinglePassSamProgram")),
    ("multi_level", re.compile(r"MultiLevelCollector|SAMRecordMultiLevelCollector")),
    # Accumulation machinery.
    ("histogram", re.compile(r"\bHistogram\s*<")),
    ("mergeable_base", re.compile(r"extends\s+MergeableMetricBase")),
    # External requirements, which are oracle-environment cost rather than porting cost.
    ("needs_r", re.compile(r"RExecutor|checkRInstallation")),
    ("needs_reference", re.compile(r"REFERENCE_SEQUENCE|ReferenceSequenceFile")),
    ("needs_intervals", re.compile(r"IntervalList|TARGET_INTERVALS")),
]

# A file is a metrics tool if it is a CommandLineProgram that writes a MetricsFile.
IS_TOOL = re.compile(r"@CommandLineProgramProperties")
WRITES_METRICS = re.compile(r"MetricsFile")


# A tool's machinery is often in a separate class it delegates to: CollectInsertSizeMetrics
# looks like a plain single-pass tool, and its InsertSizeMetricsCollector extends
# MultiLevelCollector. Scanning one file per tool therefore *understates* the machinery, which
# would put tools in the wrong stratum and flatter any delta measured within it. So references
# to collector and metrics classes are followed one level and their signals unioned in.
REFERENCED = re.compile(r"\b([A-Z][A-Za-z0-9]*(?:Collector|Metrics))\b")


def signals_of(text: str) -> set:
    return {n for n, rx in SIGNALS if rx.search(text)}


def main(root: Path) -> int:
    # Index every class in the tree so references can be resolved by name.
    by_class = {}
    for path in sorted(root.rglob("*.java")):
        by_class[path.stem] = path.read_text(errors="replace")

    tools = {}
    for path in sorted(root.rglob("*.java")):
        text = by_class[path.stem]
        if not IS_TOOL.search(text) or not WRITES_METRICS.search(text):
            continue
        name = path.stem
        present = signals_of(text)
        followed = []
        for ref in sorted(set(REFERENCED.findall(text))):
            if ref == name or ref not in by_class:
                continue
            gained = signals_of(by_class[ref]) - present
            if gained:
                followed.append(f"{ref}:{'+'.join(sorted(gained))}")
            present |= gained
        # The tool file is a thin dispatcher; the work is in the collector and bean classes it
        # references. Measuring the tool file measures the wrong thing: CollectInsertSizeMetrics
        # is 191 lines and its collector plus bean are another 335. Both numbers are kept, and
        # `footprint` is the one to size work with.
        footprint = len(text.splitlines())
        parts = []
        for ref in sorted(set(REFERENCED.findall(text))):
            if ref == name or ref not in by_class:
                continue
            n = len(by_class[ref].splitlines())
            footprint += n
            parts.append(f"{ref}({n})")
        tools[name] = {
            "path": str(path),
            "lines": len(text.splitlines()),
            "footprint": footprint,
            "footprint_parts": parts,
            "signals": sorted(present),
            "machinery_from": followed,
        }

    strata = defaultdict(list)
    for name, info in tools.items():
        strata[tuple(info["signals"])].append(name)

    print(f"{len(tools)} metrics tools in {len(strata)} strata\n")
    print(f"{'n':>3}  {'footprint':>13}  signals")
    print("-" * 78)
    for signals, names in sorted(strata.items(), key=lambda kv: -len(kv[1])):
        sizes = [tools[n]["footprint"] for n in names]
        spread = f"{min(sizes)}-{max(sizes)}"
        label = ", ".join(signals) or "(none)"
        print(f"{len(names):>3}  {spread:>13}  {label}")
        for n in sorted(names):
            print(f"                       {n} ({tools[n]['footprint']})")
        print()

    counts = Counter()
    for info in tools.values():
        for s in info["signals"]:
            counts[s] += 1
    print("signal frequency across all metrics tools:")
    for s, c in counts.most_common():
        print(f"  {c:>3}  {s}")

    out = Path("tools/stratify/strata.json")
    out.write_text(json.dumps({"tools": tools}, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1] if len(sys.argv) > 1 else "picard/src/main/java/picard")))
