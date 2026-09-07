#!/usr/bin/env python3
"""Generate the oracle jobs of `.github/workflows/ci.yml` from `manifest.json`.

The workflow used to hold one hand-written pair of steps per suite: a `docker run` that
regenerated the golden and a thirty-line Python comparator inlined in YAML. At 39 suites that was
1556 lines, the comparators had already drifted from each other, and each new tool made the file
longer than the port it was testing. At 311 tools it would have been unreadable.

Now the suites are data. This script turns them into:

  * `oracle-image`, which builds the digest-pinned image once and caches it, and
  * `oracle`, a matrix with one entry per suite, so suites run in parallel and a broken suite no
    longer hides the ones after it in a single sequential job.

Usage:
    python3 tools/conformance/generate_ci.py            # write .github/workflows/ci.yml
    python3 tools/conformance/generate_ci.py --check    # fail if the committed file is stale
"""

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TEMPLATE = HERE / "ci.template.yml"
WORKFLOW = REPO / ".github" / "workflows" / "ci.yml"
MARKER = "# {{ORACLE_JOBS}}"
COVERAGE_MARKER = "# {{COVERAGE_JOB}}"
GITIGNORE = REPO / ".gitignore"
CORPUS_BEGIN = "# {{BEGIN coverage corpora}}"
CORPUS_END = "# {{END coverage corpora}}"

CACHE_KEY = "picard-rs-oracle-${{ hashFiles('tools/oracle/**') }}"
IMAGE_TAR = "/tmp/oracle-image.tar"


def wrap(text, indent, width=96):
    """Wrap a `why` into YAML comment lines, so the reason stays next to the job that needs it."""
    words = text.split()
    lines, current = [], indent + "#"
    for word in words:
        if len(current) + 1 + len(word) > width:
            lines.append(current)
            current = indent + "#"
        current += " " + word
    lines.append(current)
    return "\n".join(lines)


def oracle_jobs(manifest):
    suites = manifest["suites"]
    probes = manifest.get("probes", [])
    oracle = manifest["oracle"]
    out = []

    out.append(
        f"""  oracle-image:
    name: Build the pinned oracle image
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Record the host CPU
        run: grep -m1 'model name' /proc/cpuinfo; grep -m1 'vendor_id' /proc/cpuinfo

      - id: cache
        uses: actions/cache@v4
        with:
          path: {IMAGE_TAR}
          key: {CACHE_KEY}

      - name: Build the oracle image
        # The probe runs during the build, so a degraded environment cannot produce an image, let
        # alone a golden. The image is built once and shared with every suite below, because
        # rebuilding it per suite would multiply a five-minute build by the number of suites.
        if: steps.cache.outputs.cache-hit != 'true'
        run: |
          docker build --platform {oracle['platform']} -t {oracle['image']} {oracle['context']}
          docker save {oracle['image']} -o {IMAGE_TAR}
"""
    )

    per_job = manifest.get("ci", {}).get("suites_per_job", 5)
    rows = []
    # Grouped by status, so a green run never lets an `unchecked` suite read as the stronger
    # claim, and chunked, because every job pays the oracle image restore once: one job per suite
    # would download the image forty-odd times per run for no extra information.
    for status in ("oracle-backed", "unchecked", "golden-pending"):
        group = [s for s in suites if s["status"] == status]
        chunks = [group[i : i + per_job] for i in range(0, len(group), per_job)]
        for n, chunk in enumerate(chunks, 1):
            ids = " ".join(s["id"] for s in chunk)
            rows.append(
                f"          - group: {status}\n"
                f"            index: \"{n}/{len(chunks)}\"\n"
                f"            suites: \"{ids}\""
            )
    matrix = "\n".join(rows)

    out.append(
        f"""  oracle:
    name: "oracle ${{{{ matrix.group }}}} ${{{{ matrix.index }}}}"
    needs: oracle-image
    runs-on: ubuntu-latest
    strategy:
      # One suite's divergence must not cancel the others: the point of the run is to learn which
      # suites diverge, not to learn that at least one does.
      fail-fast: false
      matrix:
        include:
{matrix}
    steps:
      - uses: actions/checkout@v4

      - uses: actions/cache/restore@v4
        with:
          path: {IMAGE_TAR}
          key: {CACHE_KEY}
          fail-on-cache-miss: true

      - run: docker load -i {IMAGE_TAR}

      - name: Regenerate and compare
        # Every rule that canonicalizes anything away is declared in tools/conformance/manifest.json
        # with the reason it is there; tools/conformance/compare.py is the only code that applies
        # them. Run the same thing locally with:
        #   python3 tools/conformance/run_suite.py --suites "${{{{ matrix.suites }}}}"
        run: python3 tools/conformance/run_suite.py --suites "${{{{ matrix.suites }}}}"

      - name: Publish any candidate goldens
        # A `golden-pending` suite has no golden to compare against, and one may not be committed
        # from a developer machine (decision 0008). The job publishes what the pinned container
        # produced on real x86-64; committing that artefact is what makes the suite oracle-backed.
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: candidate-goldens-${{{{ matrix.group }}}}-${{{{ strategy.job-index }}}}
          path: tools/conformance/pending/
          if-no-files-found: ignore
"""
    )

    if probes:
        probe_rows = "\n".join(f"          - probe: {p['id']}" for p in probes)
        why_lines = "\n".join(wrap(p["why"], "        ") for p in probes)
        out.append(
            f"""  oracle-probes:
    name: "probe: ${{{{ matrix.probe }}}}"
    needs: oracle-image
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        include:
{probe_rows}
    steps:
      - uses: actions/checkout@v4

      - uses: actions/cache/restore@v4
        with:
          path: {IMAGE_TAR}
          key: {CACHE_KEY}
          fail-on-cache-miss: true

      - run: docker load -i {IMAGE_TAR}

      - name: Run the probe
{why_lines}
        run: python3 tools/conformance/run_suite.py --probe ${{{{ matrix.probe }}}}
"""
        )

    return "\n".join(out)


def coverage_job(manifest):
    """The covering-array job, from the one list of tools in the manifest.

    What it emits: a build of exactly the declared port binaries, one array run per tool, one
    measurement over all of their logs, and a dirty check over the artefacts the `measured` tools
    own. How: every step is derived from the same list, so the four places that used to name a
    tool cannot disagree about which tools are covered. Why not keep the steps hand-written: they
    had to agree by hand, and this is the job whose output is the programme's coverage claim, so a
    list that has silently lost a tool is a claim that reads higher than what ran.
    """
    coverage = manifest.get("coverage")
    if not coverage:
        return ""
    tools = coverage["tools"]
    measured = [t for t in tools if t["status"] == "measured"]
    pending = [t for t in tools if t["status"] == "pending"]

    bins = " ".join(f"--bin {t['port']}" for t in tools)
    runs = "\n".join(
        f"          python3 tools/coverage/run_array.py --tool {t['tool']} --t {t['t']} \\\n"
        # A tool with no output argument is compared on its standard output; see run_array.py.
        f"            --port target/release/{t['port']}"
        f"{' --stdout' if t.get('output') == 'stdout' else ''}"
        # A tool that stamps a @PG carrying its command line: see run_array.py.
        f"{' --strip-program-records' if t.get('strip_program_records') else ''}"
        # A tool whose exit code says what it FOUND rather than that it failed.
        f"{' --exit-code-is-a-result' if t.get('exit_code_is_a_result') else ''} \\\n"
        f"            | tee /tmp/{t['tool']}.log"
        for t in tools
    )
    logs = " \\\n            ".join(f"--log /tmp/{t['tool']}.log" for t in tools)
    pending_flags = "".join(f" \\\n            --pending {t['tool']}" for t in pending)
    dirty = " \\\n            ".join(
        ["dirty=$(git status --porcelain tools/coverage/measured.json"]
        + [f"tools/coverage/corpus/{t['tool']}.t{t['t']}.txt" for t in measured]
    )

    return f"""  coverage:
    name: "Argument coverage: the covering array against both sides"
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - uses: dtolnay/rust-toolchain@master
        with:
          toolchain: "1.97.1"

      - uses: Swatinem/rust-cache@v2

      - name: Build the pinned oracle image
        run: docker build --platform {manifest['oracle']['platform']} -t {manifest['oracle']['image']} {manifest['oracle']['context']}

      - run: cargo build --release {bins}

      - name: Run the arrays, and record what the port covers
        # Every row goes through the reference and through the port, compared under the same
        # canonicalization the conformance suites use. A rejected row counts as covered when both
        # sides reject it with the same message, because reproducing a refusal is reproducing the
        # tool.
        #
        # The number that comes out is not a bug count. It is how much of the pairwise argument
        # surface the port implements, which is the figure the programme had never had until this
        # job ran: a binary written for one path answers a row that varies another path
        # differently, and that difference is the measurement, not a failure.
        run: |
{runs}
          python3 tools/coverage/measure.py {logs}{pending_flags} \\
            --out tools/coverage/measured.json \\
            --pending-out tools/coverage/pending/measured.json

      - name: The committed measurement still holds
        # The corpus and the measurement are goldens like any other: produced here, committed, and
        # re-derived on every run. A port that grows its argument surface makes this fail, which is
        # the point, and the fix is to commit the new number.
        #
        # A `pending` tool is deliberately absent from both lists: its numbers exist only in the
        # log and the published artefact until someone commits them, exactly as a golden-pending
        # suite's golden does.
        run: |
          {dirty})
          if [ -n "$dirty" ]; then
            echo "the measured argument coverage moved:"
            echo "$dirty"
            git --no-pager diff -- tools/coverage/measured.json | head -40
            echo "commit the artefacts this run produced"
            exit 1
          fi
          echo "the committed measurement matches what this run produced"

      - name: Publish the corpus and the measurement
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: covering-corpus
          path: |
            tools/coverage/corpus/
            tools/coverage/measured.json
            tools/coverage/pending/
          if-no-files-found: ignore
"""


def corpus_whitelist(manifest):
    """The .gitignore negations that keep exactly the committed corpora committed.

    Corpora are oracle output, so `tools/coverage/corpus/*` is excluded and each committed file is
    negated back in. That list is the same list as the `measured` tools, and it lived in a second
    place where it could fall behind: a corpus whose negation is missing simply stops being
    committed, and the dirty check above then fails on a file nobody deleted.
    """
    measured = [t for t in manifest.get("coverage", {}).get("tools", []) if t["status"] == "measured"]
    lines = [f"!tools/coverage/corpus/{t['tool']}.t{t['t']}.txt" for t in measured]
    return "\n".join(lines)


def render_gitignore(manifest, current):
    """The committed .gitignore with its managed corpus block replaced."""
    begin, end = current.index(CORPUS_BEGIN), current.index(CORPUS_END)
    head = current[: begin + len(CORPUS_BEGIN)]
    tail = current[end:]
    return head + "\n" + corpus_whitelist(manifest) + "\n" + tail


# `uses: owner/action@ref`, over the two shapes the workflow writes: a step that is only a `uses`,
# and one that names itself first.
USES = re.compile(r"(?P<lead>uses: )(?P<action>[\w.-]+/[\w./-]+)@(?P<ref>\S+)")


def repin(rendered, current):
    """Take every action ref from the committed workflow rather than from the template.

    The generator owns the job matrix. It does not own which version of `actions/checkout` runs,
    because `.github/workflows/` is the only place Dependabot can write, and a ref hard-coded here
    would make every grouped actions bump fail this file's own `--check`.

    An action pinned to two different refs in the committed file is left alone, so that
    disagreement reaches the guard as a diff instead of being quietly resolved to one of them.
    """
    seen = {}
    for match in USES.finditer(current):
        seen.setdefault(match.group("action"), set()).add(match.group("ref"))
    pins = {action: refs.pop() for action, refs in seen.items() if len(refs) == 1}
    return USES.sub(
        lambda m: m.group("lead") + m.group("action") + "@"
        + pins.get(m.group("action"), m.group("ref")),
        rendered,
    )


GATE = "ci"

# A job id: two spaces, a name, a colon, nothing else on the line. Job keys sit at four spaces, so
# this cannot match one of them.
JOB_ID = re.compile(r"^  (?P<id>[a-z][\w-]*):$", re.M)


def gate_job(rendered):
    """Append the one check that stands for the whole run.

    The oracle matrix names a check per suite, so a required-checks list written by hand would go
    stale the moment a suite is added, and a branch rule built on it would be guarding a set of
    jobs that no longer exists. This job needs every other one and is generated alongside them,
    which is what keeps the list honest with nobody maintaining it.
    """
    body = rendered.split("\njobs:\n", 1)[1]
    ids = [match.group("id") for match in JOB_ID.finditer(body) if match.group("id") != GATE]
    needs = "\n".join(f"      - {job}" for job in ids)
    return rendered.rstrip("\n") + f"""

  {GATE}:
    name: CI complete
    # The single context main's ruleset requires, and the single context an automatic merge waits
    # on. It runs whatever the jobs above did, because a gate that is skipped when something fails
    # reports nothing and blocks nothing.
    needs:
{needs}
    if: always()
    runs-on: ubuntu-latest
    steps:
      - name: Every job above succeeded
        # A skipped job is not a failure. A cancelled one is: a cancelled run has proved nothing,
        # and `success` on an empty result is how a green tick comes to mean less than it looks.
        env:
          RESULTS: ${{{{ toJSON(needs) }}}}
        run: |
          python3 - <<'PY'
          import json
          import os
          import sys

          results = json.loads(os.environ["RESULTS"])
          bad = sorted(
              name for name, job in results.items()
              if job["result"] not in ("success", "skipped")
          )
          for name in bad:
              print(f"{{name}}: {{results[name]['result']}}")
          print(f"{{len(results) - len(bad)}}/{{len(results)}} jobs green")
          sys.exit(1 if bad else 0)
          PY
"""


def render(manifest):
    template = TEMPLATE.read_text()
    if MARKER not in template:
        raise SystemExit(f"{TEMPLATE} has lost its {MARKER} marker")
    banner = (
        "# The jobs below this line are generated from tools/conformance/manifest.json by\n"
        "# tools/conformance/generate_ci.py. Edit the manifest, not this file: the `guard` job\n"
        "# fails if the two disagree.\n"
    )
    if COVERAGE_MARKER not in template:
        raise SystemExit(f"{TEMPLATE} has lost its {COVERAGE_MARKER} marker")
    rendered = template.replace(MARKER, banner + oracle_jobs(manifest))
    return rendered.replace(COVERAGE_MARKER, coverage_job(manifest).rstrip("\n"))


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="fail if the committed workflow is stale")
    args = ap.parse_args(argv)

    manifest = comparator.load_manifest()
    current = WORKFLOW.read_text() if WORKFLOW.exists() else ""
    rendered = gate_job(repin(render(manifest), current))
    ignore_current = GITIGNORE.read_text()
    ignore_rendered = render_gitignore(manifest, ignore_current)

    if args.check:
        if ignore_current != ignore_rendered:
            print(
                ".gitignore's corpus whitelist is stale: regenerate with "
                "tools/conformance/generate_ci.py"
            )
            return 1
        if current != rendered:
            print("ci.yml is stale: regenerate with tools/conformance/generate_ci.py")
            import difflib

            diff = difflib.unified_diff(
                current.splitlines(), rendered.splitlines(), "committed", "generated", lineterm=""
            )
            for line in list(diff)[:60]:
                print(line)
            return 1
        print(
            f"ci.yml matches the manifest: {len(manifest['suites'])} suites, "
            f"{sum(len(s['cases']) for s in manifest['suites'])} cases, "
            f"{len(manifest.get('probes', []))} probe(s), "
            f"{len(manifest.get('coverage', {}).get('tools', []))} covering array(s)"
        )
        return 0

    WORKFLOW.write_text(rendered)
    print(f"wrote {WORKFLOW} ({rendered.count(chr(10)) + 1} lines)")
    if ignore_current != ignore_rendered:
        GITIGNORE.write_text(ignore_rendered)
        print(f"wrote {GITIGNORE} (corpus whitelist)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
