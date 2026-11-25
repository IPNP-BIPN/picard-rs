#!/usr/bin/env python3
"""Time one row's two command lines inside the speed container. Run by `bench.py`, not by hand.

Ported from gatk-rs's tools/speed/timer.py. Every run is a child of this process, so `wait4`
returns exactly that child's CPU time (user + system, every thread) and peak resident set, and
each command is `exec`ed so the child IS the tool rather than a shell around it.

The two sides alternate, Picard first on even repetitions and the port first on odd ones, so a
drift in the host (a neighbour's load, a thermal step) lands on both rather than on whichever ran
last. Repetition zero warms the page cache for both and is reported but not summarised.

After the last port run the md5 of every file the port wrote is recorded, with the metrics
header's `# Started on:` line left out (it holds the wall-clock date), so a later run can say
whether an optimisation changed a byte.
"""

import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

OUT = Path("/work/out")


def prepare():
    """An empty output directory and a temporary one, as every row expects."""
    for path in (OUT, Path("/work/tmp"), Path("/work/w")):
        if path.exists():
            for stale in path.iterdir():
                if stale.is_dir() and not stale.is_symlink():
                    shutil.rmtree(stale, ignore_errors=True)
                else:
                    stale.unlink(missing_ok=True)
        path.mkdir(parents=True, exist_ok=True)


def run(command, timeout, setup=None):
    """One run: wall seconds, CPU seconds (user + system), peak RSS in KiB, exit code."""
    prepare()
    if setup:
        # Untimed: what the row needs in place before it starts (a writable copy of an input).
        subprocess.run(["bash", "-c", setup], check=True)
    start = time.perf_counter()
    child = subprocess.Popen(
        ["bash", "-c", f"cd {OUT} && exec {command}"],
        stdin=subprocess.DEVNULL,
        stdout=open(OUT.parent / "stdout.txt", "wb"),
        stderr=open(OUT.parent / "stderr.txt", "wb"),
        start_new_session=True,
    )
    timer = threading.Timer(timeout, lambda: os.killpg(child.pid, signal.SIGKILL))
    timer.start()
    _, status, usage = os.wait4(child.pid, 0)
    wall = time.perf_counter() - start
    timer.cancel()
    # Popen still holds the pid; tell it the child is gone so it does not wait again.
    child.returncode = os.waitstatus_to_exitcode(status)
    return {
        "wall": wall,
        "cpu": usage.ru_utime + usage.ru_stime,
        "rss_kib": usage.ru_maxrss,
        "exit": child.returncode,
        "timed_out": wall >= timeout,
    }


def digest(path):
    """md5 of a file, streamed, without the metrics header's date line in a small text file."""
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        head = f.read(4096)
        if path.stat().st_size < 1 << 26 and b"\0" not in head:
            for line in (head + f.read()).splitlines(keepends=True):
                if not line.startswith(b"# Started on:"):
                    md5.update(line)
            return md5.hexdigest()
        md5.update(head)
        for chunk in iter(lambda: f.read(1 << 20), b""):
            md5.update(chunk)
    return md5.hexdigest()


def outputs():
    """Every file under the output directory, with stdout, by its path relative to it."""
    found = {
        str(path.relative_to(OUT)): digest(path) for path in sorted(OUT.rglob("*")) if path.is_file()
    }
    found["<stdout>"] = digest(OUT.parent / "stdout.txt")
    return found


def tail(path, size=400):
    data = path.read_bytes()[-size:] if path.exists() else b""
    return data.decode("utf-8", "replace")


def main(spec_path):
    spec = json.loads(Path(spec_path).read_text())
    reps, timeout = spec["reps"], spec["timeout"]
    runs = {"oracle": [], "port": []}
    result = {}
    for rep in range(reps + 1):
        order = ("oracle", "port") if rep % 2 == 0 else ("port", "oracle")
        for side in order:
            if not spec.get(side):
                continue
            measured = run(spec[side], timeout, spec.get("setup"))
            runs[side].append(measured)
            if measured["exit"] != 0:
                result.setdefault("errors", {})[side] = tail(OUT.parent / "stderr.txt")
            if side == "port" and rep == reps:
                result["port_outputs"] = outputs()
    result["cold"] = runs
    result["host"] = {"cores": os.cpu_count()}
    print(json.dumps(result))


if __name__ == "__main__":
    main(sys.argv[1])
