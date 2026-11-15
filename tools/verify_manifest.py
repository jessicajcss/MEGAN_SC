#!/usr/bin/env python3
"""Check the published results against their checksums (results/MANIFEST.sha256).

Usage (from the repository root):
    python tools/verify_manifest.py            # check every file under results/ against the manifest
    python tools/verify_manifest.py --write    # maintainers only: rewrite the manifest from results/

The manifest lists the SHA-256 of every file under results/ (tables, numbers, figures, logs), so a copy of the
repository can be checked file by file. It answers "are these the published files?". The question "do I get the
same files when I recompute them?" is answered by `python reproduce.py all`, whose last step (`compare`,
tools/compare_results.py) compares the recomputed files with these published copies.

For each file the script prints one of:
    OK        identical to the published file (byte for byte), printed only with --verbose;
    CHANGED   the bytes differ from the published file;
    MISSING   the file listed in the manifest is not there;
    EXTRA     a file under results/ that the manifest does not list.
Exit status: 0 when every listed file is OK and there is no extra file, 1 otherwise.
"""
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
MANIFEST = RESULTS / "MANIFEST.sha256"


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            d.update(b)
    return d.hexdigest()


def files():
    return sorted(p for p in RESULTS.rglob("*") if p.is_file() and p != MANIFEST and not p.name.startswith("."))


def write():
    fs = files()
    lines = ["# SHA-256 of every published file under results/. Checked by tools/verify_manifest.py."]
    lines += [f"{sha256(p)}  {p.relative_to(ROOT).as_posix()}" for p in fs]
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"written {MANIFEST.relative_to(ROOT).as_posix()} ({len(fs)} files)")


def check(verbose=False):
    entries = [l.split(None, 1) for l in MANIFEST.read_text(encoding="utf-8").splitlines() if l and not l.startswith("#")]
    listed = {rel.strip() for _, rel in entries}
    counts = {"OK": 0, "CHANGED": 0, "MISSING": 0, "EXTRA": 0}
    for digest, rel in entries:
        rel = rel.strip()
        path = ROOT / rel
        status = "MISSING" if not path.is_file() else ("OK" if sha256(path) == digest else "CHANGED")
        counts[status] += 1
        if status != "OK" or verbose:
            print(f"{status:8s} {rel}")
    for p in files():
        rel = p.relative_to(ROOT).as_posix()
        if rel not in listed:
            counts["EXTRA"] += 1
            print(f"{'EXTRA':8s} {rel}")
    print(f"{len(entries)} files listed: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    return 0 if counts["OK"] == len(entries) and counts["EXTRA"] == 0 else 1


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--write" in args:
        write()
    else:
        sys.exit(check(verbose="--verbose" in args))
