#!/usr/bin/env python3
"""
Pre-commit guard: refuse to commit anything that looks like patient data.

Install as a git hook:
    ln -s ../../scripts/check_no_phi.py .git/hooks/pre-commit
    chmod +x scripts/check_no_phi.py

Or run manually over the whole tree:
    python scripts/check_no_phi.py --all
"""
import argparse
import os
import subprocess
import sys

# Extensions that are medical images or can carry them
BLOCKED_EXT = {
    ".dcm", ".dicom", ".ima", ".nii", ".gz", ".nrrd", ".nhdr", ".mha", ".mhd",
    ".img", ".hdr", ".par", ".rec", ".mnc", ".vtk", ".raw",
    ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp",
    ".zip", ".tar", ".tgz", ".7z", ".rar",
}
# Paths explicitly allowed despite the extension rules
ALLOWLIST_PREFIXES = ("docs/figures/",)

# DICOM files often have no extension at all; detect by magic bytes
DICOM_MAGIC_OFFSET = 128
DICOM_MAGIC = b"DICM"


def is_dicom(path):
    try:
        with open(path, "rb") as f:
            f.seek(DICOM_MAGIC_OFFSET)
            return f.read(4) == DICOM_MAGIC
    except Exception:
        return False


def staged_files():
    out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACM"],
                         capture_output=True, text=True)
    return [p for p in out.stdout.splitlines() if p.strip()]


def all_files():
    found = []
    for root, dirs, files in os.walk("."):
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__", ".venv", "venv"}]
        for fn in files:
            found.append(os.path.relpath(os.path.join(root, fn), "."))
    return found


def check(paths):
    problems = []
    for p in paths:
        norm = p.replace(os.sep, "/")
        if norm.startswith(ALLOWLIST_PREFIXES):
            continue
        ext = os.path.splitext(norm)[1].lower()
        if ext in BLOCKED_EXT:
            problems.append((norm, f"blocked extension '{ext}'"))
        elif os.path.isfile(norm) and is_dicom(norm):
            problems.append((norm, "file has a DICOM magic header"))
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true",
                    help="scan the whole working tree instead of the git staging area")
    a = ap.parse_args()

    paths = all_files() if a.all else staged_files()
    problems = check(paths)

    if problems:
        print("REFUSING TO PROCEED — possible patient data detected:\n", file=sys.stderr)
        for p, why in problems:
            print(f"  {p}\n      {why}", file=sys.stderr)
        print("\nMedical images are protected health information even after DICOM tags are "
              "stripped, and git history cannot be reliably scrubbed after a push.\n"
              "Remove these files, or add an explicit exception if you are certain they "
              "contain no patient data.", file=sys.stderr)
        return 1

    print(f"OK — {len(paths)} file(s) checked, no patient data detected.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
