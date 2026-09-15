#!/usr/bin/env python3
"""Runs every test_analysis_*.py in this directory and adds up the tallies."""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
passed = total = 0
failed = []
for path in sorted(HERE.glob("test_analysis_*.py")):
    result = subprocess.run([sys.executable, str(path)], capture_output=True, text=True)
    print(result.stdout, end="")
    tally = re.search(r"(\d+)/(\d+) tests passed", result.stdout)
    if tally:
        passed += int(tally.group(1))
        total += int(tally.group(2))
    if result.returncode != 0:
        failed.append(path.name)
        print(result.stderr, end="", file=sys.stderr)

print(f"\n{passed}/{total} tests passed across {len(list(HERE.glob('test_analysis_*.py')))} files")
if failed:
    print("  failing files: " + ", ".join(failed))
    raise SystemExit(1)
