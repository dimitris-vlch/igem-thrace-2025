#!/usr/bin/env python3
"""
count_ena_tsv.py

Utility script to count the number of entries in ENA-derived TSV files.
Each line is treated as an entry if the first field matches an accession-like
identifier (e.g. AY325607). Empty lines, comment lines, and headers are ignored.

Usage:
    python3 count_ena_tsv.py [directory]

By default, the current directory is used. The script scans for all files
with the .tsv extension, counts entries in each file, and reports per-file
counts as well as the total across all files.
"""

import argparse
import re
import sys
from pathlib import Path

# Regular expression to identify accession-like identifiers
ACCESSION_RE = re.compile(r'^[A-Z][A-Z0-9_.-]*\d$')

def count_file(path: Path) -> int:
    """Count valid entries in a single TSV file."""
    count = 0
    try:
        with path.open('r', encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith('#') or '\t' not in s:
                    continue
                first = s.split('\t', 1)[0]
                # Skip potential header lines
                if first.lower() in {'accession', 'id', 'accession_id'}:
                    continue
                if ACCESSION_RE.match(first):
                    count += 1
    except Exception as e:
        print(f"[Warning] Could not read {path}: {e}", file=sys.stderr)
    return count

def main():
    parser = argparse.ArgumentParser(
        description="Count the number of valid entries in TSV files containing ENA metadata."
    )
    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Directory containing the TSV files (default: current directory)."
    )
    parser.add_argument(
        "--glob",
        default="*.tsv",
        help="Glob pattern for matching TSV files (default: *.tsv)."
    )
    args = parser.parse_args()

    base = Path(args.directory)
    files = sorted(base.glob(args.glob))
    if not files:
        print(f"No files matching pattern {args.glob} found in {base}")
        sys.exit(2)

    total = 0
    name_width = max(len(p.name) for p in files)
    count_width = 0
    results = []

    for path in files:
        c = count_file(path)
        results.append((path.name, c))
        total += c
        count_width = max(count_width, len(str(c)))

    for name, c in results:
        print(f"{name.ljust(name_width)}  {str(c).rjust(count_width)}")
    print("-" * (name_width + count_width + 2))
    print(f"{'TOTAL'.ljust(name_width)}  {str(total).rjust(count_width)}")

if __name__ == "__main__":
    main()
