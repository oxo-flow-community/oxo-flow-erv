#!/usr/bin/env python3
"""Fail early when the ERVmap BED does not use the same contig names as the BAM.

`bedtools coverage -nonamecheck` turns a chromosome-name mismatch into a valid,
complete, all-zero result with exit code 0.  That is exactly how a whole
cohort's ERVmap output once became 3,220 zero rows (the official BED has no
`chr` prefix, the STAR index does), and the run looked successful.

This check runs before counting: every contig in the BED must exist in the
genome FASTA index.  A mismatch is a hard error listing both sides, so the
cause is visible instead of being buried in a QC flag.

Usage: check_bed_contigs.py ERVMAP_BED GENOME_FAI
"""
import sys

bed, fai = sys.argv[1:3]

fai_contigs = []
with open(fai) as fh:
    for line in fh:
        if line.strip():
            fai_contigs.append(line.split("\t")[0])
known = set(fai_contigs)

bed_contigs = []
with open(bed) as fh:
    for line in fh:
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        bed_contigs.append(line.split("\t")[0])

missing = sorted({c for c in bed_contigs if c not in known})
if missing:
    sys.exit(
        "ERVmap BED contig names do not match the genome index.\n"
        f"  BED contigs not in {fai}: {missing[:10]}\n"
        f"  genome index contigs (first 10): {fai_contigs[:10]}\n"
        "  The official ERVmap BED has no 'chr' prefix while a UCSC/GRCh38 STAR\n"
        "  index does; with -nonamecheck this silently yields all-zero counts.\n"
        "  Fix the BED (add the prefix) instead of changing the counting step."
    )

print(f"bed_rows={len(bed_contigs)} contigs={len(set(bed_contigs))} all_in_index=yes",
      file=sys.stderr)
