#!/usr/bin/env python3
"""Generate SQuIRE's clean_folder from the locus map (no Build/Clean run).

SQuIRE's `Count` only needs two files inside --clean_folder:
  <name>.bed          9 columns: chrom start end TE_ID milliDiv strand start end rgb
                      TE_ID = chrom|start|end|Subfamily:Family:Class|milliDiv|strand
  <name>_copies.txt   4 columns: Subfamily:Family:Class Copies Tot.Length Avg.Length

Why not `squire Build`: Build downloads the current UCSC rmsk over the network
and merges adjacent fragments with bedtools.  This workflow wants the opposite —
the same unmerged fragment coordinates the other tools see — so the folder is
written here instead.

Count uses the FULL repeat set, not just ERVs: with an ERV-only reference the
multi-mapper bookkeeping loses its anchor and counts of high-copy subfamilies
collapse.  Call this script with the full-TE map for the real run.

Usage: make_squire_clean.py LOCUS_MAP OUTDIR
"""
import csv
import os
import sys
from collections import defaultdict

map_in, outdir = sys.argv[1:3]
os.makedirs(outdir, exist_ok=True)

copies = defaultdict(int)
lengths = defaultdict(int)

with open(os.path.join(outdir, "hg38_all.bed"), "w") as fb:
    for r in csv.DictReader(open(map_in), delimiter="\t"):
        start, end = int(r["start"]), int(r["end"])
        sub = f'{r["subfamily"]}:{r["family"]}:{r["rep_class"]}'
        md = r.get("milli_div") or "0"
        te_id = f'{r["chr"]}|{start}|{end}|{sub}|{md}|{r["strand"]}'
        fb.write("\t".join([r["chr"], str(start), str(end), te_id, md,
                            r["strand"], str(start), str(end), "94,137,255"]) + "\n")
        copies[sub] += 1
        lengths[sub] += end - start

with open(os.path.join(outdir, "hg38_all_copies.txt"), "w") as fc:
    fc.write("Subfamily:Family:Class\tCopies\tTot.Length\tAvg.Length\n")
    for sub in sorted(copies):
        fc.write(f"{sub}\t{copies[sub]}\t{lengths[sub]}\t"
                 f"{lengths[sub] / copies[sub]:.2f}\n")

print(f"bed_rows={sum(copies.values())} subfamilies={len(copies)}", file=sys.stderr)
