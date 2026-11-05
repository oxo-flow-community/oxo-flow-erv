#!/usr/bin/env python3
"""Assign each ERVmap provirus a RepeatMasker ERV fragment.

ERVmap's 3,220 proviruses are curated integration sites; the other three tools
count RepeatMasker fragments.  To put ERVmap on the common annotation each
provirus is matched to the ERV fragment with the largest overlap (ties broken
by subfamily name so the table is deterministic).  Proviruses that overlap no
ERV fragment keep "." and are counted at provirus level only.

Usage: annotate_ervmap.py ERVMAP_BED LOCUS_MAP OUT_TSV

Output columns:
  ervmap_id, chr, start, end, strand, subfamily, family, region_type,
  overlap_bp, n_candidates
"""
import bisect
import csv
import sys
from collections import defaultdict

bed_in, map_in, out_tsv = sys.argv[1:4]

by_chrom = defaultdict(list)
for r in csv.DictReader(open(map_in), delimiter="\t"):
    by_chrom[r["chr"]].append((int(r["start"]), int(r["end"]), r))
for v in by_chrom.values():
    v.sort(key=lambda t: t[0])
starts = {c: [t[0] for t in v] for c, v in by_chrom.items()}

n_hit = 0
with open(out_tsv, "w") as fo:
    fo.write("ervmap_id\tchr\tstart\tend\tstrand\tsubfamily\tfamily\t"
             "region_type\toverlap_bp\tn_candidates\n")
    for line in open(bed_in):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        chrom, start, end, eid, _score, strand = f[:6]
        start, end = int(start), int(end)
        loci = by_chrom.get(chrom, [])
        # candidates are the loci that start before this provirus ends; the
        # bisect keeps the scan proportional to the local locus density
        hi = bisect.bisect_left(starts.get(chrom, []), end)
        best = None
        n_cand = 0
        for s, e, r in loci[:hi]:
            if e <= start:
                continue
            n_cand += 1
            ov = min(e, end) - max(s, start)
            key = (ov, r["subfamily"])
            if best is None or key > best[0]:
                best = (key, r, ov)
        if best is None:
            fo.write(f"{eid}\t{chrom}\t{start}\t{end}\t{strand}\t.\t.\t.\t0\t0\n")
        else:
            _, r, ov = best
            n_hit += 1
            fo.write(f"{eid}\t{chrom}\t{start}\t{end}\t{strand}\t"
                     f'{r["subfamily"]}\t{r["family"]}\t{r["region_type"]}\t{ov}\t{n_cand}\n')

print(f"proviruses_with_erv_annotation={n_hit}", file=sys.stderr)
