#!/usr/bin/env python3
"""Normalize raw ERVmap output into the common per-sample locus table.

Raw file: <sample>ERVresults.txt, produced by
  bedtools coverage -nonamecheck -a ERVmap.sorted.bed -b <bam> -counts -sorted
Columns: chr, start, end, ervmap_id, score, strand, count  (7 columns, one row
per provirus, rows in BED order).

The subfamily/family of each provirus comes from annotate_ervmap.py, which
matches it against the RepeatMasker ERV fragments once per reference build.

Usage: normalize_ervmap.py RAW_ERVresults.txt ERVMAP_ANNOTATION OUT_TSV

Output columns: locus_id, chr, start, end, strand, subfamily, family, count
locus_id is "ERVmap_<id>" so provirus rows never collide with RepeatMasker
locus ids ("<subfamily>__<chr>_<start>_<end>_<strand>").
"""
import csv
import sys

raw_in, annot_in, out_tsv = sys.argv[1:4]

annot = {r["ervmap_id"]: r for r in csv.DictReader(open(annot_in), delimiter="\t")}

n_rows = n_pos = 0
with open(raw_in) as fh, open(out_tsv, "w") as fo:
    fo.write("locus_id\tchr\tstart\tend\tstrand\tsubfamily\tfamily\tcount\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        chrom, start, end, eid, _score, strand, count = f[:7]
        a = annot.get(eid, {})
        fo.write("\t".join([
            f"ERVmap_{eid}", chrom, start, end,
            a.get("strand") or strand,
            a.get("subfamily", "."), a.get("family", "."), count,
        ]) + "\n")
        n_rows += 1
        n_pos += float(count) > 0

print(f"rows={n_rows} positive={n_pos}", file=sys.stderr)
