#!/usr/bin/env python3
"""Merge the per-sample long tables into the workflow products.

Usage: merge_expression.py PER_SAMPLE_DIR OUTDIR MIN_COUNT

PER_SAMPLE_DIR holds one directory per sample group (cohort) and one
<sample>.long.tsv.gz inside each, i.e. the source project's
{tool}/{cohort}/{sample} layout, with `group` carrying the cohort name.

Outputs (all under OUTDIR/aggregate/):
  erv_expression_long.tsv.gz     every per-sample row, concatenated
  qc_per_sample.tsv              per tool x group x sample totals, positive
                                 loci and a flag
  erv_subfamily_evidence.tsv.gz  per group x subfamily: how many of the three
                                 RepeatMasker-based tools reach MIN_COUNT in any
                                 sample, and whether ERVmap does

The QC table exists so a tool that silently produced nothing is visible instead
of looking like a clean zero.  ZERO_COUNTS is the signature of the chromosome
naming mismatch that once turned a whole cohort's ERVmap output into 3,220 zero
rows.
"""
import collections
import glob
import gzip
import os
import sys

per_dir, outdir, min_count = sys.argv[1:4]
min_count = float(min_count)

os.makedirs(f"{outdir}/aggregate", exist_ok=True)
files = sorted(glob.glob(f"{per_dir}/*/*.long.tsv.gz"))
if not files:
    sys.exit(f"no per-sample tables under {per_dir}")

long_out = f"{outdir}/aggregate/erv_expression_long.tsv.gz"
qc = collections.defaultdict(lambda: [0.0, 0])
sup = collections.defaultdict(set)          # (group, subfamily) -> tools

with gzip.open(long_out, "wt") as fo:
    fo.write("tool\tgroup\tsample\tlevel\tfeature\tregion_type\tcount\tmilli_div\n")
    for path in files:
        with gzip.open(path, "rt") as fh:
            header = next(fh)
            assert header.startswith("tool\tgroup\tsample\tlevel"), \
                f"unexpected header in {path}"
            for line in fh:
                fo.write(line)
                tool, group, sample, level, feature, _reg, cnt, _md = \
                    line.rstrip("\n").split("\t")
                cnt = float(cnt)
                if level == "locus":
                    k = (tool, group, sample)
                    qc[k][0] += cnt
                    qc[k][1] += cnt > 0
                elif cnt >= min_count:
                    sup[(group, feature)].add(tool)

with open(f"{outdir}/aggregate/qc_per_sample.tsv", "w") as fq:
    fq.write("tool\tgroup\tsample\ttotal_counts\tn_positive_loci\tflag\n")
    for (tool, group, sample), (tot, npos) in sorted(qc.items()):
        flag = "ZERO_COUNTS" if tot <= 0 else ("LOW_COUNTS" if tot < 1000 else "OK")
        fq.write(f"{tool}\t{group}\t{sample}\t{tot:.1f}\t{npos}\t{flag}\n")

rmsk_tools = {"telescope", "tetranscripts", "squire"}
with gzip.open(f"{outdir}/aggregate/erv_subfamily_evidence.tsv.gz", "wt") as fo:
    fo.write("group\tsubfamily\tregion_type\tn_tools_support\tervmap_support\n")
    for (group, sf) in sorted(sup):
        tools = sup[(group, sf)]
        reg = "internal" if sf.endswith("-int") else "ltr"
        n = len(tools & rmsk_tools)
        e = int("ervmap" in tools)
        fo.write(f"{group}\t{sf}\t{reg}\t{n}\t{e}\n")

print(f"files={len(files)} group_subfamilies={len(sup)}", file=sys.stderr)
