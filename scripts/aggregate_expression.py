#!/usr/bin/env python3
"""Put one sample's four tool outputs on one long table.

Usage: aggregate_expression.py OUTDIR GROUP SAMPLE LOCUS_MAP ERVMAP_ANNOTATION

Output: {OUTDIR}/aggregate/per_sample/{GROUP}/{SAMPLE}.long.tsv.gz
  tool  group  sample  level(locus|subfamily)  feature  region_type  count  milli_div

milli_div (RepeatMasker divergence, per mille) is a property of the locus, not
of the tool, so it is filled at locus level only.  It is left empty at
subfamily level: a subfamily spans a wide divergence range and a single summary
number would be misleading.

Parsing notes (column positions verified against real tool output):
  * Telescope report : col0 = locus id, col2 = final_count; the first line is
                       "## RunInfo", other lines starting with "#" are comments.
  * TEcount cntTable : "<locus>:<Family>:<Class>\\t<count>"; gene rows are
                       quoted ids and are dropped here (they are not ERV loci).
  * SQuIRE TEcounts  : col3 = TE_ID (chr|start|end|Subfamily:Family:Class|
                       milliDiv|strand), col15 = tot_counts.
  * SQuIRE subFcounts: col2 = Subfamily:Family:Class, col6 = tot_counts.
  * ERVmap           : already normalized to locus_id / count by
                       normalize_ervmap.py (subfamily joined from the
                       provirus -> RepeatMasker assignment).
"""
import collections
import csv
import gzip
import os
import sys

outdir, group, sample, locus_map, ervmap_annot = sys.argv[1:6]

locus2sf = {}
locus_md = {}
erv_subfamilies = set()
for r in csv.DictReader(open(locus_map), delimiter="\t"):
    locus2sf[r["locus_id"]] = (r["subfamily"], r["region_type"])
    locus_md[r["locus_id"]] = r.get("milli_div", "")
    erv_subfamilies.add(r["subfamily"])

# provirus -> subfamily, for the subfamily roll-up of ERVmap
erv2sf = {}
for r in csv.DictReader(open(ervmap_annot), delimiter="\t"):
    erv2sf[r["ervmap_id"]] = r["subfamily"]

rows = []
subsum = collections.defaultdict(float)


def add(tool, level, feature, reg, count, md=""):
    rows.append((tool, group, sample, level, feature, reg, count, md))


def squire_locus_id(te_id):
    """chr|start|end|Subfamily:Family:Class|milliDiv|strand -> locus_id"""
    p = te_id.split("|")
    if len(p) < 6:
        return None
    return f"{p[3].split(':')[0]}__{p[0]}_{p[1]}_{p[2]}_{p[5]}"


# ---- Telescope (locus) ----
rep = f"{outdir}/telescope/{group}/{sample}/{sample}.report.tsv"
if os.path.exists(rep):
    for line in open(rep):
        if line.startswith("#") or line.startswith("transcript\t"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        gid, cnt = f[0], float(f[2])
        if gid not in locus2sf:
            continue
        sf, reg = locus2sf[gid]
        add("telescope", "locus", gid, reg, cnt, locus_md.get(gid, ""))
        subsum[("telescope", sf, reg)] += cnt

# ---- TEtranscripts / TEcount (locus) ----
cnt_tbl = f"{outdir}/tetranscripts/{group}/{sample}/{sample}.cntTable"
if os.path.exists(cnt_tbl):
    with open(cnt_tbl) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            name = f[0].strip('"').split(":")[0]
            if name not in locus2sf:
                continue
            cnt = float(f[1])
            if cnt <= 0:
                continue
            sf, reg = locus2sf[name]
            add("tetranscripts", "locus", name, reg, cnt, locus_md.get(name, ""))
            subsum[("tetranscripts", sf, reg)] += cnt

# ---- SQuIRE (locus) ----
sq_te = f"{outdir}/squire/{group}/{sample}/{sample}_TEcounts.txt"
if os.path.exists(sq_te):
    with open(sq_te) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 16:
                continue
            lid = squire_locus_id(f[3])
            if lid not in locus2sf:
                continue
            sf, reg = locus2sf[lid]
            add("squire", "locus", lid, reg, float(f[15]), locus_md.get(lid, ""))

# ---- SQuIRE (subfamily, from its own subFcounts file) ----
sq_sf = f"{outdir}/squire/{group}/{sample}/{sample}_subFcounts.txt"
if os.path.exists(sq_sf):
    with open(sq_sf) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            sf = f[2].split(":")[0]
            if sf not in erv_subfamilies:
                continue          # SQuIRE counts the full repeat set; ERVs only
            reg = "internal" if sf.endswith("-int") else "ltr"
            add("squire", "subfamily", sf, reg, float(f[6]))

# ---- ERVmap (locus; subfamily from the provirus annotation) ----
ev = f"{outdir}/ervmap/{group}/{sample}/{sample}.locus.tsv"
if os.path.exists(ev):
    for r in csv.DictReader(open(ev), delimiter="\t"):
        cnt = float(r["count"])
        add("ervmap", "locus", r["locus_id"], "", cnt)
        sf = erv2sf.get(r["locus_id"].replace("ERVmap_", ""), r["subfamily"])
        if sf and sf != ".":
            reg = "internal" if sf.endswith("-int") else "ltr"
            subsum[("ervmap", sf, reg)] += cnt

for (tool, sf, reg), v in subsum.items():
    add(tool, "subfamily", sf, reg, v)

per_dir = f"{outdir}/aggregate/per_sample/{group}"
os.makedirs(per_dir, exist_ok=True)
with gzip.open(f"{per_dir}/{sample}.long.tsv.gz", "wt") as fo:
    fo.write("tool\tgroup\tsample\tlevel\tfeature\tregion_type\tcount\tmilli_div\n")
    for row in rows:
        fo.write("\t".join(map(str, row)) + "\n")

print(f"{group}/{sample}: rows={len(rows)} tools={len({r[0] for r in rows})}",
      file=sys.stderr)
