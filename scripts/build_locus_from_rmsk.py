#!/usr/bin/env python3
"""Build the locus-level ERV annotation directly from a UCSC RepeatMasker table.

One file is the source of every annotation attribute (coordinates, subfamily,
family, class, strand and milliDiv), so nothing is mixed across releases.  One
RepeatMasker fragment becomes one locus: fragments are *not* merged, because
merging rules are hand-made and cannot be reproduced, and fragments can always
be aggregated upward while merged loci cannot be split back.

Usage:
  build_locus_from_rmsk.py RMSK FAMILIES GTF MAP EXON [CHROMS]

  RMSK      UCSC rmsk table (hg38_rmsk.txt, 17 columns, 0-based coords)
  FAMILIES  comma-separated repFamily values to keep ("-" keeps every family)
  GTF/MAP/EXON  output paths, "-" to skip one of them
  CHROMS    regex of accepted chromosomes, default ^chr([0-9]+|X|Y)$

RMSK columns used (0-based): 2 milliDiv, 5 genoName, 6 genoStart, 7 genoEnd,
9 strand, 10 repName, 11 repClass, 12 repFamily.

The four ERV families (ERV1, ERVK, ERVL, ERVL-MaLR) are exactly the subfamily
universe of the ERV-SQL platform's 550-subfamily dictionary, so filtering by
repFamily here reproduces that set without shipping the dictionary itself.

Map columns:
  locus_id, subfamily, family, rep_class, region_type, chr, start, end,
  strand, milli_div
"""
import re
import sys

rmsk_in, families_in, gtf_out, map_out, exon_out = sys.argv[1:6]
chrom_re = re.compile(sys.argv[6] if len(sys.argv) > 6 else r"^chr([0-9]+|X|Y)$")

families = None if families_in == "-" else set(families_in.split(","))

recs = {}
with open(rmsk_in) as fh:
    for line in fh:
        p = line.split()
        if len(p) < 13:
            continue
        chrom, start, end, strand = p[5], int(p[6]), int(p[7]), p[9]
        rep, cls, fam, md = p[10], p[11], p[12], p[2]
        if not chrom_re.match(chrom):
            continue
        if families is not None and fam not in families:
            continue
        # duplicate rows appear in some rmsk dumps; keep the first one
        recs.setdefault((rep, chrom, start, end, strand), (fam, cls, md))

fo = open(gtf_out, "w") if gtf_out != "-" else None
fe = open(exon_out, "w") if exon_out != "-" else None
n_int = 0
with open(map_out, "w") as fm:
    fm.write("locus_id\tsubfamily\tfamily\trep_class\tregion_type\t"
             "chr\tstart\tend\tstrand\tmilli_div\n")
    for (rep, chrom, start, end, strand), (fam, cls, md) in sorted(
        recs.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[0][3])
    ):
        locus = f"{rep}__{chrom}_{start}_{end}_{strand}"
        # RepeatMasker names the internal model of an ERV "<name>-int"; every
        # other ERV fragment is an LTR (or, for ERVL-MaLR, a solo LTR).
        reg = "internal" if rep.endswith("-int") else "ltr"
        n_int += reg == "internal"
        attr = (
            f'gene_id "{locus}"; transcript_id "{locus}"; locus "{locus}"; '
            f'repName "{rep}"; repFamily "{fam}"; repClass "{cls}"; '
            f'family_id "{fam}"; class_id "{cls}"; '
            f'subfamily "{rep}"; region_type "{reg}"; milliDiv "{md}";'
        )
        line = "\t".join([chrom, "rmsk", "gene", str(start + 1), str(end),
                          ".", strand, ".", attr]) + "\n"
        if fo:
            fo.write(line)
        if fe:
            # TEcount wants one exon per locus; same coordinates, type exon
            fe.write(line.replace("\trmsk\tgene\t", "\trmsk\texon\t"))
        fm.write(f"{locus}\t{rep}\t{fam}\t{cls}\t{reg}\t{chrom}\t{start}\t{end}"
                 f"\t{strand}\t{md}\n")

if fo:
    fo.close()
if fe:
    fe.close()

print(f"loci={len(recs)} internal={n_int} ltr={len(recs) - n_int}", file=sys.stderr)
