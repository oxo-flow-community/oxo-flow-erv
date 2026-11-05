#!/usr/bin/env python3
"""Generate the synthetic fixture kit used by test/run.sh and the CI dry-run.

Everything here is fake but shaped like the real inputs, so a fresh clone can
be validated (and, with the conda envs installed, executed) without downloading
hg38 or a RepeatMasker table:

  genome.fa     two 12 kb / 8 kb contigs named chr1/chr2 so the default
                chromosome filter of the workflow applies unchanged
  genes.gtf     three protein-coding genes (featureCounts / TEcount input)
  rmsk.txt      UCSC RepeatMasker table layout, six ERV fragments from the four
                ERV families (ERV1/ERVK/ERVL/ERVL-MaLR) plus three non-ERV
                repeats that MUST be filtered out
  ERVmap.bed    four provirus intervals in the official ERVmap BED layout
  raw/*.fq.gz   paired 100 bp reads simulated from the genome

Deterministic: same seed, byte-identical output.

Usage: python3 test/fixtures/make_fixtures.py
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference")
RAW = os.path.join(HERE, "raw")

# -----------------------------------------------------------------------------
# Genome: random sequence, deterministic.
# -----------------------------------------------------------------------------
CONTIGS = {"chr1": 12000, "chr2": 8000}


def build_genome(rng):
    seqs = {}
    for name, length in CONTIGS.items():
        seqs[name] = "".join(rng.choice("ACGT") for _ in range(length))
    return seqs


# -----------------------------------------------------------------------------
# RepeatMasker-format annotation (0-based start, half-open end, the UCSC layout).
# (chrom, start, end, strand, repName, repFamily, milliDiv)
# -----------------------------------------------------------------------------
REPEATS = [
    # ERV1 provirus on chr1: 5' LTR + internal + 3' LTR
    ("chr1", 1500, 1800, "+", "MER41A", "ERV1", 35),
    ("chr1", 1800, 4200, "+", "MER41-int", "ERV1", 52),
    ("chr1", 4200, 4500, "+", "MER41A", "ERV1", 28),
    # ERVK provirus on chr1, minus strand
    ("chr1", 6000, 6300, "-", "LTR12C", "ERVK", 61),
    ("chr1", 6300, 9200, "-", "HERVK-int", "ERVK", 88),
    ("chr1", 9200, 9500, "-", "LTR12C", "ERVK", 57),
    # ERVL internal model on chr2
    ("chr2", 700, 1000, "+", "LTR57", "ERVL", 74),
    ("chr2", 1000, 2500, "+", "LTR57-int", "ERVL", 96),
    # ERVL-MaLR LTR pair on chr2
    ("chr2", 4000, 4300, "+", "MLT1A0", "ERVL-MaLR", 44),
    ("chr2", 4300, 4600, "-", "MLT1A1", "ERVL-MaLR", 39),
    # non-ERV repeats: must be excluded by the family filter
    ("chr1", 3000, 3300, "-", "AluSx", "Alu", 71),
    ("chr2", 6000, 6500, "+", "L1HS", "L1", 120),
    ("chr2", 7000, 7200, "-", "MER5A", "hAT-Charlie", 65),
]

REP_CLASS = {"ERV1": "LTR", "ERVK": "LTR", "ERVL": "LTR", "ERVL-MaLR": "LTR",
             "Alu": "SINE", "L1": "LINE", "hAT-Charlie": "DNA"}


def write_rmsk():
    with open(os.path.join(REF, "rmsk.txt"), "w") as fh:
        for i, (chrom, start, end, strand, rep, family, md) in enumerate(REPEATS, 1):
            fh.write("\t".join(map(str, [
                0,                        # bin (unused by the builder)
                2000,                     # swScore
                md, 0, 0,                 # milliDiv / milliDel / milliIns
                chrom, start, end, 0,     # genoName / start / end / genoLeft
                strand,
                rep, REP_CLASS[family], family,
                0, end - start, 0,        # repStart / repEnd / repLeft
                i,                        # id
            ])) + "\n")


# -----------------------------------------------------------------------------
# Genes: three genes, exons only (featureCounts -t exon -g gene_id).
# -----------------------------------------------------------------------------
GENES = [
    ("g1", "chr1", 500, 1400, "+", [(500, 700), (900, 1100), (1300, 1400)]),
    ("g2", "chr1", 4600, 5900, "-", [(4600, 4900), (5600, 5900)]),
    ("g3", "chr2", 2700, 3900, "+", [(2700, 3100), (3500, 3900)]),
]


def write_genes():
    with open(os.path.join(REF, "genes.gtf"), "w") as fh:
        for gid, chrom, start, end, strand, exons in GENES:
            tid = f"{gid}_T1"
            attr = f'gene_id "{gid}"; transcript_id "{tid}";'
            fh.write("\t".join([chrom, "synth", "transcript", str(start + 1), str(end),
                                ".", strand, ".", attr]) + "\n")
            for n, (s, e) in enumerate(exons, 1):
                a = attr + f' exon_number "{n}";'
                fh.write("\t".join([chrom, "synth", "exon", str(s + 1), str(e),
                                    ".", strand, ".", a]) + "\n")


# -----------------------------------------------------------------------------
# ERVmap proviruses: official BED layout chr/start/end/name/score/strand.
# -----------------------------------------------------------------------------
ERVPROVIRUS = [
    ("chr1", 1450, 4550, "5803", 500, "+"),
    ("chr1", 5950, 9550, "5804", 500, "-"),
    ("chr2", 650, 2600, "5807", 500, "+"),
    ("chr2", 3950, 4650, "5808", 500, "+"),
]


def write_ervmap():
    with open(os.path.join(REF, "ERVmap.bed"), "w") as fh:
        for row in ERVPROVIRUS:
            fh.write("\t".join(map(str, row)) + "\n")


# -----------------------------------------------------------------------------
# Reads: paired 100 bp from random fragments, both strands, no errors.
# -----------------------------------------------------------------------------
READ_LEN = 100


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_reads(seqs, sample, n_pairs, rng, start_bias=None):
    r1 = os.path.join(RAW, f"{sample}_1.fastq.gz")
    r2 = os.path.join(RAW, f"{sample}_2.fastq.gz")
    with gzip.open(r1, "wt") as f1, gzip.open(r2, "wt") as f2:
        for i in range(n_pairs):
            chrom = start_bias[i % len(start_bias)] if start_bias else rng.choice(list(seqs))
            seq = seqs[chrom]
            insert = rng.randint(2 * READ_LEN, 2 * READ_LEN + 80)
            pos = rng.randint(0, len(seq) - insert - 1)
            frag = seq[pos:pos + insert]
            a = frag[:READ_LEN]
            b = revcomp(frag[-READ_LEN:])
            name = f"{sample}_{i}"
            f1.write(f"@{name}/1\n{a}\n+\n{'I' * READ_LEN}\n")
            f2.write(f"@{name}/2\n{b}\n+\n{'I' * READ_LEN}\n")


def main():
    os.makedirs(REF, exist_ok=True)
    os.makedirs(RAW, exist_ok=True)
    rng = random.Random(20261105)
    seqs = build_genome(rng)

    with open(os.path.join(REF, "genome.fa"), "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    write_genes()
    write_rmsk()
    write_ervmap()

    # sample_02 is sequenced deeper and is deliberately biased towards the ERV
    # loci so that every quantification tool has non-zero counts on the fixture.
    write_reads(seqs, "sample_01", 400, rng)
    write_reads(seqs, "sample_02", 600, rng, start_bias=["chr1", "chr2", "chr1", "chr2", "chr1"])

    print("fixtures written to", HERE)


if __name__ == "__main__":
    main()
