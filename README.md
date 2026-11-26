# oxo-flow-erv — ERV / transposable-element quantification with four tools

[![CI](https://github.com/oxo-flow-community/oxo-flow-erv/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-erv/actions/workflows/ci.yml)

RNA-seq from FASTQ to a comparable ERV expression table. Reads are trimmed with
fastp and aligned with STAR; four published tools then quantify endogenous
retroviruses on **one** RepeatMasker-derived annotation, and their output is
aggregated into a single long table with per-sample QC and per-subfamily
evidence. Gene-level counts come from the same BAMs with `featureCounts`, so
ERV and gene counts always describe the same library.

| Tool | Level | Annotation it reads |
|---|---|---|
| [ERVmap](https://github.com/mtokuyama/ERVmap) | provirus | the official 3,220-provirus BED |
| [Telescope](https://github.com/mlbendall/telescope) | locus + subfamily | ERV fragment GTF |
| [TEtranscripts](https://github.com/mhammell-laboratory/TEtranscripts) | locus + subfamily | ERV fragment exon GTF |
| [SQuIRE](https://github.com/wyang1989/SQuIRE) | locus + subfamily | full RepeatMasker clean_folder |

Ported from [SunPast/ERV-pipeline](https://github.com/SunPast/ERV-pipeline)
(Snakemake, MIT): same tools, same parameters, same four-tool comparison. The
main change is that **ERVmap is computed here** (STAR + bedtools coverage)
instead of imported. All porting differences are listed under
[Deviations](#deviations-from-erv-pipeline).

All stages form one DAG. Run one stage only with
`oxo-flow run main.oxoflow -t ervmap` (rules are namespaced
`alignment::* / annotation::* / ervmap::* / telescope::* / tetranscripts::* /
squire::* / aggregate::*`).

**Cohorts are a directory level, not a naming convention.** Each
`[[sample_groups]]` entry becomes `{group}`, and every per-sample output keeps
the layout `{out_dir}/{tool}/{group}/{sample}/...` (the source called this level
`cohort`). Sample names may repeat across groups; the aggregation tables carry
a `group` column and are computed per group.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow ≥ 0.13.1**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-erv.git
cd oxo-flow-erv
```

### 3. Reference data (GRCh38) and runtime

- `genome_fa` — primary assembly FASTA (+ `.fai`)
- `gene_gtf` — GENCODE GTF (featureCounts + TEcount)
- `rmsk_table` — UCSC `hg38_rmsk.txt` (RepeatMasker table, 17 columns)
- `ervmap_bed` — the ERVmap 3,220-provirus BED **with the `chr` prefix**
  (the official file has none; see [ERVmap](#ervmap-how-its-run-and-pitfalls-handled))
- `star_index_dir` — a STAR 2.7.7a index; leave it unset/unbuilt and the
  workflow builds it from `genome_fa` + `gene_gtf` (use `--genomeSAindexNbases 14`
  for hg38)
- **Reads**: `reads_dir/{sample}_1.fastq.gz` (+ `_2.fastq.gz` when
  `layout = "pe"`), one entry per sample in `[[sample_groups]]`
- **Compute**: up to 10 threads / 32 GB per STAR rule; conda/mamba at runtime
  for the pinned `envs/*.yaml`

A tiny synthetic kit ships in `test/fixtures` (20 kb genome, 3 genes, 10 ERV
fragments, 4 proviruses, 2 samples in 2 groups) so a fresh clone validates and
dry-runs without downloading hg38.

## Usage

```bash
# 1. put reads at <reads_dir>/<sample>_1.fastq.gz (and _2 for paired-end),
#    and list samples in [[sample_groups]] of main.oxoflow
# 2. point the [config] keys at real references
# 3. preview
oxo-flow dry-run main.oxoflow
# 4. run
oxo-flow run main.oxoflow -j 16
# single-end data: layout = "se"
oxo-flow run main.oxoflow -j 16 layout=se
```

Primary outputs:

```text
{out_dir}/bam/{group}/{sample}/{sample}.sorted.bam(.bai)
{out_dir}/gene_counts/{group}/{sample}/{sample}_featureCounts.txt
{out_dir}/ervmap/{group}/{sample}/{sample}ERVresults.txt   raw, 3,220 rows
{out_dir}/ervmap/{group}/{sample}/{sample}.locus.tsv       normalized
{out_dir}/telescope/{group}/{sample}/{sample}.report.tsv
{out_dir}/tetranscripts/{group}/{sample}/{sample}.cntTable
{out_dir}/squire/{group}/{sample}/{sample}_TEcounts.txt
{out_dir}/squire/{group}/{sample}/{sample}_subFcounts.txt
{out_dir}/aggregate/per_sample/{group}/{sample}.long.tsv.gz
{out_dir}/aggregate/erv_expression_long.tsv.gz             tool,group,sample,...
{out_dir}/aggregate/qc_per_sample.tsv                      tool,group,sample,...
{out_dir}/aggregate/erv_subfamily_evidence.tsv.gz          group,subfamily,...
```

## Tools, versions and key parameters

Versions are pinned in `envs/*.yaml` (conda-forge + bioconda).

| Stage | Tool (version) | Key parameters |
|---|---|---|
| trim | fastp 0.23.4 | `--length_required 36 --qualified_quality_phred 20 --unqualified_percent_limit 40 --n_base_limit 5`; PE adds `--detect_adapter_for_pe` |
| align | STAR 2.7.7a | `--outSAMunmapped Within --outFilterMismatchNmax 15 --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical --outSAMstrandField intronMotif --twopassMode Basic --quantMode GeneCounts --outSAMtype BAM Unsorted`; `--limitBAMsortRAM` only when `star_limit_bam_sort_ram > 0` |
| sort/index | samtools 1.15 | `sort -m 2G` then `index` |
| gene counts | subread 2.0.1 (featureCounts) | `-p` (PE only) `-O -s 0 -t exon -g gene_id` (unstranded, to match the TCASIA matrices) |
| ERV definition | python 3.10 | `repFamily ∈ {ERV1, ERVK, ERVL, ERVL-MaLR}` — exactly the source's 550-subfamily dictionary (300 + 42 + 124 + 84); override `erv_families` to change |
| ERVmap | STAR 2.7.7a + bedtools 2.31.1 | see below |
| assign | telescope 1.0.3 (+ samtools 1.18) | `--max_iter 200 --theta_prior 200000 --ncpu 1` on a name-collated BAM |
| TEcount | TEtranscripts 2.2.3 | `--sortByPos --format BAM --mode multi --GTF <gene_gtf> --TE <te_exon_gtf>` |
| SQuIRE | SQuIRE 0.9.9.92 | `Count --read_length 36` on the pipeline BAM over the full-TE clean_folder (Map is not in the DAG) |
| aggregate | python 3.10 | per-sample long tables merged per group; source `cohort` column kept as `group` |

## ERVmap: how it's run and pitfalls handled

The source imported ERVmap results, arguing the 3,220-provirus set *is* the
method. This port computes them with the official wrapper's command chain
(STAR strict-unique align → `bedtools coverage`):

```text
STAR --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.02 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000
samtools index
bedtools coverage -nonamecheck -a ERVmap.bed -b <bam> -counts -sorted
```

Output is the same 3,220-row `<sample>ERVresults.txt`. Each pitfall below
produced a *successful-looking but wrong* result before, so each has a guard:

| Pitfall | Guard |
|---|---|
| Chromosome names: the official BED has no `chr`, the genome index does → `-nonamecheck` writes 3,220 zero rows with exit 0 | take the chr-prefixed BED **and** `check_bed_contigs.py` fails before counting when BED contigs are absent from the FASTA index |
| Contig order: BED is lexicographic, the BAM is not → `-sorted` truncates silently | BED re-sorted with `bedtools sort -faidx`; raw row count must equal the BED row count |
| Reads aligned to the genome, not a 3,220-sequence provirus multi-fasta | documented; `--genomeDir` points at the workflow's genome index |

`--limitBAMsortRAM` is passed only when `ervmap_limit_bam_sort_ram > 0` (the
project scripts pinned a host-specific value; set it to reproduce that exact
invocation).

## Deviations from ERV-pipeline

Intentional porting decisions:

- **ERVmap computed, not imported** (above).
- **Cohorts kept as `{group}`**: `config/samples.tsv` (cohort + sample +
  fq1/fq2/bam) becomes `[[sample_groups]]` + a `reads_dir` layout. One
  capability is not reproduced: the source could point a sample at an existing
  BAM instead of FASTQ; here every sample is aligned from FASTQ (no BAM
  integrity check).
- **ERV definition by family** instead of a shipped 550-subfamily dictionary;
  equivalent locus set on the standard hg38 rmsk (see Tools table). Override
  `erv_families` to change.
- **Provirus→subfamily by max-overlap** (`annotation::ervmap_annotate`) instead
  of the source's private `int_model` catalog. Sufficient for subfamily-level
  evidence; not a substitute for curated provirus models.
- **Single-end support is new** (`layout = "se"`); the official ERVmap wrapper
  is paired-end only.
- **SQuIRE Map is not in the DAG** (would rebuild a STAR index per sample).
- **Not ported**: `merge_gene_counts` / `compute_gene_factors` / external gene
  counts import (gene counts stay per sample), and `export_standard` (the
  unified long table is the deliverable). These are listed in
  `metadata.json` `excluded`.
- **conda instead of Apptainer** for SQuIRE and TEtranscripts.

Two defaults worth knowing:

- `min_count = 2` — a subfamily counts as "supported by a tool" when any sample
  reaches ≥2 counts. This is the source default and this platform's noise
  floor; it affects only `erv_subfamily_evidence.tsv`, never the per-sample
  counts. The upstream analysis compared raw counts with no threshold, so set
  `min_count = 1` to reproduce it literally.
- `squire_read_length = 36` — matches the source default and the `--read_length
  36` hard-coded in the rerun scripts. Set it to your actual read length
  otherwise.

## Verification

`oxo-flow validate` / `lint` / `dry-run` (paired- and single-end) pass on
oxo-flow 0.13.1; `test/run.sh` also executes the annotation, SQuIRE-clean and
ERVmap-annotation scripts against the fixtures.

**SQuIRE was not run live**: `squire=0.9.9.92` could not be installed in the
verification environment (conda conflict) and the container image is unavailable
there. Its command is copied from the source pipeline; `metadata.json` records
`rating: dry-run-verified`.

## Source

Ported from [SunPast/ERV-pipeline](https://github.com/SunPast/ERV-pipeline)
(Snakemake, MIT), which aligned its FASTQ stage with
`OncoHarmony-Network/TCASIA_pipeline` and took its ERV annotation from the UCSC
RepeatMasker table. Attribution in `NOTICE.md`; the upstream license is in
`LICENSE.upstream`. Method citations are in `NOTICE.md`.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community.

## Community

https://oxo-flow-community.github.io/
