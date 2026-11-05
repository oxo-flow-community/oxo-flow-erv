# oxo-flow-erv — ERV / transposable-element quantification with four tools

[![CI](https://github.com/oxo-flow-community/oxo-flow-erv/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-erv/actions/workflows/ci.yml)

> ★ Dry-run verified · ⇄ Port of [`SunPast/ERV-pipeline`](https://github.com/SunPast/ERV-pipeline) @ `02ae748`
> — same tools, same parameters, same four-tool comparison; ERVmap is now
> computed instead of imported.

RNA-seq from FASTQ to a comparable ERV expression table. Reads are trimmed with
fastp and aligned with STAR; four published tools then quantify endogenous
retroviruses on **one** annotation, and their output is put on a single long
table with per-sample QC and per-subfamily evidence:

| Tool | Level | Annotation it reads |
|---|---|---|
| [ERVmap](https://github.com/mtokuyama/ERVmap) | provirus | the official 3,220-provirus BED |
| [Telescope](https://github.com/mlbendall/telescope) | locus + subfamily | ERV fragment GTF |
| [TEtranscripts](https://github.com/mhammell-laboratory/TEtranscripts) | locus + subfamily | ERV fragment exon GTF |
| [SQuIRE](https://github.com/wyang1989/SQuIRE) | locus + subfamily | full RepeatMasker clean_folder |

Gene-level expression comes from the same BAMs with `featureCounts`, so ERV and
gene counts always describe the same library.

All stages form one DAG. Run one stage only with
`oxo-flow run main.oxoflow -t ervmap` (rules are namespaced
`alignment::* / annotation::* / ervmap::* / telescope::* / tetranscripts::* /
squire::* / aggregate::*`).

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

### 3. Requirements for a real run

- **Reference data** (GRCh38):
  - `genome_fa` — primary assembly FASTA
  - `gene_gtf` — GENCODE GTF (featureCounts + TEcount)
  - `rmsk_table` — UCSC `hg38_rmsk.txt` (RepeatMasker table, 17 columns)
  - `ervmap_bed` — the official ERVmap BED with a `chr` prefix (3,220 rows)
  - `star_index_dir` — STAR 2.7.7a index; leave it unset/unbuilt and the
    workflow builds it from `genome_fa` + `gene_gtf` (see the `[[references]]`
    block; switch `--genomeSAindexNbases` to 14 for hg38)
- **Reads**: `reads_dir/{sample}_1.fastq.gz` (+ `_2.fastq.gz` when
  `layout = "pe"`), one entry per sample in `[[sample_groups]]`
- **Compute**: up to 10 threads per rule (STAR), 32 GB for STAR, conda/mamba at
  runtime for the pinned `envs/*.yaml`

A tiny synthetic kit ships in `test/fixtures` (20 kb genome, 3 genes, 10 ERV
fragments, 4 proviruses, 2 samples) so a fresh clone validates and dry-runs
without downloading anything.

## Usage

```bash
# 1. put reads at <reads_dir>/<sample>_1.fastq.gz (and _2 for paired-end)
#    and list samples in [[sample_groups]] of main.oxoflow
# 2. point the [config] keys at real references
# 3. preview
oxo-flow dry-run main.oxoflow
# 4. run
oxo-flow run main.oxoflow -j 16
# single-end data: layout = "se" (selects fastp_se / star_align_se / ...)
oxo-flow run main.oxoflow -j 16 layout=se
```

Primary outputs:

```text
{out_dir}/bam/{sample}/{sample}.sorted.bam(.bai)
{out_dir}/gene_counts/{sample}/{sample}_featureCounts.txt
{out_dir}/ervmap/{sample}/{sample}ERVresults.txt          raw, 3,220 rows
{out_dir}/ervmap/{sample}/{sample}.locus.tsv              normalized
{out_dir}/telescope/{sample}/{sample}.report.tsv
{out_dir}/tetranscripts/{sample}/{sample}.cntTable
{out_dir}/squire/{sample}/{sample}_TEcounts.txt
{out_dir}/squire/{sample}/{sample}_subFcounts.txt
{out_dir}/aggregate/erv_expression_long.tsv.gz
{out_dir}/aggregate/qc_per_sample.tsv
{out_dir}/aggregate/erv_subfamily_evidence.tsv.gz
```

## Fidelity

| Source rule (Snakemake) | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| `fastp` | `alignment::fastp_pe` / `fastp_se` | fastp 0.23.4 | same flags; the single-end variant is new |
| `star_align` | `alignment::star_align_pe` / `star_align_se` | STAR 2.7.7a | same flags; `--limitBAMsortRAM` only when `star_limit_bam_sort_ram > 0` (source pinned ~36 GB) |
| `samtools_sort` | `alignment::samtools_sort` | samtools 1.15 | sort + index in one rule, `-m 512M` per thread |
| `gene_counts` | `alignment::featurecounts_pe` / `featurecounts_se` | subread 2.0.1 | `-p` only for paired-end; `-s 0` explicit (unstranded by design) |
| `build_locus_gtf` | `annotation::build_locus_gtf` | python 3.10 | family filter instead of the 550-name dictionary (see Deviations) |
| `build_full_locus_map` | `annotation::build_full_locus_map` | python 3.10 | identical |
| `squire_clean_folder` | `annotation::squire_clean_folder` | python 3.10 | identical |
| — (new) | `annotation::ervmap_bed_sorted` | bedtools 2.31.1 | official BED, `chr` prefix + contig-order sort for `-sorted` |
| `ervmap_provirus_map` | `annotation::ervmap_annotate` | python 3.10 | max-overlap assignment instead of the int_model catalog |
| — **(was import)** | `ervmap::ervmap_star_pe` / `ervmap_star_se` | STAR 2.7.7a | official ERVmap alignment parameters |
| — **(was import)** | `ervmap::ervmap_count` | bedtools 2.31.1 | `coverageBed -counts -sorted` |
| `ervmap_import` | `ervmap::ervmap_normalize` | python 3.10 | provirus-level table, subfamily attached |
| `telescope_assign` | `telescope::telescope_assign` | telescope 1.0.3 + samtools 1.18 | identical |
| `tetranscripts_locus` | `tetranscripts::tetranscripts_count` | TEtranscripts 2.2.3 | identical |
| `squire_count` | `squire::squire_count` | SQuIRE 0.9.9.92 | identical; conda env instead of the container |
| `aggregate_expression` | `aggregate::aggregate_sample` + `merge_aggregate` | python 3.10 | split per sample + cohort merge; `cohort` column dropped |

### ERVmap is quantified here, not imported

The Snakemake source deliberately imported ERVmap results, arguing that the
3,220-provirus set *is* the method. That reasoning does not hold once the
reference is reproducible: the official wrapper is a two-step program whose
parameters are the actual method. This port runs those steps —
STAR with `--outFilterMultimapNmax 1 --outFilterMismatchNmax 999
--outFilterMismatchNoverLmax 0.02 --alignIntronMin 20 --alignIntronMax 1000000
--alignMatesGapMax 1000000`, then
`bedtools coverage -nonamecheck -a ERVmap.bed -b <bam> -counts -sorted` —
against the same 3,220-provirus BED, so the output is the same 3,220-row
`<sample>ERVresults.txt` the old imports carried.

Two failure modes are handled explicitly because they are silent:

- **Contig order.** The official BED is in lexicographic chromosome order
  (`chr1, chr10, …, chr2`). `-sorted` compares that with the BAM header order
  and aborts mid-file, truncating the result; the BED is therefore re-sorted by
  the FASTA index before counting.
- **Chromosome names.** `-nonamecheck` turns a name mismatch into 3,220 zero
  rows and an exit code of 0. `qc_per_sample.tsv` flags that as `ZERO_COUNTS`.

### Deviations from the source pipeline

- **Cohorts → samples.** `config/samples.tsv` (cohort + sample + fq1/fq2/bam)
  becomes `[[sample_groups]]` plus the `reads_dir` layout, so every output is
  per sample and the `cohort` column is gone. Run one cohort per `out_dir` (or
  per config profile) when cohorts must stay separate.
- **ERV definition by family.** The source filtered RepeatMasker rows through a
  550-subfamily dictionary. Here the filter is `repFamily ∈ {ERV1, ERVK, ERVL,
  ERVL-MaLR}`; those four families are exactly the dictionary's contents
  (300 + 42 + 124 + 84 = 550), so the locus set is the same without shipping
  the dictionary. Override `erv_families` to change it.
- **No merged loci, exactly as before.** One RepeatMasker fragment = one locus;
  the annotation never merges fragments.
- **Provirus assignment without the private catalog.** Instead of the platform's
  `int_model` / LTR decomposition, each provirus is assigned the ERV fragment
  with the largest overlap (`annotation::ervmap_annotate`). Enough for
  subfamily-level evidence; not a substitute for the curated provirus models.
- **Single-end support is new.** The official ERVmap wrapper requires
  `--read1`/`--read2`. Under `layout = "se"` the same STAR parameters run with a
  single input file.
- **SQuIRE Map is not in the DAG.** `Count` runs on the pipeline BAM (as in the
  source); `Map` would rebuild a STAR index per sample (~32 GB) and is not
  needed. `--fetch_folder` is kept from the source call and points at an
  existing directory, since `--clean_folder` is supplied.
- **Gene-expression merge not ported.** `merge_gene_counts.py`,
  `compute_gene_factors.py` and the external-counts import belonged to the
  private warehouse; here gene counts stay per sample.
- **Database export not ported.** `export_standard.py` wrote CSVs for a private
  SQL loader; the unified long table is the deliverable here.
- **conda instead of containers.** SQuIRE and TEtranscripts ran from Apptainer
  images in the source; the port pins the same versions as conda envs.
- **Engine limits.** `min_version = 0.13.1`; the build strings in `[[references]]`
  are literal because this engine does not expand placeholders inside `build`
  (`{effective_threads}` / `{effective_memory_mb}` are not available either).

## Verification status

`oxo-flow validate` / `lint` / `dry-run` (paired- and single-end) pass on
oxo-flow 0.13.1, and `test/run.sh` also executes the annotation, SQuIRE-clean
and ERVmap-annotation scripts against the fixtures.

A live run on the synthetic fixtures was performed with STAR 2.7.7a,
samtools 1.15, fastp 0.23.4, subread 2.0.1, bedtools 2.31.1, telescope 1.0.3
and TEtranscripts 2.2.3. These stages produced non-zero, biologically
consistent output:

| Stage | Result on the fixtures |
|---|---|
| annotation | 10 ERV loci (6 internal / 4 LTR), 4 proviruses annotated |
| alignment + featureCounts | 2 samples, genes g1/g2/g3 counted |
| ERVmap | 4/4 proviruses with non-zero counts, subfamily assigned (MER41-int/ERV1, HERVK-int/ERVK) |
| Telescope | 10 features, 400 fragments, counts 59/51/44/14/… |
| TEcount | 10 ERV loci + 3 genes, counts matching the GTF locus ids |
| aggregation | long table (3 tools), QC table and subfamily evidence built on the real outputs |

The live run found four defects that dry-run could not: a STAR output name that
did not match the declared output, `--theta_prior` needing an integer, TEcount's
relative `--outdir` being deleted with the scratch directory, and tool failures
masked by a trailing `rm -rf`. All four are fixed in the rules above.

**SQuIRE was not run live**: `squire=0.9.9.92` could not be installed in the
verification environment (conda dependency conflict) and the container image is
not available there. Its command is copied from the source pipeline, where
`squire Count` was verified on real data, but this repository has not executed
it. `metadata.json` therefore records `rating: dry-run-verified`.

## Source

Ported from **[SunPast/ERV-pipeline](https://github.com/SunPast/ERV-pipeline)**
(Snakemake, commit `02ae748`, MIT). That project in turn aligned its FASTQ stage
with `OncoHarmony-Network/TCASIA_pipeline` (STAR/fastp/samtools parameters) and
took its ERV annotation from the UCSC RepeatMasker table. Attribution in
`NOTICE.md`; the upstream license is reproduced in `LICENSE.upstream`.

## Test

```bash
bash test/run.sh   # validate + lint + dry-run (pe/se) + fixture scripts
```

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community.

## Community

https://oxo-flow-community.github.io/
