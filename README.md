# ChromaFlow

**English** · [中文文档 / Chinese README](README.zh-CN.md)

A downstream analysis framework for ChIP-seq / DAP-seq / ATAC-seq built on top of nf-core pipelines. It completes the evidence chain from raw peaks to **candidate target genes**.

Organised as software: install the code once, and let each project's own `chromaflow.config` decide where data and results live. The code directory stays clean.

---

## Contents

- [Design: code / config / runtime](#design-code--config--runtime)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Help system](#help-system)
- [Core concept: peakset](#core-concept-peakset)
- [Command overview](#command-overview)
- [Pipeline in detail](#pipeline-in-detail)
- [Optional steps and how to skip them](#optional-steps-and-how-to-skip-them)
- [Sample record file and read pairing](#sample-record-file-and-read-pairing)
- [Configuration and automatic file discovery](#configuration-and-automatic-file-discovery)
- [Recipes](#recipes)
- [Output guide](#output-guide)
- [Environments](#environments)
- [Quality assessment and interpretation](#quality-assessment-and-interpretation)
- [Adapting to other species](#adapting-to-other-species)

---

## Design: code / config / runtime

| Layer | Contents | Location |
|---|---|---|
| ① Software | CLI, lib, steps, templates | **this repository** (clone once, add to PATH) |
| ② Project config | species parameters, I/O paths | one `chromaflow.config` per dataset |
| ③ Runtime I/O | fastq, references, results, logs | determined by `OUTPUT_DIR` in the config |

The repository contains **code only** — no `results/`, `reference/` or `logs/`. Switching datasets means copying a config and editing paths; the code directory never changes.

```
ChromaFlow/
├── bin/chromaflow              # unified CLI entry point
├── lib/
│   ├── common.sh               # config parsing, path resolution, peakset, logging
│   ├── help.sh                 # `chromaflow help [command]`
│   └── steps/                  # analysis steps
│       ├── 00_prepare_reference.sh   00b_normalize_gff.sh
│       ├── normalize_gff.py          make_samplesheet.py
│       ├── 01_make_samplesheet.sh    02_run_nfcore.sh
│       ├── 03_idr_analysis.sh        04_meme_motif.sh
│       ├── 05_peak_annotation.R/.sh  06_integrate.sh
│       ├── 07_gene_tracks.sh         08_go_analysis.sh
│       ├── 09_lookup.sh              plot_tracks.sh  ← plotting core
│       └── support/
├── templates/                  # project templates (config + samples_*)
└── setup_environments.sh       # conda/mamba environment installer
```

---

## Installation

```bash
git clone <repo> ~/software/ChromaFlow
export CHROMAFLOW_HOME=~/software/ChromaFlow
export PATH="$CHROMAFLOW_HOME/bin:$PATH"      # add to ~/.bashrc

bash $CHROMAFLOW_HOME/setup_environments.sh          # core conda environments
WITH_AGAT=true bash $CHROMAFLOW_HOME/setup_environments.sh   # optional AGAT env
SKIP_R=true    bash $CHROMAFLOW_HOME/setup_environments.sh   # skip R (use annotate-py)
```

---

## Quick start

```bash
# 1. Create a project
chromaflow init /data/projects/MyChIP --type chipseq
cd /data/projects/MyChIP
vim chromaflow.config        # species + paths to the files you already have
vim samples_chipseq.txt      # samples and fastq locations

# 2. Reference (optional — skip if you already have genome/GTF/chrom.sizes)
mamba activate ChromaFlow_ref
chromaflow reference          # 00  genome / GTF / chrom.sizes
chromaflow normalize-gff      # 00b optional, improves annotation accuracy
chromaflow doctor             # self-check

# 3. Main pipeline
chromaflow samplesheet        # 01  build the nf-core samplesheet
mamba activate nfcore_env
chromaflow run                # 02  run nf-core/chipseq

# 4. Downstream (the core value)
mamba activate ChromaFlow_idr        && chromaflow idr        # 03 reproducible peaks
mamba activate ChromaFlow_motif      && chromaflow motif      # 04 motif discovery
mamba activate ChromaFlow_R          && chromaflow annotate   # 05 peak annotation
mamba activate ChromaFlow_downstream && chromaflow candidates # 06 candidate targets
```

Global options (available on every command):

| Option | Effect |
|---|---|
| `-c, --config <file>` | Project config (default `$CHROMAFLOW_CONFIG` > `./chromaflow.config`) |
| `-t, --type <type>` | `chipseq` (default) / `atacseq` / `dapseq` |
| `-o, --outdir <dir>` | Redirect **this step's** output to a directory |
| `--peakset <spec>` | Which peak set downstream uses (see below) |
| `--profile <name>` | nf-core execution mode for `run`: `singularity` (default) / `apptainer` / `docker` / `podman` / `conda` / `mamba` |
| `-h, --help` | Help for the current command |

---

## Help system

```bash
chromaflow help                 # command overview
chromaflow help <command>       # details for one command
chromaflow <command> --help     # same thing
chromaflow help config          # what belongs in chromaflow.config
chromaflow help samples         # sample record format + read pairing rules
chromaflow help peakset         # peakset semantics
chromaflow help envs            # which conda env each step needs
chromaflow version
```

---

## Core concept: peakset

Every downstream step (04/05/06/07/09) answers the same question first: **which peaks?** That is controlled by `--peakset`, which is given **on the command line and never stored in the config**, so alternatives are easy to compare.

| spec | Meaning | Output dir suffix |
|---|---|---|
| `idr` (default) | reproducible peaks at IDR < 0.01 | none |
| `idr:0.05` / `idr:0.1` | looser IDR thresholds | `_idr0.05` |
| `consensus` | skip IDR, use the nf-core consensus set (signalValue back-filled) | `_consensus` |
| `<file path>` | any custom peak file | `_custom` |

Different peaksets write to **separate suffixed directories** (`annotation_consensus/`, `candidates_consensus/` …) so results coexist without overwriting.

> **Important:** to analyse one peak set, run `motif → annotate → candidates → lookup` with the **same `--peakset`**. Motifs are derived from a specific set of peak sequences; mixing peaksets between step 04 and 06 means FIMO hits will not line up with the annotated peaks.

```bash
# strict route (reproducible peaks)
chromaflow motif && chromaflow annotate && chromaflow candidates

# permissive route (when replicate agreement is poor)
chromaflow motif      --peakset consensus
chromaflow annotate   --peakset consensus
chromaflow candidates --peakset consensus
```

---

## Command overview

| Command | Step | Environment |
|---|---|---|
| `chromaflow init` | project skeleton + config template | — |
| `chromaflow reference` | 00 reference preparation *(optional)* | `ChromaFlow_ref` |
| `chromaflow normalize-gff` | 00b GFF3 normalization *(optional)* | `ChromaFlow_ref` / `ChromaFlow_agat` |
| `chromaflow samplesheet` | 01 build nf-core samplesheet | — |
| `chromaflow run` | 02 run the nf-core pipeline | `nfcore_env` |
| `chromaflow idr` | 03 IDR reproducibility | `ChromaFlow_idr` |
| `chromaflow motif` | 04 MEME-ChIP motif discovery | `ChromaFlow_motif` |
| `chromaflow annotate` | 05 peak annotation (ChIPseeker) | `ChromaFlow_R` |
| `chromaflow annotate-py` | 05 alternative: pure Python | `ChromaFlow_ref` |
| `chromaflow candidates` | 06 candidate target integration | `ChromaFlow_downstream` |
| `chromaflow lookup` | 09 prior-gene binding lookup | `ChromaFlow_downstream` |
| `chromaflow tracks` | 07 single-gene browser tracks | `ChromaFlow_downstream` |
| `chromaflow signal` | *(optional)* deepTools aggregate profiles | `ChromaFlow_downstream` |
| `chromaflow go-prep` / `go-plot` | 08 GO export / visualization | `ChromaFlow_ref` |
| `chromaflow spacing` | *(support)* motif spacing analysis | `ChromaFlow_ref` |
| `chromaflow doctor` | environment and file self-check | any |

---

## Pipeline in detail

### 00 Reference preparation *(optional)*
Assembles genome FASTA, GTF and chrom.sizes, and can merge organelle genomes (chloroplast / mitochondrion) into the nuclear genome so organelle-derived false peaks are identifiable.

Source files may be compressed or not — **compression is detected from the file's magic bytes, not the extension**. Plain `.fa` / `.gff3` files are used in place with no copy; gzip/bgzip, bzip2, xz, zstd and zip are decompressed once into the reference directory and reused afterwards; a file that is gzipped despite being named `.gff3` (a common download artefact) is still handled correctly, and directory auto-detection never picks a compressed file for a step that needs plain text.

Everything this step produces is **registered in a manifest** (`<OUTPUT_DIR>/.chromaflow/reference.manifest`) and picked up automatically later — you never copy generated paths back into the config. If you already have a genome, GTF and chrom.sizes, just place them in `REFERENCE_DIR` and skip the step entirely.

### 00b GFF normalization *(optional)*
GFF3 files from Phytozome and similar sources may **lack `exon` rows** or use inconsistent UTR naming, which makes the gene model incomplete and breaks ChIPseeker's UTR/exon/intron classification (typical symptom: 0% UTR in the annotation pie chart).

```bash
chromaflow normalize-gff                    # auto: AGAT if present, else built-in
chromaflow normalize-gff --engine builtin   # force the dependency-free engine
chromaflow normalize-gff --engine agat      # force AGAT
```

AGAT is powerful but its Perl dependency tree frequently conflicts during environment solving, so it is **installed separately and optionally**. Without it, the built-in Python normalizer adds exons from CDS, infers 5'/3' UTRs, synthesizes missing gene features and sorts the output. Skipping this step entirely is also fine — annotation falls back to the plain GFF3 and only warns.

### Measured values: genome size, read length, GC

These three used to be hand-entered constants in the config. They are now measured, because a hand-entered value is a value nobody re-checks when the genome or the sequencing run changes.

| Value | Measured by | How | Actually used for |
|---|---|---|---|
| `MACS_GSIZE` | `reference` (00) | khmer `unique-kmers.py -k <read length>` — the command the old logs asked you to run manually | MACS's Poisson background, λ = N·d/gsize |
| `READ_LENGTH` | `samplesheet` (01) | samples 2000 reads per fastq, takes the maximum | only used by nf-core to *derive* the genome size when `--macs_gsize` is absent |
| `GENOME_GC` | `reference` (00) | one pass over the FASTA | informational only — see the note below |

Results are written to `<OUTPUT_DIR>/.chromaflow/derived.env`, **not back into your config**:

```
config      = what you declare   (always wins; setting a value disables its measurement)
derived.env = what was measured  (regenerated when the input changes)
```

This is the same separation as the reference manifest. The config stays a hand-written declaration of intent that can live in git without a script rewriting it under you, and the measured values are scoped to the outputs they describe, so deleting `outputs/` invalidates them correctly. If you nevertheless want them archived in the config, `chromaflow reference --write-config` appends them as a delimited, idempotent block (with a backup).

Measurement is cached against a signature of the genome file and `k`, so it runs once; `--force` re-measures and `--no-measure` skips it.

**If khmer is unavailable**, the fallback is a non-N base count (an upper bound, ~5–10% high), and step 02 can instead omit `--macs_gsize` and pass `--read_length` so that nf-core computes the size internally — <cite index="5-1">nf-core/chipseq runs the same khmer tool in its own `KHMER_UNIQUEKMERS` process when the size is not supplied</cite>. Either way nothing has to be entered by hand.

**Do these numbers actually matter?**

- `MACS_GSIZE` is real but forgiving. It scales the Poisson background, so an error changes p-values slightly. A few percent — for example a placeholder 100,000,000 against a measured 103,854,640 for *Chlamydomonas*, a 3.7% difference — will not visibly move the peak list. Being wrong by a factor of two would.
- `READ_LENGTH` is inert whenever `MACS_GSIZE` is supplied: <cite index="7-1">nf-core only uses it to calculate the MACS genome size when `--macs_gsize` isn't provided</cite>. It is measured because it is nearly free and it keeps the k-mer size honest.
- `GC_CONTENT` has been **removed**. It was never used for anything: step 04 built its background model with `fasta-get-markov`, which estimates base and dinucleotide frequencies from the actual sequence, and the config value only appeared in a log line — which meant the log confidently printed whatever number you had typed, regardless of the genome. Step 04 now reports the GC read back from the model it actually used.

### 01–02 Samplesheet and nf-core
`samplesheet` builds the nf-core input table from your sample record file (see [below](#sample-record-file-and-read-pairing)); `run` invokes nf-core/chipseq (or atacseq), producing BAM, bigWig, MACS peaks and consensus peaks.

**Execution mode is controlled by `NFCORE_PROFILE`** (command line `--profile` wins over the config). It sets nextflow's `-profile` *and* selects which execution path `chromaflow run` takes:

| Value | Mode | Notes |
|---|---|---|
| `singularity` (default) | container | images pulled automatically; best reproducibility |
| `apptainer` | container | for Apptainer 1.x installations |
| `docker` / `podman` | container | those engines instead |
| `conda` / `mamba` | conda | one env per process, ~30–60 min on first run; for networks where image pulls fail |

```bash
chromaflow run --profile conda      # one-off
NFCORE_PROFILE=conda                # or make it the project default in the config
```

Combinations such as `NFCORE_PROFILE="singularity,test"` are passed through as-is. Image cache: `$OUTPUT_DIR/.singularity_cache` (exported as both `NXF_SINGULARITY_CACHEDIR` and `NXF_APPTAINER_CACHEDIR`). Conda env cache: `$OUTPUT_DIR/.conda_envs`. Conda mode rewrites channels to a Tsinghua mirror by default — set `NFCORE_CONDA_MIRROR=none` to keep your own channel configuration.

### 03 IDR reproducibility
Runs IDR on every replicate pair and keeps peaks that are **statistically reproducible**. Output has 8 columns: `chr start end name signalValue strand globalIDR n_pairs`.

```bash
chromaflow idr                          # default: union — any passing pair keeps the peak
chromaflow idr --min-pairs 2            # keep only peaks supported by ≥2 pairs
chromaflow idr --exclude-rep REP1       # rerun without a replicate
chromaflow idr --exclude-rep REP1,REP3
```

- **`--min-pairs N`**: `n_pairs` records how many pairs support each peak. With ≥3 replicates, `--min-pairs 2` is a safer way to tighten than dropping a replicate — no data is discarded, the evidence bar is simply raised.
- **`--exclude-rep`**: use it to test the impact of a suspicious replicate empirically (**read the pairwise table in the log first**).

### 04 MEME-ChIP motif discovery
De novo motif discovery on summit ±100 bp, Tomtom comparison against a known motif database, CentriMo central enrichment and a genome-wide FIMO scan.

- Uses **all peak sequences** (no top-N truncation). `MEME_SEARCHSIZE=0` means unlimited; runtime grows roughly quadratically, so set e.g. `MEME_SEARCHSIZE=200000` when sequences are numerous.
- Output directories carry the peakset suffix: `motif/`, `motif_consensus/` …
- Gaps in MEME numbering (MEME-1, 3, 4) are normal: MEME-ChIP clusters and filters redundant motifs, and filtered ones produce no FIMO output.

### 05 Peak annotation
`annotatePeak` assigns each peak to a gene and classifies its genomic feature.

- `PROMOTER_UP` / `PROMOTER_DOWN` map directly to ChIPseeker's `tssRegion=c(-UP,+DOWN)` and thus **define the "Promoter" slice of the pie chart**.
- ChIPseeker classifies by priority `Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic`. If the promoter window extends downstream of the TSS, the front of the 5'UTR is absorbed by "Promoter" — that is by design.
- Each peak is assigned to **one** gene (nearest TSS); downstream steps always use the `geneId` column.

Produces `peak_annotation_full.tsv`, `genes.bed`, `background_universe.txt`, the annotation pie chart, TSS distance plots and more. `annotate-py` is an R-free alternative with the same output contract.

### 06 Candidate target integration (04 motifs × 05 annotation)
Maps FIMO motif hits (genomic coordinates) back onto peaks and intersects them with the annotated positional features:

> **candidate = acceptable position (`TARGET_FEATURES`) + FIMO motif inside the peak**

```bash
chromaflow candidates
chromaflow candidates --features "Promoter,1st Intron,Distal Intergenic" --pval 1e-4
```

Produces `candidate_targets.tsv` (per peak: coordinates, per-motif counts, verdict), per-motif `motif_<name>_hits.tsv` / `_genes.txt`, `motif_summary.tsv` (with CentriMo), `high_confidence_genes.txt` and `go/`.

### 07 Single-gene browser tracks
```bash
chromaflow tracks --genes GENE1,GENE2
chromaflow tracks --gene-list mylist.txt [--outdir DIR]
```

- **Gene column auto-detection**: a header containing `gene_id` / `Gene` / `Genes` / `gene_name` / `id` is used and skipped, so `query_binding.tsv`, `candidate_targets.tsv` and `candidate_genes.txt` can be passed directly; without a header the first column is used.
- **Peak labels carry signalValue** (e.g. `IDR_peak_1 (sig=6299)`) for quantitative comparison.
- **A log2(ChIP/Input) difference track is always overlaid** (deepTools: replicates merged with `bigwigAverage`, then `bigwigCompare`). It has its own y-axis; positive = enriched (red), negative = below background (blue), flat = indistinguishable from background.

> Coverage track y-values (e.g. 0–6) are normalized sequencing depth, whereas signalValue is MACS fold enrichment (can reach thousands). They are **different scales** — do not judge enrichment strength by peak height.

Plotting lives in `lib/steps/plot_tracks.sh` (the STYLE block at the top controls colours, heights, fonts, DPI). Both `tracks` and `lookup --plot` call it, so **visual tweaks happen in one place**.

### 08 GO enrichment
`go-prep` exports foreground and background gene lists for external tools (agriGO, g:Profiler …); `go-plot` turns the returned results into a bubble chart.

### 09 Prior-gene lookup — independent evidence breakdown
Takes a list of genes you care about (for example phosphate-starvation genes) and evaluates them **independently** from the 05 annotation and 04 FIMO output — it does **not** depend on the candidate table, so it still works when `candidates` did not surface what you expected.

```bash
chromaflow lookup --gene-list phosphate_genes.txt
chromaflow lookup --gene-list phosphate_genes.txt --peakset consensus --plot
chromaflow lookup --gene-list phosphate_genes.txt --min-signal 50 --level candidate
```

`query_binding.tsv` is an **evidence breakdown**; the `tier` column explains why each gene is (or is not) a candidate:

| tier | Meaning |
|---|---|
| `candidate` | right position + motif (highest confidence) |
| `in_target` | peak at an acceptable position, but **no motif support** |
| `motif_only` | motif present, but **position not acceptable** (e.g. distal) |
| `peak_only` | peak present, neither position nor motif criteria met |
| `none` | no peak |

It also reports three-tier statistics and an **ID sanity check**, distinguishing "gene exists but has no peak" (a real biological result) from "ID not present in the annotation" (usually an ID-format or species mismatch), the latter written to `id_not_found.txt`.

---

## Optional steps and how to skip them

`reference`, `normalize-gff`, `signal` and `spacing` are optional. Skipping them never breaks the main path.

| Situation | What to do | Consequence |
|---|---|---|
| No suitable organelle genome | `chromaflow reference --skip-organelle`, or set `ORGANELLE_SOURCES=()` / `SKIP_ORGANELLE=true` | Nuclear genome only; everything else identical |
| Some organelle files missing | nothing — they are skipped with a warning | Add `--strict` if you would rather fail |
| Genome/GTF/chrom.sizes already prepared | skip `reference` entirely, put the files in `REFERENCE_DIR` | They are auto-detected |
| AGAT will not install | run `normalize-gff` anyway | Falls back to the built-in Python engine |
| No need for GFF normalization at all | skip `normalize-gff` | `annotate` uses the plain GFF3 and warns if `exon` rows are absent |
| No GTF needed (downstream only) | `chromaflow reference --skip-gtf` | Only step 02 requires a GTF |
| No R available | use `chromaflow annotate-py` | Same output contract; 06–09 unaffected |
| No motif database | leave `MEME_DB` empty | Tomtom comparison skipped, de novo discovery unaffected |
| No annotation_info table | leave `ANNOTATION_INFO` empty | GO export limited, everything else unaffected |

`chromaflow doctor` reports optional features as `[INFO]` / `[SKIP]` rather than failures, so its summary reflects what actually blocks you.

---

## Sample record file and read pairing

Located automatically: `--samples <file>` > `SAMPLES_FILE` in the config > `<PROJECT_DIR|INPUT_DIR>/samples_<type>.{txt,tsv,csv}` > `samples.{txt,tsv,csv}`.

Three accepted layouts (TAB, comma or whitespace separated):

```
# 1) header-based (recommended) — any column order, optional columns
sample_id   antibody   control_id   replicate   fastq_path
Dof_ChIP    crDOF      Input        1           /data/raw/Dof_rep1
Input       -          -            1           /data/raw/input

# 2) positional (legacy, still supported)
Dof_ChIP    crDOF      Input        1           /data/raw/Dof_rep1

# 3) explicit read files (skips pairing)
sample_id   replicate  fastq_1                  fastq_2
Dof_ChIP    1          /data/raw/A_1.fq.gz      /data/raw/A_2.fq.gz
```

Column aliases (case-insensitive): `sample|name`, `ab|target`, `control|input`, `rep`, `fastq|path|dir`, `r1|read1`, `r2|read2`.

**Read pairing is inferred from the right-most mate token in the filename**, so there is no fixed suffix list to maintain. Covered out of the box:

```
_1/_2   _R1/_R2   .1/.2   -1/-2   _R1_001/_R2_001   _good_1/_good_2
_clean_1/_clean_2   _val_1/_val_2   _1P/_2P   sampleR1/sampleR2
```

Multi-lane files produce multiple rows for the same sample. Files whose mate is absent are treated as single-end with a warning. Names such as `sample_rep1` / `sample_rep2` are **not** mistaken for mates.

For unusual naming schemes, override the detection (`{R}` is replaced by 1 and 2):

```bash
# in the samples file
#@ pair_pattern: _good_{R}.fq.gz
#@ single_end: true

# in chromaflow.config
FASTQ_PAIR_PATTERN="_R{R}_001.fastq.gz"
FASTQ_SINGLE_END=true

# on the command line
chromaflow samplesheet --pair-pattern '_R{R}_001.fastq.gz'
chromaflow samplesheet --single-end
chromaflow samplesheet --list-fastq        # preview pairing, write nothing
```

---

## Configuration and automatic file discovery

The config declares **only what you already have**. Files produced by the pipeline are never written back into it.

```bash
SPECIES_SHORT="Cre"
PROJECT_DIR="/data/projects/CrDof_ChIP"
REFERENCE_DIR="${PROJECT_DIR}/inputs/reference"

GENOME_FASTA_SOURCE=".../genome.fa"         # plain or .gz/.bz2/.xz/.zst — auto-detected
GENOME_GFF3_SOURCE=".../genes.gff3"         # same
ANNOTATION_INFO=""                          # optional
MEME_DB=""                                  # optional
ORGANELLE_SOURCES=("chloroplast:NC_005353.1" "mitochondria:NC_001638.1")   # optional
```

Derived files (`*_combined.fa`, `*.fai`, `chrom.sizes`, `*.gtf`, `*.agat.gff3`) are resolved in this order:

```
explicit value in config  >  manifest registered by step 00/00b
                          >  scan of REFERENCE_DIR  >  default name
```

The manifest lives at `<OUTPUT_DIR>/.chromaflow/reference.manifest`. Directory scanning prefers normalized annotations (`*.agat.gff3` over a plain `*.gff3`) and indexed genomes, so after running `normalize-gff` the annotation step uses the new file with no config edit. Explicit assignments in the config always win, so existing configs keep working unchanged.

---

## Recipes

### A. Standard analysis: find direct targets
```bash
mamba activate ChromaFlow_idr        && chromaflow idr
mamba activate ChromaFlow_motif      && chromaflow motif
mamba activate ChromaFlow_R          && chromaflow annotate
mamba activate ChromaFlow_downstream && chromaflow candidates
chromaflow tracks --gene-list results/downstream/chipseq/candidates/high_confidence_genes.txt
```

### B. Poor reproducibility: fall back to consensus
When IDR leaves only a few dozen peaks, use the nf-core consensus set (higher sensitivity, lower specificity):
```bash
for step in motif annotate candidates; do chromaflow $step --peakset consensus; done
chromaflow lookup --peakset consensus --gene-list my_genes.txt
```

### C. Check whether prior genes are bound
```bash
chromaflow lookup --gene-list my_genes.txt                      # IDR
chromaflow lookup --gene-list my_genes.txt --peakset consensus  # consensus
chromaflow tracks --gene-list results/downstream/chipseq/lookup/candidate_genes.txt \
                  --outdir figures/my_targets
```

### D. Tighten to a high-confidence core set
```bash
chromaflow idr --min-pairs 2
# then rerun motif → annotate → candidates with the same peakset
```

### E. Test the influence of one replicate
```bash
chromaflow idr                                        # read the pairwise table first
chromaflow idr --exclude-rep REP1 -o results_noREP1   # compare peak counts
```

### F. Widen the genomic regions of interest
```bash
chromaflow candidates --features "Promoter,5' UTR,1st Exon,1st Intron,Distal Intergenic"
```

---

## Output guide

```
results/downstream/<type>/
├── idr/
│   ├── idr_peaks_threshold_0.01.bed   # reproducible peaks (signalValue / globalIDR / n_pairs)
│   └── recommended_peaks.bed          # automatically chosen set
├── motif[_peakset]/
│   └── meme-chip_output/              # meme-chip.html, fimo_out_*/, centrimo_out/
├── annotation[_peakset]/
│   ├── peak_annotation_full.tsv       # per peak: gene, feature, signalValue
│   ├── genes.bed                      # gene coordinates (used by tracks / lookup)
│   ├── background_universe.txt        # GO background
│   └── anno_pie.png / dist_to_tss.png / tss_avgprof.png
├── candidates[_peakset]/
│   ├── candidate_targets.tsv          # per-peak verdict (with per-motif counts)
│   ├── motif_<name>_hits.tsv / _genes.txt
│   ├── motif_summary.tsv              # hit counts + CentriMo central enrichment
│   ├── high_confidence_genes.txt      # candidate genes, ranked by strength
│   └── go/                            # input for agriGO / g:Profiler
├── lookup[_peakset]/
│   ├── query_binding.tsv              # evidence breakdown (tier column)
│   └── candidate_genes.txt / bound_genes.txt / id_not_found.txt
└── tracks[_peakset]/                  # per-gene browser shots + tracks.ini
```

---

## Environments

| Environment | Purpose | Key tools |
|---|---|---|
| `nfcore_env` | main pipeline | nextflow, nf-core |
| `ChromaFlow_ref` | reference prep, GFF normalization, GO | samtools, gffread, bedtools |
| `ChromaFlow_agat` | *(optional)* AGAT-based normalization | AGAT |
| `ChromaFlow_idr` | IDR | idr, numpy<1.24 |
| `ChromaFlow_motif` | motif discovery | MEME Suite, bedtools |
| `ChromaFlow_R` | ChIPseeker annotation | R, ChIPseeker, clusterProfiler |
| `ChromaFlow_downstream` | integration, lookup, plotting | deepTools, pyGenomeTracks, bedtools |

`ChromaFlow_agat` is deliberately separate and not installed by default: AGAT's Perl dependency tree often conflicts with the other packages. Install it with `WITH_AGAT=true bash setup_environments.sh` only if you want AGAT specifically.

Run `chromaflow doctor` to check tools and reference files.

---

## Quality assessment and interpretation

### Reading the key metrics

| Metric | Source | Meaning |
|---|---|---|
| NSC / RSC | phantompeakqualtools (printed in 03) | strand cross-correlation of enrichment; NSC > 1.05, RSC > 1.0 are ideal |
| pairwise pass rate | the pair table in 03 | agreement between replicates |
| `n_pairs` | column 8 of the 03 output | how many pairs support that peak |
| signalValue | MACS fold enrichment | **quantitative** enrichment (shown in peak labels) |
| globalIDR | IDR | reproducibility (-log10 scale; larger is better) |
| CentriMo E-value | 04 | whether the motif is **centrally enriched** at summits — positional evidence of direct binding |

### The right order for judging replicate quality

1. **Read the pairwise table first**, not the per-replicate peak count or RSC. A replicate with more peaks is not necessarily more consistent — greater depth or background produces many irreproducible marginal peaks.
2. **Compare NSC against Input.** If each ChIP's NSC is close to the Input's, overall enrichment is limited and an RSC "PASS" means little.
3. **Tighten with `--min-pairs`, not by deleting replicates.** Removing a replicate destroys pairs; `--min-pairs 2` raises the evidence bar while keeping all data.

### Layers of evidence

Each layer answers a different question — **do not conflate them**:

- **MACS / IDR** → is there reproducible enrichment here (statistical test vs Input)
- **ChIPseeker** → which gene, and which part of it, does the peak fall in
- **FIMO** → does the peak sequence contain the motif (**match quality, not binding strength**)
- **CentriMo** → is the motif centrally enriched (supporting evidence for direct binding)
- **Candidate call** → the intersection of positional and sequence evidence

So a gene can be "candidate grade" yet show a weak track signal — candidacy is about **position and sequence**, while strength must be read from signalValue and the difference track. Moving from "candidate" to "direct functional target" still requires independent validation such as RNA-seq differential expression, EMSA or ChIP-qPCR.

---

## Adapting to other species

1. Edit `SPECIES_NAME` / `SPECIES_SHORT`, the source file paths, `MACS_GSIZE` and `ORGANELLE_SOURCES` in the config. Species without a suitable organelle genome simply leave `ORGANELLE_SOURCES=()` or pass `--skip-organelle`.
2. Run `chromaflow normalize-gff` once and confirm that UTR proportions look sane in the annotation pie chart.
3. Point `MEME_DB` at a database for the right clade (JASPAR plants / vertebrates …), or leave it empty to skip Tomtom.
4. Set `KEGG_ORGANISM` to the KEGG species code; if gene IDs do not match KEGG's, that part is skipped automatically.
5. **Gene ID format**: `lookup` strips only trailing version suffixes (`.1`, `.v6.1`), so `AT1G01010.1`, `ENSG00000141510.17`, `LOC_Os01g01010.1` and `Cre02.g075050_4532` all work; if an ID still maps to multiple genes after normalization it is flagged as ambiguous rather than guessed. Watch the **ID check** line in the lookup output to catch format mismatches early.

---

## Dependency versions

nf-core/chipseq 2.1.0 · nf-core/atacseq 2.1.2 · MACS3 · IDR · MEME Suite 5.x · ChIPseeker 1.42 · deepTools 3.5 · pyGenomeTracks · AGAT *(optional)* · bedtools 2.31
