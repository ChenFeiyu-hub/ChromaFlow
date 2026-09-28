# ChromaFlow 1.2.7

[中文指南](README.zh-CN.md) · [Offline deployment](docs/OFFLINE.md) · [Downstream commands](docs/DOWNSTREAM.md) · [Motif selection and tracks](docs/SELECTION.md) · [Release notes](CHANGELOG.md)

ChromaFlow wraps nf-core/chipseq 2.1.0 (ChIP-seq and DAP-seq) and nf-core/atacseq 2.1.2, then provides IDR, motif discovery, peak annotation, candidate-gene integration, signal plots, gene tracks and GO exports. Source code, project configuration and analysis outputs live separately.

Version 1.2.7 shares motif filters between candidates and lookup, passes verified
selections to tracks, fixes CentriMo local/central labels, and sizes dense gene/motif
tracks automatically. Exact gene IDs and the 1.2.6 annotation policy are preserved.
See [selection and filtering commands](docs/SELECTION.md) and the
[downstream guide](docs/DOWNSTREAM.md).

## Install the software

Use Linux or WSL2 Ubuntu, Bash, Python 3 and standard GNU command-line tools. Prepare/download Linux containers on Linux/WSL, not in Windows PowerShell. In WSL, `D:\biosoft` is `/mnt/d/biosoft`.

Install the current repository version:

```bash
mkdir -p "$HOME/biosoft"
git clone https://github.com/ChenFeiyu-hub/ChromaFlow.git "$HOME/biosoft/ChromaFlow"
export CHROMAFLOW_HOME="$HOME/biosoft/ChromaFlow"
export PATH="$CHROMAFLOW_HOME/bin:$PATH"
hash -r
chromaflow version
```

Record `git -C "$CHROMAFLOW_HOME" rev-parse HEAD` to identify the exact revision.
For an offline server, transfer a checked-out copy or an available release package.
If you already have the versioned ChromaFlow archive and its checksum, install it with:

```bash
cd ~/nfcore-upload   # directory containing the transferred release
sha256sum -c ChromaFlow-1.2.7.tar.gz.sha256
mkdir -p "$HOME/biosoft/ChromaFlow-1.2.7"
tar -xzf ChromaFlow-1.2.7.tar.gz \
  -C "$HOME/biosoft/ChromaFlow-1.2.7" --strip-components=1
export CHROMAFLOW_HOME="$HOME/biosoft/ChromaFlow-1.2.7"
export PATH="$CHROMAFLOW_HOME/bin:$PATH"
chromaflow version
```

Install into a new version directory. Add the two `export` lines to your shell startup file after checking the release. Existing projects continue to use their own config; do not replace a working project config with the template.

## Choose an environment plan

```bash
# Inspect first; no packages or environments are changed.
chromaflow setup --mode offline --components runtime --dry-run

# Choose the runtime appropriate to this computer:
chromaflow setup --mode offline --components runtime
# OR, on a connected Linux/WSL bundle preparation computer:
chromaflow setup --mode prepare --components runtime
# OR, for online execution using nf-core's Conda profile:
chromaflow setup --mode conda --components runtime

# Optional, independent ChromaFlow analysis environments:
chromaflow setup --components downstream
# Options: --skip-r, --with-agat, --solver mamba|micromamba|conda
```

| Mode | Runtime environment | Included packages | Used for |
|---|---|---|---|
| `offline` (default) | `ChromaFlow_runtime` | Python, Java | Run transferred container bundles |
| `prepare` | `nfcore_env` | Python, Java, nf-core/tools 3.3.2 | Download and package offline bundles |
| `conda` | `nfcore_env` | Python, Java, Nextflow 25.10.4 | Online nf-core Conda execution |

`--components all` is the default; use `runtime` for a minimal installation. Existing named environments are kept, not updated or deleted. If reusing `nfcore_env` for a different purpose, check that the required tools are present. Apptainer/Singularity is a system prerequisite and is not installed by this script.

**Offline mode describes nf-core execution. The Conda environment installer still needs package access.** A workflow bundle includes Nextflow, plugins, workflow code/configuration and container images. It excludes reference genomes, reads, motif databases and ChromaFlow's downstream environments. Reuse working downstream environments on an offline server, or provision them separately on a compatible Linux system.

## Prepare and transfer both nf-core bundles

No separate manual `nf-core pipelines download` or plugin-install commands are needed:

```bash
conda activate nfcore_env  # created with --mode prepare
chromaflow offline prepare --type chipseq --bundle "$HOME/nfcore-offline/chipseq-2.1.0"
chromaflow offline prepare --type atacseq --bundle "$HOME/nfcore-offline/atacseq-2.1.2"
chromaflow offline check --bundle "$HOME/nfcore-offline/chipseq-2.1.0" --verify --runtime-check
chromaflow offline check --bundle "$HOME/nfcore-offline/atacseq-2.1.2" --verify --runtime-check
chromaflow offline pack --bundle "$HOME/nfcore-offline/chipseq-2.1.0" --archive "$HOME/chipseq-2.1.0.tar"
chromaflow offline pack --bundle "$HOME/nfcore-offline/atacseq-2.1.2" --archive "$HOME/atacseq-2.1.2.tar"
```

Transfer the two archives, their `.sha256` files and the ChromaFlow release. Follow [the deployment guide](docs/OFFLINE.md) for shared `/opt` installation, target-side validation, permissions, project configuration, troubleshooting and version upgrades. DAP-seq shares the chipseq bundle. Existing valid bundles for these versions can be reused after upgrading ChromaFlow.

## Start a project

```bash
chromaflow init "$HOME/projects/MyChIP" --type chipseq
cd "$HOME/projects/MyChIP"
# Edit chromaflow.config and samples_chipseq.txt before continuing.
conda activate ChromaFlow_ref
chromaflow samplesheet --type chipseq --strict
chromaflow reference
chromaflow normalize-gff --engine builtin  # when annotation needs normalization
chromaflow reference                     # refresh GTF after a changed annotation
chromaflow doctor --type chipseq

conda activate ChromaFlow_runtime
chromaflow run --type chipseq --offline \
  --bundle /opt/chromaflow/bundles/chipseq-2.1.0 --dry-run
chromaflow run --type chipseq --offline \
  --bundle /opt/chromaflow/bundles/chipseq-2.1.0 --yes
```

The template requires your species prefix and source FASTA/GFF3. The species label never determines gene-ID parsing. Leave `MACS_GSIZE` and `READ_LENGTH` unset to use measured values. `samplesheet` measures read length; `reference` can estimate effective genome size with unique k-mers. Derived values and reference paths are recorded under `<OUTPUT_DIR>/.chromaflow/`. Explicit config values take precedence, including old explicit GFF/GTF paths.

The default output anchor is `<PROJECT_DIR>/outputs`; set `OUTPUT_DIR="${PROJECT_DIR}"` to place `results/`, `logs/`, `samplesheets/` and `work/` directly under the project. Reference preparation and normalization are unnecessary when the project already has validated references. Normalization validation and downstream GTF/TxDb acceptance are separate checks; see [GFF normalization](GFF_NORMALIZATION.md).

## Compare IDR and consensus downstream analyses

`--peakset idr` selects the IDR 0.01 set. `--peakset consensus` uses nf-core consensus peaks without IDR filtering. Neither requires rerunning nf-core. Run IDR once for the former, then use the same peakset throughout each downstream branch:

```bash
conda activate ChromaFlow_idr
chromaflow idr --type chipseq --min-pairs 1

for PEAKSET in idr consensus; do
  conda activate ChromaFlow_R
  chromaflow annotate --type chipseq --peakset "$PEAKSET"
  conda activate ChromaFlow_motif
  chromaflow motif --type chipseq --peakset "$PEAKSET"
  conda activate ChromaFlow_downstream
  chromaflow candidates --type chipseq --peakset "$PEAKSET"
  chromaflow signal --type chipseq --peakset "$PEAKSET"
  conda activate ChromaFlow_ref
  chromaflow go-prep --type chipseq --peakset "$PEAKSET" --set promoter
done
```

This IDR wrapper accepts a single biological condition with at least two replicate peak files. It rejects mixed conditions instead of comparing unrelated samples. Missing/empty requested IDR sets fail explicitly. Use `--peakset idr:0.05` to choose a different computed threshold deliberately. Full commands and interpretation limits are in [DOWNSTREAM.md](docs/DOWNSTREAM.md).

IDR is optional. For a consensus-only analysis, skip `idr` and run the downstream
commands with `--peakset consensus`. Set `DEFAULT_PEAKSET=consensus` in the project
config to omit that repeated option; explicit `--peakset` always takes precedence.
The template default remains `idr`.

## Annotate, select candidate genes and plot the same evidence

The project config controls annotation and downstream location selection separately:

```bash
# Example project settings; choose the promoter window for your analysis.
DEFAULT_PEAKSET=consensus
PROMOTER_UP=3000
PROMOTER_DOWN=0
TARGET_FEATURES="Promoter,5' UTR,1st Exon"
TARGET_MAX_TSS_DIST=0
FIMO_PVAL=1e-4
```

`PROMOTER_UP/DOWN` define the strand-relative inclusive TSS window. ChIPseeker's
priority is Promoter > 5′ UTR > 3′ UTR > Exon > Intron > Downstream > Intergenic;
feature classification and nearest-TSS gene assignment are different fields.
`TARGET_FEATURES` selects accepted categories afterwards. The optional upstream
distance rescue applies to the already assigned gene and is disabled with zero.
Changing the promoter window requires rerunning annotation.

The built-in promoter audit also reports **every peak overlapping the whole
configured promoter**, including peaks assigned to other genes. Inspect
`promoter_peak_overlaps.tsv`, `promoter_peak_summary.tsv`, and `annotation_policy.json`.
This audit adds regional evidence without reassigning peaks or promoting candidates.

After `annotate` and `motif` have completed for the same branch:

```bash
conda activate ChromaFlow_downstream
chromaflow candidates --type chipseq --peakset consensus
chromaflow lookup --type chipseq --peakset consensus \
  --gene-list genes_of_interest.txt --level candidate \
  --outdir results/downstream/chipseq/lookup_selected
chromaflow tracks --type chipseq --peakset consensus \
  --selection results/downstream/chipseq/lookup_selected \
  --outdir results/downstream/chipseq/tracks_selected

# Optional: select a specific model and query/plot in one command.
chromaflow lookup --type chipseq --peakset consensus \
  --gene-list genes_of_interest.txt --level candidate --motif MEME-4 --plot \
  --outdir results/downstream/chipseq/lookup_MEME4
```

Lookup previews exact query/annotation IDs and asks **y/N**. It does not strip
suffixes or infer aliases; `--yes` explicitly accepts exact matching for unattended
runs. `--level all` (default) retains all query rows and plots exactly matched genes;
`in_target` selects accepted-location evidence; `candidate` also requires a passing
motif hit. `--plot` follows that level. Unmatched IDs remain in `id_check.tsv`.

Candidates and lookup share these optional filters:

| Option | Config | Meaning |
|---|---|---|
| `--motif MEME-4,MEME-7` | `MOTIF_IDS` | Exact run-local model labels or raw IDs |
| `--pval 1e-5` | `FIMO_PVAL` | Site-level FIMO p-value maximum; default 1e-4 |
| `--qval 0.05` | `FIMO_QVAL` | Site-level FIMO q-value maximum |
| `--motif-evalue 0.05` | `MOTIF_EVALUE_MAX` | Discovery model E-value maximum |
| `--centrimo-evalue 0.05` | `CENTRIMO_EVALUE_MAX` | Positional-enrichment E-value maximum |
| `--min-signal 2` | CLI default 0 | Minimum peak signal score |

Optional config filters are empty/disabled by default. Filters are applied to all
peaks before selecting each gene's best evidence. MEME-4 is not necessarily
`fimo_out_4`; selecting an unscanned model or requesting missing statistics fails
explicitly. These commands filter existing FIMO output, not unscanned DNA.

`selected_genes.txt`, `selected_evidence.tsv`, `motif_catalog.tsv`, and `selection.json`
record the selection and its provenance. `tracks --selection` checks file hashes,
branch and promoter settings, then uses the saved genes/motif hits while retaining
all regional peaks and neighboring genes. `candidates --plot` also uses its selected
genes. The full `candidate_targets.tsv` includes failing peaks: use the selected
list for candidate-only plots. The legacy `high_confidence_genes.txt` filename is an
alias for computational candidates, not a claim of experimental validation.

CentriMo reports local/central/unknown mode and passes-threshold/above-threshold/
not-reported status; its display threshold is E≤0.05 and does not filter candidates
unless explicitly requested. Significant local enrichment is not automatically
central enrichment. Track heights adapt to transcript/motif counts; set
`TRACK_GENE_HEIGHT` and `TRACK_MOTIF_HEIGHT` to `auto` or a positive height in cm.
Peak labels preserve three significant digits. Full semantics and examples are in
[Selection and tracks](docs/SELECTION.md).

## Outputs and logs

```text
<OUTPUT_DIR>/
  results/chipseq/                         nf-core results
  results/downstream/chipseq/
    idr/                                  pairwise IDR and threshold sets
    annotation/   annotation_consensus/   separate branches
    motif/        motif_consensus/
    signal/       signal_consensus/       created when signal is requested
    go/           go_consensus/
    candidates/   candidates_consensus/
    lookup/       lookup_consensus/
    tracks/       tracks_consensus/
  logs/downstream/chipseq/
    <timestamp>_<command>_<peakset>_<unique>.log
```

Every downstream command through `chromaflow` automatically records combined terminal stdout/stderr, command, version, config path and final exit code. You do not need to add `tee`. Each invocation gets a new log; `-o` changes that step's results directory, not the project log directory. Main nf-core and reference logs retain their existing locations under `logs/`. Argument/config-resolution errors that occur before a project can be identified are printed to the terminal. User interruption can leave a partial log without an end marker.

Software-generated messages, help and plot labels are English. User paths, gene names and raw biological data are preserved. Existing logs are not rewritten. Calling internal scripts directly bypasses the unified CLI logger.

`annotation*/dist_to_tss.png` now shows explicit upstream/downstream distance bins, percentages and counts. `dist_to_tss.tsv` records the bin boundaries and denominator. The plot summarizes annotation distances; it does not measure enrichment against a genomic background. It can be redrawn from an existing ChIPseeker table without repeating annotation:

```bash
Rscript "$CHROMAFLOW_HOME/lib/steps/plot_tss.R" \
  results/downstream/chipseq/annotation/peak_annotation_full.tsv \
  results/downstream/chipseq/annotation/tss_plot_v1.2.7
```

## Command reference

| Command | Purpose |
|---|---|
| `setup`, `init`, `doctor` | Environments, project creation, checks |
| `offline prepare/check/pack` | Build, verify and package local nf-core bundles |
| `samplesheet`, `reference`, `normalize-gff`, `run` | Inputs, reference preparation, nf-core |
| `idr`, `motif`, `annotate`, `annotate-py` | Reproducibility, motifs, annotation |
| `candidates`, `lookup`, `tracks` | Candidate integration, prior-gene queries, tracks |
| `signal`, `go-prep`, `go-plot`, `spacing` | Signal plots, GO export/plots, motif spacing |

Use `chromaflow help <command>`, `chromaflow help config`, `chromaflow help samples` and `chromaflow help envs`. `annotate-py` is a simpler alternative and does not reproduce every ChIPseeker field or R plot; use the R path for the full annotation workflow.

## Validation

```bash
python3 -m unittest discover -s tests -v
```

The v1.2.7 validation run passed 113 tests, covering offline bundles, CLI execution,
logs, peakset routing, setup plans, GFF preservation, TxDb identity aliases,
interval-score recovery, shared motif filters, CentriMo labels and selection
provenance. IDR postprocessing integration tests use real bedtools with deterministic
IDR output; they skip when bedtools is unavailable. Linux/Bash-specific checks also
skip on unsupported platforms. Real pyGenomeTracks rendering was checked separately
with synthetic neighboring genes and dense transcript models. These checks do not
replace a complete run on your own samples.
