# Downstream analysis and logging

For shared motif filters, `lookup --level ... --plot`, and verified
`tracks --selection` handoff, see [Selection and filtering](SELECTION.md).

Run after nf-core finishes. Commands assume the current directory contains the intended `chromaflow.config`; add `-c /absolute/path/to/project.config` when using another file. Replace `chipseq` with `dapseq` or `atacseq` to select that project's nf-core results.

## Two separate peak sets

| Selection | Meaning | Result-directory suffix |
|---|---|---|
| `--peakset idr` / `idr:0.01` | IDR threshold 0.01 | none |
| `--peakset idr:0.05` | IDR threshold 0.05 | `_idr0.05` |
| `--peakset consensus` | nf-core consensus, without IDR filtering | `_consensus` |
| `--peakset /path/to/peaks.bed` | Explicit custom intervals | `_custom` |

The template default is `idr` at 0.01. Set `DEFAULT_PEAKSET=consensus` in the project config to change the default; an explicit `--peakset` takes precedence. `IDR_THRESHOLD` in the legacy config only controls the optional support summary. The recommendation symlink from IDR is informational. An explicit `idr` request does not silently use that symlink or fall back to consensus.

Use the same peakset for motif, annotation, candidates, lookup, tracks, signal, spacing and GO commands. Custom files share the `_custom` suffix, so use separate projects/output anchors when comparing several custom sets. `-o` redirects only one step's outputs; it does not rewire subsequent steps' inputs.

## IDR once, then both branches

```bash
# Save this block as a Bash script and run it with bash; stop on any failed step.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
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

This uses `TARGET_FEATURES` and `FIMO_PVAL` from the project config. The new
template selects `Promoter,5' UTR,1st Exon`; an existing explicit setting is
preserved. Add `--features` only when intentionally overriding that setting.

IDR uses replicate peak files in the default BWA/MACS3 merged-library layout. The wrapper checks nf-core `_REP<number>` filenames and requires exactly one biological condition. It is not a grouping interface for multiple conditions. There are `n*(n-1)/2` pairs; `--min-pairs 1` retains merged regions supported by any pair. For three replicates, `--min-pairs 2` or `3` is stricter, but pair support is counted over the merged region and does not mean every base is supported by every pair. Different min-pairs runs overwrite the same IDR directory; use separate output anchors to compare them.

IDR pair failure aborts the command. `idr.status` records `running`, `failed` or `complete`; downstream commands reject a marked incomplete run. Empty threshold sets are allowed as an IDR result, but cannot be consumed by downstream steps requiring peaks. Explicitly choose another threshold or consensus after reviewing the results.

The exported BED has eight columns:

```text
chrom start end name signalValue strand globalIDR n_pairs
```

`signalValue` is the maximum MACS column-7 value among retained replicate input peaks overlapping the merged region. The IDR fit still ranks by p-value; its p-value output is not relabeled as enrichment. `globalIDR` is the maximum passing `-log10(global IDR)` in the merged region, not a new fitted IDR for the whole region. NarrowPeak and broadPeak IDR output columns differ, as documented by the [IDR project](https://github.com/nboley/idr#output). These scores help rank peaks within an analysis; they are not a differential-binding test across conditions.

## What annotation and motif analysis measure

R annotation uses the validated GFF3 to construct TxDb, assigns nearest genes/TSS and exports `peak_annotation_full.tsv`, `genes.bed`, target lists and plots. Feature overlap and nearest-gene assignment can refer to different genes; inspect the corresponding fields before interpreting a target relationship. `annotate-py` is a simplified alternative and does not provide every R field/plot.

The TSS chart counts rows with finite `distanceToTSS` from the annotation table. Negative and positive distances are upstream and downstream relative to the gene strand. Zero is its own bin; all other bin endpoints are documented in `dist_to_tss.tsv` (`min_abs_bp`/`max_abs_bp`, inclusive). Missing distances are excluded from the percentage denominator and reported. Peaks dropped earlier because their sequence is absent from TxDb are reported in the annotation log, not counted as chart missing values. Distance fractions are descriptive, not genomic-background enrichment p-values.

```bash
# Redraw an existing R annotation table into a new directory. Base R suffices.
Rscript "$CHROMAFLOW_HOME/lib/steps/plot_tss.R" \
  results/downstream/chipseq/annotation/peak_annotation_full.tsv \
  results/downstream/chipseq/annotation/tss_plot_v1.2.7
```

MEME-ChIP uses 200 bp windows: summit-centered for a narrowPeak input, midpoint-centered for the exported IDR/consensus BED inputs. FIMO results describe those extracted sequences; this is not a whole-genome FIMO scan. `MEME_DB` can point to a local motif database; without it, known-motif matching is unavailable. `MEME_SEARCHSIZE=0` uses all sequence and can be expensive for large consensus sets.

Candidate selection combines an accepted location feature with a motif hit. Use `--max-tss-dist N` only when you intentionally want the optional upstream-distance rescue rule. A candidate is a computational hypothesis, not proof of transcriptional regulation.

## Additional commands

```bash
conda activate ChromaFlow_downstream
chromaflow lookup --type chipseq --peakset consensus --gene-list genes_of_interest.txt --plot
chromaflow tracks --type chipseq --peakset consensus --gene-list genes_of_interest.txt
chromaflow tracks --type chipseq --peakset idr --genes GENE1,GENE2

conda activate ChromaFlow_ref
chromaflow go-prep --type chipseq --peakset consensus --set all
chromaflow go-prep --type chipseq --peakset idr --set dist
chromaflow go-plot --type chipseq --peakset consensus --input enrichment_results.tsv --top 20
chromaflow spacing --type chipseq --peakset consensus
```

`lookup` can inspect location-only evidence when motif results are unavailable. `tracks` requires a gene list and appropriate bigWigs/GTF. `signal` produces `peak_midpoint_heatmap.png`, `peak_midpoint_profile.png`, TSS profiles when gene coordinates exist, and correlation/PCA when at least two bigWigs exist. It currently discovers all bigWigs under the default pipeline location; inspect the printed sample list before interpreting replicate correlations. Alternative aligner layouts may require adapting paths and are not certified by these commands.

## Review every gene in a prior list

```bash
chromaflow lookup --type chipseq --peakset consensus --gene-list genes_of_interest.txt --level all
chromaflow tracks --type chipseq --peakset consensus --gene-list genes_of_interest.txt --motifs --flank 5000
# Repeat with --peakset idr after that branch's IDR, annotation and motif steps succeed.
```

`lookup` reads the full annotation and FIMO results directly, groups peaks by the
assigned nearest-TSS gene, and reports one best peak per query gene. It prioritizes
evidence tier, then signal. `has_peak=no` does not establish absence of regional
signal. IDs match exactly, including punctuation, suffixes and case. The command
shows examples from the query and annotation namespaces and asks y/N before
generating results. `--yes` explicitly accepts this exact comparison in unattended
runs; it never permits automatic conversion. Zero exact matches abort. Missing IDs
remain `id_status=not_found` in the full query table and `id_check.tsv`.
`lookup --plot` follows `--level`: the default `all` plots every exactly matched
query gene, including genes without assigned peaks; `in_target` and `candidate`
plot only their selected genes. Unmatched IDs stay in the ID reports. The motif
overlay uses the stored shared filters; see [Selection and tracks](SELECTION.md).

`tracks --motifs` reuses the selected branch's FIMO files at `FIMO_PVAL`, converts
1-based inclusive hits to BED coordinates, and exports `motif_hits.bed` and
`motif_hits.tsv`. Gray input-window intervals are recovered from `peak_sequences.fa`
when present, rather than from extraction requests that may have been skipped.
No-hit regions outside those windows were not tested. The default pipeline extracts
forward genomic sequence; reverse-oriented window headers require explicit
coordinate conversion and are rejected. `plot_status.tsv` records each requested
gene's plotted/missing_coordinates/render_failed status. Review missing IDs; a
successful invocation with some missing coordinates does not mean every gene plotted.

The gene track contains neighboring genes, and peak/coverage tracks are independent
of candidate feature filtering. Thus a peak assigned to gene B can still appear in
gene A's plotted region. ChIPseeker's feature label and its nearest-TSS assignment
need not refer to the same gene. Overlap with another gene body does not automatically
exclude a promoter label: promoter has annotation priority. Conversely, a peak
classified outside `TARGET_FEATURES` may fail the candidate rule, and a peak assigned
to B will not become an assigned peak for A just because it lies upstream of A.

`TARGET_MAX_TSS_DIST=N` adds an upstream-distance rule for the **already assigned**
gene, consistently in candidates and lookup. `candidates --max-tss-dist N` overrides
it for that invocation only; use config when the two commands should share a rule.
It does not enumerate all nearby promoters. Keep it unset/zero unless intentionally
using rescue. The shared feature parser recognizes `Exon (... exon 1 of N)` as
`1st Exon`, and likewise intron rank 1; rank follows the annotation, not ID syntax
or left-to-right genomic order. First exons already labeled Promoter/UTR by
ChIPseeker retain those higher-priority categories.

## Complete promoter-region audit (built in)

`annotate` explicitly uses gene-level nearest TSS, overlap=TSS, both peak strands,
and priority Promoter > 5UTR > 3UTR > Exon > Intron > Downstream > Intergenic. It
preserves raw labels, appends `feature_class`, and verifies Promoter labels against
the configured distance window. `annotation_policy.json` records settings and a
table hash. Candidates, lookup and tracks reject a changed window or altered table;
legacy tables without this record warn and should be regenerated for acceptance.

Annotate audits **all genes**; lookup and tracks audit their query genes using the
same code and `PROMOTER_UP/DOWN`. It checks any positive intersection of the selected
peakset with the entire strand-aware promoter, independently of ChIPseeker's assigned
gene or TARGET_FEATURES. Gene-level TSS comes from genes.bed (the same gene model
used for annotation); this does not enumerate every transcript-specific alternative
TSS. Inclusive `[-UP,+DOWN]`, including TSS, is exported as 0-based half-open BED:
positive strand `[tss0-UP,tss0+DOWN+1)`; negative strand `[tss0-DOWN,tss0+UP+1)`.
Coordinates are clipped to zero and, when available, reference sequence lengths.

- `promoter_regions.bed`: exact windows, with gene IDs and strand.
- `promoter_peak_overlaps.tsv`: every gene/peak overlap and its length, original
  assigned gene IDs/labels, self/other/unannotated relation, and overlapping gene
  bodies. `promoter_overlap_gene_body_ids` limits body overlap to the peak's
  intersection with the query promoter; `overlapped_gene_body_ids` uses the full peak.
- `promoter_peak_summary.tsv`: counts per gene, including genes with no overlaps
  and missing query IDs. Other-assigned counts do not mean every peak is misassigned.
- `promoter_audit.json`: settings, coordinate convention, input hashes and scope.

The audit is supplementary regional evidence, not peak reassignment or automatic
candidate promotion. One peak can overlap several promoters. IDR audits show only
the selected IDR peaks; consensus is the comparison branch, and both reuse the same
bigWig signal. Tracks mark the promoter window and expand context when --flank is
shorter than its configured extents, so the full window remains visible.

`MAKE_PLOTS=false` disables R annotation figures while retaining tables/audits.
`TRACK_AUTO_TOPN` and `TRACK_GENES_TOPN` are legacy settings and do not limit an
explicit gene list; they are absent from the new template. `--topn-tracks` is retired.

IDR and consensus motif discovery are independent: a label such as MEME-1 need not
denote the same motif across the two branches. Check motif matrices/logos or database
matches before comparing motif identities. Neither motif occurrence nor promoter
proximity alone establishes direct regulation.

`spacing` requires at least two FIMO result files. GO preparation exports both raw IDs and a legacy `.stripped.txt` variant with a trailing numeric suffix removed. **Prefer raw IDs unless you have verified the external database's identifier namespace.** Suffix removal is not a universal cross-species ID conversion. Use matching foreground/background namespaces. `go-prep` does not itself submit anything to a website. `KEGG_ORGANISM=""` disables the R annotation step's online KEGG query.

## Automatic logs

```text
<OUTPUT_DIR>/logs/downstream/<type>/<timestamp>_<command>_<peakset>_<unique>.log
```

All downstream CLI commands automatically tee stdout/stderr. Start time, software version, config path, shell-escaped command and final exit status are recorded. Tool diagnostics are visible rather than discarded for required IDR, deepTools, bedtools and track-rendering work. Optional R analyses may still emit warnings and continue; inspect them even when the exit code is zero.

The logger starts after resolving a usable project config. Malformed CLI/config errors before this point are terminal-only. Killing the process can leave a log without an end marker. The `Rscript plot_tss.R` standalone helper is outside the CLI logger; its generated TSV makes the plotted counts inspectable. Logs do not retroactively recover output from older runs.

New directories are created only when their corresponding steps are requested; old empty directories from earlier versions are not required inputs. Existing generated files can coexist with new ones (for example old `peak_summit_*.png` and new `peak_midpoint_*.png`); inspect filenames and timestamps when assembling figures.
