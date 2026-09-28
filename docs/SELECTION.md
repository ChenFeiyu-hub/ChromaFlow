# Shared motif selection and tracks (1.2.7)

`candidates` and `lookup` use the same motif filters and location policy. All
filters operate on existing MEME-ChIP/FIMO results; they do not run motif discovery
again. Gene IDs remain exact, case-sensitive IDs. No suffix conversion is applied.

## One project, one selected branch

Set `DEFAULT_PEAKSET=consensus` in the project's existing `chromaflow.config` to use
consensus for downstream commands without repeating `--peakset consensus`.
An explicit `--peakset` overrides this default. Projects without this setting
retain the previous `idr` default. This choice does not remove other results.

Activate `ChromaFlow_downstream`, then run in the project directory:

```bash
chromaflow candidates --type chipseq
chromaflow lookup --type chipseq --gene-list genes.txt --level candidate --plot
```

Lookup shows query and annotation ID examples and asks y/N. `--yes` explicitly
accepts exact matching in unattended runs. `--level candidate` selects candidate
genes for both the table and plots; `--level in_target` selects genes with accepted
location evidence. The default `--level all` keeps all queries in the table and
plots all exactly matched genes. Unmatched IDs remain in `id_check.tsv`.

## Optional filters

| CLI option (candidates and lookup) | Config default | Meaning |
|---|---|---|
| `--motif MEME-4,MEME-7` | `MOTIF_IDS=""` | Exact motif labels or raw motif IDs from this run; comma separated |
| `--pval 1e-5` | `FIMO_PVAL=1e-4` | Maximum site-level FIMO p-value |
| `--qval 0.05` | `FIMO_QVAL=""` | Maximum site-level FIMO q-value |
| `--motif-evalue 0.05` | `MOTIF_EVALUE_MAX=""` | Maximum discovery E-value recorded in MEME-ChIP's HTML payload |
| `--centrimo-evalue 0.05` | `CENTRIMO_EVALUE_MAX=""` | Maximum CentriMo positional-enrichment E-value |
| `--min-signal 2` | CLI default 0 | Minimum peak signal score |

Empty optional settings disable that filter. Filters combine with AND; a list of
motif IDs means any of those models can support a candidate. Position remains
controlled by `TARGET_FEATURES` and `TARGET_MAX_TSS_DIST`. These filters are applied
to every peak before lookup selects a best peak; a weaker peak with the requested
motif is not hidden by a stronger peak with another motif.

MEME-4 is a run-local label, **not** `fimo_out_4`. MEME-ChIP normally scans cluster
representatives, so some discovered models have no FIMO output. Explicitly selecting
an unscanned model fails with an explanation; it is not reported as zero motif
occurrences. Missing q-values or discovery E-values also fail when the corresponding
filter is requested. `motif_catalog.tsv` records raw IDs, labels, discovery E-values,
available scan paths, PWM hashes when present, and CentriMo status.

Tightening a threshold filters stored hits. Relaxing it cannot recover hits that
FIMO did not export; a new scan is required. The software prints this limitation.
Models missing from an HTML display may still be scanned: absent discovery statistics
are left unavailable, never invented from site-level p-values.

Use separate result directories to retain several selections:

```bash
chromaflow lookup --type chipseq --gene-list genes.txt --level candidate \
  --motif MEME-4 --outdir results/downstream/chipseq/lookup_MEME4 --plot

chromaflow lookup --type chipseq --gene-list genes.txt --level candidate \
  --qval 0.05 --outdir results/downstream/chipseq/lookup_q005
```

## Pass the selection to tracks

```bash
chromaflow tracks --type chipseq \
  --selection results/downstream/chipseq/lookup_MEME4 \
  --outdir results/downstream/chipseq/tracks_MEME4
```

`--selection` accepts a directory or its `selection.json`. It loads the selected
genes and the stored filtered motif-hit BED, checks pipeline/peakset, promoter
parameters, and source/output hashes. Run it with the original project's config
and branch. Do not combine it with `--genes`, `--gene-list` or `--motifs`.

Automatic selection outputs:

- `selected_genes.txt`: genes chosen by the command and requested level.
- `selected_evidence.tsv`: passing motif occurrences associated with selected
  genes/peaks (all relevant peaks, not only lookup's one-line best-peak summary).
- `motif_hits.tsv/.bed`: all hits passing the same motif filters, including
  neighboring context hits. Nearby motif hits are not automatically assigned to
  the query gene; use `selected_evidence.tsv` for those associations.
- `motif_catalog.tsv`: model identities, statistics, scan availability and status.
- `selection.json`: branch, parameters, filters, source/output hashes and version.

`candidates --plot` uses its selected candidate genes. `candidate_genes.txt` is also
available. `high_confidence_genes.txt` remains a compatibility alias for the same
computational candidates; its historical name is not evidence of experimental
validation. `candidate_targets.tsv` contains both passing and failing peaks: do not
treat its entire gene column as the selected candidate list.

Tracks always keep **all peaks in the selected branch** and neighboring gene
structures in the plotted region. Only the gene list and motif overlay are selected.
The existing `tracks --genes G1,G2 --motifs` and `--gene-list FILE` modes remain
available; those use the project's FIMO p-value cutoff, without additional selection
filters. Use `--selection` when tables and plots must share an advanced filter.

Gene and motif track heights now grow with the region's transcript/hit count.
`TRACK_GENE_HEIGHT=auto` and `TRACK_MOTIF_HEIGHT=auto` preserve every isoform; either
may be set to a positive height in cm. Per-region layouts are in `.layout/`.
Peak signal labels use three significant digits. Coverage and the log2 ratio are
visual evidence; the latter uses pseudocount 1 and 50 bp bins and is not a significance
test. It should not be equated to MACS peak signalValue.

## CentriMo interpretation

CentriMo in local mode searches across positions. A significant local result is
not automatically a centrally enriched result. The report uses a **0.05 E-value
threshold for labels**; this does not filter candidate genes unless
`--centrimo-evalue` is explicitly set.

- `centrimo_mode`: `local`, `central`, or `unknown` if command metadata is unavailable.
- `positional_enrichment`: `passes_threshold`, `above_threshold`, or `not_reported`.
- `centrally_enriched`: `yes/no` only for a reported result in known central mode;
  otherwise `not_assessed`.

An E-value of 1.8 is above the report threshold. A missing row may have been omitted
by CentriMo's output threshold and is not proof of absent enrichment. E-values are
not p-values; neither is a gene-specific significance test. See the official
[CentriMo documentation](https://meme-suite.org/meme/doc/centrimo.html) and
[FIMO documentation](https://meme-suite.org/meme/doc/fimo.html).

No new dependencies are required for this maintenance release. Existing reference,
nf-core, annotation and motif results can be reused. Rerun candidates/lookup/tracks
to regenerate their corrected reports. Their complete English terminal output is
still logged automatically under `logs/downstream/<type>/`.
