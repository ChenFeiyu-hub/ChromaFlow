# TSS boundaries and non-IDR metadata correction for ChromaFlow 1.2.8

This local update follows the five previously installed 1.2.8 hotfixes and
records `tss-idr-metadata-v1`. No new environment or container is needed.

## What changes

TSS plots now use the prepared reference's `CHROM_SIZES`, rather than a coverage
vector ending at the last observed peak. ChIPseeker still defines the unique
gene TSS windows. A bundled overlap calculation counts each overlapping peak
once per base, reverses negative-strand windows and keeps full-width matrices
even when only one window remains. The display window is still -3000 to +3000
bp, independently of the configured promoter annotation window. The mean is
over **full-width windows with at least one peak overlap**; it is not a mean
over all genes, a sequencing-coverage profile or a normalized TSS enrichment
score. True chromosome-edge windows are excluded, not padded or truncated.
Identical TSS/strand windows are counted once, following `getPromoters`.

Every report exports `tss_windows.tsv` and `tss_window_summary.tsv`. Reasons are
`retained`, `no_peak_overlap` and `outside_reference`. The
`legacy_end_bound_exclusion` flag identifies valid retained windows that the old
last-peak boundary would have removed. Invalid sequence names, peak coordinates
or missing chromosome sizes produce an explicit plot failure.

`peak_globalIDR` is populated only for an explicitly selected ChromaFlow IDR
BED8 set. It retains the existing -log10 scale. Consensus/single-replicate
narrowPeak and custom inputs have `NA`; column 7 is never guessed to be IDR.
The same rule protects lookup, including older annotation tables. For custom
BED paths, IDR semantics are unknown even if the file happens to have 8 columns.
Use the native `--peakset idr[:threshold]` selector for IDR evidence.

For existing batch results, `batch reports` validates and corrects only the
known narrowPeak signalValue-as-IDR defect. It requires explicit non-IDR input
provenance, matching source peak IDs/coordinates/signal values, unchanged saved
outputs and the original annotation policy identity. It backs up the original
table, policy, stage record and (when present) candidate selection manifest under
`.chromaflow/batch/annotation_metadata_backup/`. Only the IDR column becomes NA;
all other table cells, gene IDs, coordinates, TSS assignments, motif evidence,
candidate rows and thresholds are preserved. Candidate `selection.json` updates
only the two affected annotation source identities so tracks validation still
works. The normal output integrity checks and strictly verified dependency
receipt keep expensive downstream work reusable. Keep these backups and batch
state with the project. Unrelated edits are not accepted as metadata repair.

## Install and refresh

Wait for existing ChromaFlow commands to finish, upload this archive and run:

```bash
cd /home/cfy/nfcore-upload
tar -xzf ChromaFlow-1.2.8-tss-idr-hotfix.tar.gz
KIT=/home/cfy/nfcore-upload/ChromaFlow-1.2.8-tss-idr-hotfix
python3 "$KIT/install.py" /home/cfy/biosoft/ChromaFlow-1.2.8

CF=/home/cfy/biosoft/ChromaFlow-1.2.8/bin/chromaflow
cd /home/cfy/Data/At_DAPseq/Dof
"$CF" batch reports -c chromaflow.config --only AT1G29160_col,AT3G45610_col &&
"$CF" batch reports -c chromaflow.config
```

The two small reports previously lost every window. Expect a written
`tss_avgprof.png`, a nonzero retained-window count and an all-field annotation
comparison pass after the metadata correction. The full refresh covers the
other completed annotations too. It does not run nf-core, MACS, MEME or candidate
integration. Repeating it verifies and reuses corrected reports. All new plots
remain under `results/analyses/<id>/reports/annotation/`; original annotation
images remain historical, so use the refreshed report directory for review.
No separate SHA256 file or manual transfer check is required.

## Resume test and feedback

After the full refresh returns zero, use the updated acceptance helper bundled
in this same archive. It snapshots the corrected baseline before a real resume:

```bash
KIT=/home/cfy/nfcore-upload/ChromaFlow-1.2.8-tss-idr-hotfix
python3 "$KIT/acceptance.py" run
TEST_RC=$?
printf '[TEST_SHELL_EXIT] %s\n' "$TEST_RC"
```

It checks the same full-selection native/downstream reuse as the prior test,
plus reference-bound TSS audits and non-IDR metadata. Expect 18 annotate,
14 motif and 14 candidate stages to resume, with unchanged peak/gene counts.
The known AT1G51700 model failure still permits BATCH_EXIT=1; the acceptance
test should return 0 when reuse and report-file checks pass. No MACS modeling
parameters or small-input policies are changed. This test retries the partial
native group using Nextflow resume; it can take time. Visual review remains
necessary even when the file/numerical checks pass.

Run feedback packaging separately, including after a failed test:

```bash
python3 "$KIT/acceptance.py" pack
```

Return the printed `At_Dof-tss-idr-*.tar.gz` archive. It includes the existing
bounded diagnostics plus all new TSS window audits, corrected annotation tables,
repair backups and chromosome-size metadata. FASTQ, BAM, whole reference FASTA
and containers are excluded. There is no GitHub publication in this update.

## Validation limits

Local tests compare the actual R overlap calculation against independent
per-base enumeration, cover both strands, factor sequence names, one/no windows,
true reference edges and the prior false boundary. Python tests verify the
field-only repair, strict refusal of mismatched inputs, rollback, downstream
reuse and tracks selection provenance. Full TxDb/ChIPseeker reconstruction and
the new complete figures must still be accepted on the server; local fixture
tests do not claim to replace that check.
