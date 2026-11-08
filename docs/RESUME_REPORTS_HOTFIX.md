# ChromaFlow 1.2.8 resume and reports update (2026-11-07)

Before testing report refresh, apply the subsequent
[report-table guard correction](REPORT_GUARD_HOTFIX.md). The initial guard
omitted the `feature_class` postprocessing step when comparing annotations.

This update follows `fimo-zero-hits-v1`, `read-length-v1` and `peak-isolation-v1`.
Install it only after the current batch has finished. It requires no downloads,
new Conda environments or separate checksum file. The installer backs up the
changed files, verifies the existing installation internally and refuses
unrecognized local modifications. The CLI version remains 1.2.8; installation
provenance records `resume-reports-v1`.

## What changes

- Keep the known ignored Open MPI `mca_btl_openib` / `librdmacm.so.1` loader
  warning in diagnostic logs, but exclude that exact line from the motif tool
  identity. Its process ID no longer forces MEME-ChIP and candidates to rerun.
  Tool paths, versions, unknown diagnostics, parameters, inputs and output
  checks remain significant. A nonzero version-probe exit is still an error.
- Migrate old motif cache records only when a saved environment record proves
  that the ignored line is the only difference and all outputs still match.
  Preserve the original dependency identity so candidates can also resume.
  Print a reason before each stage that actually runs.
- Render annotation intersections from ChIPseeker's detailed memberships using
  bundled base R code; `ggupset` is no longer required for this plot.
  A peak can have several detailed feature memberships; the exact intersection
  bars partition the peaks. These memberships do not change the exclusive
  annotation category or promoter/nearest-TSS assignment.
- Render the TSS profile from the unweighted tag matrix as the mean peak overlap
  over retained windows. Its display window stays -3000 to +3000 bp, independently
  of the configured promoter window. Export the numerical values. Empty/all-zero
  matrices are explicitly skipped; invalid matrices are reported as failures.
  `plot_status.tsv` distinguishes written, skipped and failed plots.
- Read the actual installed CLI version for new `selection.json` files.
  `batch reports` corrects the known 1.2.7 label in saved 1.2.8 candidate results
  only after verifying every recorded candidate output. It preserves all
  filters and scientific fields, backs up the original metadata, and updates
  the recorded output hash so resume remains valid.

## Install the single uploaded archive

```bash
cd /home/cfy/nfcore-upload
tar -xzf ChromaFlow-1.2.8-resume-reports-hotfix.tar.gz
python3 ChromaFlow-1.2.8-resume-reports-hotfix/install.py \
  /home/cfy/biosoft/ChromaFlow-1.2.8
```

Do not reinstall the original 1.2.8 release over the patched directory. The new
installer requires the three predecessor installation manifests and is
idempotent. Backups live under the installation's `.chromaflow-hotfixes/`.

## Test 1: refresh reports without rerunning scientific stages

```bash
CF=/home/cfy/biosoft/ChromaFlow-1.2.8/bin/chromaflow
cd /home/cfy/Data/At_DAPseq/Dof
"$CF" batch reports -c chromaflow.config
```

This command reads the saved batch plan and requires the unchanged config. It
does not read FASTQs, run nf-core/MEME-ChIP, or change the batch status summary.
It rebuilds the annotation in memory with the same parameters and requires
exact equality with the existing annotation table before publishing plots.
Reports are separate from cached scientific outputs:

```text
results/analyses/<analysis_id>/reports/annotation/
  anno_pie.png
  dist_to_tss.png + dist_to_tss.tsv
  anno_upset.png + annotation_intersections.tsv
  tss_avgprof.png + tss_avgprof.tsv   # only with valid nonzero data
  genomic_distribution.tsv
  plot_status.tsv
  report.json
```

The original annotation directory is retained for provenance. Inspect the new
`reports/annotation` figures after a refresh. Repeating `batch reports` reuses
verified refreshed reports. `--only AT3G45610_col` selects an exact saved analysis
without narrowing or overwriting `batch_status.tsv`. `MAKE_PLOTS=false` is
respected. A failed render leaves prior reports and scientific results intact.
Logs are under `logs/batch/*_reports.log` and `logs/analyses/<id>/*_reports.log`.

For the supplied run, expect 18 annotation reports and skips for the one failed
model plus two no-peak analyses. A TSS curve may legitimately be `skipped` for a
small peak set. A `failed` plot or annotation-equality failure needs review.
Local validation covers the base R renderer and controlled workflow producers;
the complete ChIPseeker reconstruction still needs this server acceptance run.

## Test 2: full-selection resume

Run this after Test 1 succeeds, using the same complete sample table:

```bash
"$CF" batch run -c chromaflow.config --samples samples_dof.tsv \
  --resume --until candidates
BATCH_RC=$?
printf '[BATCH_EXIT] %s\n' "$BATCH_RC"
"$CF" batch status -c chromaflow.config
```

The partially completed native group is retried with Nextflow resume. Successful
native tasks remain eligible for its cache. Expect the 14 complete analyses'
motif/candidate stages to show `[RESUME]`, including any one-time MPI-diagnostic
migration messages. The 18 successful annotations should also resume. A change
in scientific inputs, environment, group membership or recorded files still
causes a justified rerun; see `[RUN]` reasons.

With the original data and MACS settings, AtDof1.7 is expected to remain
`model_failed`, and the batch may still exit 1. Expected counts remain:
14 complete, 4 insufficient_peaks, 2 no_peaks, 1 model_failed. This update does
not convert low-evidence samples into successful analyses.

## Test 3: isolated AtDof1.7 MFOLD comparison

The archive includes a project-specific test helper outside the software
runtime. It reads the current failed task and checks its original command. It
uses that task's cleaned treatment/control BAMs and the installed MACS3 3.0.1
image. It never executes the saved shell script. Review the three commands:

```bash
KIT=/home/cfy/nfcore-upload/ChromaFlow-1.2.8-resume-reports-hotfix
python3 "$KIT/atdof_model_check.py" --plan-only
python3 "$KIT/atdof_model_check.py"
MODEL_RC=$?
printf '[MODEL_EXIT] %s\n' "$MODEL_RC"
```

Comparisons: MFOLD 5,50 (baseline), 2,50, and 2,100. Everything else stays
equivalent to the reviewed command: BAM mode, keep-dup all, the recorded
effective genome size, the same control, q=0.05 and bandwidth=300. The broader
ranges test model eligibility; they do not establish an optimal parameter.
There is no `--nomodel`, fixed fragment-length assumption or automatic adoption.

Each invocation creates a new `experiments/AT1G51700_mfold_<timestamp>/` containing
`experiment.json`, `mfold_summary.tsv`, logs and any MACS outputs. Input BAMs are
mounted read-only. All three comparisons are attempted even if one model fails.
Exit 1 means at least one comparison failed; a baseline model failure is expected.
Batch results/status are left intact. Successful modeling still requires review
of fragment estimates, warnings and peaks, especially with the previously
observed high duplicate rate and shallow control. Do not select parameters only
to increase the candidate count.

## Feedback to return

Return the new `*_reports.log`, latest `*_run.log`, `batch_status.tsv`,
`.chromaflow/batch/reports.json`, each refreshed `plot_status.tsv`, and the new
experiment's `experiment.json`, `mfold_summary.tsv` and `macs3.log` files.
No FASTQs, BAMs, reference FASTA or container images are needed.

Implementation references:
[ChIPseeker detailed memberships](https://github.com/YuLab-SMU/ChIPseeker/blob/master/R/upsetplot.R)
and [MACS3 3.0.1 modeling and peak calling](https://github.com/macs3-project/MACS/blob/v3.0.1/MACS3/Commands/callpeak_cmd.py).
