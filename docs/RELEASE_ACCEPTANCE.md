# ChromaFlow 1.2.9 release acceptance

This release consolidates all six validated 1.2.8 fixes in a complete installation.
Install into a new directory and retain the existing working project, config,
sample table, offline bundles, Conda environments, batch state and work cache.
There is no separate hotfix installation or manual transfer-checksum step.

The reference server acceptance was completed on 2026-11-08. See the
[final validation record](RELEASE_VALIDATION.md) for evidence and scope.
The commands below remain available to validate another installation or project.

## Scope and gates

| Gate | Required evidence |
|---|---|
| Clean package | CLI version, complete runtime manifest, LF scripts, archive layout/modes, valid docs |
| Python/Linux regression | All packaged tests pass without skipped cases; real bedtools for IDR interval postprocessing |
| R regression | Window/strand/boundary tests, rendered plots, real R-to-Python annotation guard |
| Upgrade compatibility | Saved and proposed scientific plan fields agree before any analysis is started |
| Real batch resume | Completed annotation/motif/candidate stages and successful native tasks reused; outcomes preserved |
| Real reports | Full ChIPseeker reconstruction passes; all expected report files and numerical audits valid |
| Repeated reports | Second report command reuses identical records and files |
| Final inspection | Review returned logs, numerical evidence and PNGs; file presence alone does not close this gate |

Tests of expensive producers in the synthetic suite use fixtures. They are named
as such and do not claim to run real nf-core, MEME discovery or IDR model fitting.
Real server evidence from the existing dataset supplies the integration check.
The DAP-seq batch scope still ends at candidates; it does not add batch IDR,
lookup, tracks, ATAC or scheduler support. Single-project commands retain their
existing regression coverage. The validated GFF normalization engine still has
its own 1.2.3 engine version; this is independent of the CLI release version.

## Install and run

Wait for any active project commands to finish. Extract the full archive into a
new directory; do not overlay either an old installation or project outputs:

```bash
cd "$HOME/nfcore-upload"
CF_HOME="$HOME/biosoft/ChromaFlow-1.2.9"
test ! -e "$CF_HOME" && mkdir -p "$CF_HOME" &&
tar -xzf ChromaFlow-1.2.9.tar.gz -C "$CF_HOME" --strip-components=1
"$CF_HOME/bin/chromaflow" version
```

If the directory already exists, stop and inspect it rather than extracting
over it. Reusing an already correctly extracted installation is fine.

Run the single packaged acceptance helper from any working directory:

```bash
CF_HOME="$HOME/biosoft/ChromaFlow-1.2.9"
PROJECT=/path/to/existing/DAPseq-project
python3 "$CF_HOME/tools/acceptance.py" run --project "$PROJECT"
TEST_RC=$?
printf '[TEST_SHELL_EXIT] %s\n' "$TEST_RC"
```

The helper uses the **full saved sample selection**, config and configured
Conda environments. It checks the proposed plan, runs the full Python suite in
the reference environment and three real R tests in the annotation environment,
then performs one `batch run --resume --until candidates` and two `batch reports`
commands. It records all output in English under
`<project>/diagnostics/release-acceptance-<timestamp>-<pid>/`.
It does not change thresholds, model settings, controls or environment packages.
Report refresh uses the separate `results/analyses/<id>/reports/annotation/`
directories. It may archive the previous report directory once when rendering
the updated layout. Completed scientific outputs remain subject to strict checks.

The helper checks the complete release runtime against its included compatibility
metadata. This is an internal cache check, not a manual upload/checksum workflow.
An unexpected installation, changed plan, missing dependency, failed test or
skipped Python test stops acceptance before starting the resumed analysis.
It intentionally does not narrow the full selection: changing native group
membership would invalidate the resume test. The known failed native task can
therefore be retried while successful tasks remain cached.

## Expected At_Dof outcome

The accepted 1.2.8 baseline contains 21 analyses: 14 complete, 4 with fewer than
10 peaks, 2 with no peaks and 1 with MACS3 model failure (40 model pairs, 100
required). Expect 18 annotation, 14 motif and 14 candidates stages to resume.
Expect 72 report PNGs for 18 annotations. The other three analyses have no
applicable annotation report. The updated layout should separate long y-axis
titles from three-decimal tick labels.

`BATCH_EXIT=1` remains expected while the known model failure persists.
`ACCEPTANCE_EXIT=0` means automated regressions, upgrade reuse, report generation,
idempotence and numerical/file audits passed. It does not mean every analysis
succeeded or that low/no-peak outcomes establish biological absence.
`ACCEPTANCE_EXIT=1` indicates evidence needs review; `2` means preflight/setup or
execution did not finish; `130` means interruption. Do not delete work/state to
make a test pass. Correct the cause and preserve the failed test directory.

Cached candidate selections may correctly say version 1.2.8 after an upgrade;
new 1.2.9 reports identify their own renderer version. Original scientific
provenance is retained. TSS profiles show mean peak overlap per retained window,
not sequencing coverage or a standard TSS enrichment score.

## Collect feedback even after failure

Run this separately from the test, so a nonzero test exit does not prevent it:

```bash
python3 "$CF_HOME/tools/acceptance.py" pack --project "$PROJECT" \
  --outdir "$HOME/nfcore-upload"
```

Return the printed `ChromaFlow-1.2.9-acceptance-*.tar.gz`. The archive contains
regression summaries/logs, before/after state, plan comparison, report repetition
evidence, all annotation PNGs/TSVs and bounded native failure diagnostics. It
excludes FASTQ, BAM, whole reference FASTA, containers and unrelated work files.
No additional checksum file is generated. Packing success is separate from
test success. For an earlier run, pass its exact `--test-dir`.

## Local regression only

```bash
cd "$CF_HOME"
python3 tests/run_suite.py --json /tmp/chromaflow-regression.json
Rscript tests/test_tss_windows.R "$CF_HOME"
Rscript tests/test_annotation_plots.R "$CF_HOME" /tmp/chromaflow-r-plots
Rscript tests/test_report_annotation_guard.R "$CF_HOME" python3 /tmp/chromaflow-r-guard
```

If bedtools is not on PATH, set `CHROMAFLOW_TEST_BEDTOOLS_IMAGE` to an existing
local SIF/IMG and make Apptainer available. No image download occurs in tests.
The local Python suite plus R fixtures do not replace the server acceptance
command or final visual inspection. `release/validation.json` records the local
and completed reference-server gates. The earlier test archives retain their
build-time pending status; the final validation record closes those gates for
the identical runtime. A different installation or project requires its own
checks; the reference acceptance is not a claim that a future run has passed.
