# ChromaFlow 1.2.9 completed release validation

The reference server acceptance and final evidence review passed on 2026-11-08.
No release-blocking defect was found within the agreed scope. The tested runtime
is recorded in [`release/compatibility.json`](../release/compatibility.json);
machine-readable validation is in [`release/validation.json`](../release/validation.json).

| Gate | Result |
|---|---|
| Python regression on the target server | 193 tests; no failures, errors or skips |
| R regression in the configured annotation environment | TSS windows, annotation plots and real R/Python annotation guard passed |
| Upgrade from the validated 1.2.8 fixes | Scientific plan fields and stage identities preserved |
| Completed downstream reuse | All 46 stages reused: 18 annotation, 14 motif, 14 candidates |
| Native resume in the partial single-end group | 583 cached tasks; samplesheet validation and MultiQC completed; the one known model failure retried |
| Actual ChIPseeker report reconstruction | All 18 analyses passed, covering 27,240 annotation rows |
| Repeated report generation | All 18 report sets reused with unchanged records and file metadata |
| Report figures | All 72 PNGs decode; all 18 updated TSS profiles visually inspected; the other 54 PNGs unchanged |
| Independent numerical review | All 6,001 points in each TSS profile reconstructed from intervals and reference bounds; maximum rounding difference 1.39e-16 |

The plan check, report commands and overall acceptance exited 0. The raw batch
command exited 1 because the existing MACS3 failure remained explicitly reported.
Final analysis outcomes were unchanged: 14 complete, four with fewer than ten
peaks, two with no peaks, and one model failure. That sample still had 40 paired
model peaks against a requirement of 100. No calling thresholds were changed and
no failed or low-peak result was promoted to success.

Scientific preservation was checked against the preceding accepted server round:
72 included annotation files, 18 selected-peak files, 14 candidate selection
manifests and 126 report numerical/status tables were identical. Cached candidate
provenance correctly retains version 1.2.8; new report metadata identifies 1.2.9.
The TSS axis layout now separates three-decimal ticks from the y-axis title.

The local release checks also passed Python/Linux regression, actual R fixtures,
four IDR interval postprocessing tests using real bedtools, offline workflow and
plugin startup for chipseq 2.1.0 and atacseq 2.1.2, source/document checks and fresh
package installation checks. Runtime code is unchanged by the final publication
cleanup, which removes the unused root `.gitignore` and `conf/` and updates docs
and validation metadata. Existing test archives retain their original bytes.

## Interpretation and scope

- Real-data integration covers DAP-seq batches through candidates, upgrade/resume,
  native failure isolation and report reconstruction. Successful scientific
  producers were reused in the final acceptance; this was not a fresh raw-data
  recomputation of every analysis.
- The bounded feedback contains stage records, captured output metadata,
  selected scientific files and reports. FASTQ, BAM, whole reference FASTA and
  complete candidate table bodies are excluded.
- IDR statistical model fitting, batch ATAC and scheduler execution are not
  additional claims of this acceptance. Existing single-project commands retain
  their regression coverage.
- Sparse/no-peak outputs do not establish biological absence. TSS profiles show
  mean peak overlap per retained window, not read coverage or a standard TSS
  enrichment score.
- Fontconfig initialization, undefined `monochromeLogs`, and explicitly ignored
  OpenMPI `openib` warnings remained in the evidence without new task failures.

Use the [acceptance guide](RELEASE_ACCEPTANCE.md) to validate another installation
or project. A successful reference acceptance does not override local input,
environment or output checks.
