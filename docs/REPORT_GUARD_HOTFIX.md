# Report-table comparison correction for ChromaFlow 1.2.8

This patch follows `resume-reports-v1`. It fixes the report refresh rejection
caused by comparing an intermediate R table with the finished annotation table.

Normal `annotate` runs R and then `annotation_policy.py validate`. The Python
step appends `feature_class` and serializes TSV cells with CSV quoting where
needed. The initial report refresh compared raw R text before that step, so
unchanged annotations were rejected. This is a report-guard implementation bug;
the rejection alone does not mean that annotation, genome data or peak calls
changed. TxDb construction and annotation can succeed before this check fails.

The corrected guard writes a temporary R table, applies the same existing
feature policy in memory, and compares every parsed field and row against the
saved table. It still checks IDs, coordinates, annotations, distances, scores,
row order and `feature_class`. It introduces no numeric tolerance or ID rewriting.
Neither comparison input is edited. Real differences now print the first
affected rows/columns and both values in the log.

## Install

After current ChromaFlow commands finish, upload the single archive to the
server and run:

```bash
cd /home/cfy/nfcore-upload
tar -xzf ChromaFlow-1.2.8-report-guard-hotfix.tar.gz
python3 ChromaFlow-1.2.8-report-guard-hotfix/install.py \
  /home/cfy/biosoft/ChromaFlow-1.2.8
```

The installer requires the four existing patch manifests (FIMO zero hits,
read length, peak isolation, resume/reports), backs up changed files and records
`report-table-guard-v1`. Reinstalling is idempotent. No separate checksum file,
container download or environment installation is required. Cache compatibility
is preserved across the verified patch chain; normal input/output checks remain.

## Retest one small analysis, then all reports

```bash
CF=/home/cfy/biosoft/ChromaFlow-1.2.8/bin/chromaflow
cd /home/cfy/Data/At_DAPseq/Dof

"$CF" batch reports -c chromaflow.config --only AT5G60850_col &&
"$CF" batch reports -c chromaflow.config
```

The uploaded example is AT5G60850_col with four annotated peaks. The pasted
AT1G29160_col error with two peaks comes from the same comparison step. A
successful comparison will print:

```text
[VERIFIED] Reconstructed annotation matches all ... rows and ... columns after the same feature_class policy; saved files unchanged
```

Reports remain under `results/analyses/<id>/reports/annotation/`. Empty valid
TSS profiles may be `skipped`; `failed` or a field-level mismatch still requires
review. `--only` does not overwrite the full batch status summary. These commands
do not rerun alignment, peak calling, MEME-ChIP or candidate integration.

The previous failed guard ran before report publication or selection-version
correction. Existing annotation/candidate results therefore remain in place;
there is no need to clear results, work directories or batch state.

## Validation scope

Tests reproduce the old failure using actual R serialization followed by the
real Python annotation policy. They then execute the production R comparison
function with the corrected Python comparator. Unchanged tables pass; altered
gene IDs, coordinates, scores, feature classes, ordering and promoter policy are
rejected. Saved files remain unchanged. Full Bioconductor reconstruction is
still checked by the server retest; the local test does not claim to run it.

After reports pass, continue the resume and isolated MFOLD checks described in
[the previous guide](RESUME_REPORTS_HOTFIX.md). If another comparison fails,
return its new log: it will now contain the specific field differences.
