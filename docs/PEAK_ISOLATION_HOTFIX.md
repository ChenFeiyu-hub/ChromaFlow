# Peak-task failure isolation for ChromaFlow 1.2.8

Patch identifier: `peak-isolation-v1`. Prerequisites: the installed
`fimo-zero-hits-v1` and `read-length-v1` patches. The CLI version remains 1.2.8;
the installer and batch `[HOTFIX]` line identify this update.

## What changes

One failed MACS3 task previously stopped its shared nf-core group and marked all
members `nfcore_failed`. The batch now allows that exact MACS3 process to fail
without terminating independent samples. Every executed group writes a fresh
Nextflow task trace. ChromaFlow requires an unambiguous completed/cached or failed
MACS3 task for every requested biological replicate and checks expected output
files before allowing an analysis into downstream processing.

The MACS3 command, thresholds, effective genome size, read length, sample/control
assignments and annotation/motif policies are unchanged. There is no automatic
`--nomodel`, `--extsize` or species-specific rescue. Shared-reference or other
fatal nf-core errors retain the normal group failure boundary.

| Status | Meaning |
|---|---|
| `model_failed` | That sample's task stderr reports too few paired peaks to build the MACS3 model. |
| `peak_calling_failed` | A MACS3 task failed for another/unclassified reason. |
| `nfcore_incomplete` | Peak tasks succeeded but an expected replicate/consensus output is absent. |
| `nfcore_failed` | Native workflow failed, or its current task trace cannot be accepted; group outputs are not promoted. |
| `complete` | That analysis completed candidate integration. |
| `no_peaks` / `insufficient_peaks` | No peaks, or fewer than 10 peaks; candidates were not evaluated. |

If any biological replicate fails, the whole analysis is withheld from downstream
processing. No empty peak file is substituted. Failure/unevaluated candidate
counts remain NA. A partially successful group is stored as `partial`, never
`complete`; the overall ChromaFlow batch returns exit code 1 if a failure remains.

The native Nextflow run can return 0 after finishing the surviving tasks under
its scoped `ignore` strategy. ChromaFlow then reports `[PARTIAL]` and the individual
failures and returns 1. A native nonzero exit is not accepted as a partial success.
No unrelated process receives an error-strategy override. Existing optional-task
behavior defined by nf-core is retained.

Task traces are written to `logs/batch/<group>/<run-id>_tasks.tsv`. The generated
`.chromaflow/batch/runs/<group>/peak_tasks.json` records parsed task outcomes and
failure diagnostics; `outputs.json` contains only accepted analyses. The final
`[SUMMARY]` line counts actual analysis statuses.

## Install

After the active batch has ended, upload only
`ChromaFlow-1.2.8-peak-isolation-hotfix.tar.gz` to `/home/cfy/nfcore-upload`.
No separate SHA256 file or manual checksum step is needed.

```bash
cd /home/cfy/nfcore-upload &&
tar -xzf ChromaFlow-1.2.8-peak-isolation-hotfix.tar.gz &&
python3 ChromaFlow-1.2.8-peak-isolation-hotfix/install.py \
    /home/cfy/biosoft/ChromaFlow-1.2.8
```

The installer saves replaced files under
`.chromaflow-hotfixes/peak-isolation-v1/backup_<time>/`, verifies compatibility with
the installed prerequisite patches and records the update. This preserves the
existing internal resume identities; unexpected local runtime edits are not
silently assigned to a previously successful cache. No environments or containers
need to be rebuilt.

## Next server test: all 21 Dof analyses

Use the existing `chromaflow.config` and renamed local FASTQ sample table. Preserve
`reference/`, `work/`, `.chromaflow/` and all existing results. Do not delete caches
or split the sample table for this test.

```bash
CF=/home/cfy/biosoft/ChromaFlow-1.2.8/bin/chromaflow
cd /home/cfy/Data/At_DAPseq/Dof || exit 1

if "$CF" batch plan -c chromaflow.config --samples samples_dof.tsv --check-env; then
    "$CF" batch run -c chromaflow.config --samples samples_dof.tsv \
        --resume --until candidates
    batch_rc=$?
    printf '[BATCH_EXIT] %s\n' "$batch_rc"
    "$CF" batch status -c chromaflow.config
fi
```

This command intentionally checks status even when the batch exits 1. Expected
plan: 21 analyses in `se_k100` (18), `se_k130` (2) and `pe_k100` (1).

Acceptance criteria:

1. The startup line includes `peak-isolation-v1`.
2. Existing successful PE and 130 bp stages remain eligible for normal verified
   reuse. The failed SE 100 bp group is retried using its existing Nextflow cache.
3. If AtDof1.7 (`AT1G51700_col`) has the same model error, only that analysis is
   marked `model_failed`, with the model-pair count and task log in its reason.
   Other targets are evaluated independently; further failures remain visible.
4. AtDof1.5 (`AT1G29160_col`) can retain `insufficient_peaks` with two peaks and NA
   candidates. This is not a crash and is not evidence of zero candidate genes.
5. All 21 rows receive a status. Check the final `[SUMMARY]` and
   `results/summary/batch_status.tsv`; a remaining `model_failed` correctly keeps
   the batch exit code at 1.

For follow-up, return the final batch status table and any new failed task's
reported log. If the same group still reports `nfcore_failed`, provide its newest
batch nf-core log and task trace. Do not assume this is the same model error.

Unchanged partial groups retry the failed native tasks on the next `--resume`;
successful native tasks and completed downstream analyses remain eligible for
reuse. The internal input/output checks remain active. A successful-looking file
from an old run is not enough to admit a failed analysis.

## Validation scope

Tests cover task-trace validation, single and multiple failures, failed biological
replicates, missing outputs, unrelated workflow errors, partial-group recovery,
real candidate integration, and chaining the update onto previous hotfixes.
The real Nextflow 25.10.4 engine was also exercised with tiny local tasks: one
task failed, other tasks completed, cached successes were republished after an
output-directory archive, and unrelated task errors remained fatal. These are
software checks; the next server run validates the complete nf-core workflow on
the actual data.

Upstream semantics: [Nextflow errorStrategy](https://docs.seqera.io/nextflow/reference/process/directives/error-strategy),
[task trace fields](https://docs.seqera.io/nextflow/reference/config/trace).
