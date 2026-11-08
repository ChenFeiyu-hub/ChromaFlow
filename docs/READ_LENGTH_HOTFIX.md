# Read-length compatibility hotfix for ChromaFlow 1.2.8

Patch identifier: `read-length-v1` (2026-11-06). Apply this after the
`fimo-zero-hits-v1` patch. The CLI version remains 1.2.8; installation manifests
and the batch `[HOTFIX]` message identify both patches.

## Behavior

nf-core/chipseq 2.1.0 accepts only 50, 75, 100, 150 or 200 for `read_length`.
That parameter is used to estimate effective genome size when `macs_gsize` is
absent. Previously ChromaFlow supplied both, so a batch with a real read length
of 130 failed validation even after ChromaFlow had measured its effective size.

Both container and Conda launchers now supply `--macs_gsize` alone when that value
is available. Otherwise they retain the upstream `--read_length` estimation path.
Batch still measures the real read length, uses that k for `unique-kmers.py`, and
records it in the group config and size metadata. No rounding to 150, FASTQ
editing or schema-validation bypass is involved. Single-project users with a
nonstandard read length should first measure size with, for example,
`chromaflow reference -c chromaflow.config --kmer 130` in their reference
environment, then inspect `chromaflow run --dry-run`.

The pinned ATAC-seq launcher uses the same alternative size inputs; its behavior
is covered by the launcher regression tests. Batch remains DAP-seq only.

Upstream contracts:

- [chipseq 2.1.0 schema](https://github.com/nf-core/chipseq/blob/2.1.0/nextflow_schema.json)
- [chipseq 2.1.0 size estimation](https://github.com/nf-core/chipseq/blob/2.1.0/workflows/chipseq.nf)
- [atacseq 2.1.2 size estimation](https://github.com/nf-core/atacseq/blob/2.1.2/subworkflows/local/prepare_genome.nf)

## Install into the existing installation

Wait for any active ChromaFlow job to finish, then upload the archive and checksum
to `/home/cfy/nfcore-upload`. Run as the installation owner:

```bash
cd /home/cfy/nfcore-upload &&
sha256sum -c ChromaFlow-1.2.8-read-length-hotfix.tar.gz.sha256 &&
tar -xzf ChromaFlow-1.2.8-read-length-hotfix.tar.gz &&
(cd ChromaFlow-1.2.8-read-length-hotfix && sha256sum -c MANIFEST.sha256) &&
python3 ChromaFlow-1.2.8-read-length-hotfix/install.py \
    /home/cfy/biosoft/ChromaFlow-1.2.8
```

The installer requires the recorded FIMO hotfix, verifies the complete expected
runtime tree and each payload hash, and saves backups under
`.chromaflow-hotfixes/read-length-v1/`. Repeating the installer verifies the
already-installed patch. Unexpected local edits or a broken predecessor manifest
stop installation before files are replaced. No environment rebuild or container
download is required.

This narrow update preserves the predecessor's stage cache identities because
batch always supplies `macs_gsize`; the removed parameter did not affect any
previously successful computation. All input, environment, parameter and output
hash checks remain active. Failed stages are never reused. Later unrelated code
edits disable this compatibility mapping.

## Continue the existing Arabidopsis Dof project

Keep the current `chromaflow.config` and `samples_dof.tsv`. They already contain
the copied/renamed local FASTQs and the sequence-name-corrected FASTA. Do not
replace them with an older deployment template. Keep `work/`, `.chromaflow/`,
`reference/` and the three completed smoke results for verified resume.

```bash
CF=/home/cfy/biosoft/ChromaFlow-1.2.8/bin/chromaflow
cd /home/cfy/Data/At_DAPseq/Dof || exit 1

"$CF" batch plan -c chromaflow.config --samples samples_dof.tsv --check-env &&
"$CF" batch run -c chromaflow.config --samples samples_dof.tsv \
    --only AT1G29160_col,AT3G50410_col --resume --until nfcore &&
"$CF" batch run -c chromaflow.config --samples samples_dof.tsv \
    --resume --until candidates

"$CF" batch status -c chromaflow.config
```

Expected full plan: **21 analyses, 3 groups**:

| Group | Targets | Effective-size k |
|---|---:|---:|
| `se_k100` | 18 | 100 |
| `se_k130` | 2 | 130 |
| `pe_k100` | 1 | 100 |

The first run tests the two 130 bp analyses, `AT1G29160_col` (AtDof1.5/COG1) and
`AT3G50410_col` (OBP1), through nf-core. Its `stopped_after_nfcore` status and NA
candidate counts are expected. On success, the `&&` chain continues all 21
analyses through candidates. On failure it stops before the full run; inspect the
reported log rather than changing MACS thresholds or forcing a successful status.

The expanded `se_k100` sample set changes its group signature. Nextflow is called
with `-resume` and may reuse matching tasks; downstream stages of its two smoke
analyses are recomputed. The unchanged `pe_k100` group and the successfully tested
`se_k130` group remain eligible for reuse. Reference preparation and the k=100
estimate remain eligible too. Eligibility always depends on the normal checks.

The last `batch status` reports the most recent run's selection. If the 130 bp
test failed, it contains only those two targets, not a full-project outcome.
After a full run, inspect `results/summary/batch_status.tsv` for all 21 analyses.
`complete` means candidates finished; `no_peaks` and `insufficient_peaks` describe
limited inputs, with unassessed candidate counts recorded as NA. A failed nf-core
group affects its members; a downstream failure affects its analysis, while other
independent work continues. The batch exits nonzero if a failure was recorded.

The patch addresses parameter validation. Actual alignment, MACS model fitting
and the remaining samples' downstream results still need the server run.

## Local validation (2026-11-06)

- Full regression suite: 146 tests, 142 passed and 4 skipped because bedtools was
  unavailable in the local test environment.
- Nextflow 25.10.4 and the cached `nf-validation@1.1.3` plugin reproduced the
  rejection of 130 for both pinned pipeline schemas. Both accepted the explicit
  size without `read_length`. These tests defined no analysis processes.
- Actual release-to-FIMO-to-read-length installations reused all eligible fixture
  stages, including candidates, and preserved existing analysis-output hashes.
  Repeated installation, corrupt payload/predecessor refusal, unrelated-code
  refusal and rollback after an injected write failure were checked.
- Miniature FASTQs and a controlled estimator verified k=130 propagation through
  batch grouping, the estimator call and saved metadata. Scientific producers
  were test fixtures; this does not replace the real server run.
