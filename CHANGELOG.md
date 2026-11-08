# Changelog

## v1.2.9 — consolidated batch release and acceptance (2026-11-08)

- Complete target-server acceptance: 193 Python tests without skips, three R
  suites, all 46 completed downstream stages reused, and 18 report sets checked
  numerically and visually. Repeated report generation preserves all outputs.
- Remove the obsolete root `.gitignore` and `conf/nextflow.config`; current
  runtime configuration remains in project settings and the nf-core launcher.
- Include all six validated 1.2.8 fixes in one complete installation. Preserve
  reviewed scientific stage identities across a side-by-side upgrade, subject
  to exact runtime matching and the existing input/environment/output checks.
- Keep cached 1.2.8 candidate provenance unchanged; newly rendered reports and
  newly computed results identify the new release. Retain the narrow historical
  reporter-typo repair for saved 1.2.8 plans.
- Separate three-decimal TSS tick labels from the axis title; profile values,
  reference bounds, windows, gene/peak assignments and thresholds are unchanged.
- Package one acceptance helper for Python/Linux and actual R regressions,
  full-plan upgrade/resume checks, repeated report refresh and bounded feedback.
  Report skipped tests and pending visual/server verification explicitly.
- Consolidate English/Chinese installation and testing documentation. Deliver
  full archives without separate local-transfer checksum files or extra hotfix
  kits. Offline download integrity and normal cache validation are retained.
- Keep the DAP batch endpoint at candidates and preserve the known biological/
  data sufficiency statuses. Do not change MACS modeling or promote failed cases.

## v1.2.8 TSS and IDR metadata hotfix (2026-11-07)

- Use reference chromosome lengths for full-width TSS profiles; export retained,
  zero-overlap and true-edge window audits with strand-correct peak counts.
- Interpret globalIDR only for an explicit ChromaFlow IDR BED8 selection; use
  NA for non-IDR/custom annotations and lookup, including legacy lookup input.
- Repair verified old batch IDR labels with backups, strict output checks and
  preserved dependency identities. Refresh candidate source provenance without
  changing candidate evidence. Keep all biological filters and MACS settings.
- Bundle updated server resume/report acceptance and feedback collection in one
  local archive. No new environment, manual checksum or GitHub update is needed.

## v1.2.8 report-table guard correction (2026-11-07)

- Compare reconstructed annotations after the same Python feature_class policy
  used by normal annotate; compare parsed TSV cells instead of intermediate R
  text against the postprocessed table. Keep all scientific fields and row
  order strict and report concrete differences without editing either input.
- Add actual R-to-Python guard integration tests, including the former false
  rejection and rejection of changed scientific evidence. Retain verified cache
  compatibility through the existing four-patch chain.

## v1.2.8 resume-reports hotfix (2026-11-07)

- Exclude only the known ignored Open MPI loader warning from motif environment
  identity. Verify historical signatures and outputs before migrating cached
  results; keep downstream dependency identities stable and log rerun reasons.
- Replace the optional ggupset dependency with bundled base R intersection
  plots, and validate/export TSS profiles with explicit skipped/failed status.
- Add `batch reports` for separately staged diagnostic figures. Require exact
  reconstructed annotation equality and preserve scientific results/status.
- Derive selection metadata from the installed CLI version. Repair the known
  old label in verified batch results with a backup and updated cache metadata.
- Provide a single local update archive and a separate, non-promoting AtDof1.7
  MFOLD sensitivity helper. Do not change MACS defaults or biological filters.

## v1.2.8 peak-isolation hotfix (2026-11-06)

- Isolate the pinned batch MACS3 process with a scoped Nextflow error strategy;
  validate a fresh task trace and required replicate/consensus files before
  admitting each analysis to downstream processing. Keep other native failures fatal.
- Report `model_failed`, `peak_calling_failed` and `nfcore_incomplete` separately.
  Store partial groups as partial and retain a nonzero batch exit while any analysis
  fails. Do not create empty peaks or change MACS thresholds/model parameters.
- Preserve completed-stage cache eligibility across the checked hotfix chain.
  Retry failed native tasks while retaining successful downstream reuse.
- Add a final status-count summary and describe the fewer-than-ten-peak cutoff
  accurately as ChromaFlow's small-input policy. Close subprocess output streams.
- Deliver a single local installation archive and normal batch continuation
  commands without separate checksum files or manual SHA256 steps.

## v1.2.8 read-length hotfix (2026-11-06)

- Send only `--macs_gsize` when effective genome size is already available, in
  both the container and Conda launchers. Keep upstream `--read_length` fallback
  when size is absent. This removes the unnecessary enum rejection of 130 bp
  batches without changing the actual read length or the k-mer estimate.
- Chain verified installation provenance with the FIMO zero-hit patch so the
  same selected inputs remain eligible for resume, including candidates. Changes
  to sample-group membership still invalidate that group's dependent stages.
- Cover both launchers, ChIP/DAP and ATAC arguments, real-k batch metadata,
  chained upgrade/resume and malformed/tampered patch records with regressions.
- Add installation and full 21-target Dof continuation commands to the guide.

## v1.2.8 FIMO zero-hit hotfix (2026-11-06)

- Accept FIMO 5.5.7's complete three-comment, headerless zero-hit TSV in the
  shared candidates/lookup reader and direct motif tracks export. Retain motif
  identity and a zero hit count; leave input results and thresholds unchanged.
- Continue rejecting blank, incomplete or malformed tables and footer-only TSVs
  whose companion GFF contains hits. Record a checked GFF in selection provenance.
- A checked hotfix installer records before/after runtime hashes. Only this
  reader-only update preserves reference/nf-core/annotation/motif cache identities;
  candidate integration is recomputed, and all normal output checks remain active.
- Add zero-hit integration and upgrade/resume regression tests. This targeted
  hotfix does not address the separate 130 bp nf-core read-length compatibility issue.

## v1.2.8 — manifest-driven DAP-seq batches

- Add batch init, plan, run, status, exact analysis selection and staged execution
  through candidates. Keep one user config and TSV with explicit controls,
  technical runs, biological replicates and SE/PE layout.
- Share reference normalization, per-k effective size and compatible nf-core
  executions. Reuse saved BWA indexes; bind downstream inputs by exact output
  manifests. Handle singleton narrowPeak output without claiming reproducibility.
- Isolate annotation, motif discovery and candidate results per analysis; retain
  the existing promoter/ID rules. Scope motif names in cross-analysis summaries.
- Switch environments explicitly, record executable paths, constrain local
  executor resources and MEME threads, and save English terminal logs.
- Resume only matching stages with verified outputs. Continue independent units
  after failure; keep NA distinct from a successful zero-candidate result.
- Add synthetic workflow/candidate integration and recovery regression coverage.
  Package a separate curated Dof test manifest; no dataset IDs enter core code.


## v1.2.7 — shared motif selection and auditable tracks

Documentation review (2026-09-28): expand both READMEs with repository installation,
consensus-only operation, annotation policy, exact-ID lookup, shared filter options
and selection-to-tracks commands. Align current deployment/architecture guides and
replace historical personal deployment paths with generic public examples. Runtime code
is unchanged from the tested v1.2.7 package.

- Share exact-model, FIMO p/q, discovery E-value and optional CentriMo filters
  between candidates and lookup; filter every peak before ranking. Reject absent
  requested statistics/unscanned models; record model identities and PWM hashes.
- Export selected genes/evidence and a hashed selection manifest. Add tracks
  --selection and candidates --plot; lookup --level controls its plots too.
  Preserve all context peaks and neighboring isoforms. Reject changed source
  files, cross-branch selections and changed promoter parameters.
- Separate local/central CentriMo modes, E-value threshold status and unavailable
  results. Defaults retain candidates without imposing extra motif filters.
- Add opt-in DEFAULT_PEAKSET project setting; CLI branch selection wins.
- Scale gene/motif track heights to regional content and preserve three significant
  digits in peak labels; identify ratio inputs/pseudocount in logs.
- Correct SPP quality tags (-2/-1/0/1/2 = very low/low/medium/high/very high).
- Keep exact gene IDs, annotation rules, GFF/TxDb code and offline bundles unchanged.

## v1.2.6 — explicit annotation policy, promoter audits and exact gene IDs

- Make ChIPseeker's gene-level nearest-TSS assignment and feature priority explicit.
  Verify that Promoter labels agree with the configured inclusive TSS-distance
  window; record parameters and the annotation hash, and reject stale parameter
  combinations in candidates, lookup and tracks. Honor MAKE_PLOTS for R figures.
- Share feature parsing and location rules between candidates and lookup. Recognize
  exon/intron rank 1 in ChIPseeker text, including negative-strand transcripts,
  without species-specific naming rules. Preserve raw annotations and append
  feature_class. The project template accepts Promoter, 5' UTR and 1st Exon;
  existing explicit TARGET_FEATURES values remain authoritative.
- Embed promoter-region audits in annotate (all genes), lookup and tracks (query
  genes). Use strand-aware gene-level windows from the same PROMOTER_UP/DOWN,
  include TSS and endpoints, clip to reference boundaries, and report every
  overlapping selected peak, including other-assigned and unannotated peaks.
  These supplementary associations never replace ChIPseeker's assigned gene.
- Remove lookup's automatic ID suffix normalization. Show exact-ID examples and
  match counts, require y/N confirmation, and support explicit --yes for unattended
  runs. Zero matches fail; partial unmatched IDs are separately reported. Use exact
  gene IDs in tracks too. lookup --plot now draws the complete query list.
- Show all selected peaks and neighboring genes in the track region, mark promoter
  windows, and expand short flanks to contain the promoter. Improve label margins;
  verify both strands with real pyGenomeTracks using synthetic data.
- Retire unused automatic top-N options from the template. Use standard ChromaFlow
  commands for the entire retest; no project-specific execution kit is needed.
- Preserve GFF/GTF normalization, offline bundles and successful nf-core outputs.

## v1.2.5 — IDR interval ordering and gene-list motif tracks

- Re-sort bedtools annotate output before recovering MACS signal with map.
  A real-bedtools regression spanning interval bins and sequence names reproduces
  the 1.2.4 failure and verifies the fix. Preserve IDR fitting and thresholds.
- Reject GO preparation from a marked failed/incomplete IDR run, including when
  older annotation files remain. Rerun annotation after rebuilding IDR.
- Respect the R annotation consistency flag instead of counting missing `NA`
  feature-gene values as conflicting gene IDs. Correct the dropped-peak log's
  units. Candidate location/motif selection is unchanged.
- Add opt-in `tracks --motifs --flank BP`: existing FIMO hits, input windows
  recovered from the extracted FASTA, and unfiltered regional signal/peaks/genes.
  Export hit coordinates and per-query plotting status. No new genome scan.
- Document one-best-peak lookup semantics, all-gene plotting, feature/nearest-TSS
  ambiguity, and the difference between candidates' distance rescue and lookup.
- Validated GFF/GTF, successful nf-core results and existing environments/bundles
  can be reused. Actual full-data rendering remains a server-side validation.

## v1.2.4 — downstream logs, clearer plots and deployment modes

- Capture each downstream CLI invocation's stdout/stderr in a unique project
  `logs/downstream/<type>/` log, with command, version, config and exit status.
  Translate owned runtime/help/plot messages to English; retain biological IDs.
- Create feature output directories on demand. Respect peakset selection for
  signal, GO preparation, Python annotation and motif spacing; avoid silent IDR
  fallback. Honor output overrides in IDR and motif commands too.
- Replace ChIPseeker's diverging percentage TSS chart with explicit directional
  bins, counts, percentages and a TSV containing boundaries/denominators. Add a
  standalone redraw helper. Label deepTools center plots as peak midpoints.
- Surface failures from IDR, bedtools, deepTools, integration and track rendering.
  Reject mixed-condition IDR input and mark incomplete IDR runs. Recover actual
  MACS signalValue from replicate peaks rather than reusing the p-value rank
  column; handle broadPeak global-IDR column placement. Correct the Python
  consensus overlap search for long/nested peaks without a fixed-length cutoff.
- Add `chromaflow setup` with offline/prepare/conda modes and runtime/downstream
  components. Offline runtime needs only Python/Java; nf-core process Conda envs
  are unnecessary with the local container bundle. Keep existing environments.
- Document both bundles, shared /opt installation, target validation, config,
  upgrades and downstream branches. Use a generic species template; leave
  measured values unset, optional organelles empty and online KEGG disabled.
- Include the 1.2.3 TxDb identity hotfix. Preserve the GFF normalization and GTF
  conversion implementations; validated references and existing bundles can be
  reused. This update does not require repeating a successful nf-core run.
- New archive checksums use LF. Tests cover logger errors/uniqueness, branch
  isolation, environment plans and real-bedtools IDR output processing. Full
  ChIPseeker/MEME/deepTools runs and fresh Conda solving remain target checks.

### 1.2.3 TxDb identity hotfix (2026-09-23)

- Fix the acceptance check's assumption that TxDb `tx_name`/`gene_id` must
  equal GFF `ID`/`Parent`. Match only explicit GFF ID, Name, transcript_id and
  gene_id relationships, with separate transcript and gene namespaces.
- Retain exact gene-edge, sequence, strand, TSS and exported-span checks.
  Reject ambiguous names, missing/duplicate transcripts and split identities;
  do not strip species-specific version suffixes or match by coordinates alone.
- Export lossless gene edges and auditable GFF-to-TxDb identity maps. The
  existing ChIPseeker importer and its public naming policy remain unchanged.
- Normalized GFF/GTF and offline bundles need no regeneration for this fix.
  The original 1.2.3 archives remain immutable; distribute a separate hotfix.

## v1.2.3 — conservative GFF3 normalization and validated GTF regeneration

- Preserve every source feature record, including existing TE/unknown gene types,
  IDs, multiple Parent edges, coordinates, phase and attributes. No species or
  identifier-prefix branches. Unknown types are not missing nodes.
- Check syntax, incompatible ID collisions, dangling parents and cycles; allow
  compatible discontinuous features. Recover missing genes only with explicit
  `--repair-missing-parents` and consistent evidence for the original Parent ID.
- Add only evidence-supported exons and terminal UTRs. Do not invent exons from
  transcript spans, shrink coordinates or classify internal CDS gaps as UTRs.
- Fix `auto` to builtin; validate source-record preservation for external engines
  too. Default filename is `.normalized.gff3`; explicit legacy names still work.
- Atomically publish validated files; record hashes, actual engine/version, policy,
  offline role map, warnings and a per-addition TSV audit. Cache reuse requires
  matching content and rules; unknown/stale outputs are regenerated.
- Regenerate GTF when its source/rules/converter change; verify transcript-to-gene
  assignments before publishing. Reference preparation retains the registered
  normalized GFF instead of replacing it with its original source.
- Expose TxDb construction warnings. Add a server-side TxDb/TSS acceptance tool.
- Remove fixed MACS_GSIZE/READ_LENGTH placeholder examples from the main config;
  explicit overrides remain supported. Existing nf-core bundles remain usable.
- See `GFF_NORMALIZATION.md` for policy, compatibility boundaries and migration.

## v1.2.2 — fixed ATAC offline plugin version

- Pin the unversioned nf-validation declaration in the downloaded copy of
  nf-core/atacseq 2.1.2 to 1.1.3. Explicit upstream versions remain unchanged;
  unknown releases still require a fixed plugin declaration.
- Verify actual ATAC 2.1.2 help startup with Nextflow 25.10.4 in offline mode
  using the locally installed nf-validation 1.1.3; all 21 WSL tests pass.
- Document preparing, transferring and selecting separate ChIP/DAP and ATAC bundles.
  Full ATAC container download and target-server sample analysis remain deployment checks.

## v1.2.1 — offline plugin directory fix

- Recognize Nextflow's `id-version` plugin installation directories while retaining
  `id@version` in plugin declarations and install commands. Version 1.2.0 incorrectly
  rejected a successfully installed `nf-validation@1.1.3` plugin.
- Resume the same bundle after this failure by reusing its existing workflow,
  containers and installed plugin; no cache deletion or directory renaming is needed.
- Correct the test fixtures to match Nextflow's actual layout and cover recovery
  without additional downloads, wrong versions and empty JAR files.

## v1.2.0 — local nf-core deployment

- Add offline prepare/check/pack commands for fixed-release Linux container bundles.
- Bundle pipeline code, local nf-core configs, Nextflow runtime, pinned plugins and
  Singularity images; include a SHA256 inventory and relocatable tar export.
- Add run --offline --bundle, --dry-run and --yes; check bundle/type/version before launch.
- Keep online Conda/container runs; full offline Conda environment migration is not included.
- Write Nextflow debug logs directly under OUTPUT_DIR/logs and report the path on failure.
- Honor the CLI --profile in container runs and check the executable actually selected.
- See OFFLINE.zh-CN.md for deployment, verification and upgrades. Linux data analysis
  must be validated on the target host; bundle checks are not biological validation.

## v1.1.0 — maintenance release

### 1. Sample records and read pairing (`make_samplesheet.py`, `01`)
- Sample files can now use a **header row** with any column order, optional
  columns and aliases (`sample|name`, `ab|target`, `control|input`, `rep`,
  `fastq|path|dir`, `r1/r2`). Positional 4–5 column files still work.
- New option to give `fastq_1` / `fastq_2` explicitly and skip pairing.
- **Pairing is no longer a fixed suffix list.** The right-most mate token in
  the filename is detected, covering `_1/_2`, `_R1/_R2`, `.1/.2`, `-1/-2`,
  `_R1_001`, `_good_1`, `_clean_1`, `_val_1`, `_1P/_2P`, `sampleR1`, …
  Replicate-style names (`sample_rep1`/`sample_rep2`) are not mistaken for mates.
- Overrides: `#@ pair_pattern: _good_{R}.fq.gz` / `#@ single_end: true` in the
  sample file, `FASTQ_PAIR_PATTERN` / `FASTQ_SINGLE_END` in the config, or
  `--pair-pattern` / `--single-end` on the command line.
- `--list-fastq` previews pairing without writing; `--strict` fails on issues.
- The sample file is located automatically (`--samples` > `SAMPLES_FILE` >
  `samples_<type>.{txt,tsv,csv}` > `samples.{txt,tsv,csv}`).

### 2. Config no longer records script outputs
- Config declares **sources only**: `GENOME_FASTA_SOURCE`, `GENOME_GFF3_SOURCE`,
  `ANNOTATION_INFO`, `MEME_DB`, `ORGANELLE_SOURCES`.
- Derived files are resolved as
  `config explicit > manifest > directory scan > default name`.
- Steps 00/00b register their products in
  `<OUTPUT_DIR>/.chromaflow/reference.manifest`; directory scanning prefers
  `*.agat.gff3` over plain `*.gff3` and indexed genomes over unindexed ones.
- `ensure_chrom_sizes` regenerates chrom.sizes from the `.fai` on demand.
- Old configs with explicit paths keep working unchanged.

### 3. AGAT decoupled
- New `lib/steps/normalize_gff.py`: dependency-free normalizer (adds exons from
  CDS, infers 5'/3' UTRs, synthesizes gene features, sorts output).
- `normalize-gff --engine auto|agat|builtin|gffread`; `auto` falls back to the
  built-in engine when AGAT is absent.
- AGAT moved out of `ChromaFlow_ref` into an optional `ChromaFlow_agat` env
  (`WITH_AGAT=true bash setup_environments.sh`).
- Step 05 no longer hard-fails on an un-normalized GFF3; it warns instead.

### 3b. Compression handling (`detect_compression`, `ensure_plain_file`)
- Source files are classified by **magic bytes, not extension**: gzip/bgzip,
  bzip2, xz, zstd and zip are decompressed automatically.
- Already-decompressed `.fa` / `.gff3` are used **in place, with no copy**
  (previously a plain `.fa` source was copied to `<SPECIES>_nuclear.fa`).
- Files that are gzipped despite a plain name (e.g. `x.gene.gff3` containing
  gzip data) are detected and written to `x.gene.plain.gff3` without touching
  the source.
- Decompressed products are cached and reused unless the source is newer or
  `--force` is given; organelle FASTAs and the `normalize-gff` input go
  through the same path.
- Directory auto-detection now skips compressed files for formats that must be
  plain text (GFF3/GTF/chrom.sizes/FASTA).

### 4. Optional steps are genuinely skippable
- `reference --skip-organelle | --skip-gtf | --force | --strict`.
- Missing organelle genomes now warn and continue instead of exiting.
- Undefined/empty `ORGANELLE_SOURCES` is safe under `set -u`.
- `READ_LENGTH` has a default (previously an unbound variable in steps 00/02).
- `doctor` reports optional features as `[INFO]`/`[SKIP]`, not failures.

### 4b. nf-core execution mode (`NFCORE_PROFILE`) — bug fix
- **Bug:** `chromaflow run` always dispatched to `02_run_nfcore.sh`, which
  hard-requires Apptainer/Singularity *before* reading the profile. Setting
  `NFCORE_PROFILE=conda` (as `setup_environments.sh` advised) therefore exited
  with "未检测到 Apptainer/Singularity", and `02_run_nfcore_conda.sh` was
  reachable only by invoking it manually.
- `run` now routes on the profile: `conda`/`mamba` → the conda script,
  everything else → the container script.
- New `--profile` flag; resolution order `--profile > NFCORE_PROFILE > singularity`.
  Comma combinations (`"singularity,test"`) are supported.
- The container script checks the engine implied by the profile
  (singularity/apptainer/docker/podman) instead of always demanding
  Singularity, and suggests `--profile conda` on failure.
- `02_run_nfcore_conda.sh` no longer hardcodes `-profile conda`; the Tsinghua
  channel mirror can be disabled with `NFCORE_CONDA_MIRROR=none`.
- Defaults added for `NFCORE_PROFILE`, the pipeline versions and `MAX_CPUS` /
  `MAX_MEMORY` / `MAX_TIME` (previously unbound-variable crashes if a config
  omitted them); missing `MACS_GSIZE` now gives a readable error.
- `NXF_APPTAINER_CACHEDIR` / `APPTAINER_CACHEDIR` are exported alongside
  `NXF_SINGULARITY_CACHEDIR` so images are not downloaded twice.
- `doctor` prints the resolved profile, the path it selects, and checks the
  matching engine. Documented in `chromaflow help run`, `help config`, the
  config template and both READMEs.

### 4c. Measured values replace hand-entered constants
- `MACS_GSIZE` is now measured by step 00 with khmer `unique-kmers.py` (the
  command the old logs printed as a manual instruction), cached against a
  genome+k signature; `--force` re-measures, `--no-measure` skips.
  Fallback without khmer: non-N base count (upper bound).
- `READ_LENGTH` is measured by step 01 from the fastq files in the samplesheet
  (2000 reads sampled per file, maximum taken).
- Step 02 passes `--macs_gsize` when known; otherwise it omits it and passes
  `--read_length`, letting nf-core's own KHMER_UNIQUEKMERS process compute the
  value. Neither path requires manual entry.
- Values go to `<OUTPUT_DIR>/.chromaflow/derived.env`, **not** back into the
  user's config. Precedence: config explicit > derived.env > measure now.
  Opt-in `chromaflow reference --write-config` appends a delimited, idempotent
  block (values quoted, config backed up) for those who want them archived.
- `GC_CONTENT` **removed**: it was only ever echoed in a log line, so the log
  reported whatever the user had typed. Step 04's background model comes from
  `fasta-get-markov`, and it now reports the GC read back from that model.
- `doctor` shows each measured value with its origin (config / measured / not
  yet measured).

### 4d. Species parameters clarified + early chromosome-name check
- Documented that only `SPECIES_SHORT` is functional (prefix for derived
  files); `SPECIES_NAME` and `GENOME_VERSION` are log labels. Nothing parses
  them to infer organism or gene-ID format.
- `doctor` now compares chromosome names between the genome FASTA and the
  GFF3 and FAILs on a complete mismatch (the classic `Chr1` vs `chr01` vs `1`
  problem when switching to rice). Previously this surfaced only at step 05,
  after the whole nf-core run.
- Added a worked rice example to both READMEs.

### 4e. conda channel order fix (caused a real run failure)
- **Bug:** the generated `.condarc` listed `bioconda` before `conda-forge`.
  nf-core checks this and emits only a WARN, then the run fails much later
  inside a process — observed as `PHANTOMPEAKQUALTOOLS` failing with
  `/usr/bin/which: no run_spp.R in (...)` after ~all alignment work was done.
  run_spp.R comes from the process's conda env, which the solver had built
  with mismatched cross-channel builds.
- Channel order corrected to `conda-forge > bioconda > defaults` with
  `channel_priority: strict`. A stale `.condarc` written by an older
  ChromaFlow is detected (bioconda first) and rebuilt automatically.
- Added `--skip-spp` / `SKIP_SPP=true` to skip phantompeakqualtools, which is
  the most conda-fragile module and only produces the NSC/RSC QC metric.
- Added `-- <args>` passthrough on `chromaflow run` for arbitrary
  nextflow/nf-core arguments (`-- --skip_preseq -resume`).

### 4j. MULTIQC / kaleido pin
- **Bug:** with unpinned dependencies the solver installs kaleido >= 1.0 and
  python 3.14 into the MultiQC environment. Kaleido 1.0 removed its bundled
  chromium and requires a system Google Chrome, so MULTIQC fails with
  `RuntimeError: Kaleido requires Google Chrome to be installed` while
  exporting static plots — after every analysis process has succeeded.
- The generated conda config now pins `python=3.12` and
  `python-kaleido=0.2.1` for the MULTIQC process, with the multiqc version
  read from the pipeline assets. The `withName` selectors were verified not to
  over-match MULTIQC_CUSTOM_PEAKS / MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.
- Added `--skip-multiqc`.

### 4h. conda downloads: mamba actually enabled + retry/mirror options
- **Bug:** the script exported `NXF_CONDA_USE_MAMBA=true`, but no such
  Nextflow environment variable exists — the real setting is `conda.useMamba`
  in a config file. Plain conda had been doing every environment build.
  conda's downloader aborts the whole environment when a transfer drops
  (`curl: (18) end of response with N bytes missing`), which is what killed
  the phantompeakqualtools build at ~12 KB/s. `useMamba` (or
  `useMicromamba`) is now written into the generated config when available.
- `.condarc` gained download resilience: `remote_max_retries: 5`,
  `remote_read_timeout_secs: 600`, `remote_connect_timeout_secs: 30`,
  `remote_backoff_factor: 2`.
- `NFCORE_CONDA_MIRROR` now accepts `tuna|bfsu|nju|sjtu|none` or a custom
  mirror base URL, so a throttling mirror can be swapped without editing code.
- **Fixed a heredoc bug found in testing:** interpolating a newline-terminated
  variable immediately before `EOF` meant the terminator was never on its own
  line, so the first heredoc swallowed the else-branch and no `.condarc` was
  written at all. The retry keys are now inlined.

### 4g. conda env build timeout — the actual root cause of 4e/4f symptoms
- **Bug:** nextflow's `conda.createTimeout` defaults to 20 minutes.
  Building the phantompeakqualtools env (r-base + Bioconductor + boost,
  >1 GB) exceeds that on a slow link, so conda is killed and the run fails
  with `Failed to create Conda environment ... status : 120`
  (plus a Groovy `IOException: Stream closed` from the killed reader thread).
- **The partial env directory is left behind.** Nextflow decides an
  environment exists from the directory alone, so the next run reuses the
  broken one — which is where the earlier `which: no run_spp.R` failure
  actually came from. The channel-order fix in 4e was correct practice but
  was *not* the cause of that symptom.
- Fixes: `conda.createTimeout` set to 3 h (`NFCORE_CONDA_TIMEOUT` to
  override) in the generated config, and `clean_partial_conda_envs` deletes
  env dirs lacking `conda-meta/history` before launching.
- `--no-conda-patch` now only drops the dependency pins; the timeout always
  applies.
- Corrected the framing of nf-core's "problem with your Conda configuration"
  warning: `conda config --show channels` merges system, user and CONDARC
  files, so it can persist even with a correct project `.condarc`, and module
  environment.yml files declare their own channels. It is usually noise.

### 4f. conda dependency pins (deeptools / numpy 2.x)
- **Bug:** under `-profile conda`, DEEPTOOLS_COMPUTEMATRIX fails with
  `AttributeError: module 'numpy' has no attribute 'NAN'`. nf-core's conda
  specs pin the tool but not its dependencies, so the solver installs numpy
  2.x while deeptools 3.5.x still calls `np.NAN` (removed in numpy 2.0).
  Container mode is unaffected. The old setup_environments.sh header mentioned
  this as a manual `pip install 'numpy<2.0'` note; it is now handled
  automatically.
- Conda mode generates `<OUTPUT_DIR>/.chromaflow/nfcore_conda_patch.config`
  and passes it with `-c`: pins `numpy=1.26` for `DEEPTOOLS_*` processes and
  sets `conda.createTimeout` (default 90 min, `NFCORE_CONDA_TIMEOUT` to
  override). The deeptools version is read from the downloaded pipeline
  assets so the pin tracks nf-core rather than being hardcoded.
- Added `--skip-plot-profile` and `--no-conda-patch` escape hatches.
- **Fixed a bug introduced in 4e:** the stale-`.condarc` detector matched
  `  - bioconda` anywhere in the first five lines, which is also true of the
  *correct* file, so it rebuilt and warned on every run. It now compares which
  channel comes first.

### 4i. QC fallback when phantompeakqualtools is skipped
- Step 03 already gated its NSC/RSC block on the output directory existing, so
  `--skip-spp` never broke it — but the log went silent. It now falls back to
  listing FRiP values found under the results tree and points at
  plotFingerprint / plotProfile / MultiQC / the IDR pairwise table.
- Documented in both READMEs what `--skip-spp` actually costs (only NSC/RSC,
  the QualityTag and the cross-correlation plot; MACS3 estimates its own
  fragment size so peaks are unchanged) and why IDR is not a substitute:
  IDR measures replicate agreement, not signal-to-noise, and uniformly poor
  replicates can agree with each other.

### 4k. Documentation: conda activate instead of mamba activate
- All docs and hint messages now say `conda activate <env>`. mamba 2.x is a
  standalone binary whose shell function is installed by `mamba shell init`,
  separate from `conda init`, so `mamba activate` fails with "running as a
  subprocess and can't modify the parent shell" on setups that only ran
  `conda init`. No script ever called activate, so nothing was broken — but
  the instructions did not work as written for those users.

### 4l. Scientific-notation genome coordinates (step 05 → 06/07)
- **Bug:** `chromaflow candidates` failed with
  `ValueError: invalid literal for int() with base 10: '8e+06'`.
  In `05_peak_annotation.R`, `start(g) - 1` promotes GRanges integer
  coordinates to double; `write.table` renders doubles with `as.character()`,
  which uses scientific notation for round numbers, so a gene starting at
  8,000,000 was written to genes.bed as `8e+06`. Only affects genomes with
  genes at such coordinates, which is why it surfaced late.
- Writer fixed at the source: `as.integer(start(g) - 1L)` plus a global
  `options(scipen = 999)` covering every table the script writes
  (peak_annotation_full.tsv had the same exposure).
- Readers made tolerant so existing outputs do not require rerunning step 05:
  `06_integrate.sh` parses via `int(float(...))` and skips unparseable lines
  with a warning; `plot_tracks.sh` (both blocks) does the same;
  `06_signal_deeptools.sh` normalizes with `printf "%d"` in its awk.

### 4m. --peakset consensus: peak_score was silently all zero
- **Bug:** with `--peakset consensus`, every `peak_score` in the annotation
  output was 0. nf-core builds `consensus_peaks.bed` with
  `awk 'FNR>1 {print $1,$2,$3,$4,"0","+"}'` — column 5 is a hardcoded `"0"`.
  `get_consensus_with_signal` exists to back-fill signalValue from the
  replicate narrowPeaks, but it required bedtools and fell back to the raw
  consensus file **silently** when bedtools was missing. `ChromaFlow_R`, the
  documented environment for `annotate`, does not include bedtools — so the
  documented workflow hit the silent path.
- Added `lib/steps/support/consensus_signal.py`: a dependency-free Python
  back-fill (bisect-based interval overlap, max signalValue per consensus
  interval), used automatically when bedtools is absent. Verified to match
  the bedtools result, including half-open boundary semantics.
- Every fallback path now warns explicitly instead of degrading silently, and
  a failed back-fill no longer leaves a partial file in the cache location
  (which would have been reused as valid on the next run).
- `bedtools` added to the ChromaFlow_R environment recipe; `doctor` now
  reports which back-fill path is available.

### 4n. annotation gene vs nearest-TSS gene are now distinguished
- ChIPseeker's `annotation` (feature the peak overlaps, possibly a
  neighbouring gene) and `geneId` (nearest TSS) can name different genes.
  Step 06 filtered on the former while reporting the latter, so a peak in
  gene A's intron that sits 1.4 kb upstream of gene B's TSS was recorded
  against B but rejected as "Intron" — a real upstream site for B was lost.
- Step 05 now emits `feature_simple`, `feature_geneId` and
  `gene_consistent`, and logs how many peaks disagree plus the PROMOTER_UP
  hint. Step 06 emits `feature_gene` and `pass_reason` in
  candidate_targets.tsv and reports the mismatch count.
- New `--max-tss-dist N` / `TARGET_MAX_TSS_DIST`: a peak also passes the
  position filter when it lies within N bp upstream of geneId's TSS,
  regardless of what it overlaps. Default 0 (off) so existing results are
  unchanged. Documented alongside the alternative (widen PROMOTER_UP).

### 5. Help and documentation
- `lib/help.sh`: `chromaflow help`, `chromaflow help <command>`,
  `chromaflow <command> --help`, plus topics `config`, `samples`, `peakset`,
  `envs`; `chromaflow version`.
- `README.md` is now English and links to `README.zh-CN.md` (Chinese), which
  has been updated with the same new sections.
