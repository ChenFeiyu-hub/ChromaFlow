# Changelog

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

### 5. Help and documentation
- `lib/help.sh`: `chromaflow help`, `chromaflow help <command>`,
  `chromaflow <command> --help`, plus topics `config`, `samples`, `peakset`,
  `envs`; `chromaflow version`.
- `README.md` is now English and links to `README.zh-CN.md` (Chinese), which
  has been updated with the same new sections.
