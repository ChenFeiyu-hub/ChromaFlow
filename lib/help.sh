#!/usr/bin/env bash
# ============================================================================
# lib/help.sh — `chromaflow help [command]` 的帮助文本
# ============================================================================
# 顶层帮助保持一屏可读；每个子命令的细节用 `chromaflow help <command>` 展开。
# ============================================================================

CHROMAFLOW_VERSION="${CHROMAFLOW_VERSION:-1.2.9}"

help_main() {
cat <<'EOF'
ChromaFlow — ChIP-seq / DAP-seq / ATAC-seq downstream analysis toolkit
Peaks → reproducibility → motifs → annotation → candidate target genes.

USAGE
  chromaflow <command> [options]
  chromaflow help <command>            Detailed help for one command
  chromaflow doctor                    Check tools, references and outputs

SETUP
  setup --mode offline|prepare|conda   Install the selected environment set
  init <dir> [--type T]                Create a project skeleton + config
  doctor [-c cfg] [--type T]           Pre-flight self-check
  offline prepare|check|pack           Prepare and verify a local nf-core bundle
  version                              Print version

REFERENCE (optional, run once per species)
  reference                            00  Genome / GTF / chrom.sizes
  normalize-gff                        00b Normalize GFF3 (AGAT or built-in)

MAIN PIPELINE
  samplesheet [--type T]               01  Build the nf-core samplesheet
  run [--type T]                       02  Run nf-core/chipseq | atacseq
  batch init|plan|run|status|reports    DAP-seq batch workflow and diagnostic reports

DOWNSTREAM
  idr [--min-pairs N] [--exclude-rep R]   03  Reproducible peaks (IDR)
  motif                                   04  MEME-ChIP de novo motifs
  annotate                                05  Peak annotation (ChIPseeker/R)
  annotate-py                             05  Same, pure Python (no R needed)
  candidates [--features F] [--pval P]    06  motif x annotation -> candidates
  lookup --gene-list F [--plot]           09  Query prior genes for binding
  tracks [--genes G | --gene-list F]      07  Single-gene browser tracks
  go-prep [--set promoter|dist|all]       08a Export gene lists for GO tools
  go-plot --input FILE [--top N]          08b Plot GO enrichment results
  signal                                  (opt) deepTools aggregate profiles
  spacing                                 (opt) motif spacing analysis

GLOBAL OPTIONS
  -c, --config <file>   Project config (default: $CHROMAFLOW_CONFIG, ./chromaflow.config)
  -t, --type <type>     chipseq (default) | atacseq | dapseq
      --profile <name>  nf-core execution mode for `run`: singularity
                        (default) | apptainer | docker | podman | conda | mamba
      --offline         Run from a prepared container bundle without online lookups
      --bundle <dir>    Bundle directory (offline prepare/check/pack and run)
      --dry-run         Validate and print the run command without starting analysis
      --yes             Skip the run confirmation prompt
  -o, --outdir <dir>    Redirect THIS step's output to <dir>
      --peakset <spec>  Which peaks downstream uses:
                          idr (default) | idr:0.05 | consensus | <file path>
                        Non-default peaksets write to suffixed dirs, e.g.
                        annotation_consensus/, so results can coexist.
  -h, --help            Show this help

OPTIONAL STEPS
  reference / normalize-gff / signal / spacing are optional. Skipping them
  does not break the main path — for example, a species without an organelle
  genome can run `chromaflow reference --skip-organelle`, and annotation
  works with a plain GFF3 if you never run normalize-gff.

ENVIRONMENTS (activate suitable tools; see `chromaflow help envs`)
  conda activate ChromaFlow_idr && chromaflow idr

DOCS
  README.md (English) · README.zh-CN.md (Chinese) · ARCHITECTURE.md
EOF
}

help_envs() {
cat <<'EOF'
ENVIRONMENTS
  ChromaFlow_runtime     Python, Java                    -> offline run
  nfcore_env             nf-core CLI (prepare mode) or    -> offline prepare or
                         Nextflow (conda mode)               online Conda run
  ChromaFlow_ref         samtools, gffread, bedtools,     -> reference, samplesheet,
                         python, matplotlib                  annotate-py, go-*, spacing
  ChromaFlow_idr         idr (python 3.10, numpy<1.24)    -> idr
  ChromaFlow_motif       MEME Suite, bedtools, samtools   -> motif
  ChromaFlow_R           ChIPseeker, clusterProfiler      -> annotate
  ChromaFlow_downstream  deepTools, pyGenomeTracks        -> candidates, lookup,
                                                             tracks, signal
  ChromaFlow_agat        AGAT (optional, separate env)    -> normalize-gff --engine agat

  Activate with `conda activate <env>`. `mamba activate` only works if the
  mamba shell function is installed (mamba 2.x needs `mamba shell init`,
  separate from `conda init`); the environments themselves are identical
  either way, so conda activate is the portable choice.

  Install:  chromaflow setup --mode offline --components runtime
            chromaflow setup --mode prepare --components runtime
            chromaflow setup --mode conda --components runtime
            chromaflow setup --components downstream [--skip-r] [--with-agat]
            chromaflow setup --dry-run              # inspect without installing

  Existing environments are kept. Installation needs package access; an offline
  workflow bundle does not include these separate downstream environments.

  AGAT is optional and installed in its own environment on purpose: its Perl
  dependency tree often conflicts with other packages. Without AGAT,
  `chromaflow normalize-gff` defaults to a conservative Python normalizer.
EOF
}

help_config() {
cat <<'EOF'
PROJECT CONFIG (chromaflow.config)

  Declare only what you already have; everything the pipeline produces is
  discovered automatically.

  Required
    SPECIES_SHORT="Cre"                 Prefix for derived files
                                        (<SHORT>_combined.fa, <SHORT>_genes.gtf).
                                        The only species field with an effect —
                                        change it when you change species or the
                                        derived files will collide.
    PROJECT_DIR="/data/projects/X"      Anchor for inputs/ and outputs/

  Labels only (log headers; wrong values do not affect analysis)
    SPECIES_NAME="Chlamydomonas_reinhardtii"
    GENOME_VERSION="v6.1"

  Species identity comes from the genome and GFF3 files, not from these
  strings — nothing parses them to infer gene ID format or organism.

  Source files (things you downloaded; compressed or not, either is fine)
    GENOME_FASTA_SOURCE=".../genome.fa"        # or .fa.gz/.bz2/.xz/.zst
    GENOME_GFF3_SOURCE=".../genes.gff3"        # or .gff3.gz ...
    ANNOTATION_INFO=".../annotation_info.txt"    # optional, for GO
    MEME_DB=".../JASPAR2024_plants.meme"         # optional, for tomtom
    ORGANELLE_SOURCES=("chloroplast:NC_005353.1" "mitochondria:NC_001638.1")
                                                 # optional, leave empty to skip

  Auto-resolved (do NOT put script outputs in config)
    GENOME_COMBINED, GENOME_NUCLEAR, CHROM_SIZES, GENOME_GTF, GENOME_GFF3
    Resolution order:
      explicit value in config  >  manifest registered by step 00/00b
                                >  scan of REFERENCE_DIR  >  default name
    Manifest: <OUTPUT_DIR>/.chromaflow/reference.manifest

  Execution mode (step 02)
    NFCORE_PROFILE=singularity          singularity|apptainer|docker|podman
                                        |conda|mamba — also selects whether
                                        `run` uses the container or conda path
    NFCORE_CONDA_MIRROR=tuna            set to `none` to keep your channels
    NFCORE_CHIPSEQ_VERSION, NFCORE_ATACSEQ_VERSION

  Measured automatically — normally absent from the config
    MACS_GSIZE      measured by `reference` (khmer unique-kmers)
    READ_LENGTH     measured by `samplesheet` (sampled from the fastq files)
    GENOME_GC       measured by `reference` (informational only)
    Setting any of them in the config disables measurement for that value.
    Stored in <OUTPUT_DIR>/.chromaflow/derived.env.

  Common parameters
    MAX_CPUS, MAX_MEMORY, MAX_TIME
    PROMOTER_UP / PROMOTER_DOWN, TARGET_DIST_KB, TARGET_FEATURES
    FIMO_PVAL, MEME_SEARCHSIZE, KEGG_ORGANISM
    SAMPLES_FILE, FASTQ_PAIR_PATTERN, FASTQ_SINGLE_END
    SKIP_ORGANELLE=true                 Never merge organelle genomes

  Template: templates/chromaflow.config
EOF
}

help_samples() {
cat <<'EOF'
SAMPLE RECORD FILE

  Located automatically: --samples <file> > SAMPLES_FILE in config >
  <PROJECT_DIR|INPUT_DIR>/samples_<type>.{txt,tsv,csv} > samples.{txt,tsv,csv}

  Three accepted layouts (TAB, comma or whitespace separated):

    1) Header-based (recommended) - any column order, optional columns
         sample_id  antibody  control_id  replicate  fastq_path
         Dof_ChIP   crDOF     Input       1          /data/raw/Dof_rep1

    2) Positional (legacy, still supported)
         sample_id <TAB> antibody <TAB> control_id <TAB> replicate <TAB> fastq_path

    3) Explicit read files (skips pairing)
         sample_id  replicate  fastq_1  fastq_2

  Column aliases: sample|name, ab|target, control|input, rep,
                  fastq|path|dir, r1|read1, r2|read2

  READ PAIRING
    Detected automatically from the right-most mate token in the filename,
    so no suffix list to maintain. Covered patterns include:
      _1/_2   _R1/_R2   .1/.2   -1/-2   _R1_001/_R2_001   _good_1/_good_2
      _clean_1/_clean_2   _val_1/_val_2   _1P/_2P   sampleR1/sampleR2
    Files whose mate is missing are treated as single-end with a warning.
    Names like sample_rep1 / sample_rep2 are NOT mistaken for mates.

    Override when a naming scheme is unusual ({R} is replaced by 1 and 2):
      In the samples file:   #@ pair_pattern: _good_{R}.fq.gz
                             #@ single_end: true
      In the config:         FASTQ_PAIR_PATTERN="_R{R}_001.fastq.gz"
                             FASTQ_SINGLE_END=true
      On the command line:   chromaflow samplesheet --pair-pattern '_R{R}_001.fastq.gz'
                             chromaflow samplesheet --single-end

  Preview pairing without writing anything:
      chromaflow samplesheet --list-fastq

  Templates: templates/samples_chipseq.txt, samples_atacseq.txt, samples_dapseq.txt
EOF
}

help_command() {
    case "$1" in
    batch) python3 "${CHROMAFLOW_HOME}/lib/batch.py" --help ;;
    setup) bash "${CHROMAFLOW_HOME}/setup_environments.sh" --help ;;
    offline) cat <<'EOF'
chromaflow offline prepare|check|pack --bundle DIR [options]

  Prepare on a connected Linux machine (Java + nf-core CLI required):
    chromaflow offline prepare --type dapseq --bundle /data/offline/chipseq-2.1.0
    chromaflow offline check --bundle /data/offline/chipseq-2.1.0 --verify --runtime-check
    chromaflow offline pack --bundle /data/offline/chipseq-2.1.0 --archive /data/chipseq-2.1.0.tar

  Transfer the tar and its .sha256 to the offline Linux server; verify and extract:
    sha256sum -c chipseq-2.1.0.tar.sha256
    tar -xf chipseq-2.1.0.tar -C /data/offline

  Run there (Java and Apptainer/Singularity must already be installed):
    chromaflow offline check --bundle /data/offline/chipseq-2.1.0 --verify --runtime-check
    chromaflow run -c chromaflow.config --type dapseq --profile apptainer \
      --offline --bundle /data/offline/chipseq-2.1.0

  PREPARE OPTIONS
    --revision VERSION          Pipeline release (default from type/project config)
    --nextflow-version VERSION  Fixed Nextflow release (default 25.10.4)
    --parallel-downloads N      Concurrent container downloads (default 4)
    --mode singularity         Default: pipeline, shared configs, runtime, plugins, images
    --mode code-only           Downloads no analysis environments; cannot run --offline

  CHECK OPTIONS
    --verify                   Verify SHA256 of all files, including large images
    --runtime-check            Linux: start bundled Nextflow with --help in offline mode
                               Does not execute analysis or replace a small-data test

  This release supports full offline runs with the singularity/apptainer profile.
  Conda online mode remains available. Copying .conda_envs is not a portable
  offline deployment and is not automatically certified by this command.
  Bundle excludes FASTQ, references and ChromaFlow downstream environments.
  Keep all those files locally available on the target server.
  For version upgrades prepare a NEW bundle directory, validate it, then change
  NFCORE_OFFLINE_BUNDLE and the matching NFCORE_*_VERSION in the project config.
  See docs/OFFLINE.md for current deployment instructions and troubleshooting.
EOF
;;
    init) cat <<'EOF'
chromaflow init <project_dir> [--type chipseq|atacseq|dapseq]

  Creates inputs/{raw,reference}, outputs/, a chromaflow.config with
  PROJECT_DIR filled in, and a matching sample record template.
EOF
;;
    reference) cat <<'EOF'
chromaflow reference [--skip-organelle] [--skip-gtf] [--force] [--strict]
  Env: ChromaFlow_ref                                        (OPTIONAL STEP)

  Prepares the source genome/GFF3, optionally merges organelle genomes,
  builds the faidx index, chrom.sizes and a GTF.

  Compression is detected from the file's magic bytes, not its extension:
  gzip/bgzip, bzip2, xz, zstd and zip are decompressed automatically, a file
  that is already plain text is used in place (no copy), and a file that is
  gzipped despite being named .gff3 is still handled correctly.

  --no-measure       Skip measuring effective genome size and GC
  --kmer N           k for the genome-size measurement (default: measured
                     read length, or 150 if not known yet)
  --write-config     Also write measured values back into chromaflow.config
                     as a delimited block (off by default — see below)

  --skip-organelle   Nuclear genome only. Use this when the species has no
                     suitable chloroplast/mitochondrion assembly. Missing
                     organelle files are skipped with a warning by default.
  --skip-gtf         Skip GFF3 -> GTF (only needed by step 02).
  --force            Rebuild products that already exist.
  --strict           Fail instead of warning when an organelle is missing.

  MEASURED VALUES
    The effective genome size (MACS -g) is measured with khmer's
    unique-kmers.py — the command the old logs asked you to run by hand — and
    the genome GC is measured in the same pass. Results go to
    <OUTPUT_DIR>/.chromaflow/derived.env, NOT into your config:
      config   = what you declare      (wins if set)
      derived  = what the scripts measured
    Re-measured only when the genome or k changes, or with --force.
    Without khmer, the fallback is a non-N base count (an upper bound), and
    step 02 can instead pass --read_length so nf-core measures it internally.

  Products are registered in <OUTPUT_DIR>/.chromaflow/reference.manifest;
  you do not need to copy their paths back into chromaflow.config.

  Already have a genome, GTF and chrom.sizes? Put them in REFERENCE_DIR and
  skip this step entirely — they are auto-detected.
EOF
;;
    normalize-gff) cat <<'EOF'
chromaflow normalize-gff [--input F] [--output F] [--engine E] [--force]
  Env: ChromaFlow_ref (or ChromaFlow_agat)                   (OPTIONAL STEP)

  Preserve source records and validate IDs/Parent relationships. Add exons
  from CDS/UTR evidence and terminal UTRs from existing exons/CDS/strand.

  --engine auto      Always use builtin (reproducible conservative default)
          builtin    Python 3; preserve original types, coordinates and edges
          agat       External AGAT; same preservation checks apply
          gffread    External gffread; same preservation checks apply
  --repair-missing-parents  Explicitly recover missing shared gene Parent IDs
                            only with consistent transcript evidence
  --type-map FILE    Local JSON mapping of additional feature types to roles
  --force            Rebuild; does not permit destructive structural changes

  Outputs: *.normalized.gff3, *.report.json and *.changes.tsv.
  Cache requires matching input/rules/engine/output hashes. Missing Parent,
  cycles, malformed rows or destructive external-engine changes fail.

  Validated output is registered in the manifest. Explicit GENOME_GFF3 in
  config takes precedence. Run reference afterwards to regenerate GTF;
  GTF/TxDb compatibility must be checked separately. See GFF_NORMALIZATION.md.
EOF
;;
    samplesheet) cat <<'EOF'
chromaflow samplesheet [--type T] [--samples F] [--pair-pattern P]
                       [--single-end] [--list-fastq] [--strict]

  Builds the nf-core samplesheet from your sample record file.
  See `chromaflow help samples` for the file format and pairing rules.

  --list-fastq   Print detected pairs and exit (nothing written)
  --strict       Exit non-zero if validation finds problems
EOF
;;
    run) cat <<'EOF'
chromaflow run [--type chipseq|atacseq|dapseq] [--profile PROFILE]
  Env: ChromaFlow_runtime (offline) or nfcore_env (online Conda)

  Runs nf-core/chipseq (chipseq, dapseq) or nf-core/atacseq with -resume.
  Requires: samplesheet, genome FASTA, GTF and MACS_GSIZE or measured read length.

  OFFLINE MODE
    --offline --bundle DIR      Use a prepared, matching container bundle
    --profile apptainer         Or singularity; full offline Conda is not supported
    --dry-run                   Validate and print command; no analysis tasks
    --yes                       Skip the confirmation prompt
    Config equivalents: NFCORE_OFFLINE=true, NFCORE_OFFLINE_BUNDLE=/path/to/bundle
    Nextflow debug logs are saved directly to OUTPUT_DIR/logs/nextflow_*.log.
    See chromaflow help offline.

  EXECUTION MODE — controlled by NFCORE_PROFILE
    The profile sets nextflow's -profile AND selects the execution path:

      singularity  (default) container mode, images pulled automatically
      apptainer              same, for Apptainer 1.x installations
      docker / podman        container mode with those engines
      conda / mamba          conda mode: nextflow builds one env per process
                             (slow first run, ~30-60 min; for restricted
                              networks where image pulls fail)

    Resolution order: --profile > NFCORE_PROFILE in the config > singularity
    Combinations work too, e.g. NFCORE_PROFILE="singularity,test".

      chromaflow run --profile conda        # one-off conda run
      NFCORE_PROFILE=conda                  # or make it the project default

    Image cache: $OUTPUT_DIR/.singularity_cache (NXF_SINGULARITY_CACHEDIR,
    also exported as NXF_APPTAINER_CACHEDIR).
    Conda env cache: $OUTPUT_DIR/.conda_envs (NXF_CONDA_CACHEDIR).
    Conda mode rewrites channels to a Tsinghua mirror by default; set
    NFCORE_CONDA_MIRROR=none to keep your own channel configuration.

  PASSING EXTRA ARGUMENTS
    --skip-spp            skip phantompeakqualtools (also SKIP_SPP=true in
                          the config). That module is the most fragile part of
                          the pipeline under conda and contributes only a QC
                          metric (NSC/RSC) — skipping it does not affect peaks.
    --skip-plot-profile   skip the deepTools profile/heatmap QC plots
    --skip-multiqc        skip the aggregate MultiQC report (runs last;
                          skipping it costs no analysis output)
    --no-conda-patch      do not inject the dependency-pin config (below)
    -- <args...>          everything after `--` is passed through verbatim,
                          e.g. `chromaflow run -- --skip_preseq -resume`

  CONDA MODE FIXES (written to .chromaflow/nfcore_conda_patch.config, -c)
    1. createTimeout (default 3 h). Nextflow allows only 20 min to build an
       environment. A big R/Bioconductor env (phantompeakqualtools pulls
       r-base + bioconductor-* + boost, >1 GB) routinely exceeds that; conda
       is then killed and the run fails with
         Failed to create Conda environment ... status : 120
       Worse, the half-built env directory is left behind, and the next run
       treats it as ready — which surfaces later as a bare "command not
       found" (e.g. run_spp.R). ChromaFlow also deletes such partial
       environments before launching (they lack conda-meta/history).
       Override with NFCORE_CONDA_TIMEOUT.
    2. useMamba. Nextflow only honours `conda.useMamba` in a config file —
       there is no NXF_CONDA_USE_MAMBA environment variable, so exporting one
       (as older versions of this script did) silently left plain conda in
       charge. conda's downloader gives up on a dropped connection
       (`curl: (18) end of response with N bytes missing`); mamba retries and
       parallelises. Enabled automatically when mamba (or micromamba) exists.
    3. Dependency pins for processes whose deps nf-core leaves unpinned:
       DEEPTOOLS_* gets numpy=1.26 (deeptools 3.5.x calls np.NAN, removed in
       numpy 2.0); MULTIQC gets python=3.12 + python-kaleido=0.2.1 (kaleido
       >=1.0 dropped its bundled chromium and demands a system Google Chrome,
       so MultiQC dies exporting static plots). Tool versions are read from
       the pipeline assets so the pins track nf-core.
    --no-conda-patch drops pin 3; timeout and useMamba always apply.

    The generated .condarc also sets download resilience
    (remote_max_retries=5, remote_read_timeout_secs=600) and supports
    switching mirror: NFCORE_CONDA_MIRROR=tuna|bfsu|nju|sjtu|none, or any
    base URL of your own mirror. Switch when one mirror throttles.
    Container mode needs none of this: image dependencies are fixed at build
    time.

  ABOUT THE nf-core CHANNEL WARNING
    "There is a problem with your Conda configuration" comes from
    `conda config --show channels`, which merges the system, user (~/.condarc)
    and CONDARC files. It can therefore appear even when the project .condarc
    is correct, because your ~/.condarc lists bioconda first. Module
    environments declare their own channels, so this warning is usually noise
    — check for "status : 120" before chasing it.

  CONDA CHANNEL ORDER
    Conda mode writes a .condarc with conda-forge > bioconda > defaults and
    strict priority. That order is required by nf-core: with bioconda first,
    the solver mixes incompatible builds and a process fails later with a
    "command not found" (classically run_spp.R). A stale .condarc from an
    older ChromaFlow is detected and rebuilt automatically.

  Force a clean rerun:
    rm -rf "$OUTPUT_DIR"/work "$OUTPUT_DIR"/.nextflow*
EOF
;;
    idr) cat <<'EOF'
chromaflow idr [--min-pairs N] [--exclude-rep REP1[,REP3]]
  Env: ChromaFlow_idr

  Pairwise IDR for a single biological condition with >=2 replicates; columns:
    chr start end name signalValue strand globalIDR n_pairs

  signalValue: maximum overlapping MACS signal from retained replicates.
  globalIDR: -log10 scale; merged regions retain the maximum passing value.
  Failed or empty pair outputs abort. Downstream --peakset idr selects 0.01;
  it does not fall back to consensus. See docs/DOWNSTREAM.md.

  --min-pairs N     Keep peaks supported by at least N replicate pairs.
                    Preferred over dropping a replicate: no data is lost.
  --exclude-rep R   Rerun without a suspicious replicate to measure impact.

  Read the pairwise table printed in the log before judging replicates.
EOF
;;
    motif) cat <<'EOF'
chromaflow motif [--peakset SPEC]
  Env: ChromaFlow_motif

  MEME-ChIP on summit +/-100bp for narrowPeak, or midpoint +/-100bp for BED:
  de novo motifs, Tomtom vs MEME_DB, CentriMo, and FIMO in extracted sequences.
  MEME_SEARCHSIZE=0 uses all sequences (runtime grows ~quadratically).

  Use the SAME --peakset for motif, annotate, candidates and lookup, or
  FIMO hits will not line up with the annotated peaks.
EOF
;;
    annotate|annotate-py) cat <<'EOF'
chromaflow annotate      [--peakset SPEC]   Env: ChromaFlow_R
chromaflow annotate-py   [--peakset SPEC]   Env: ChromaFlow_ref

  Assigns each peak to one gene (nearest TSS) and classifies its feature.
  PROMOTER_UP/PROMOTER_DOWN map to ChIPseeker tssRegion=c(-UP,+DOWN).
  Priority: Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic.

  annotate-py is a simpler R-free alternative (bedtools + Python); fields
  and plots differ from ChIPseeker. Use annotate for the full R workflow.

  Outputs: peak_annotation_full.tsv, genes.bed, background_universe.txt, plots.
  R annotate also validates/records annotation_policy.json and audits every
  gene's full configured promoter against all selected peaks. Regional overlaps
  are in promoter_peak_overlaps.tsv; assignments/candidate rules are unchanged.
EOF
;;
    candidates) cat <<'EOF'
chromaflow candidates [--features "Promoter,1st Intron"] [--pval 1e-4]
                      [--max-tss-dist 2000] [--min-signal N] [--plot]
                      [--motif MEME-4,MEME-7] [--qval 0.05]
                      [--motif-evalue 0.05] [--centrimo-evalue 0.05]
  Env: ChromaFlow_downstream

  candidate = peak at an accepted feature (TARGET_FEATURES) AND containing
  a FIMO motif hit.

  Feature values: Promoter | 5' UTR | 3' UTR | 1st Exon | Other Exon |
                  1st Intron | Other Intron | Downstream | Distal Intergenic

  Outputs: candidate_targets.tsv, motif_<name>_hits.tsv/_genes.txt,
           motif_summary.tsv, candidate_genes.txt, go/
  Shared motif filters: --motif (exact motif label/raw ID within this run),
  --pval (FIMO site p), --qval (FIMO site q), --motif-evalue (discovery E),
  --centrimo-evalue (positional enrichment E). Optional filters default off.
  Missing requested statistics or an explicitly unscanned model cause an error.
  No rescanning; a relaxed threshold cannot recover hits absent from FIMO output.
  CentriMo local/central mode and report threshold are recorded separately.
  --plot draws selected candidates. Or use tracks --selection OUTPUT_DIR.
  selected_genes.txt, selected_evidence.tsv, motif_catalog.tsv and selection.json
  preserve filter settings and input hashes. high_confidence_genes.txt is a
  legacy alias for computational candidates, not experimental validation.
EOF
;;
    lookup) cat <<'EOF'
chromaflow lookup --gene-list FILE [--plot] [--min-signal N] [--level L] [--yes]
  Env: ChromaFlow_downstream

  Independent evidence breakdown for genes you care about (does not depend
  on the candidate table). query_binding.tsv tier column:
    candidate  = right position + motif      in_target = position, no motif
    motif_only = motif, wrong position       peak_only = peak only
    none       = no peak

  Exact IDs only: no version/suffix removal, aliases or regex matching.
  Shows query/annotation ID examples and asks y/N before analysis. --yes is
  explicit unattended acceptance, not permission to change IDs. Zero matches fail.
  id_check.tsv/id_not_found.txt distinguish unmatched IDs from no assigned peak.
  Reports one best assigned peak per query gene; separately reports ALL regional
  promoter overlaps, including peaks assigned to other genes. No reassignment.
  Uses the same TARGET_FEATURES/TARGET_MAX_TSS_DIST policy as candidates.
  --level all|in_target|candidate controls both the table and --plot selection.
  Default all plots all exactly matched query genes; unmatched IDs stay in the
  ID report. candidate plots only candidates. No further ID transformation.
  Shared filters: --motif, --pval, --qval, --motif-evalue, --centrimo-evalue.
  Writes selected_genes.txt, selected_evidence.tsv and selection.json.
  Filter every peak first, then select the best tier/signal peak per gene.
EOF
;;
    tracks) cat <<'EOF'
chromaflow tracks [--genes G1,G2 | --gene-list FILE | --selection DIR] [--outdir DIR]
                  [--motifs] [--flank BP]
  Env: ChromaFlow_downstream

  pyGenomeTracks browser shots per gene, with a log2(ChIP/Input) difference
  track and signalValue in the peak label.
  Gene lists auto-detect a gene column (gene_id/gene/gene_name/id header),
  so query_binding.tsv can be passed directly. candidate_targets.tsv contains
  both passing and failing peaks: prefer candidates/selected_genes.txt or
  --selection candidates_output for the filtered set.
  --selection accepts a result directory or selection.json from candidates/lookup.
  It uses stored genes and filtered motif hits; verifies peakset, promoter settings
  and file hashes. Use the same --peakset/config. Do not combine --selection with
  --genes, --gene-list or --motifs. Context peaks/neighbor genes remain complete.
  --motifs overlays this peakset's existing FIMO hits and extracted input
  windows. It does not rescan untested DNA or filter tracks by TARGET_FEATURES.
  --flank sets context on both sides of the whole gene (default: TRACK_FLANK).
  Context is expanded if needed to include the full configured promoter.
  Shows every selected peak and neighboring gene in the plot region, without
  filtering by the query gene's assigned/candidate peaks. Exact IDs only.
  Exports promoter_regions.bed, promoter_peak_overlaps.tsv and a per-gene summary.
  plot_status.tsv records plotted, missing_coordinates and render_failed genes.

  Gene/motif track heights scale with regional transcript/hit counts, retaining
  all isoforms. Config TRACK_GENE_HEIGHT / TRACK_MOTIF_HEIGHT: auto or cm.
  Styling lives in lib/steps/plot_tracks.sh (STYLE block at the top).
EOF
;;
    go-prep|go-plot) cat <<'EOF'
chromaflow go-prep [--set promoter|dist|all]     Env: ChromaFlow_ref
chromaflow go-plot --input RESULT [--top N]      Env: ChromaFlow_ref

  go-prep exports foreground/background gene lists for agriGO, g:Profiler
  and similar tools; go-plot turns the returned table into a bubble chart.
EOF
;;
    doctor) cat <<'EOF'
chromaflow doctor [-c config] [--type T]

  Checks software, per-step conda tools, reference files and which outputs
  already exist. Optional features (organelles, AGAT, GO annotation) are
  reported as [SKIP]/[INFO], never as failures.
  Exit code: 0 when there is no FAIL.
EOF
;;
    peakset) cat <<'EOF'
--peakset SPEC   Which peak set downstream steps consume

  idr           (default) reproducible peaks at IDR < 0.01
  idr:0.05      other IDR thresholds -> output dirs suffixed _idr0.05
  consensus     skip IDR, use the nf-core consensus set. nf-core writes
                column 5 of consensus_peaks.bed as a hardcoded "0", so
                ChromaFlow back-fills signalValue by intersecting with the
                replicate narrowPeak files (bedtools, or a bundled Python
                fallback). Without that step every peak_score would be 0.
  <file>        any BED you provide -> output dirs suffixed _custom

  Each peakset writes to its own directory, so runs can be compared side by
  side. Keep the peakset constant across motif -> annotate -> candidates.
EOF
;;
    envs)    help_envs ;;
    config)  help_config ;;
    samples) help_samples ;;
    *) log_error "No help topic: $1"
       echo "Topics: init reference normalize-gff samplesheet run idr motif annotate"
       echo "        candidates lookup tracks go-prep go-plot doctor peakset"
       echo "        config samples envs"
       return 1 ;;
    esac
}
