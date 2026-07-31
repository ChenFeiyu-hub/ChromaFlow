#!/usr/bin/env bash
# ============================================================================
# lib/help.sh — `chromaflow help [command]` 的帮助文本
# ============================================================================
# 顶层帮助保持一屏可读；每个子命令的细节用 `chromaflow help <command>` 展开。
# ============================================================================

CHROMAFLOW_VERSION="${CHROMAFLOW_VERSION:-1.1.0}"

help_main() {
cat <<'EOF'
ChromaFlow — ChIP-seq / DAP-seq / ATAC-seq downstream analysis toolkit
Peaks → reproducibility → motifs → annotation → candidate target genes.

USAGE
  chromaflow <command> [options]
  chromaflow help <command>            Detailed help for one command
  chromaflow doctor                    Check tools, references and outputs

SETUP
  init <dir> [--type T]                Create a project skeleton + config
  doctor [-c cfg] [--type T]           Pre-flight self-check
  version                              Print version

REFERENCE (optional, run once per species)
  reference                            00  Genome / GTF / chrom.sizes
  normalize-gff                        00b Normalize GFF3 (AGAT or built-in)

MAIN PIPELINE
  samplesheet [--type T]               01  Build the nf-core samplesheet
  run [--type T]                       02  Run nf-core/chipseq | atacseq

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

ENVIRONMENTS (each step needs its own conda env; see `chromaflow help envs`)
  mamba activate ChromaFlow_idr && chromaflow idr

DOCS
  README.md (English) · README.zh-CN.md (中文) · ARCHITECTURE.md
EOF
}

help_envs() {
cat <<'EOF'
ENVIRONMENTS
  nfcore_env             nextflow, nf-core CLI            -> run
  ChromaFlow_ref         samtools, gffread, bedtools,     -> reference, samplesheet,
                         python, matplotlib                  annotate-py, go-*, spacing
  ChromaFlow_idr         idr (python 3.10, numpy<1.24)    -> idr
  ChromaFlow_motif       MEME Suite, bedtools, samtools   -> motif
  ChromaFlow_R           ChIPseeker, clusterProfiler      -> annotate
  ChromaFlow_downstream  deepTools, pyGenomeTracks        -> candidates, lookup,
                                                             tracks, signal
  ChromaFlow_agat        AGAT (optional, separate env)    -> normalize-gff --engine agat

  Install:  bash setup_environments.sh              # core envs
            WITH_AGAT=true bash setup_environments.sh   # + AGAT env
            SKIP_R=true    bash setup_environments.sh   # no R (use annotate-py)

  AGAT is optional and installed in its own environment on purpose: its Perl
  dependency tree often conflicts with other packages. Without AGAT,
  `chromaflow normalize-gff` falls back to a built-in Python normalizer.
EOF
}

help_config() {
cat <<'EOF'
PROJECT CONFIG (chromaflow.config)

  Declare only what you already have; everything the pipeline produces is
  discovered automatically.

  Required
    SPECIES_SHORT="Cre"                 Short species tag
    PROJECT_DIR="/data/projects/X"      Anchor for inputs/ and outputs/

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

  Rebuilds a complete GFF3 so ChIPseeker/TxDb reports UTR/exon/intron
  correctly (symptom of skipping it: 0% UTR in the annotation pie chart).

  --engine auto      Use AGAT if available, else the built-in normalizer
          agat       Force AGAT (most thorough)
          builtin    Pure Python, zero extra dependencies:
                     adds exons from CDS, infers 5'/3' UTRs, synthesizes
                     missing gene features, sorts output
          gffread    Tidy hierarchy with gffread

  Skipping this step does not break annotation: step 05 falls back to the
  plain GFF3 and only warns that exon-level features may be missing.

  The output is registered as the current GENOME_GFF3, so annotate picks it
  up without editing the config.
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
  Env: nfcore_env

  Runs nf-core/chipseq (chipseq, dapseq) or nf-core/atacseq with -resume.
  Requires: samplesheet, genome FASTA, GTF and MACS_GSIZE.

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

  Force a clean rerun:
    rm -rf "$OUTPUT_DIR"/work "$OUTPUT_DIR"/.nextflow*
EOF
;;
    idr) cat <<'EOF'
chromaflow idr [--min-pairs N] [--exclude-rep REP1[,REP3]]
  Env: ChromaFlow_idr

  Pairwise IDR across replicates; output columns:
    chr start end name signalValue strand globalIDR n_pairs

  --min-pairs N     Keep peaks supported by at least N replicate pairs.
                    Preferred over dropping a replicate: no data is lost.
  --exclude-rep R   Rerun without a suspicious replicate to measure impact.

  Read the pairwise table printed in the log before judging replicates.
EOF
;;
    motif) cat <<'EOF'
chromaflow motif [--peakset SPEC]
  Env: ChromaFlow_motif

  MEME-ChIP on summit +/-100bp: de novo motifs, Tomtom vs MEME_DB,
  CentriMo central enrichment, FIMO genome-wide scan.
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

  annotate-py is the R-free alternative (bedtools + Python) with the same
  output contract, so steps 06-09 work either way.

  Outputs: peak_annotation_full.tsv, genes.bed, background_universe.txt, plots.
EOF
;;
    candidates) cat <<'EOF'
chromaflow candidates [--features "Promoter,1st Intron"] [--pval 1e-4]
  Env: ChromaFlow_downstream

  candidate = peak at an accepted feature (TARGET_FEATURES) AND containing
  a FIMO motif hit.

  Feature values: Promoter | 5' UTR | 3' UTR | 1st Exon | Other Exon |
                  1st Intron | Other Intron | Downstream | Distal Intergenic

  Outputs: candidate_targets.tsv, motif_<name>_hits.tsv/_genes.txt,
           motif_summary.tsv, high_confidence_genes.txt, go/
EOF
;;
    lookup) cat <<'EOF'
chromaflow lookup --gene-list FILE [--plot] [--min-signal N] [--level L]
  Env: ChromaFlow_downstream

  Independent evidence breakdown for genes you care about (does not depend
  on the candidate table). query_binding.tsv tier column:
    candidate  = right position + motif      in_target = position, no motif
    motif_only = motif, wrong position       peak_only = peak only
    none       = no peak

  Also validates IDs, separating "gene has no peak" from "ID not in the
  annotation" (written to id_not_found.txt).
EOF
;;
    tracks) cat <<'EOF'
chromaflow tracks [--genes G1,G2 | --gene-list FILE] [--outdir DIR]
  Env: ChromaFlow_downstream

  pyGenomeTracks browser shots per gene, with a log2(ChIP/Input) difference
  track and signalValue in the peak label.
  Gene lists auto-detect a gene column (gene_id/gene/gene_name/id header),
  so candidate_targets.tsv or query_binding.tsv can be passed directly.

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
  consensus     skip IDR, use the nf-core consensus set (signalValue is
                back-filled from the replicate narrowPeak files)
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
