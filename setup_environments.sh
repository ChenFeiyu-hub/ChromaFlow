#!/usr/bin/env bash
# Install only the environments required for the selected execution mode.
set -euo pipefail
export LANGUAGE=en LC_MESSAGES=C LC_NUMERIC=C LC_TIME=C
unset LC_ALL
MODE=offline
COMPONENTS=all
DRY_RUN="${CHROMAFLOW_DRY_RUN:-false}"
SKIP_R="${SKIP_R:-false}"
WITH_AGAT="${WITH_AGAT:-false}"
SOLVER=""
HERE="$(cd "$(dirname "$0")" && pwd)"
usage() {
    cat <<'EOF'
Usage: chromaflow setup [options]
       bash setup_environments.sh [options]

  --mode offline|prepare|conda    Default: offline (singularity/apptainer aliases)
  --components runtime|downstream|all   Default: all
  --skip-r                       Use annotate-py instead of the R annotation step
  --with-agat                    Also install the optional AGAT environment
  --solver mamba|micromamba|conda Select an installed solver
  --dry-run                      Print the installation plan; change nothing

offline: ChromaFlow_runtime (Python + Java). Nextflow, plugins and workflow
         tools come from the transferred container bundle; no nf-core CLI.
prepare: nfcore_env (Python + Java + nf-core/tools) for a connected Linux/WSL
         machine. ChromaFlow downloads its pinned Nextflow into each bundle.
conda:   nfcore_env (Python + Java + Nextflow) for the online Conda profile.
         Nextflow builds workflow tool environments on demand.

Downstream environments are independent of these three modes. Existing named
environments are kept. This installer does not silently update or delete them.
Package installation needs network access or a configured local package mirror.
The offline mode describes workflow execution, not an offline package installer.
Apptainer/Singularity must be installed separately by your system administrator.
EOF
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode|--components|--solver)
            [[ $# -ge 2 ]] || { echo "[ERROR] $1 requires a value" >&2; exit 2; }
            case "$1" in --mode) MODE="$2" ;; --components) COMPONENTS="$2" ;; --solver) SOLVER="$2" ;; esac
            shift 2 ;;
        --skip-r) SKIP_R=true; shift ;;
        --with-agat) WITH_AGAT=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$MODE" in singularity|apptainer) MODE=offline ;; offline|prepare|conda) ;; *) echo '[ERROR] Invalid --mode' >&2; exit 2 ;; esac
case "$COMPONENTS" in runtime|downstream|all) ;; *) echo '[ERROR] Invalid --components' >&2; exit 2 ;; esac
case "$SOLVER" in ''|mamba|micromamba|conda) ;; *) echo '[ERROR] Invalid --solver' >&2; exit 2 ;; esac
if [[ -z "$SOLVER" ]]; then
    for candidate in mamba micromamba conda; do
        if command -v "$candidate" >/dev/null; then SOLVER="$candidate"; break; fi
    done
fi
if [[ "$DRY_RUN" == true ]]; then
    SOLVER="${SOLVER:-mamba}"
elif [[ -z "$SOLVER" ]] || ! command -v "$SOLVER" >/dev/null; then
    echo '[ERROR] Install Miniforge/Conda first, or select an available solver.' >&2
    exit 1
fi
echo "[PLAN] mode=$MODE components=$COMPONENTS solver=$SOLVER dry_run=$DRY_RUN"
echo '[INFO] Conda package installation may require network access.'
print_command() { printf '[CMD] '; printf '%q ' "$@"; printf '\n'; }
env_exists() { "$SOLVER" env list | awk -v name="$1" '$1==name {found=1} END {exit !found}'; }
create_env() {
    local name="$1"; shift
    if [[ "$DRY_RUN" != true ]] && env_exists "$name"; then
        echo "[SKIP] Existing environment: $name (kept unchanged)"
        return 0
    fi
    local args=(create -y -n "$name" --override-channels -c conda-forge -c bioconda --strict-channel-priority "$@")
    print_command "$SOLVER" "${args[@]}"
    [[ "$DRY_RUN" == true ]] || "$SOLVER" "${args[@]}"
}
if [[ "$COMPONENTS" != downstream ]]; then
    case "$MODE" in
        offline) create_env ChromaFlow_runtime python=3.11 openjdk=17 ;;
        prepare) create_env nfcore_env python=3.11 openjdk=17 nf-core=3.3.2 ;;
        conda) create_env nfcore_env python=3.11 openjdk=17 nextflow=25.10.4 ;;
    esac
    if [[ "$MODE" != conda ]]; then
        if command -v apptainer >/dev/null || command -v singularity >/dev/null; then
            echo '[OK] Apptainer/Singularity is available on PATH.'
        else
            echo '[INFO] Install/load Apptainer or Singularity before preparing/running container bundles.'
        fi
    fi
fi
if [[ "$COMPONENTS" != runtime ]]; then
    create_env ChromaFlow_ref python=3.11 samtools=1.21 gffread bedtools=2.31 khmer pigz matplotlib-base
    create_env ChromaFlow_downstream python=3.12 deeptools=3.5.6 'numpy<2' pygenometracks bedtools=2.31
    create_env ChromaFlow_idr python=3.10 idr=2.0.4.2 numpy=1.23.5 bedtools=2.31
    create_env ChromaFlow_motif python=3.12 meme=5.5.7 bedtools=2.31 samtools=1.21
    if [[ "$WITH_AGAT" == true ]]; then create_env ChromaFlow_agat agat; fi
    if [[ "$SKIP_R" == true ]]; then
        echo '[SKIP] ChromaFlow_R: use annotate-py with ChromaFlow_ref.'
    elif [[ "$DRY_RUN" != true ]] && env_exists ChromaFlow_R; then
        echo '[SKIP] Existing environment: ChromaFlow_R (kept unchanged)'
    else
        print_command bash "$HERE/env/install_chipseeker.sh" --method conda --solver "$SOLVER"
        [[ "$DRY_RUN" == true ]] || bash "$HERE/env/install_chipseeker.sh" --method conda --solver "$SOLVER"
    fi
fi
echo '[DONE] Environment plan completed. See docs/OFFLINE.md and docs/DOWNSTREAM.md.'
