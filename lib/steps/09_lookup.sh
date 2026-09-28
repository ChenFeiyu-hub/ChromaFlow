#!/usr/bin/env bash
# Exact-ID query, annotation evidence and complete promoter overlap audit.
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs
PIPELINE="${1:-chipseq}"; shift || true
source "${CHROMAFLOW_HOME}/lib/selection.sh"
selection_defaults
FPVAL="${FIMO_PVAL}";
GENE_LIST=""; DO_PLOT=false; MIN_SIGNAL=0; LEVEL=all
while [[ $# -gt 0 ]]; do
    case "$1" in
        --motif|--qval|--motif-evalue|--centrimo-evalue) selection_option "$1" "${2:-}"; shift 2 ;;
        --pval) FPVAL="$2"; shift 2 ;;
        --gene-list) GENE_LIST="$2"; shift 2 ;;
        --plot) DO_PLOT=true; shift ;;
        --min-signal) MIN_SIGNAL="$2"; shift 2 ;;
        --level) LEVEL="$2"; shift 2 ;;
        *) die "Unknown lookup option: $1" ;;
    esac
done
[[ -s "${GENE_LIST}" ]] || die '--gene-list must point to a query gene list'
init_pipeline_dirs "${PIPELINE}"
open_log 09_lookup "${PIPELINE}"
log_section "Step 09: Exact-ID gene lookup (${PIPELINE})"
resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
ANNO_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}"
MOTIF_DIR="${RESULTS_DIR}/downstream/${PIPELINE}/motif${PEAKSET_SUFFIX}"
OUT="$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/lookup${PEAKSET_SUFFIX}")"
[[ -s "${ANNO_OUT}/peak_annotation_full.tsv" && -s "${ANNO_OUT}/genes.bed" ]] || die 'Run annotate first: full annotation and genes.bed are required'
mkdir -p "${OUT}"
rm -f -- "${OUT}/selection.json"
mapfile -t FIMO_TSVS < <(find "${MOTIF_DIR}" -path '*/fimo_out_*/fimo.tsv' -type f 2>/dev/null | sort)
echo "  Annotation: ${ANNO_OUT}/peak_annotation_full.tsv"
echo "  Gene list: ${GENE_LIST}"
echo "  TARGET_FEATURES: ${TARGET_FEATURES}; FIMO p-value threshold: ${FPVAL}"
echo "  Promoter: -${PROMOTER_UP}/+${PROMOTER_DOWN}; upstream rescue: ${TARGET_MAX_TSS_DIST:-0}"
[[ ${#FIMO_TSVS[@]} -gt 0 ]] || echo '  FIMO absent: location-only evidence; unscanned regions are not motif negatives.'
MC_DIR=$(find "${MOTIF_DIR}" -type d -name meme-chip_output -print -quit 2>/dev/null || true)
export CFSEL_QUERY="${GENE_LIST}"
export CFSEL_ROOT="${MC_DIR:-${MOTIF_DIR}}" CFSEL_PIPELINE="${PIPELINE}"
if python3 "${CHROMAFLOW_HOME}/lib/steps/support/lookup.py" \
    "${ANNO_OUT}/peak_annotation_full.tsv" "${GENE_LIST}" "${OUT}" \
    "${TARGET_FEATURES}" "${FPVAL}" "${MIN_SIGNAL}" "${LEVEL}" \
    "${PROMOTER_UP}" "${PROMOTER_DOWN}" "${TARGET_MAX_TSS_DIST:-0}" "${CHROMAFLOW_YES:-false}" \
    "${PEAKSET_FILE}" "${CHROM_SIZES}" "${FIMO_TSVS[@]}"; then
    :
else
    rc=$?
    [[ "$(cat "${OUT}/lookup.status" 2>/dev/null)" == cancelled ]] || printf 'failed\n' > "${OUT}/lookup.status"
    exit "${rc}"
fi
if [[ "${DO_PLOT}" == true ]]; then
    if [[ -s "${OUT}/selected_genes.txt" ]]; then
        bash "${CHROMAFLOW_HOME}/lib/steps/07_gene_tracks.sh" "${PIPELINE}" \
            --selection "${OUT}" --outdir "${OUT}/tracks"
    else
        echo '[SKIP] No genes passed the requested level; no tracks requested.'
    fi
fi
echo "[REPORT] ${OUT}/query_binding.tsv"
echo "[REPORT] ${OUT}/promoter_peak_overlaps.tsv"
echo "[DONE] $(date)"
