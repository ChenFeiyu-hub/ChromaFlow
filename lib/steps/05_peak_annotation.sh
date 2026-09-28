#!/usr/bin/env bash
# ============================================================================
# 05_peak_annotation.sh — Peak 注释（ChIPseeker/R 包装器）
# ============================================================================
# 环境: conda activate ChromaFlow_R   (ChIPseeker + clusterProfiler)
# 默认调用 05_peak_annotation.R；无 R 环境时可改用 05_peak_annotation_py.sh。
# 输出契约与旧版一致，06/07/08/09 无需改动。
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
init_pipeline_dirs "${PIPELINE}"
open_log "05_peak_annotation" "${PIPELINE}"
log_section "Step 05: Peak Annotation — ChIPseeker (${PIPELINE})"
[[ "${PROMOTER_UP}" =~ ^[0-9]+$ && "${PROMOTER_DOWN}" =~ ^[0-9]+$ ]] || die 'PROMOTER_UP and PROMOTER_DOWN must be nonnegative integers'

check_tool Rscript || { echo "  Requires R: conda activate ChromaFlow_R (or use chromaflow annotate-py)"; exit 1; }
if [[ ! -s "${GENOME_GFF3}" ]]; then
    log_error "Missing GFF3: ${GENOME_GFF3}(required by ChIPseeker to build TxDb)"
    echo "  Options:"
    echo "    · Set GENOME_GFF3_SOURCE=<your .gff3/.gff3.gz> in config, then run chromaflow reference"
    echo "    · Place an existing GFF3 in ${REFERENCE_DIR}/ for automatic detection"
    echo "    · Or explicitly set GENOME_GFF3=... in config"
    exit 1
fi
# 未规范化的 GFF（无 exon 行）只警告，不阻断——normalize-gff 是可选步骤
if ! grep -qP '\texon\t' "${GENOME_GFF3}" 2>/dev/null; then
    log_warn "This GFF3 lacks exon records; UTR/exon/intron classification may be unreliable"
    echo "         Run first: chromaflow normalize-gff   (the builtin engine does not require AGAT)"
fi

# peakset：决定用哪套 peak + 输出目录后缀（--peakset: idr|consensus|<路径>）
resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
INPUT_PEAKS="${PEAKSET_FILE}"; SRC="${CHROMAFLOW_PEAKSET:-idr}"
[[ -s "${INPUT_PEAKS}" ]] || die "Peak-set file is empty: ${INPUT_PEAKS}"
ANNO_OUT="$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}")"
mkdir -p "${ANNO_OUT}"

echo "  Input peaks: ${INPUT_PEAKS} (${SRC})"
echo "  GFF3:       ${GENOME_GFF3}"
echo "  Promoter window: -${PROMOTER_UP}/+${PROMOTER_DOWN}   Target-gene distance: ${TARGET_DIST_KB}kb"
echo "  MAKE_PLOTS: ${MAKE_PLOTS}; assignment: gene-level nearest TSS; priority: Promoter > UTR > Exon > Intron"
echo "  KEGG Species:  ${KEGG_ORGANISM:-(unset; skipping KEGG)}"
echo ""

Rscript "${CHROMAFLOW_HOME}/lib/steps/05_peak_annotation.R" \
    "${GENOME_GFF3}" "${INPUT_PEAKS}" "${ANNO_OUT}" \
    "${PROMOTER_UP}" "${PROMOTER_DOWN}" "${ANNOTATION_INFO:-}" \
    "${TARGET_DIST_KB}" "${PIPELINE}" "${KEGG_ORGANISM:-}" "${MAKE_PLOTS}"

python3 "${CHROMAFLOW_HOME}/lib/steps/support/annotation_policy.py" validate \
    "${ANNO_OUT}/peak_annotation_full.tsv" --up "${PROMOTER_UP}" --down "${PROMOTER_DOWN}"
python3 "${CHROMAFLOW_HOME}/lib/steps/support/promoter_audit.py" \
    --genes-bed "${ANNO_OUT}/genes.bed" --peaks "${INPUT_PEAKS}" \
    --annotation "${ANNO_OUT}/peak_annotation_full.tsv" --outdir "${ANNO_OUT}" \
    --up "${PROMOTER_UP}" --down "${PROMOTER_DOWN}" --chrom-sizes "${CHROM_SIZES}"

echo ""
echo "============================================"
echo " Annotation completed (${PIPELINE})"
echo "============================================"
echo "  Output directory: ${ANNO_OUT}"
echo "  Next step: chromaflow candidates --type ${PIPELINE}   (integrate step 04 motifs into candidate tables)"
echo "[DONE] $(date)"
