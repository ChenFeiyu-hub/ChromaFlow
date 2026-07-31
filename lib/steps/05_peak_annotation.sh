#!/usr/bin/env bash
# ============================================================================
# 05_peak_annotation.sh — Peak 注释（ChIPseeker/R 包装器）
# ============================================================================
# 环境: mamba activate ChromaFlow_R   (ChIPseeker + clusterProfiler)
# 默认调用 05_peak_annotation.R；无 R 环境时可改用 05_peak_annotation_py.sh。
# 输出契约与旧版一致，06/07/08/09 无需改动。
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
init_pipeline_dirs "${PIPELINE}"
open_log "05_peak_annotation" "${PIPELINE}"
log_section "Step 05: Peak Annotation — ChIPseeker (${PIPELINE})"

check_tool Rscript || { echo "  需 R：mamba activate ChromaFlow_R（或用 chromaflow annotate-py）"; exit 1; }
if [[ ! -s "${GENOME_GFF3}" ]]; then
    log_error "缺 GFF3: ${GENOME_GFF3}（ChIPseeker 建 TxDb 需要）"
    echo "  处理办法（任选其一）："
    echo "    · 在 config 设 GENOME_GFF3_SOURCE=<你的 .gff3/.gff3.gz>，再跑 chromaflow reference"
    echo "    · 已有 gff3 时直接放进 ${REFERENCE_DIR}/ ，会被自动识别"
    echo "    · 或显式指定: GENOME_GFF3=... 写进 config"
    exit 1
fi
# 未规范化的 GFF（无 exon 行）只警告，不阻断——normalize-gff 是可选步骤
if ! grep -qP '\texon\t' "${GENOME_GFF3}" 2>/dev/null; then
    log_warn "该 GFF3 无 exon 行，UTR/exon/intron 判定会失真"
    echo "         建议先跑: chromaflow normalize-gff   （无 AGAT 时会自动用内置引擎）"
fi

# peakset：决定用哪套 peak + 输出目录后缀（--peakset: idr|consensus|<路径>）
resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
INPUT_PEAKS="${PEAKSET_FILE}"; SRC="${CHROMAFLOW_PEAKSET:-idr}"
[[ -s "${INPUT_PEAKS}" ]] || die "peakset peak 文件为空: ${INPUT_PEAKS}"
ANNO_OUT="$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}")"
mkdir -p "${ANNO_OUT}"

echo "  输入 peaks: ${INPUT_PEAKS} (${SRC})"
echo "  GFF3:       ${GENOME_GFF3}"
echo "  启动子窗口: -${PROMOTER_UP}/+${PROMOTER_DOWN}   靶基因距离: ${TARGET_DIST_KB}kb"
echo "  KEGG 物种:  ${KEGG_ORGANISM:-（未设，跳过 KEGG）}"
echo ""

Rscript "${CHROMAFLOW_HOME}/lib/steps/05_peak_annotation.R" \
    "${GENOME_GFF3}" "${INPUT_PEAKS}" "${ANNO_OUT}" \
    "${PROMOTER_UP}" "${PROMOTER_DOWN}" "${ANNOTATION_INFO:-}" \
    "${TARGET_DIST_KB}" "${PIPELINE}" "${KEGG_ORGANISM:-}"

echo ""
echo "============================================"
echo " 注释完成 (${PIPELINE})"
echo "============================================"
echo "  输出目录: ${ANNO_OUT}"
echo "  下一步: chromaflow candidates --type ${PIPELINE}   (整合 04 motif → 候选表)"
echo "[DONE] $(date)"
