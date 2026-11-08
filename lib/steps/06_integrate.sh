#!/usr/bin/env bash
# ============================================================================
# 06_integrate.sh — motif×注释候选表与 GO 输入导出
# ============================================================================
# 环境: conda activate ChromaFlow_downstream  (含 pyGenomeTracks + python)
# 用法: chromaflow candidates --type chipseq [--pval 1e-4]
#       [--features "Promoter,1st Intron,Distal Intergenic"]
#       [--max-tss-dist 2000]   TSS 上游该距离内也算合适位置（救回邻近基因体内的上游 peak）
# ----------------------------------------------------------------------------
# 直接读 04(MEME-ChIP) 的 fimo_out_*/fimo.tsv（不重跑 FIMO）与 centrimo.tsv。
#   1) 每个 motif 的 FIMO 命中(基因组坐标) → 映射到 05 注释的 peak/基因
#   2) 按「合适位置(TARGET_FEATURES) + 含 motif」筛出候选
#   3) 每个 motif 独立输出：命中明细(位置/链/p/距TSS) + 基因列表
#   4) 汇总表整合 CentriMo(中心富集)信息
#   5) 导出候选基因列表(agriGO/g:Profiler)；使用 tracks --gene-list 显式绘图
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true
source "${CHROMAFLOW_HOME}/lib/selection.sh"
selection_defaults
FPVAL_CLI=""; FEAT_CLI=""; MIN_SIGNAL=0; DO_PLOT=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --motif|--qval|--motif-evalue|--centrimo-evalue) selection_option "$1" "${2:-}"; shift 2 ;;
        --min-signal) MIN_SIGNAL="$2"; shift 2 ;;
        --plot) DO_PLOT=true; shift ;;
        --pval) FPVAL_CLI="$2"; shift 2 ;;
        --topn-tracks) die '--topn-tracks is retired; use tracks --gene-list to choose every gene to plot' ;;
        --features) FEAT_CLI="$2"; shift 2 ;;
        --max-tss-dist) MAXTSS_CLI="$2"; shift 2 ;;
        *) die "Unknown candidates option: $1" ;;
    esac
done
FPVAL="${FPVAL_CLI:-${FIMO_PVAL:-1e-4}}"
FEATURES="${FEAT_CLI:-${TARGET_FEATURES:-Promoter}}"
# 距离补救阈值（bp）：peak 在 geneId 的 TSS 上游该距离内也算合适位置，
# 用于救回「压在邻近基因体内、但确实位于目标基因上游」的 peak。0=关闭。
MAX_TSS_DIST="${MAXTSS_CLI:-${TARGET_MAX_TSS_DIST:-0}}"

init_pipeline_dirs "${PIPELINE}"
open_log "06_integrate" "${PIPELINE}"
log_section "Step 06: Candidate target integration (04 motifs x 05 annotation) — ${PIPELINE}"

resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
ANNO_OUT="$(downstream_dir "${PIPELINE}")/annotation${PEAKSET_SUFFIX}"
MOTIF_DIR="$(downstream_dir "${PIPELINE}")/motif${PEAKSET_SUFFIX}"
OUT="$(step_outdir "$(downstream_dir "${PIPELINE}")/candidates${PEAKSET_SUFFIX}")"
GENES_BED="${ANNO_OUT}/genes.bed"
mkdir -p "${OUT}"
rm -f -- "${OUT}/selection.json"

ANNO_TSV="${ANNO_OUT}/peak_annotation_full.tsv"
[[ -s "${ANNO_TSV}" ]] || die "Missing step 05 annotation table: ${ANNO_TSV}(run chromaflow annotate first)"
python3 "${CHROMAFLOW_HOME}/lib/steps/support/annotation_policy.py" check \
    "${ANNO_TSV}" --up "${PROMOTER_UP}" --down "${PROMOTER_DOWN}"

MC_DIR=""
if [[ -d "${MOTIF_DIR}" ]]; then
    MC_DIR=$(find "${MOTIF_DIR}" -type d -name "meme-chip_output" -print -quit)
fi
[[ -z "${MC_DIR}" ]] && MC_DIR="${MOTIF_DIR}"
mapfile -t FIMO_TSVS < <(find "${MC_DIR}" -path "*/fimo_out_*/fimo.tsv" 2>/dev/null | sort)
[[ ${#FIMO_TSVS[@]} -gt 0 ]] || die "No step 04 fimo_out_*/fimo.tsv found (run chromaflow motif first): ${MC_DIR}"
CENTRIMO="${MC_DIR}/centrimo_out/centrimo.tsv"

echo "  Step 05 annotation:      ${ANNO_TSV}"
echo "  FIMO files:    ${#FIMO_TSVS[@]} files"
echo "  CentriMo:     $([[ -s "${CENTRIMO}" ]] && echo present || echo absent)"
echo "  FIMO p-value threshold:    ${FPVAL}"
echo "  Accepted location features: ${FEATURES}"
echo ""

export CFSEL_ROOT="${MC_DIR}" CFSEL_PIPELINE="${PIPELINE}" CFSEL_PEAKS="${PEAKSET_FILE}"
export CFSEL_UP="${PROMOTER_UP}" CFSEL_DOWN="${PROMOTER_DOWN}" CFSEL_MIN_SIGNAL="${MIN_SIGNAL}"
# Invalidate a previous selection before this run can fail.
python3 "${CHROMAFLOW_HOME}/lib/steps/support/candidates.py" "${ANNO_TSV}" "${OUT}" "${FPVAL}" "${CENTRIMO}" "${GENES_BED}" "${FEATURES}" "${MAX_TSS_DIST}" "${FIMO_TSVS[@]}"


CAND="${OUT}/candidate_targets.tsv"
HICONF="${OUT}/high_confidence_genes.txt"
[[ -s "${CAND}" ]] || die "Candidate table generation failed"

# ─────────────── 自动: 导出 agriGO/g:Profiler 基因列表 ───────────────
echo ""
echo "[GO export] Generate external enrichment inputs"
GO_OUT="${OUT}/go"; mkdir -p "${GO_OUT}"
if [[ -f "${HICONF}" ]]; then
    sort -u "${HICONF}" > "${GO_OUT}/target_candidates.raw.txt"
    sed -E 's/\.[0-9]+$//' "${GO_OUT}/target_candidates.raw.txt" | sort -u > "${GO_OUT}/target_candidates.stripped.txt"
fi
if [[ -s "${ANNO_OUT}/background_universe.txt" ]]; then
    cp "${ANNO_OUT}/background_universe.txt" "${GO_OUT}/background.raw.txt"
    sed -E 's/\.[0-9]+$//' "${ANNO_OUT}/background_universe.txt" | sort -u > "${GO_OUT}/background.stripped.txt"
fi
cat > "${GO_OUT}/README_GO.txt" <<EOF
External GO enrichment: candidate genes (accepted location + motif)
Foreground: target_candidates.raw.txt
Background: background.raw.txt (use as the custom background)
Choose a service that supports your species and identifier namespace.
After obtaining enrichment results: chromaflow go-plot --input <results> --type ${PIPELINE}
EOF
echo "  → ${GO_OUT}/target_candidates.stripped.txt"

echo ""
echo "============================================"
echo " Integration completed (${PIPELINE})"
echo "============================================"
echo "  Candidate table (per peak):  ${CAND}"
echo "  Motif summary:      ${OUT}/motif_summary.tsv"
echo "  Per-motif details:   ${OUT}/motif_<name>_hits.tsv / _genes.txt"
echo "  Candidate genes:      ${HICONF}"
echo "  GO export:         ${GO_OUT}/"
echo "  (Generate gene tracks separately: chromaflow tracks --gene-list ${HICONF})"
echo "[DONE] $(date)"

if [[ "${DO_PLOT}" == true ]]; then
    if [[ -s "${OUT}/selected_genes.txt" ]]; then
        bash "${CHROMAFLOW_HOME}/lib/steps/07_gene_tracks.sh" "${PIPELINE}" --selection "${OUT}" --outdir "${OUT}/tracks"
    else
        echo '[SKIP] No genes passed selection; no tracks requested.'
    fi
fi
