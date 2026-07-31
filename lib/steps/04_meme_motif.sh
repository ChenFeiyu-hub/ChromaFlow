#!/usr/bin/env bash
# ============================================================================
# 04_meme_motif.sh — MEME-ChIP de novo motif 发现
# ============================================================================
# 支持: chipseq / atacseq / dapseq
# 环境: mamba activate ChromaFlow_motif
# 用法: chromaflow motif --type [chipseq|atacseq|dapseq]
# ============================================================================

set -euo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
init_pipeline_dirs "${PIPELINE}"

open_log "04_meme_motif" "${PIPELINE}"

log_section "Step 04: MEME-ChIP Motif Discovery (${PIPELINE})"

# ─────────────────────────────────────────────
# 检查工具
# ─────────────────────────────────────────────
MISSING=0
for tool in meme-chip bedtools; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "[ERROR] ${tool} 未安装"; MISSING=1
    fi
done
if [[ ${MISSING} -eq 1 ]]; then
    echo "  请: mamba activate ChromaFlow_motif"; exit 1
fi

# ─────────────────────────────────────────────
# 路径配置
# ─────────────────────────────────────────────
# peakset：决定用哪套 peak 做 motif + 输出目录后缀（--peakset: idr|consensus|<路径>）
# 使 04 与 05/06/09 用同一套 peak；不同 peakset 的 motif 结果并存在 motif<后缀>/
resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
INPUT_PEAKS="${PEAKSET_FILE}"
PEAK_SOURCE="${CHROMAFLOW_PEAKSET:-idr}"
MOTIF_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/motif${PEAKSET_SUFFIX}"
mkdir -p "${MOTIF_OUT}"

[[ -s "${INPUT_PEAKS}" ]] || { echo "[ERROR] peakset peak 文件为空: ${INPUT_PEAKS}"; exit 1; }

echo "  Pipeline:    ${PIPELINE}"
echo "  peakset:     ${PEAK_SOURCE}  → $(basename "${INPUT_PEAKS}")"
echo "  输入 peaks:  ${INPUT_PEAKS}"
echo "  输出目录:    ${MOTIF_OUT}"
N_PEAKS=$(wc -l < "${INPUT_PEAKS}")
echo "  Peak 数量:   ${N_PEAKS}"
echo ""

if [[ ${N_PEAKS} -lt 10 ]]; then
    echo "[ERROR] peak 数量太少 (<10)"; exit 1
fi

# ─────────────────────────────────────────────
# STEP 1: 基因组背景模型
# ─────────────────────────────────────────────
echo "[STEP 1] 生成基因组背景模型"
# 背景 GC 由 fasta-get-markov 直接从真实序列估计（含 1 阶二核苷酸频率），
# 不需要在 config 里手填 GC 值——手填的数既不会被用到，也可能与基因组不符。
BG_MODEL="${MOTIF_OUT}/genome_background.txt"
if [[ -f "${BG_MODEL}" ]]; then
    echo "  [SKIP] 已存在"
else
    fasta-get-markov -m 1 "${GENOME_COMBINED}" "${BG_MODEL}"
    echo "  [DONE]"
fi
# 从背景模型里读回实测 GC，仅作展示
BG_GC="$(awk '$1=="C"||$1=="G" {s+=$2} END {if (s>0) printf "%.1f", 100*s}' "${BG_MODEL}" 2>/dev/null || true)"
[[ -n "${BG_GC}" ]] && echo "  实测背景 GC: ${BG_GC}%"
echo ""

# ─────────────────────────────────────────────
# STEP 2: 提取 peak 序列
# ─────────────────────────────────────────────
echo "[STEP 2] 提取 peak 序列"

SUMMIT_BED="${MOTIF_OUT}/peak_regions_200bp.bed"
PEAK_FASTA="${MOTIF_OUT}/peak_sequences.fa"

n_cols=$(head -1 "${INPUT_PEAKS}" | awk '{print NF}')

if [[ ${n_cols} -ge 10 ]]; then
    awk 'BEGIN{OFS="\t"} {
        summit = $2 + $10
        start = summit - 100; if(start<0) start=0
        end = summit + 100
        print $1, start, end, $4, $5, "+"
    }' "${INPUT_PEAKS}" > "${SUMMIT_BED}"
    echo "  格式: narrowPeak, summit ± 100bp"
else
    awk 'BEGIN{OFS="\t"} {
        center = int(($2+$3)/2)
        start = center - 100; if(start<0) start=0
        end = center + 100
        print $1, start, end, ($4!=""?$4:"peak_"NR), ($5!=""?$5:0), "+"
    }' "${INPUT_PEAKS}" > "${SUMMIT_BED}"
    echo "  格式: BED, peak 中心 ± 100bp"
fi

bedtools getfasta -fi "${GENOME_COMBINED}" -bed "${SUMMIT_BED}" -fo "${PEAK_FASTA}" -name+
N_SEQS=$(grep -c "^>" "${PEAK_FASTA}")
echo "  提取序列: ${N_SEQS} 条 × 200bp"
echo ""

# ─────────────────────────────────────────────
# 用全部 peak 序列做 MEME-ChIP
# ─────────────────────────────────────────────
# 输入 peak 已经过 IDR/consensus 严格筛选，不再二次截断到 top-N，
# 直接用全部序列，忠实反映当前 peakset 的筛选结果。
INPUT_FASTA="${PEAK_FASTA}"
echo "  使用全部 ${N_SEQS} 条 peak 序列（不再截断 top-N）"
echo ""

# ─────────────────────────────────────────────
# STEP 4: MEME-ChIP
# ─────────────────────────────────────────────
echo "[STEP 3] 运行 MEME-ChIP"
echo ""

MEMECHIP_OUT="${MOTIF_OUT}/meme-chip_output"

DB_ARGS=""
if [[ -f "${MEME_DB}" ]]; then
    DB_ARGS="-db ${MEME_DB}"
    echo "  motif 数据库: $(basename ${MEME_DB})"
else
    echo "  ⚠  未找到 motif 数据库，仅做 de novo 发现"
fi

# 关键：MEME 内部用 -meme-searchsize 限制参与 motif 搜索的序列字符数（默认会截断）。
# 设 0 = 不限制，用全部 peak 序列，忠实反映筛选结果。
# 注意：MEME 运行时间随 searchsize 约平方增长，序列很多时(如 consensus 1161)会较慢；
# 需要限制时把 MEME_SEARCHSIZE 设成一个具体字符数即可。
MEME_SEARCHSIZE="${MEME_SEARCHSIZE:-0}"
echo "  -meme-searchsize=${MEME_SEARCHSIZE}（0=用全部 ${N_SEQS} 条序列，不截断）"

meme-chip \
    -oc "${MEMECHIP_OUT}" \
    -bfile "${BG_MODEL}" \
    ${DB_ARGS} \
    -minw 6 \
    -maxw 15 \
    -meme-nmotifs 10 \
    -meme-searchsize "${MEME_SEARCHSIZE}" \
    -meme-p 8 \
    -centrimo-local \
    "${INPUT_FASTA}"

echo "  [DONE] MEME-ChIP 完成"

# ─────────────────────────────────────────────
# 总结
# ─────────────────────────────────────────────
echo ""
echo "============================================"
echo " Motif 分析总结 (${PIPELINE})"
echo "============================================"
echo ""
echo "  输入 peaks:     ${N_PEAKS} (${PEAK_SOURCE})"
echo "  GC 背景模型:    由基因组序列实测 (${BG_GC:-?}% GC)"
echo "  输出目录:       ${MOTIF_OUT}"
echo "  MEME-ChIP 报告: ${MEMECHIP_OUT}/meme-chip.html"
echo ""
echo "[DONE] $(date)"
echo ""
echo "下一步(用同一 peakset): chromaflow annotate --type ${PIPELINE} ${CHROMAFLOW_PEAKSET:+--peakset ${CHROMAFLOW_PEAKSET}}"