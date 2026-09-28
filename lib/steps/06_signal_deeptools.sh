#!/usr/bin/env bash
# ============================================================================
# 06_signal_deeptools.sh — deepTools 信号可视化
# ============================================================================
# 支持: chipseq / atacseq / dapseq
# 环境: conda activate ChromaFlow_downstream   (含 deepTools)
# 用法: chromaflow signal --type [chipseq|atacseq|dapseq]
# ----------------------------------------------------------------------------
# 补 nf-core 之不足：以「本 TF 的 peak summit」为中心的结合强度热图，
# 以及以 TSS 为中心的 metaprofile —— DAP-seq/ChIP-seq 文章的标配图。
# 产出:
#   peak_summit_heatmap.png / _profile.png   结合是否集中在 summit
#   tss_metaprofile.png / tss_heatmap.png     信号相对基因 TSS 的分布
#   replicate_correlation.png / pca.png       重复一致性（≥2 个 bigWig 时）
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
init_pipeline_dirs "${PIPELINE}"
open_log "06_signal_deeptools" "${PIPELINE}"
log_section "Step 06: deepTools Signal Visualization (${PIPELINE})"

resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
SIGNAL_OUT="$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/signal${PEAKSET_SUFFIX}")"
ANNO_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}"
mkdir -p "${SIGNAL_OUT}"

# ─────────────── 工具检查 ───────────────
for tool in computeMatrix plotHeatmap plotProfile; do
    check_tool "${tool}" || { echo "  Run: conda activate ChromaFlow_downstream"; exit 1; }
done

# ─────────────── 定位 bigWig ───────────────
BW_DIR="${RESULTS_DIR}/${PIPELINE}/bwa/merged_library/bigwig"
mapfile -t BIGWIGS < <(find "${BW_DIR}" -maxdepth 1 \( -name "*.bigWig" -o -name "*.bw" \) 2>/dev/null | sort)
if [[ ${#BIGWIGS[@]} -eq 0 ]]; then
    # 兜底：在 pipeline 结果里递归找
    mapfile -t BIGWIGS < <(find "${RESULTS_DIR}/${PIPELINE}" \( -name "*.bigWig" -o -name "*.bw" \) 2>/dev/null | sort)
fi
[[ ${#BIGWIGS[@]} -gt 0 ]] || die "No bigWig found (expected ${BW_DIR}); confirm that nf-core completed."

LABELS=(); for bw in "${BIGWIGS[@]}"; do LABELS+=("$(basename "${bw}" | sed 's/\.\(bigWig\|bw\)$//')"); done
echo "  Found ${#BIGWIGS[@]} bigWig files:"; printf "    %s\n" "${LABELS[@]}"; echo ""

# ─────────────── 选择 peaks ───────────────
PEAK_BED="${PEAKSET_FILE}"; PEAK_SOURCE="${CHROMAFLOW_PEAKSET:-idr}"
[[ -s "${PEAK_BED}" ]] || die "Peak set is missing or empty: ${PEAK_BED}"
# 标准化为 BED（computeMatrix 需要）
PEAK_BED_STD="${SIGNAL_OUT}/peaks.bed"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,(NF>=4?$4:"peak_"NR)}' "${PEAK_BED}" | sort -k1,1 -k2,2n > "${PEAK_BED_STD}"
echo "  Peaks: $(wc -l < "${PEAK_BED_STD}") (${PEAK_SOURCE})"

# ─────────────── 基因 BED（TSS 参考）───────────────
GENES_BED6="${SIGNAL_OUT}/genes.bed6"
if [[ -s "${ANNO_OUT}/genes.bed" ]]; then
    # 复用 step 05 产物: chr,start,end,gid,0,strand,name → BED6
    # %d 强制定点输出：兼容旧版 05 可能写出的 "8e+06" 科学计数坐标
    awk 'BEGIN{OFS="\t"} NF>=6 {printf "%s\t%d\t%d\t%s\t0\t%s\n", $1,$2,$3,$4,$6}' \
        "${ANNO_OUT}/genes.bed" > "${GENES_BED6}"
elif [[ -f "${GENOME_GTF}" ]]; then
    GF="gene"; grep -qP '\tgene\t' "${GENOME_GTF}" || GF="transcript"
    python3 - "${GENOME_GTF}" "${GF}" "${GENES_BED6}" <<'PYEOF'
import sys, re
gtf, feat, out = sys.argv[1], sys.argv[2], sys.argv[3]
g = {}
for line in open(gtf):
    if line.startswith("#"): continue
    c = line.rstrip("\n").split("\t")
    if len(c) < 9 or c[2] != feat: continue
    m = re.search(r'gene_id "([^"]+)"', c[8])
    if not m: continue
    gid = m.group(1); st, en, sd = int(c[3])-1, int(c[4]), c[6]
    d = g.setdefault(gid, {"c":c[0],"s":st,"e":en,"d":sd})
    d["s"] = min(d["s"], st); d["e"] = max(d["e"], en)
with open(out,"w") as w:
    for gid,d in sorted(g.items()):
        w.write(f'{d["c"]}\t{d["s"]}\t{d["e"]}\t{gid}\t0\t{d["d"]}\n')
PYEOF
else
    log_warn "No gene BED available; skipping TSS metaprofile"
fi
echo ""

# ─────────────── A. peak summit 中心热图 ───────────────
echo "[A] Peak midpoint +/-${SIGNAL_PEAK_WIN} bp signal heatmap"
MAT_PEAK="${SIGNAL_OUT}/matrix_peaks.gz"
computeMatrix reference-point \
    --referencePoint center \
    -R "${PEAK_BED_STD}" \
    -S "${BIGWIGS[@]}" \
    --samplesLabel "${LABELS[@]}" \
    -a "${SIGNAL_PEAK_WIN}" -b "${SIGNAL_PEAK_WIN}" \
    --binSize "${DEEPT_BIN_SIZE}" \
    --skipZeros --missingDataAsZero \
    -p "${THREADS}" \
    -o "${MAT_PEAK}"

plotHeatmap -m "${MAT_PEAK}" -o "${SIGNAL_OUT}/peak_midpoint_heatmap.png" \
    --refPointLabel 'Peak midpoint' --colorMap Blues \
    --heatmapHeight 12 --whatToShow 'heatmap and colorbar'
plotProfile -m "${MAT_PEAK}" -o "${SIGNAL_OUT}/peak_midpoint_profile.png" \
    --refPointLabel 'Peak midpoint' --perGroup
echo "  → peak_midpoint_heatmap.png, peak_midpoint_profile.png"
echo ""

# ─────────────── B. TSS metaprofile ───────────────
if [[ -s "${GENES_BED6}" ]]; then
    echo "[B] TSS ±${SIGNAL_TSS_WIN}bp metaprofile"
    MAT_TSS="${SIGNAL_OUT}/matrix_tss.gz"
    computeMatrix reference-point \
        --referencePoint TSS \
        -R "${GENES_BED6}" \
        -S "${BIGWIGS[@]}" \
        --samplesLabel "${LABELS[@]}" \
        -a "${SIGNAL_TSS_WIN}" -b "${SIGNAL_TSS_WIN}" \
        --binSize "${DEEPT_BIN_SIZE}" \
        --skipZeros --missingDataAsZero \
        -p "${THREADS}" \
        -o "${MAT_TSS}"
    plotProfile -m "${MAT_TSS}" -o "${SIGNAL_OUT}/tss_metaprofile.png" \
        --refPointLabel TSS --perGroup
    plotHeatmap -m "${MAT_TSS}" -o "${SIGNAL_OUT}/tss_heatmap.png" \
        --refPointLabel TSS --colorMap Blues \
        --heatmapHeight 12 --whatToShow 'heatmap and colorbar'
    echo "  → tss_metaprofile.png, tss_heatmap.png"
    echo ""
fi

# ─────────────── C. 重复相关性（≥2 bigWig）───────────────
if [[ ${#BIGWIGS[@]} -ge 2 ]]; then
    for tool in multiBigwigSummary plotCorrelation plotPCA; do check_tool "$tool" || exit 1; done
    echo "[C] Replicate concordance (Spearman correlation + PCA)"
    SUMMARY="${SIGNAL_OUT}/bigwig_summary.npz"
    multiBigwigSummary bins -b "${BIGWIGS[@]}" --labels "${LABELS[@]}" \
        --binSize 1000 -p "${THREADS}" -o "${SUMMARY}"
    plotCorrelation -in "${SUMMARY}" \
        --corMethod spearman --whatToPlot heatmap --skipZeros \
        --plotNumbers --colorMap RdBu \
        -o "${SIGNAL_OUT}/replicate_correlation.png"
    plotPCA -in "${SUMMARY}" -o "${SIGNAL_OUT}/pca.png"
    echo "  → replicate_correlation.png, pca.png"
    echo ""
fi

# ─────────────── 总结 ───────────────
echo "============================================"
echo " Signal visualization summary (${PIPELINE})"
echo "============================================"
echo "  bigWig count:  ${#BIGWIGS[@]}"
echo "  Peaks:      $(wc -l < "${PEAK_BED_STD}") (${PEAK_SOURCE})"
echo "  Output directory:   ${SIGNAL_OUT}"
echo ""
echo "  Interpretation notes:"
echo "    · peak_midpoint_heatmap: Inspect whether signal is concentrated near peak centers"
echo "    · tss_metaprofile:     Proximal-promoter binding may produce enrichment near TSS"
echo "    · replicate_correlation: Inspect replicate correlation alongside other QC metrics"
echo ""
echo "[DONE] $(date)"
