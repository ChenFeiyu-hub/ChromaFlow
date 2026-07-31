#!/usr/bin/env bash
# ============================================================================
# support01_motif_spacing_analysis.sh — Motif 共现与间距分析
# ============================================================================
# 分析多个 motif 在 ChIP/ATAC peaks 内的空间关系
#
# 环境: mamba activate ChromaFlow_ref
#
# 用法:
#   # 自动发现 MEME-ChIP 的所有 fimo 结果
#   chromaflow spacing --type chipseq
#
#   # 手动指定 2 个或多个 fimo.tsv 文件
#   chromaflow spacing --type chipseq \
#       results/.../fimo_out_1/fimo.tsv \
#       results/.../fimo_out_2/fimo.tsv
#
#   # 也可以指定自定义 FIMO 结果（如先验 motif 扫描）
#   chromaflow spacing --type chipseq \
#       results/.../fimo_out_1/fimo.tsv \
#       results/.../fimo_dof/fimo.tsv \
#       results/.../fimo_custom/fimo.tsv
# ============================================================================

set -euo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
shift || true
init_pipeline_dirs "${PIPELINE}"

# ─────────────────────────────────────────────
# 收集 FIMO 输入文件
# ─────────────────────────────────────────────
FIMO_FILES=()

if [[ $# -gt 0 ]]; then
    # 手动指定模式
    for f in "$@"; do
        if [[ -f "$f" ]]; then
            FIMO_FILES+=("$f")
        else
            echo "[WARN] 文件不存在: $f"
        fi
    done
else
    # 自动发现 MEME-ChIP 的所有 fimo 输出
    MEMECHIP_DIR="${RESULTS_DIR}/downstream/${PIPELINE}/motif/meme-chip_output"
    for d in "${MEMECHIP_DIR}"/fimo_out_*/fimo.tsv; do
        [[ -f "$d" ]] && FIMO_FILES+=("$d")
    done
fi

if [[ ${#FIMO_FILES[@]} -lt 2 ]]; then
    echo "[ERROR] 至少需要 2 个 fimo.tsv 文件进行共现分析"
    echo ""
    echo "用法:"
    echo "  bash $0 chipseq                          # 自动发现"
    echo "  bash $0 chipseq fimo1.tsv fimo2.tsv      # 手动指定"
    exit 1
fi

# ─────────────────────────────────────────────
# 路径配置
# ─────────────────────────────────────────────
MOTIF_DIR="${RESULTS_DIR}/downstream/${PIPELINE}/motif"
SPACING_OUT="${MOTIF_DIR}/motif_spacing"
mkdir -p "${SPACING_OUT}"

# 获取 peak 文件
REC_PEAKS="${RESULTS_DIR}/downstream/${PIPELINE}/idr/recommended_peaks.bed"
CONSENSUS=$(get_consensus_peaks "${PIPELINE}")

if [[ -f "${REC_PEAKS}" ]] && [[ -s "${REC_PEAKS}" ]]; then
    PEAK_BED="${REC_PEAKS}"
elif [[ -f "${CONSENSUS}" ]]; then
    PEAK_BED="${CONSENSUS}"
else
    echo "[ERROR] 未找到 peak 文件"; exit 1
fi

open_log "support_motif_spacing" "${PIPELINE}"

log_section "Motif Spacing Analysis (${PIPELINE})"

echo "  FIMO 输入 (${#FIMO_FILES[@]} 个):"
for f in "${FIMO_FILES[@]}"; do
    motif_name=$(grep -v "^#" "$f" | grep -v "^$" | awk 'NR==2{print $2}')
    n_hits=$(grep -v "^#" "$f" | grep -v "^$" | tail -n +2 | wc -l)
    echo "    ${motif_name}: ${n_hits} hits — $(basename $(dirname $f))"
done
echo ""
echo "  Peak 文件: ${PEAK_BED} ($(wc -l < "${PEAK_BED}") peaks)"
echo "  输出目录:  ${SPACING_OUT}"
echo ""

# ─────────────────────────────────────────────
# STEP 1: FIMO hits → BED → intersect peaks
# ─────────────────────────────────────────────
echo "[STEP 1] 将 FIMO hits 映射到 peaks"
echo ""

HITS_IN_PEAKS="${SPACING_OUT}/hits_in_peaks.tsv"
echo -e "motif_id\tmotif_alt\tpeak_chr\tpeak_start\tpeak_end\thit_start\thit_stop\tstrand\tscore\tp_value\tmatched_seq" \
    > "${HITS_IN_PEAKS}"

# 标准化 peak BED 为 3 列
PEAK_BED3="${SPACING_OUT}/_peaks.bed3"
awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' "${PEAK_BED}" | sort -k1,1 -k2,2n > "${PEAK_BED3}"

for fimo_tsv in "${FIMO_FILES[@]}"; do
    # FIMO hits → BED（基因组坐标）
    FIMO_BED="${SPACING_OUT}/_fimo_tmp.bed"
    grep -v "^#" "${fimo_tsv}" | grep -v "^$" | awk 'NR>1{OFS="\t"; print $3, $4-1, $5, $1"|"$2"|"$6"|"$7"|"$8"|"$10}' \
        | sort -k1,1 -k2,2n > "${FIMO_BED}"

    # intersect: 只保留落在 peak 内的 hits，同时记录 peak 坐标
    bedtools intersect -a "${FIMO_BED}" -b "${PEAK_BED3}" -wa -wb 2>/dev/null \
        | awk 'BEGIN{OFS="\t"}{
            split($4, info, "|")
            # info: motif_id, motif_alt, strand, score, pvalue, matched_seq
            print info[1], info[2], $5, $6, $7, $2+1, $3, info[3], info[4], info[5], info[6]
        }' >> "${HITS_IN_PEAKS}"
done

rm -f "${SPACING_OUT}/_fimo_tmp.bed" "${SPACING_OUT}/_peaks.bed3"

n_mapped=$(tail -n +2 "${HITS_IN_PEAKS}" | wc -l)
echo "  落在 peaks 内的 hits: ${n_mapped}"

tail -n +2 "${HITS_IN_PEAKS}" | awk '{print $2}' | sort | uniq -c | sort -rn | while read cnt malt; do
    echo "    ${malt}: ${cnt} hits in peaks"
done
echo ""

# ─────────────────────────────────────────────
# STEP 2: Python 共现与间距分析
# ─────────────────────────────────────────────
echo "[STEP 2] 共现与间距分析"
echo ""

python3 - "${HITS_IN_PEAKS}" "${SPACING_OUT}" << 'PYEOF'
import sys, os
from collections import defaultdict
from itertools import combinations

hits_file = sys.argv[1]
out_dir = sys.argv[2]

# ── 读取 hits ──
hits = []
with open(hits_file) as f:
    header = f.readline()
    for line in f:
        p = line.strip().split("\t")
        if len(p) < 11:
            continue
        hits.append({
            "motif_id": p[0], "motif_alt": p[1],
            "peak_key": f"{p[2]}:{p[3]}-{p[4]}",
            "peak_chr": p[2], "peak_start": int(p[3]), "peak_end": int(p[4]),
            "hit_start": int(p[5]), "hit_stop": int(p[6]),
            "strand": p[7], "score": p[8], "pvalue": p[9],
            "seq": p[10],
        })

if not hits:
    print("  [WARN] 无落在 peak 内的 motif hits")
    sys.exit(0)

# ── 按 peak 分组 ──
peak_hits = defaultdict(list)
for h in hits:
    peak_hits[h["peak_key"]].append(h)

n_peaks = len(peak_hits)
motif_names = sorted(set(h["motif_alt"] for h in hits))

print(f"  含 motif 的 peaks: {n_peaks}")
print(f"  Motif 类型 ({len(motif_names)}): {', '.join(motif_names)}")
print()

# ── 各 motif peak 出现率 ──
print("  === 各 motif 在 peaks 中的出现率 ===")
motif_peak_sets = {}
for m in motif_names:
    pks = set(h["peak_key"] for h in hits if h["motif_alt"] == m)
    motif_peak_sets[m] = pks
    print(f"    {m}: {len(pks)}/{n_peaks} peaks ({len(pks)/n_peaks*100:.1f}%)")
print()

# ── 两两共现分析 ──
print("  === 共现分析 ===")
pair_cooc = {}

for m1, m2 in combinations(motif_names, 2):
    shared = motif_peak_sets[m1] & motif_peak_sets[m2]
    pair_cooc[(m1,m2)] = shared
    print(f"    {m1} vs {m2}:")
    print(f"      共现: {len(shared)} peaks")
    if motif_peak_sets[m1]:
        print(f"      共现率 (相对 {m1}): {len(shared)/len(motif_peak_sets[m1])*100:.1f}%")
    if motif_peak_sets[m2]:
        print(f"      共现率 (相对 {m2}): {len(shared)/len(motif_peak_sets[m2])*100:.1f}%")
    print()

# ── 间距分析 ──
print("  === Peak 内间距分析 ===")

for (m1, m2), shared_peaks in pair_cooc.items():
    if not shared_peaks:
        print(f"    {m1} vs {m2}: 无共现 peak\n")
        continue

    details = []
    for pk in shared_peaks:
        h1_list = [h for h in peak_hits[pk] if h["motif_alt"] == m1]
        h2_list = [h for h in peak_hits[pk] if h["motif_alt"] == m2]

        for h1 in h1_list:
            for h2 in h2_list:
                c1 = (h1["hit_start"] + h1["hit_stop"]) / 2
                c2 = (h2["hit_start"] + h2["hit_stop"]) / 2
                center_dist = abs(c2 - c1)

                if h1["hit_stop"] <= h2["hit_start"]:
                    gap = h2["hit_start"] - h1["hit_stop"]
                elif h2["hit_stop"] <= h1["hit_start"]:
                    gap = h1["hit_start"] - h2["hit_stop"]
                else:
                    gap = 0

                orient = f"{h1['strand']}/{h2['strand']}"
                details.append({
                    "peak": pk, "peak_chr": h1["peak_chr"],
                    "m1_start": h1["hit_start"], "m1_stop": h1["hit_stop"], "m1_strand": h1["strand"],
                    "m2_start": h2["hit_start"], "m2_stop": h2["hit_stop"], "m2_strand": h2["strand"],
                    "center_dist": center_dist, "gap": gap, "orientation": orient,
                    "m1_seq": h1["seq"], "m2_seq": h2["seq"],
                })

    if not details:
        print(f"    {m1} vs {m2}: 无有效配对\n")
        continue

    distances = [d["center_dist"] for d in details]
    gaps = [d["gap"] for d in details]
    n = len(distances)
    ds = sorted(distances)
    gs = sorted(gaps)

    print(f"    {m1} vs {m2}:")
    print(f"      共现 peaks: {len(shared_peaks)}")
    print(f"      motif 配对数: {n}")
    print()
    print(f"      Center-to-center 距离 (bp):")
    print(f"        最小: {ds[0]:.0f}   Q25: {ds[n//4]:.0f}   中位数: {ds[n//2]:.0f}   Q75: {ds[3*n//4]:.0f}   最大: {ds[-1]:.0f}")
    print()
    print(f"      Gap 距离 (bp, 边到边):")
    print(f"        最小: {gs[0]}   Q25: {gs[n//4]}   中位数: {gs[n//2]}   Q75: {gs[3*n//4]}   最大: {gs[-1]}")
    print()

    # 间距分布
    print(f"      间距分布 (center-to-center):")
    bins = [(0,10),(10,20),(20,30),(30,50),(50,75),(75,100),(100,150),(150,200)]
    for lo, hi in bins:
        cnt = sum(1 for d in distances if lo <= d < hi)
        if cnt > 0:
            bar = "#" * min(cnt, 80)
            print(f"        {lo:3d}-{hi:3d}bp: {cnt:4d} {bar}")
    cnt_200 = sum(1 for d in distances if d >= 200)
    if cnt_200 > 0:
        print(f"        >=200bp:  {cnt_200:4d} {'#'*min(cnt_200,80)}")
    print()

    # 方向统计
    orient_counts = defaultdict(int)
    for d in details:
        orient_counts[d["orientation"]] += 1
    print(f"      链方向:")
    for o, c in sorted(orient_counts.items(), key=lambda x: -x[1]):
        print(f"        {o}: {c} ({c/n*100:.1f}%)")
    print()

    # ── 输出详细间距表 ──
    key = f"{m1}_vs_{m2}"
    out_file = os.path.join(out_dir, f"spacing_{key}.tsv")
    with open(out_file, "w") as f:
        f.write("peak\tpeak_chr\tm1_start\tm1_stop\tm1_strand\tm2_start\tm2_stop\tm2_strand\t"
                "center_distance\tgap_distance\torientation\tm1_sequence\tm2_sequence\n")
        for d in sorted(details, key=lambda x: x["center_dist"]):
            f.write(f"{d['peak']}\t{d['peak_chr']}\t{d['m1_start']}\t{d['m1_stop']}\t{d['m1_strand']}\t"
                    f"{d['m2_start']}\t{d['m2_stop']}\t{d['m2_strand']}\t"
                    f"{d['center_dist']:.0f}\t{d['gap']}\t{d['orientation']}\t"
                    f"{d['m1_seq']}\t{d['m2_seq']}\n")
    print(f"    详细间距表: {out_file}")

    # ── 对接设计建议 ──
    print()
    print(f"    --- 分子对接 DNA 设计建议 ---")
    median_gap = gs[n//2]
    q25_gap = gs[n//4]
    q75_gap = gs[3*n//4]
    m1_width = len(details[0]["m1_seq"])
    m2_width = len(details[0]["m2_seq"])
    suggested_len = m1_width + median_gap + m2_width + 10  # 两侧各加 5bp flanking
    print(f"      推荐 DNA 长度: ~{suggested_len}bp")
    print(f"      = {m1_width}bp (motif1) + {median_gap}bp (中位 gap) + {m2_width}bp (motif2) + 10bp (flanking)")
    print(f"      Gap 范围 (Q25-Q75): {q25_gap}-{q75_gap}bp")

    # 典型示例
    typical = [d for d in details if q25_gap <= d["gap"] <= q75_gap]
    if typical:
        ex = min(typical, key=lambda x: abs(x["gap"] - median_gap))
    else:
        ex = min(details, key=lambda x: abs(x["gap"] - median_gap))
    print(f"      典型示例 (gap={ex['gap']}bp):")
    print(f"        Peak: {ex['peak']}")
    print(f"        {m1}: {ex['m1_seq']} ({ex['m1_strand']}, {ex['m1_start']}-{ex['m1_stop']})")
    print(f"        {m2}: {ex['m2_seq']} ({ex['m2_strand']}, {ex['m2_start']}-{ex['m2_stop']})")
    print()

    # ── 输出共现区域 BED ──
    cooc_bed = os.path.join(out_dir, f"cooccurrence_{key}.bed")
    regions = set()
    for d in details:
        start = min(d["m1_start"], d["m2_start"]) - 10
        end = max(d["m1_stop"], d["m2_stop"]) + 10
        regions.add((d["peak_chr"], max(0, start), end))
    with open(cooc_bed, "w") as f:
        for chrom, s, e in sorted(regions):
            f.write(f"{chrom}\t{s}\t{e}\n")
    print(f"    共现区域 BED: {cooc_bed}")

    # ── 提取序列 ──
    cooc_fa = os.path.join(out_dir, f"cooccurrence_{key}.fa")
    print(f"    共现区域 FASTA: {cooc_fa}")
    print(f"    (运行后由 bedtools getfasta 生成)")

# ── 汇总 ──
print()
print("  分析完成")

PYEOF

# ─────────────────────────────────────────────
# STEP 3: 提取共现区域序列
# ─────────────────────────────────────────────
echo ""
echo "[STEP 3] 提取共现区域 DNA 序列"

for cooc_bed in "${SPACING_OUT}"/cooccurrence_*.bed; do
    [[ -f "${cooc_bed}" ]] || continue
    cooc_fa="${cooc_bed%.bed}.fa"
    if [[ -s "${cooc_bed}" ]]; then
        bedtools getfasta -fi "${GENOME_COMBINED}" -bed "${cooc_bed}" -fo "${cooc_fa}" 2>/dev/null || true
        n_seqs=$(grep -c "^>" "${cooc_fa}" 2>/dev/null || echo 0)
        echo "  $(basename ${cooc_bed}): ${n_seqs} 条序列 → $(basename ${cooc_fa})"
    fi
done

# ─────────────────────────────────────────────
# 总结
# ─────────────────────────────────────────────
echo ""
echo "============================================"
echo " Motif 间距分析完成 (${PIPELINE})"
echo "============================================"
echo ""
echo "  输出目录: ${SPACING_OUT}"
echo "  关键文件:"
echo "    - hits_in_peaks.tsv:        落在 peaks 内的 motif hits"
echo "    - spacing_*.tsv:            两两 motif 间距详表"
echo "    - cooccurrence_*.bed/fa:    共现区域坐标与序列"
echo ""
echo "[DONE] $(date)"