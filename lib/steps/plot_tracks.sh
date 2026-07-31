#!/usr/bin/env bash
# ============================================================================
# plot_tracks.sh — 纯绘图核心：单基因浏览器波峰图 (pyGenomeTracks)
# ============================================================================
# 只负责「画图」，不含任何业务逻辑（不知道基因从哪来、用哪套 peak）。
# 主流程(07 tracks)、反查(lookup) 等都通过调用本脚本出图，便于统一做美术优化。
#
# 用法（全部为显式参数）:
#   plot_tracks.sh --gene-list F --genes-bed F --gtf F --bigwig-dir DIR --outdir DIR
#                  [--peaks F] [--flank N] [--ymax N] [--title-suffix S]
# ----------------------------------------------------------------------------
# 美术相关参数集中在本文件顶部的 STYLE 段，改这里即可统一调整外观。
# ============================================================================
set -uo pipefail
source "${CHROMAFLOW_HOME:?需 CHROMAFLOW_HOME}/lib/common.sh"

# ─────────────── STYLE（美术：改这里，不动逻辑）───────────────
STYLE_BW_COLORS=("#08519c" "#3182bd" "#6baed6" "#9ecae1" "#c6dbef")
STYLE_BW_HEIGHT=2
STYLE_BW_BINS=500
STYLE_PEAK_COLOR="#e6550d"
STYLE_PEAK_HEIGHT=0.6
STYLE_GENE_COLOR="#333333"
STYLE_GENE_HEIGHT=2
STYLE_FONTSIZE=10
STYLE_DPI=200
# 差值(log2 ChIP/Input)轨道：常驻，独立 y 轴，正=富集(红)/负=低于背景(蓝)
STYLE_RATIO_COLOR="#d73027"
STYLE_RATIO_NEG_COLOR="#4575b4"
STYLE_RATIO_HEIGHT=2.5
STYLE_RATIO_BINSIZE=50
STYLE_RATIO_PSEUDOCOUNT=1
STYLE_INPUT_REGEX="input|control|igg|mock"   # 用于从 bigWig 里识别对照(不分大小写)

# ─────────────── 解析参数 ───────────────
GENE_LIST=""; GENES_BED=""; GTF=""; BW_DIR=""; OUT=""; PEAKS=""
FLANK=3000; YMAX=""; TITLE_SUFFIX=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gene-list)    GENE_LIST="$2"; shift 2 ;;
        --genes-bed)    GENES_BED="$2"; shift 2 ;;
        --gtf)          GTF="$2"; shift 2 ;;
        --bigwig-dir)   BW_DIR="$2"; shift 2 ;;
        --outdir)       OUT="$2"; shift 2 ;;
        --peaks)        PEAKS="$2"; shift 2 ;;
        --flank)        FLANK="$2"; shift 2 ;;
        --ymax)         YMAX="$2"; shift 2 ;;
        --title-suffix) TITLE_SUFFIX="$2"; shift 2 ;;
        *) shift ;;
    esac
done

[[ -s "${GENE_LIST}" ]] || die "plot_tracks: 缺 --gene-list 或为空"
[[ -s "${GENES_BED}" ]] || die "plot_tracks: 缺 --genes-bed（基因坐标）"
[[ -f "${GTF}" ]]       || die "plot_tracks: 缺 --gtf"
[[ -n "${OUT}" ]]       || die "plot_tracks: 缺 --outdir"
mkdir -p "${OUT}"
check_tool pyGenomeTracks || { echo "  需 pyGenomeTracks（mamba activate ChromaFlow_downstream）"; exit 1; }

# ─────────────── bigWig ───────────────
mapfile -t BIGWIGS < <(find "${BW_DIR}" -maxdepth 1 \( -name "*.bigWig" -o -name "*.bw" \) 2>/dev/null | sort)
[[ ${#BIGWIGS[@]} -gt 0 ]] || die "plot_tracks: ${BW_DIR} 下无 bigWig"
LABELS=(); for bw in "${BIGWIGS[@]}"; do LABELS+=("$(basename "${bw}" | sed 's/\.\(bigWig\|bw\)$//')"); done

# 区分 ChIP 与 Input（按文件名）：用于差值图
CHIP_BW=(); INPUT_BW=()
for bw in "${BIGWIGS[@]}"; do
    bn=$(basename "${bw}")
    if echo "${bn}" | grep -qiE "${STYLE_INPUT_REGEX}"; then INPUT_BW+=("${bw}"); else CHIP_BW+=("${bw}"); fi
done

# ─────────────── 差值轨道：log2(ChIP/Input)，常驻 ───────────────
# 用 deepTools：多个 ChIP 重复先 bigwigAverage 合并，再 bigwigCompare 对 Input 求 log2 比值。
# 只算一次（本次调用的所有基因共用），缓存于输出目录。
RATIO_BW=""
if [[ ${#CHIP_BW[@]} -gt 0 && ${#INPUT_BW[@]} -gt 0 ]] && command -v bigwigCompare &>/dev/null; then
    RDIR="${OUT}/.ratio"; mkdir -p "${RDIR}"
    CHIP_MEAN="${RDIR}/chip_mean.bw"
    if [[ ${#CHIP_BW[@]} -eq 1 ]]; then
        CHIP_MEAN="${CHIP_BW[0]}"
    elif command -v bigwigAverage &>/dev/null; then
        bigwigAverage -b "${CHIP_BW[@]}" --binSize "${STYLE_RATIO_BINSIZE}" -o "${CHIP_MEAN}" >/dev/null 2>&1 \
            || CHIP_MEAN="${CHIP_BW[0]}"
    else
        CHIP_MEAN="${CHIP_BW[0]}"   # 无 bigwigAverage 时退用第一个 ChIP
    fi
    RATIO_BW="${RDIR}/log2ratio.bw"
    if ! bigwigCompare -b1 "${CHIP_MEAN}" -b2 "${INPUT_BW[0]}" \
            --operation log2 --pseudocount "${STYLE_RATIO_PSEUDOCOUNT}" \
            --binSize "${STYLE_RATIO_BINSIZE}" -o "${RATIO_BW}" >/dev/null 2>&1; then
        RATIO_BW=""; echo "  ⚠ 差值轨道生成失败，仅画覆盖轨道"
    else
        echo "  差值轨道: log2(ChIP合并/Input) 已生成"
    fi
elif [[ ${#INPUT_BW[@]} -eq 0 ]]; then
    echo "  （未识别到 Input bigWig，跳过差值轨道；对照名需含 ${STYLE_INPUT_REGEX}）"
fi

# ─────────────── peaks 标准化（把 signalValue 富集值写进标签，量化"高于背景"）───────────────
PEAKS_STD="${OUT}/peaks.bed"
if [[ -s "${PEAKS}" ]]; then
    # peakset 文件第5列=signalValue(MACS 富集倍数)。标签显示为 "name (sig=NNN)"，
    # 让图上直接看到定量富集，而不是只靠波峰形状（形状受各图 y 轴自适应影响会误导）。
    awk 'BEGIN{OFS="\t"}{
        nm=(NF>=4?$4:"peak_"NR);
        if(NF>=5 && ($5+0)>0) nm=nm" (sig="sprintf("%.0f",$5)")";
        print $1,$2,$3,nm }' "${PEAKS}" | sort -k1,1 -k2,2n > "${PEAKS_STD}"
else
    : > "${PEAKS_STD}"
fi

# ─────────────── 统一 y 上限（未给则按区域自动算）───────────────
if [[ -z "${YMAX}" ]] && python3 -c "import pyBigWig" 2>/dev/null; then
    YMAX=$(python3 - "${GENES_BED}" "${GENE_LIST}" "${FLANK}" "${BIGWIGS[@]}" <<'PYEOF'
import sys, pyBigWig
genes_bed, glist, flank = sys.argv[1], sys.argv[2], int(sys.argv[3])
bws = sys.argv[4:]
coords={}
for line in open(genes_bed):
    c=line.rstrip("\n").split("\t"); coords[c[3]]=(c[0],int(c[1]),int(c[2]))
    if len(c)>=7: coords[c[6]]=(c[0],int(c[1]),int(c[2]))
regions=[]
for g in open(glist):
    g=g.strip()
    if g in coords:
        ch,s,e=coords[g]; regions.append((ch,max(0,s-flank),e+flank))
mx=0.0
for p in bws:
    bw=pyBigWig.open(p)
    for ch,s,e in regions:
        if ch in bw.chroms():
            e2=min(e, bw.chroms()[ch])
            try:
                v=bw.stats(ch,s,e2,type="max")[0]
                if v: mx=max(mx,v)
            except: pass
    bw.close()
print(int(mx*1.1)+1 if mx>0 else 0)
PYEOF
)
    [[ "${YMAX}" == "0" ]] && YMAX=""
fi
[[ -n "${YMAX}" ]] && echo "  统一 y 上限: ${YMAX}" || echo "  y 上限: 自动(每图自适应)"

# 差值轨道的独立 y 范围（基于区域数据，与覆盖轨道 0–6 量纲分开，避免误读）
RATIO_MIN=""; RATIO_MAX=""
if [[ -n "${RATIO_BW}" ]] && python3 -c "import pyBigWig" 2>/dev/null; then
    read -r RATIO_MIN RATIO_MAX < <(python3 - "${GENES_BED}" "${GENE_LIST}" "${FLANK}" "${RATIO_BW}" <<'PYEOF'
import sys, pyBigWig
genes_bed, glist, flank, rbw = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
coords={}
for line in open(genes_bed):
    c=line.rstrip("\n").split("\t"); coords[c[3]]=(c[0],int(c[1]),int(c[2]))
    if len(c)>=7: coords[c[6]]=(c[0],int(c[1]),int(c[2]))
regions=[]
for g in open(glist):
    g=g.strip()
    if g in coords:
        ch,s,e=coords[g]; regions.append((ch,max(0,s-flank),e+flank))
bw=pyBigWig.open(rbw); lo=0.0; hi=0.0
for ch,s,e in regions:
    if ch in bw.chroms():
        e2=min(e,bw.chroms()[ch])
        try:
            mx=bw.stats(ch,s,e2,type="max")[0]; mn=bw.stats(ch,s,e2,type="min")[0]
            if mx is not None: hi=max(hi,mx)
            if mn is not None: lo=min(lo,mn)
        except: pass
bw.close()
import math
hi=math.ceil(hi+0.5); lo=math.floor(min(lo,-0.5))
print(lo, hi)
PYEOF
)
    echo "  差值轨道 y 范围: ${RATIO_MIN} ~ ${RATIO_MAX}"
fi

# ─────────────── 生成 tracks.ini ───────────────
INI="${OUT}/tracks.ini"
{
    i=0
    for bw in "${BIGWIGS[@]}"; do
        col="${STYLE_BW_COLORS[$(( i % ${#STYLE_BW_COLORS[@]} ))]}"
        echo "[bw_${i}]"
        echo "file = ${bw}"
        echo "title = ${LABELS[$i]}"
        echo "height = ${STYLE_BW_HEIGHT}"
        echo "color = ${col}"
        echo "min_value = 0"
        [[ -n "${YMAX}" ]] && echo "max_value = ${YMAX}"
        echo "file_type = bigwig"
        echo "number_of_bins = ${STYLE_BW_BINS}"
        echo ""
        i=$((i+1))
    done
    # 差值轨道（常驻）：独立 y 轴 + 正红/负蓝，直观显示相对 Input 的富集
    if [[ -n "${RATIO_BW}" ]]; then
        echo "[spacer]"; echo "height = 0.4"; echo ""
        echo "[ratio]"
        echo "file = ${RATIO_BW}"
        echo "title = log2(ChIP/Input)"
        echo "height = ${STYLE_RATIO_HEIGHT}"
        echo "color = ${STYLE_RATIO_COLOR}"
        echo "negative_color = ${STYLE_RATIO_NEG_COLOR}"
        [[ -n "${RATIO_MIN}" ]] && echo "min_value = ${RATIO_MIN}"
        [[ -n "${RATIO_MAX}" ]] && echo "max_value = ${RATIO_MAX}"
        echo "file_type = bigwig"
        echo "number_of_bins = ${STYLE_BW_BINS}"
        echo ""
    fi
    if [[ -s "${PEAKS_STD}" ]]; then
        echo "[spacer]"; echo "height = 0.3"; echo ""
        echo "[peaks]"
        echo "file = ${PEAKS_STD}"
        echo "title = peaks"
        echo "color = ${STYLE_PEAK_COLOR}"
        echo "height = ${STYLE_PEAK_HEIGHT}"
        echo "display = interleaved"
        echo "labels = true"
        echo "fontsize = ${STYLE_FONTSIZE}"
        echo "file_type = bed"
        echo ""
    fi
    echo "[spacer]"; echo "height = 0.3"; echo ""
    echo "[genes]"
    echo "file = ${GTF}"
    echo "title = genes"
    echo "prefered_name = gene_name"
    echo "merge_transcripts = false"
    echo "height = ${STYLE_GENE_HEIGHT}"
    echo "color = ${STYLE_GENE_COLOR}"
    echo "file_type = gtf"
    echo "fontsize = ${STYLE_FONTSIZE}"
    echo ""
    echo "[x-axis]"
    echo "fontsize = ${STYLE_FONTSIZE}"
} > "${INI}"

# ─────────────── 逐基因渲染 ───────────────
echo "[渲染]"
declare -A COORD
while IFS=$'\t' read -r ch s e gid sc st nm; do
    COORD["${gid}"]="${ch}:${s}:${e}:${st}"
    [[ -n "${nm}" ]] && COORD["${nm}"]="${ch}:${s}:${e}:${st}"
done < "${GENES_BED}"

N=$(grep -cvE '^\s*$' "${GENE_LIST}"); ok_n=0
while read -r g; do
    [[ -z "${g}" ]] && continue
    info="${COORD[${g}]:-}"
    if [[ -z "${info}" ]]; then echo "  [跳过] 找不到坐标: ${g}"; continue; fi
    IFS=':' read -r ch s e st <<< "${info}"
    rs=$(( s - FLANK )); [[ ${rs} -lt 0 ]] && rs=0
    re=$(( e + FLANK ))
    region="${ch}:${rs}-${re}"
    out="${OUT}/${g}.png"
    if pyGenomeTracks --tracks "${INI}" --region "${region}" \
        --title "${g} (${region}, ${st})${TITLE_SUFFIX}" --dpi "${STYLE_DPI}" \
        -o "${out}" >/dev/null 2>&1; then
        echo "  [OK] ${g}  ${region}  → $(basename "${out}")"
        ok_n=$((ok_n+1))
    else
        echo "  [FAIL] ${g}  ${region}"
    fi
done < "${GENE_LIST}"

echo "  成功: ${ok_n}/${N}   输出: ${OUT}"
echo "  ini:  ${INI}  (改 STYLE 段或此文件可调美术)"
