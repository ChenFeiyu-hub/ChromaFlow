#!/usr/bin/env bash
# ============================================================================
# 07_gene_tracks.sh — 单基因波峰图入口（薄封装，绘图交给 plot_tracks.sh）
# ============================================================================
# 环境: mamba activate ChromaFlow_downstream
# 用法:
#   chromaflow tracks --type chipseq --genes GENE1,GENE2
#   chromaflow tracks --type chipseq --gene-list mylist.txt [--outdir DIR]
#   chromaflow tracks --type chipseq --gene-list mylist.txt --peakset consensus
# 本步只负责：解析基因来源 + 用 --peakset 选 peak + 定位 bigWig/GTF/genes.bed，
# 然后调用 plot_tracks.sh 出图。美术在 plot_tracks.sh 里调。
# ============================================================================
set -uo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true
init_pipeline_dirs "${PIPELINE}"

GENES_ARG=""; GENE_LIST_FILE=""; OUTDIR_ARG=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genes)     GENES_ARG="$2"; shift 2 ;;
        --gene-list) GENE_LIST_FILE="$2"; shift 2 ;;
        --outdir)    OUTDIR_ARG="$2"; shift 2 ;;
        *) shift ;;
    esac
done

open_log "07_gene_tracks" "${PIPELINE}"
log_section "Step 07: Per-gene Browser Tracks (${PIPELINE})"

# peakset：选 peak 文件 + 输出后缀（--peakset 全局参数经 CHROMAFLOW_PEAKSET 传入）
resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
ANNO_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}"
TRACK_OUT="${OUTDIR_ARG:-$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/tracks${PEAKSET_SUFFIX}")}"
mkdir -p "${TRACK_OUT}"
echo "  peakset:  ${CHROMAFLOW_PEAKSET:-idr}  → $(basename "${PEAKSET_FILE}")"
echo "  输出目录: ${TRACK_OUT}"

BW_DIR="${RESULTS_DIR}/${PIPELINE}/bwa/merged_library/bigwig"
[[ -d "${BW_DIR}" ]] || BW_DIR="${RESULTS_DIR}/${PIPELINE}"

# genes.bed：优先 05 产物，否则从 GTF 现算
GENES_BED="${ANNO_OUT}/genes.bed"
if [[ ! -s "${GENES_BED}" ]]; then
    echo "  (未见 05 的 genes.bed，从 GTF 现算)"
    GENES_BED="${TRACK_OUT}/genes.bed"
    GF="gene"; grep -qP '\tgene\t' "${GENOME_GTF}" || GF="transcript"
    python3 - "${GENOME_GTF}" "${GF}" "${GENES_BED}" <<'PYEOF'
import sys, re
gtf, feat, out = sys.argv[1], sys.argv[2], sys.argv[3]
g={}
for line in open(gtf):
    if line.startswith("#"): continue
    c=line.rstrip("\n").split("\t")
    if len(c)<9 or c[2]!=feat: continue
    m=re.search(r'gene_id "([^"]+)"',c[8])
    if not m: continue
    gid=m.group(1); mn=re.search(r'gene_name "([^"]+)"',c[8]); nm=mn.group(1) if mn else gid
    st,en,sd=int(c[3])-1,int(c[4]),c[6]
    d=g.setdefault(gid,{"c":c[0],"s":st,"e":en,"d":sd,"n":nm})
    d["s"]=min(d["s"],st); d["e"]=max(d["e"],en)
with open(out,"w") as w:
    for gid,d in sorted(g.items()): w.write(f'{d["c"]}\t{d["s"]}\t{d["e"]}\t{gid}\t0\t{d["d"]}\t{d["n"]}\n')
PYEOF
fi

# 确定要画的基因（必须 --genes 或 --gene-list；不再默认取 top-N）
GENES_TO_PLOT="${TRACK_OUT}/_genes_to_plot.txt"
if [[ -n "${GENE_LIST_FILE}" && -f "${GENE_LIST_FILE}" ]]; then
    normalize_gene_list "${GENE_LIST_FILE}" "${GENES_TO_PLOT}"
    echo "  基因来源: 列表文件 ${GENE_LIST_FILE} （识别到 $(grep -cvE '^\s*$' "${GENES_TO_PLOT}") 个基因）"
elif [[ -n "${GENES_ARG}" ]]; then
    echo "${GENES_ARG}" | tr ',' '\n' | sed 's/ //g' | grep -v '^$' > "${GENES_TO_PLOT}"
    echo "  基因来源: --genes 参数"
else
    die "需指定基因：--genes GENE1,GENE2  或  --gene-list 文件"
fi
[[ -s "${GENES_TO_PLOT}" ]] || die "基因列表为空"
echo ""

# 交给纯绘图脚本
bash "${CHROMAFLOW_HOME}/lib/steps/plot_tracks.sh" \
    --gene-list "${GENES_TO_PLOT}" \
    --genes-bed "${GENES_BED}" \
    --gtf "${GENOME_GTF}" \
    --bigwig-dir "${BW_DIR}" \
    --peaks "${PEAKSET_FILE}" \
    --outdir "${TRACK_OUT}" \
    --flank "${TRACK_FLANK}" \
    ${TRACK_MAX:+--ymax "${TRACK_MAX}"}

echo "[DONE] $(date)"
