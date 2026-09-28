#!/usr/bin/env bash
# ============================================================================
# 08_go_analysis.sh — GO 富集辅助（面向外部 agriGO / g:Profiler 工作流）
# ============================================================================
# 环境: conda activate ChromaFlow_ref   (需 python3 + matplotlib)
# 两个子命令：
#   chromaflow go-prep --type chipseq [--set promoter|dist|all]
#       → 导出可直接上传 agriGO/g:Profiler 的基因 ID 列表 + 背景全集
#   chromaflow go-plot --input <agriGO结果.txt> [--type chipseq] [--top 20]
#       → 把 agriGO/g:Profiler 结果表画成气泡/条形图
# ----------------------------------------------------------------------------
# 思路：富集在外部网站跑（衣藻等物种 agriGO 不支持时用 g:Profiler），
# 本脚本只负责「出干净的输入」和「把结果画成图」，与在线库解耦。
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

MODE="${1:?Usage: go-prep or go-plot}"; shift || true
PIPELINE="chipseq"; SET="promoter"; INPUT=""; TOPN=20
while [[ $# -gt 0 ]]; do
    case "$1" in
        --type)  PIPELINE="$2"; shift 2 ;;
        --set)   SET="$2";      shift 2 ;;
        --input) INPUT="$2";    shift 2 ;;
        --top)   TOPN="$2";     shift 2 ;;
        *) shift ;;
    esac
done
init_pipeline_dirs "${PIPELINE}"
PEAKSET_SUFFIX=""
case "${CHROMAFLOW_PEAKSET:-idr}" in
    idr|idr:0.01) ;;
    idr:*) PEAKSET_SUFFIX="_idr${CHROMAFLOW_PEAKSET#idr:}" ;;
    consensus) PEAKSET_SUFFIX=_consensus ;;
    *) PEAKSET_SUFFIX=_custom ;;
esac
ANNO_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}"
GO_OUT="$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/go${PEAKSET_SUFFIX}")"
mkdir -p "${GO_OUT}"

# ============================================================================
# 子命令 1: go-prep — 导出基因列表
# ============================================================================
if [[ "${MODE}" == "prep" ]]; then
    # A failed upstream IDR run must not silently export an older annotation.
    resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
    open_log "08_go_prep" "${PIPELINE}"
    log_section "Step 08a: Export GO inputs (${PIPELINE}, set=${SET})"

    case "${SET}" in
        promoter) SRC="${ANNO_OUT}/promoter_target_genes.txt" ;;
        dist)     SRC=$(find "${ANNO_OUT}" -name "target_genes_within_*kb.txt" | head -1) ;;
        all)      SRC="${ANNO_OUT}/all_nearest_genes.txt" ;;
        *) die "--set must be promoter|dist|all" ;;
    esac
    [[ -s "${SRC}" ]] || die "Target gene list not found: ${SRC}(run chromaflow annotate first)"
    UNIV="${ANNO_OUT}/background_universe.txt"

    # 原始 ID
    cut -f1 "${SRC}" | grep -v '^$' | sort -u > "${GO_OUT}/target_${SET}.raw.txt"
    # 去转录本/版本后缀（agriGO/g:Profiler 常需裸基因 ID）
    sed -E 's/\.[0-9]+$//' "${GO_OUT}/target_${SET}.raw.txt" | sort -u > "${GO_OUT}/target_${SET}.stripped.txt"
    n_raw=$(wc -l < "${GO_OUT}/target_${SET}.raw.txt")
    n_str=$(wc -l < "${GO_OUT}/target_${SET}.stripped.txt")

    if [[ -s "${UNIV}" ]]; then
        cp "${UNIV}" "${GO_OUT}/background.raw.txt"
        sed -E 's/\.[0-9]+$//' "${UNIV}" | sort -u > "${GO_OUT}/background.stripped.txt"
        n_bg=$(wc -l < "${GO_OUT}/background.stripped.txt")
    else
        n_bg=0; log_warn "No background_universe.txt; supply an appropriate custom background for enrichment"
    fi

    cat > "${GO_OUT}/README_GO.txt" <<EOF
GO enrichment: external analysis instructions
=====================================
Target gene set (${SET}): ${n_raw} genes (after suffix removal: ${n_str})
Background universe: ${n_bg} genes

Input files:
  target_${SET}.stripped.txt   ← Foreground with suffix removed; verify ID namespace before use
  background.stripped.txt       ← Custom background; use the same ID namespace as the foreground

Choose a compatible enrichment tool:
  · For supported plant species, agriGO can use target_${SET}.stripped.txt
  · For other species, use a supported service or a custom GMT
  · If IDs do not match, check target_${SET}.raw.txt and the reference namespace

After obtaining enrichment results:
  chromaflow go-plot --input <downloaded-results-table> --type ${PIPELINE} --top 20
EOF

    echo "  Target genes(${SET}): ${n_raw} (suffix removed: ${n_str})"
    echo "  Background: ${n_bg}"
    echo ""
    echo "  Exported to: ${GO_OUT}/"
    echo "    · target_${SET}.stripped.txt   (foreground)"
    echo "    · background.stripped.txt       (background)"
    echo "    · README_GO.txt                 (tool selection notes)"
    echo ""
    echo "  → After external enrichment analysis:"
    echo "    chromaflow go-plot --input results.tsv --type ${PIPELINE}"
    echo "[DONE] $(date)"
    exit 0
fi

# ============================================================================
# 子命令 2: go-plot — 可视化外部结果
# ============================================================================
if [[ "${MODE}" == "plot" ]]; then
    open_log "08_go_plot" "${PIPELINE}"
    log_section "Step 08b: Plot GO results (${PIPELINE})"
    [[ -n "${INPUT}" && -f "${INPUT}" ]] || die "--input must point to an enrichment result table"
    python3 -c "import matplotlib" 2>/dev/null || die "Requires matplotlib (conda activate ChromaFlow_ref)"

    python3 - "${INPUT}" "${GO_OUT}" "${TOPN}" <<'PYEOF'
import sys, csv, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

inp, outdir, topn = sys.argv[1], sys.argv[2], int(sys.argv[3])

# 读表：自动判定分隔符
with open(inp, encoding="utf-8", errors="ignore") as f:
    sample = f.read(4096); f.seek(0)
    delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    rows = list(csv.reader(f, delimiter=delim))
rows = [r for r in rows if any(c.strip() for c in r)]
header = rows[0]
hl = [h.strip().lower() for h in header]

def find(*names):
    # 先精确匹配列名（避免 "term" 误配 "term_type"）
    for n in names:
        for i,h in enumerate(hl):
            if h == n: return i
    # 再退子串匹配
    for n in names:
        for i,h in enumerate(hl):
            if n in h: return i
    return None

# 兼容 agriGO / g:Profiler / clusterProfiler 常见列名
i_term = find("term_name","description","go_name","name","term")
i_id   = find("term_id","go_acc","go_id","native","goid","acc")
i_fdr  = find("fdr","adjusted","adj_p","q_value","qvalue","p.adjust","padj")
i_p    = find("p_value","pvalue","p-value","p_val")
i_cat  = find("term_type","source","category","ontology","namespace")
i_cnt  = find("intersection_size","queryitem","count","num","genes")

if i_term is None: i_term = 0
pcol = i_fdr if i_fdr is not None else i_p

def num(x):
    try: return float(x)
    except: return None

data=[]
for r in rows[1:]:
    if len(r) <= max(x for x in [i_term,pcol] if x is not None): continue
    p = num(r[pcol]) if pcol is not None else None
    if p is None or p<=0: continue
    term = r[i_term].strip()[:55]
    cat  = r[i_cat].strip() if i_cat is not None and i_cat<len(r) else ""
    cnt  = num(r[i_cnt]) if (i_cnt is not None and i_cnt<len(r)) else None
    data.append({"term":term,"p":p,"cat":cat,"cnt":cnt})

if not data:
    sys.exit("No usable rows found; check enrichment result columns and delimiters")

import math
for d in data: d["score"] = -math.log10(d["p"])
data.sort(key=lambda d:d["score"], reverse=True)
data = data[:topn]
data.reverse()  # 画图从下往上

terms=[d["term"] for d in data]
scores=[d["score"] for d in data]
cnts=[d["cnt"] if d["cnt"] else 20 for d in data]
cats=[d["cat"] for d in data]

# 类别配色
uniq=sorted(set(c for c in cats if c))
palette=["#2b8cbe","#31a354","#e6550d","#756bb1","#c51b8a"]
cmap={c:palette[i%len(palette)] for i,c in enumerate(uniq)}
colors=[cmap.get(c,"#888888") for c in cats]

# 气泡图：x=-log10(FDR), y=term, size=gene count
fig,ax=plt.subplots(figsize=(8, max(3,0.4*len(terms)+1)))
smax=max(cnts) if cnts else 1
sizes=[40+300*(c/smax) for c in cnts]
ax.scatter(scores, range(len(terms)), s=sizes, c=colors, alpha=0.8, edgecolors="k", linewidths=0.4)
ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=8)
ax.set_xlabel("-log10(FDR)"); ax.set_title(f"GO enrichment (top {len(terms)})")
ax.grid(axis="x", ls=":", alpha=0.5)
if uniq:
    from matplotlib.lines import Line2D
    handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=cmap[c],markersize=8,label=c) for c in uniq]
    ax.legend(handles=handles, title="category", fontsize=8, loc="lower right")
plt.tight_layout()
plt.savefig(f"{outdir}/go_enrichment_bubble.png", dpi=200); plt.close()

# 附带条形图
fig,ax=plt.subplots(figsize=(8, max(3,0.4*len(terms)+1)))
ax.barh(range(len(terms)), scores, color=colors, edgecolor="k", linewidth=0.4)
ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=8)
ax.set_xlabel("-log10(FDR)"); ax.set_title(f"GO enrichment (top {len(terms)})")
plt.tight_layout()
plt.savefig(f"{outdir}/go_enrichment_bar.png", dpi=200); plt.close()

print(f"  Parsed {len(rows)-1} rows; plotting top {len(terms)} terms")
print(f"  Column mapping: term={i_term} p={pcol} cat={i_cat} count={i_cnt}")
print(f"  → go_enrichment_bubble.png, go_enrichment_bar.png")
PYEOF
    echo "  Output: ${GO_OUT}"
    echo "[DONE] $(date)"
    exit 0
fi

die "Unknown MODE: ${MODE} (use prep or plot)"
