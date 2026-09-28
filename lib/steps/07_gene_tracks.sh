#!/usr/bin/env bash
# ============================================================================
# 07_gene_tracks.sh — 单基因波峰图入口（薄封装，绘图交给 plot_tracks.sh）
# ============================================================================
# 环境: conda activate ChromaFlow_downstream
# 用法:
#   chromaflow tracks --type chipseq --genes GENE1,GENE2
#   chromaflow tracks --type chipseq --gene-list mylist.txt [--outdir DIR]
#   chromaflow tracks --type chipseq --gene-list mylist.txt --peakset consensus
# 本步只负责：解析基因来源 + 用 --peakset 选 peak + 定位 bigWig/GTF/genes.bed，
# 然后调用 plot_tracks.sh 出图。美术在 plot_tracks.sh 里调。
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true
init_pipeline_dirs "${PIPELINE}"

GENES_ARG=""; GENE_LIST_FILE=""; OUTDIR_ARG=""; WITH_MOTIFS=false; FLANK_ARG=""; SELECTION=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genes)     GENES_ARG="$2"; shift 2 ;;
        --gene-list) GENE_LIST_FILE="$2"; shift 2 ;;
        --outdir)    OUTDIR_ARG="$2"; shift 2 ;;
        --motifs)    WITH_MOTIFS=true; shift ;;
        --selection) SELECTION="$2"; shift 2 ;;
        --flank)     FLANK_ARG="$2"; shift 2 ;;
        *) die "Unknown tracks option: $1" ;;
    esac
done

open_log "07_gene_tracks" "${PIPELINE}"
log_section "Step 07: Per-gene Browser Tracks (${PIPELINE})"
if [[ -n "${TRACK_AUTO_TOPN:-}${TRACK_GENES_TOPN:-}" ]]; then
    echo '[INFO] Legacy TRACK_AUTO_TOPN/TRACK_GENES_TOPN are unused; the explicit gene list controls plotting.'
fi

# peakset：选 peak 文件 + 输出后缀（--peakset 全局参数经 CHROMAFLOW_PEAKSET 传入）
resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
if [[ -n "${SELECTION}" ]]; then
    [[ -z "${GENES_ARG}${GENE_LIST_FILE}" && "${WITH_MOTIFS}" == false ]] || die '--selection supplies genes and motif hits; do not combine with --genes, --gene-list or --motifs'
    SELECTION=$(python3 "${CHROMAFLOW_HOME}/lib/steps/support/motif_selection.py" \
        "${SELECTION}" "${CHROMAFLOW_PEAKSET}" "${PIPELINE}" "${PEAKSET_FILE}" "${PROMOTER_UP}" "${PROMOTER_DOWN}")
    GENE_LIST_FILE="${SELECTION}/selected_genes.txt"
    if [[ ! -s "${GENE_LIST_FILE}" ]]; then
        echo '[SKIP] Selection is valid but has no selected genes; no tracks generated.'
        exit 0
    fi
    echo "[SELECTION] Validated ${SELECTION}/selection.json; using its stored gene list and filtered context hits."
fi
ANNO_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}"
TRACK_OUT="${OUTDIR_ARG:-$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/tracks${PEAKSET_SUFFIX}")}"
mkdir -p "${TRACK_OUT}"
echo "  peakset:  ${CHROMAFLOW_PEAKSET:-idr}  → $(basename "${PEAKSET_FILE}")"
echo "  Output directory: ${TRACK_OUT}"

BW_DIR="${RESULTS_DIR}/${PIPELINE}/bwa/merged_library/bigwig"
[[ -d "${BW_DIR}" ]] || BW_DIR="${RESULTS_DIR}/${PIPELINE}"

# genes.bed：优先 05 产物，否则从 GTF 现算
GENES_BED="${ANNO_OUT}/genes.bed"
if [[ ! -s "${GENES_BED}" ]]; then
    echo "  (Step 05 genes.bed not found; deriving gene intervals from GTF)"
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
    python3 - "${CHROMAFLOW_HOME}" "${GENE_LIST_FILE}" "${GENES_TO_PLOT}" <<'PYEOF'
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.argv[1])/'lib/steps/support'))
from promoter_audit import read_query
Path(sys.argv[3]).write_text(''.join(g+'\n' for g in read_query(sys.argv[2])), encoding='utf-8')
PYEOF
    echo "  Gene source: list file ${GENE_LIST_FILE} (detected $(grep -cvE '^\s*$' "${GENES_TO_PLOT}") genes)"
elif [[ -n "${GENES_ARG}" ]]; then
    echo "${GENES_ARG}" | tr ',' '\n' | sed 's/ //g' | grep -v '^$' > "${GENES_TO_PLOT}"
    echo "  Gene source: --genes argument"
else
    die "Specify genes with --genes GENE1,GENE2 or --gene-list FILE"
fi
[[ -s "${GENES_TO_PLOT}" ]] || die "Gene list is empty"
echo ""

ANNO_ARGS=()
if [[ -s "${ANNO_OUT}/peak_annotation_full.tsv" ]]; then
    python3 "${CHROMAFLOW_HOME}/lib/steps/support/annotation_policy.py" check \
        "${ANNO_OUT}/peak_annotation_full.tsv" --up "${PROMOTER_UP}" --down "${PROMOTER_DOWN}"
    ANNO_ARGS=(--annotation "${ANNO_OUT}/peak_annotation_full.tsv")
else
    log_warn 'Full annotation not available; regional peaks will be marked unannotated.'
fi
python3 "${CHROMAFLOW_HOME}/lib/steps/support/promoter_audit.py" \
    --genes-bed "${GENES_BED}" --peaks "${PEAKSET_FILE}" --gene-list "${GENES_TO_PLOT}" \
    --outdir "${TRACK_OUT}" --up "${PROMOTER_UP}" --down "${PROMOTER_DOWN}" \
    --chrom-sizes "${CHROM_SIZES}" "${ANNO_ARGS[@]}"
FLANK="${FLANK_ARG:-${TRACK_FLANK}}"
[[ "${FLANK}" =~ ^[0-9]+$ ]] || die '--flank must be a nonnegative integer'
# Include the entire promoter even when requested context is shorter.
for extent in "${PROMOTER_UP}" "${PROMOTER_DOWN}"; do
    if (( extent > FLANK )); then FLANK="$extent"; fi
done
echo "  Plot flank: ${FLANK} bp (at least the configured promoter extent)"
EXTRA_TRACKS=(--promoter-bed "${TRACK_OUT}/promoter_regions.bed" --chrom-sizes "${CHROM_SIZES}")
if [[ -n "${SELECTION}" ]]; then
    EXTRA_TRACKS+=(--motif-bed "${SELECTION}/motif_hits.bed")
    if [[ -f "${SELECTION}/motif_input_windows.bed" ]]; then
        EXTRA_TRACKS+=(--scan-windows "${SELECTION}/motif_input_windows.bed")
    fi
    printf '%s\n' "${SELECTION}/selection.json" > "${TRACK_OUT}/selection_source.txt"
elif [[ "${WITH_MOTIFS}" == true ]]; then
    MOTIF_DIR="${RESULTS_DIR}/downstream/${PIPELINE}/motif${PEAKSET_SUFFIX}"
    mapfile -t FIMO_FILES < <(find "${MOTIF_DIR}" -path '*/fimo_out_*/fimo.tsv' -type f 2>/dev/null | sort)
    [[ ${#FIMO_FILES[@]} -gt 0 ]] || die "No FIMO results for this peakset: ${MOTIF_DIR}; run motif first"
    SCAN_ARGS=()
    if [[ -s "${MOTIF_DIR}/peak_sequences.fa" ]]; then
        SCAN_ARGS=(--scan-fasta "${MOTIF_DIR}/peak_sequences.fa" --windows-output "${TRACK_OUT}/motif_input_windows.bed")
    fi
    python3 "${CHROMAFLOW_HOME}/lib/steps/support/fimo_tracks.py" \
        --output "${TRACK_OUT}/motif_hits.bed" --audit "${TRACK_OUT}/motif_hits.tsv" \
        --chrom-sizes "${CHROM_SIZES}" --pval "${FIMO_PVAL}" "${SCAN_ARGS[@]}" "${FIMO_FILES[@]}"
    EXTRA_TRACKS+=(--motif-bed "${TRACK_OUT}/motif_hits.bed")
    if [[ ${#SCAN_ARGS[@]} -gt 0 ]]; then
        EXTRA_TRACKS+=(--scan-windows "${TRACK_OUT}/motif_input_windows.bed")
    else
        log_warn 'Motif input FASTA not found; only the available FIMO hits will be shown.'
    fi
fi

# 交给纯绘图脚本
bash "${CHROMAFLOW_HOME}/lib/steps/plot_tracks.sh" \
    --gene-list "${GENES_TO_PLOT}" \
    --genes-bed "${GENES_BED}" \
    --gtf "${GENOME_GTF}" \
    --bigwig-dir "${BW_DIR}" \
    --peaks "${PEAKSET_FILE}" \
    --outdir "${TRACK_OUT}" \
    --flank "${FLANK}" \
    ${TRACK_MAX:+--ymax "${TRACK_MAX}"} "${EXTRA_TRACKS[@]}"

echo "[DONE] $(date)"
