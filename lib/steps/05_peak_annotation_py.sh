#!/usr/bin/env bash
# ============================================================================
# [注] 这是 05 的 Python/bedtools 备选实现（无 R 环境时用）。默认走 ChIPseeker 版 05_peak_annotation.R
# 05_peak_annotation.sh — Peak 注释 + 靶基因（Python + bedtools，无需 R）
# ============================================================================
# 相比旧版的改动：
#   [架构] 所有路径从项目配置派生（OUTPUT_DIR 等），不再写回代码目录
#   [修复] TSS 距离改用 `bedtools closest -D b`（链感知，上游恒为负）
#   [新增] 启动子窗口可配（PROMOTER_UP/DOWN），并防越界到相邻基因
#   [新增] 距离过滤靶基因集（TARGET_DIST_KB），三档置信度输出
#   [新增] 本地 GO 富集输入（gene2go + 背景全集，衣藻可用；不再依赖 agriGO）
#   [新增] matplotlib 出图：基因组分布 + TSS 距离直方图
#   [健壮] N_PEAKS=0 除零防护
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
init_pipeline_dirs "${PIPELINE}"
open_log "05_peak_annotation" "${PIPELINE}"
log_section "Step 05: Peak Annotation (${PIPELINE})"

# ─────────────── 路径与输入选择 ───────────────
resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
INPUT_PEAKS="${PEAKSET_FILE}"; PEAK_SOURCE="${CHROMAFLOW_PEAKSET:-idr}"
[[ -s "${INPUT_PEAKS}" ]] || die "Peak set is missing or empty: ${INPUT_PEAKS}"
ANNO_OUT="$(step_outdir "$(downstream_dir "${PIPELINE}")/annotation${PEAKSET_SUFFIX}")"
mkdir -p "${ANNO_OUT}"

echo "  Pipeline:      ${PIPELINE}"
echo "  Input peaks:    ${INPUT_PEAKS} (${PEAK_SOURCE})"
echo "  GTF:           ${GENOME_GTF}"
echo "  Promoter window:    -${PROMOTER_UP} / +${PROMOTER_DOWN} bp"
echo "  Target-gene distance threshold:  ${TARGET_DIST_KB} kb"
echo "  Output directory:      ${ANNO_OUT}"
echo ""

for tool in bedtools python3; do check_tool "${tool}" || exit 1; done

# ─────────────── STEP 1: 从 GTF 提取注释区间 ───────────────
echo "[STEP 1] Extract gene intervals from GTF"
GENES_BED="${ANNO_OUT}/genes.bed"
TSS_BED="${ANNO_OUT}/tss.bed"
PROMOTER_BED="${ANNO_OUT}/promoters.bed"
EXON_BED="${ANNO_OUT}/exons.bed"

GENE_FEATURE="gene"
grep -qP '\tgene\t' "${GENOME_GTF}" || { GENE_FEATURE="transcript"; echo "  GTF has no gene records; aggregate transcript intervals by gene"; }

python3 - "${GENOME_GTF}" "${GENE_FEATURE}" "${GENES_BED}" <<'PYEOF'
import sys, re
gtf, feat, out = sys.argv[1], sys.argv[2], sys.argv[3]
genes = {}
with open(gtf) as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != feat: continue
        m = re.search(r'gene_id "([^"]+)"', c[8])
        if not m: continue
        gid = m.group(1)
        mn = re.search(r'gene_name "([^"]+)"', c[8])
        name = mn.group(1) if mn else gid
        st, en, strand = int(c[3]) - 1, int(c[4]), c[6]
        g = genes.setdefault(gid, {"chr": c[0], "start": st, "end": en, "strand": strand, "name": name})
        g["start"] = min(g["start"], st); g["end"] = max(g["end"], en)
rows = sorted((g["chr"], g["start"], g["end"], gid, "0", g["strand"], g["name"]) for gid, g in genes.items())
with open(out, "w") as w:
    for r in rows: w.write("\t".join(map(str, r)) + "\n")
print(f"  Gene count: {len(rows)}")
PYEOF

# TSS
awk 'BEGIN{OFS="\t"}{ tss=($6=="+")?$2:$3-1; print $1,tss,tss+1,$4,$5,$6,$7 }' \
    "${GENES_BED}" | sort -k1,1 -k2,2n > "${TSS_BED}"
echo "  TSS: $(wc -l < "${TSS_BED}") records"

# Promoter：窗口可配；用 chrom.sizes 防止越界（右端也裁）
awk -v up="${PROMOTER_UP}" -v down="${PROMOTER_DOWN}" 'BEGIN{OFS="\t"}{
    if($6=="+"){s=$2-up; e=$2+down} else {s=$3-down; e=$3+up}
    if(s<0)s=0; print $1,s,e,$4,$5,$6,$7
}' "${TSS_BED}" | sort -k1,1 -k2,2n > "${PROMOTER_BED}.tmp"
if [[ -f "${CHROM_SIZES}" ]]; then
    bedtools slop -i "${PROMOTER_BED}.tmp" -g "${CHROM_SIZES}" -b 0 \
        | awk 'BEGIN{OFS="\t"}{print}' > "${PROMOTER_BED}"
    rm -f "${PROMOTER_BED}.tmp"
else
    mv "${PROMOTER_BED}.tmp" "${PROMOTER_BED}"
fi
echo "  Promoter: $(wc -l < "${PROMOTER_BED}") regions"

# Exon（合并）
python3 - "${GENOME_GTF}" "${EXON_BED}" <<'PYEOF'
import sys
gtf, out = sys.argv[1], sys.argv[2]
ex = []
with open(gtf) as f:
    for line in f:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "exon": continue
        ex.append((c[0], int(c[3]) - 1, int(c[4])))
ex.sort()
with open(out, "w") as w:
    for e in ex: w.write(f"{e[0]}\t{e[1]}\t{e[2]}\n")
PYEOF
bedtools merge -i "${EXON_BED}" > "${EXON_BED}.m" && mv "${EXON_BED}.m" "${EXON_BED}"
echo "  Exon (merged): $(wc -l < "${EXON_BED}") regions"
echo ""

# ─────────────── STEP 2: 基因组分布注释（分层互斥）───────────────
echo "[STEP 2] Annotate peaks (Promoter -> Exon -> Intron -> Intergenic)"
N_PEAKS=$(wc -l < "${INPUT_PEAKS}")
[[ "${N_PEAKS}" -gt 0 ]] || die "Input peaks are empty"

INPUT_BED="${ANNO_OUT}/input_peaks.bed"
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,(NF>=4?$4:"peak_"NR),(NF>=5?$5:0)}' \
    "${INPUT_PEAKS}" | sort -k1,1 -k2,2n > "${INPUT_BED}"

pct() { awk "BEGIN{printf \"%.1f\", ($1)/($2)*100}"; }

bedtools intersect -a "${INPUT_BED}" -b "${PROMOTER_BED}" -u > "${ANNO_OUT}/peaks_in_promoter.bed"
bedtools intersect -a "${INPUT_BED}" -b "${PROMOTER_BED}" -v \
    | bedtools intersect -a - -b "${EXON_BED}" -u > "${ANNO_OUT}/peaks_in_exon.bed"
bedtools intersect -a "${INPUT_BED}" -b "${PROMOTER_BED}" -v \
    | bedtools intersect -a - -b "${EXON_BED}" -v \
    | bedtools intersect -a - -b "${GENES_BED}" -u > "${ANNO_OUT}/peaks_in_intron.bed"
bedtools intersect -a "${INPUT_BED}" -b "${PROMOTER_BED}" -v \
    | bedtools intersect -a - -b "${GENES_BED}" -v > "${ANNO_OUT}/peaks_intergenic.bed"

n_promoter=$(wc -l < "${ANNO_OUT}/peaks_in_promoter.bed")
n_exon=$(wc -l < "${ANNO_OUT}/peaks_in_exon.bed")
n_intron=$(wc -l < "${ANNO_OUT}/peaks_in_intron.bed")
n_intergenic=$(wc -l < "${ANNO_OUT}/peaks_intergenic.bed")

LABEL="Annotation distribution"; [[ "${PIPELINE}" == "atacseq" ]] && LABEL="Accessible region distribution"
echo "  ${LABEL}:"
printf "    %-12s %6d (%s%%)\n" "Promoter" "${n_promoter}"   "$(pct "${n_promoter}"   "${N_PEAKS}")"
printf "    %-12s %6d (%s%%)\n" "Exon"     "${n_exon}"       "$(pct "${n_exon}"       "${N_PEAKS}")"
printf "    %-12s %6d (%s%%)\n" "Intron"   "${n_intron}"     "$(pct "${n_intron}"     "${N_PEAKS}")"
printf "    %-12s %6d (%s%%)\n" "Intergenic" "${n_intergenic}" "$(pct "${n_intergenic}" "${N_PEAKS}")"
echo ""

# ─────────────── STEP 3: 最近 TSS + 链感知距离（修复 -D b）───────────────
echo "[STEP 3] Assign the nearest TSS and compute strand-relative distance (-D b)"
NEAREST="${ANNO_OUT}/peaks_nearest_tss.tsv"
FULL_ANNO="${ANNO_OUT}/peak_annotation_full.tsv"

CLOSEST_ARGS=(-a "${INPUT_BED}" -b "${TSS_BED}" -D b -t first)
[[ -f "${CHROM_SIZES}" ]] && CLOSEST_ARGS+=(-g "${CHROM_SIZES}")
bedtools closest "${CLOSEST_ARGS[@]}" > "${NEAREST}"

if [[ -s "${NEAREST}" ]]; then
    echo -e "peak_chr\tpeak_start\tpeak_end\tpeak_id\tpeak_score\tgene_chr\ttss_start\ttss_end\tgene_id\tgene_score\tstrand\tgene_name\tdist_to_TSS" > "${FULL_ANNO}"
    cat "${NEAREST}" >> "${FULL_ANNO}"
    echo "  Annotation table: $(wc -l < "${NEAREST}") lines → $(basename "${FULL_ANNO}")"
    echo "  TSS distance distribution (negative=upstream, positive=downstream):"
    awk '{ d=($NF<0)?-$NF:$NF
        if(d<=500)b="0-500bp"; else if(d<=1000)b="500bp-1kb"; else if(d<=3000)b="1-3kb";
        else if(d<=5000)b="3-5kb"; else if(d<=10000)b="5-10kb"; else b=">10kb"
        up=($NF<0)?"upstream":"downstream"; c[b]++ } END{
        split("0-500bp,500bp-1kb,1-3kb,3-5kb,5-10kb,>10kb",o,",")
        for(i=1;i<=6;i++){printf "    %-12s %d\n",o[i],c[o[i]]+0}
    }' "${NEAREST}"
fi
echo ""

# ─────────────── STEP 4: 三档靶基因 + 结合强度排序 ───────────────
GENE_LABEL="Target genes"; [[ "${PIPELINE}" == "atacseq" ]] && GENE_LABEL="Associated genes"
echo "[STEP 4] Export ${GENE_LABEL} (three selection rules)"

# (a) 启动子靶基因（最高置信）
PROMOTER_GENES="${ANNO_OUT}/promoter_target_genes.txt"
bedtools intersect -a "${PROMOTER_BED}" -b "${INPUT_BED}" -u \
    | awk '{print $4}' | sort -u > "${PROMOTER_GENES}"
PROMOTER_GENES_NAMED="${ANNO_OUT}/promoter_target_genes_named.tsv"
echo -e "gene_id\tgene_name" > "${PROMOTER_GENES_NAMED}"
bedtools intersect -a "${PROMOTER_BED}" -b "${INPUT_BED}" -u \
    | awk 'BEGIN{OFS="\t"} !s[$4]++{print $4,$7}' >> "${PROMOTER_GENES_NAMED}"

# (b) 距离过滤靶基因（最近 TSS 在阈值内）+ 结合强度排序
DIST_BP=$(( TARGET_DIST_KB * 1000 ))
DIST_GENES="${ANNO_OUT}/target_genes_within_${TARGET_DIST_KB}kb.txt"
DIST_RANKED="${ANNO_OUT}/target_genes_within_${TARGET_DIST_KB}kb_ranked.tsv"
if [[ -s "${NEAREST}" ]]; then
    # 每基因取其关联 peak 的最大 score 作为结合强度，按强度降序
    awk -v cut="${DIST_BP}" 'BEGIN{OFS="\t"}{
        d=($NF<0)?-$NF:$NF; if(d>cut) next
        gid=$9; if(gid==".") next
        if($5+0>best[gid]){best[gid]=$5+0; name[gid]=$12; dist[gid]=$NF}
    } END{ for(g in best) print g,name[g],best[g],dist[g] }' "${NEAREST}" \
        | sort -k3,3nr > "${DIST_RANKED}.body"
    echo -e "gene_id\tgene_name\tpeak_score\tdist_to_TSS" > "${DIST_RANKED}"
    cat "${DIST_RANKED}.body" >> "${DIST_RANKED}"; rm -f "${DIST_RANKED}.body"
    cut -f1 "${DIST_RANKED}" | tail -n +2 | sort -u > "${DIST_GENES}"
else
    touch "${DIST_GENES}"; echo -e "gene_id\tgene_name\tpeak_score\tdist_to_TSS" > "${DIST_RANKED}"
fi

# (c) 全部最近基因（仅供参考）
ALL_GENES="${ANNO_OUT}/all_nearest_genes.txt"
if [[ -s "${NEAREST}" ]]; then
    awk '$9!="."{print $9}' "${NEAREST}" | sort -u > "${ALL_GENES}"
else touch "${ALL_GENES}"; fi

echo "  (a) Promoter ${GENE_LABEL}:              $(wc -l < "${PROMOTER_GENES}")"
echo "  (b) ${GENE_LABEL} within ${TARGET_DIST_KB} kb of the nearest TSS: $(wc -l < "${DIST_GENES}")  (includes score-ranked table)"
echo "  (c) All nearest genes (reference):           $(wc -l < "${ALL_GENES}")"
echo ""

# ─────────────── STEP 5: 本地 GO 富集输入（不依赖 agriGO）───────────────
echo "[STEP 5] Generate GO enrichment inputs"
GENE2GO="${ANNO_OUT}/gene2go.tsv"
UNIVERSE="${ANNO_OUT}/background_universe.txt"
awk '{print $4}' "${GENES_BED}" | sort -u > "${UNIVERSE}"   # 背景 = 全基因组基因

if [[ -n "${ANNOTATION_INFO:-}" && -f "${ANNOTATION_INFO}" ]]; then
    # Phytozome annotation_info.txt：抽取 gene/transcript → GO 列
    python3 - "${ANNOTATION_INFO}" "${GENE2GO}" <<'PYEOF'
import sys, re
info, out = sys.argv[1], sys.argv[2]
with open(info) as f:
    header = f.readline().rstrip("\n").split("\t")
def col(*names):
    for n in names:
        for i,h in enumerate(header):
            if h.strip().lower()==n: return i
    return None
gi = col("locusname","gene","genename","locus")
ti = col("transcriptname","transcript","pacid")
goi = col("go")
pairs=set()
if goi is not None:
    with open(info) as f:
        next(f)
        for line in f:
            c=line.rstrip("\n").split("\t")
            if len(c)<=goi: continue
            gid = c[gi] if gi is not None and gi<len(c) and c[gi] else (c[ti] if ti is not None and ti<len(c) else None)
            if not gid: continue
            gid = re.sub(r'\.\d+$','',gid)   # 去transcripts后缀
            for go in re.findall(r'GO:\d+', c[goi]):
                pairs.add((gid,go))
with open(out,"w") as w:
    w.write("gene_id\tGO\n")
    for g,go in sorted(pairs): w.write(f"{g}\t{go}\n")
print(f"  gene2go mappings: {len(pairs)} mappings (genes: {len(set(p[0] for p in pairs))} genes)")
PYEOF
    echo "  Background universe:   $(($(wc -l < "${UNIVERSE}"))) genes → $(basename "${UNIVERSE}")"
    echo "  GO mappings:    $(basename "${GENE2GO}")"
    echo ""
    echo "  Downstream enrichment options (check species and identifier support):"
    echo "    · Local: hypergeometric test / R clusterProfiler(enricher, TERM2GENE=gene2go)"
    echo "    · Online: g:Profiler (select your species or use a custom GMT)"
    echo "    · Check species coverage before choosing an online enrichment service"
else
    log_warn "ANNOTATION_INFO is unset; cannot generate gene2go mappings"
    echo "  Set ANNOTATION_INFO to your Phytozome annotation_info.txt"
    echo "  Or use ${PROMOTER_GENES} with a compatible enrichment tool and matching background"
fi
echo ""

# ─────────────── STEP 6: 出图（matplotlib，可关）───────────────
if [[ "${MAKE_PLOTS}" == "true" ]] && python3 -c "import matplotlib" 2>/dev/null; then
    echo "[STEP 6] Generate plots"
    python3 - "${ANNO_OUT}" "${n_promoter}" "${n_exon}" "${n_intron}" "${n_intergenic}" "${NEAREST}" <<'PYEOF'
import sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out, npr, nex, nin, nig, nearest = sys.argv[1], *map(int, sys.argv[2:6]), sys.argv[6]
# 分布条形
cats=["Promoter","Exon","Intron","Intergenic"]; vals=[npr,nex,nin,nig]
plt.figure(figsize=(5,3.2))
plt.bar(cats, vals, color=["#2b8cbe","#74a9cf","#bdc9e1","#d9d9d9"])
plt.ylabel("peaks"); plt.title("Genomic distribution"); plt.tight_layout()
plt.savefig(f"{out}/genomic_distribution.png", dpi=200); plt.close()
# TSS 距离直方
try:
    d=[]
    with open(nearest) as f:
        for line in f:
            v=line.rstrip("\n").split("\t")[-1]
            try: x=int(v)
            except: continue
            if abs(x)<=10000: d.append(x)
    if d:
        plt.figure(figsize=(5,3.2))
        plt.hist(d, bins=40, color="#2b8cbe")
        plt.axvline(0, color="k", lw=0.8)
        plt.xlabel("distance to TSS (bp, - upstream)"); plt.ylabel("peaks")
        plt.title("Peak distance to nearest TSS"); plt.tight_layout()
        plt.savefig(f"{out}/tss_distance_hist.png", dpi=200); plt.close()
except FileNotFoundError: pass
print("  Generated: genomic_distribution.png, tss_distance_hist.png")
PYEOF
    echo ""
fi

# ─────────────── 总结 ───────────────
echo "============================================"
echo " Peak annotation summary (${PIPELINE})"
echo "============================================"
echo "  Input peaks:            ${N_PEAKS} (${PEAK_SOURCE})"
echo "  Promoter peaks:        ${n_promoter} ($(pct "${n_promoter}" "${N_PEAKS}")%)"
echo "  Promoter ${GENE_LABEL}: $(wc -l < "${PROMOTER_GENES}")"
echo "  ≤${TARGET_DIST_KB}kb ${GENE_LABEL}:    $(wc -l < "${DIST_GENES}")"
echo "  Output directory:              ${ANNO_OUT}"
echo "[DONE] $(date)"
