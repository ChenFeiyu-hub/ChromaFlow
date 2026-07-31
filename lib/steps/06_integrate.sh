#!/usr/bin/env bash
# ============================================================================
# 06_integrate.sh — motif×注释 候选表 + 自动波峰图 + 自动 GO 导出
# ============================================================================
# 环境: mamba activate ChromaFlow_downstream  (含 pyGenomeTracks + python)
# 用法: chromaflow candidates --type chipseq [--pval 1e-4] [--topn-tracks 15]
#       [--features "Promoter,1st Intron,Distal Intergenic"]
# ----------------------------------------------------------------------------
# 直接读 04(MEME-ChIP) 的 fimo_out_*/fimo.tsv（不重跑 FIMO）与 centrimo.tsv。
#   1) 每个 motif 的 FIMO 命中(基因组坐标) → 映射到 05 注释的 peak/基因
#   2) 按「合适位置(TARGET_FEATURES) + 含 motif」筛出候选
#   3) 每个 motif 独立输出：命中明细(位置/链/p/距TSS) + 基因列表
#   4) 汇总表整合 CentriMo(中心富集)信息
#   5) 自动导出候选基因列表(agriGO/g:Profiler) + 自动跑波峰图(07)
# ============================================================================
set -uo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true
FPVAL_CLI=""; TOPN_CLI=""; FEAT_CLI=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pval) FPVAL_CLI="$2"; shift 2 ;;
        --topn-tracks) TOPN_CLI="$2"; shift 2 ;;
        --features) FEAT_CLI="$2"; shift 2 ;;
        *) shift ;;
    esac
done
FPVAL="${FPVAL_CLI:-${FIMO_PVAL:-1e-4}}"
TOPN="${TOPN_CLI:-${TRACK_AUTO_TOPN:-15}}"
FEATURES="${FEAT_CLI:-${TARGET_FEATURES:-Promoter}}"

init_pipeline_dirs "${PIPELINE}"
open_log "06_integrate" "${PIPELINE}"
log_section "Step 06: 候选靶基因整合 (04 motif × 05 注释) — ${PIPELINE}"

resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
ANNO_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}"
MOTIF_DIR="${RESULTS_DIR}/downstream/${PIPELINE}/motif${PEAKSET_SUFFIX}"
OUT="$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/candidates${PEAKSET_SUFFIX}")"
GENES_BED="${ANNO_OUT}/genes.bed"
mkdir -p "${OUT}"

ANNO_TSV="${ANNO_OUT}/peak_annotation_full.tsv"
[[ -s "${ANNO_TSV}" ]] || die "缺 05 注释表: ${ANNO_TSV}（先跑 chromaflow annotate）"

MC_DIR=$(find "${MOTIF_DIR}" -type d -name "meme-chip_output" 2>/dev/null | head -1)
[[ -z "${MC_DIR}" ]] && MC_DIR="${MOTIF_DIR}"
mapfile -t FIMO_TSVS < <(find "${MC_DIR}" -path "*/fimo_out_*/fimo.tsv" 2>/dev/null | sort)
[[ ${#FIMO_TSVS[@]} -gt 0 ]] || die "未找到 04 的 fimo_out_*/fimo.tsv（先跑 chromaflow motif）: ${MC_DIR}"
CENTRIMO="${MC_DIR}/centrimo_out/centrimo.tsv"

echo "  05 注释:      ${ANNO_TSV}"
echo "  FIMO 文件:    ${#FIMO_TSVS[@]} 个"
echo "  CentriMo:     $([[ -s "${CENTRIMO}" ]] && echo 有 || echo 无)"
echo "  FIMO p 阈:    ${FPVAL}"
echo "  合适位置特征: ${FEATURES}"
echo ""

python3 - "${ANNO_TSV}" "${OUT}" "${FPVAL}" "${CENTRIMO}" "${GENES_BED}" "${FEATURES}" "${FIMO_TSVS[@]}" <<'PYEOF'
import sys, csv, re, collections, os
anno_f, outdir, pmax_s, centrimo, genes_bed, feat_s = sys.argv[1:7]
fimo_files = sys.argv[7:]
pmax = float(pmax_s)
target_features = [x.strip() for x in feat_s.split(",") if x.strip()]

def norm_feat(a):
    a = a or ""
    if a.startswith("Promoter"): return "Promoter"
    if "5' UTR" in a: return "5' UTR"
    if "3' UTR" in a: return "3' UTR"
    if a.startswith("1st Exon"): return "1st Exon"
    if a.startswith("Other Exon") or a.startswith("Exon"): return "Other Exon"
    if a.startswith("1st Intron"): return "1st Intron"
    if a.startswith("Other Intron") or a.startswith("Intron"): return "Other Intron"
    if a.startswith("Downstream"): return "Downstream"
    return "Distal Intergenic"

# ---- 基因 TSS（算 motif 到 TSS 距离）----
gene_tss = {}
if os.path.exists(genes_bed):
    for line in open(genes_bed):
        c = line.rstrip("\n").split("\t")
        if len(c) < 6: continue
        strand = c[5]; tss = int(c[1]) if strand == "+" else int(c[2]) - 1
        gene_tss[c[3]] = (c[0], tss, strand)
def dist_tss(gid, pos):
    if gid not in gene_tss: return ""
    _, tss, strand = gene_tss[gid]
    d = pos - tss
    return d if strand == "+" else -d

# ---- 05 注释 → peaks ----
anno_rows = list(csv.DictReader(open(anno_f), delimiter="\t"))
peaks = []; by_chr = collections.defaultdict(list)
for i, r in enumerate(anno_rows):
    ch = r.get("seqnames") or ""
    try: s = int(float(r.get("start"))) - 1; e = int(float(r.get("end")))
    except: continue
    pid = r.get("peak_id") or f"peak_{i}"
    rec = {"chr":ch,"s":s,"e":e,"pid":pid,"gene":r.get("geneId") or "",
           "feat":norm_feat(r.get("annotation","")),
           "dist":r.get("distanceToTSS",""),"score":r.get("peak_score","")}
    peaks.append(rec); by_chr[ch].append(rec)
for ch in by_chr: by_chr[ch].sort(key=lambda x:x["s"])

# ---- 读 FIMO 各 motif 命中 → 映射到 peak/基因 ----
motif_hits = collections.defaultdict(list)
motif_labels = []; sample = []; mode = "?"
for ff in fimo_files:
    rows = [l for l in open(ff) if l.strip() and not l.startswith("#")]
    if not rows: continue
    seen = set()
    for r in csv.DictReader(rows, delimiter="\t"):
        seq = (r.get("sequence_name") or "").strip()
        if not seq: continue
        try: p = float(r.get("p-value") or 1)
        except: p = 1
        if p > pmax: continue
        try: fs = int(float(r.get("start"))); fe = int(float(r.get("stop")))
        except: continue
        strand = r.get("strand","")
        m = re.match(r'^(.+):(\d+)-(\d+)', seq)
        if m: chrom=m.group(1); win=int(m.group(2)); gs=win+fs-1; ge=win+fe; mode="A(名字带坐标)"
        else: chrom=seq; gs=fs-1; ge=fe; mode="B(基因组坐标在start/stop)"
        label = (r.get("motif_alt_id") or r.get("motif_id") or "motif").strip()
        seen.add(label); mid = (gs+ge)//2
        for pk in by_chr.get(chrom, []):
            if pk["s"] > ge: break
            if gs < pk["e"] and pk["s"] < ge:
                motif_hits[label].append({"chr":chrom,"s":gs,"e":ge,"strand":strand,"p":p,
                    "gene":pk["gene"],"pid":pk["pid"],"feat":pk["feat"],
                    "dist":dist_tss(pk["gene"], mid)})
        if len(sample) < 3: sample.append(f"{seq}[{fs}-{fe}]")
    motif_labels += sorted(seen)
motif_labels = sorted(set(motif_labels))
total = sum(len(v) for v in motif_hits.values())
print(f"  motif 标签: {motif_labels}")
print(f"  FIMO 样例: {sample}")
print(f"  命中总数(p<{pmax}): {total}；坐标模式: {mode}")

# ---- CentriMo：每 motif 是否中心富集 ----
centrimo_info = {}; centrimo_top = ""
if centrimo and os.path.exists(centrimo):
    crows = [l for l in open(centrimo) if l.strip() and not l.startswith("#")]
    if crows:
        for j, r in enumerate(csv.DictReader(crows, delimiter="\t")):
            mid = (r.get("motif_alt_id") or r.get("motif_id") or "").strip()
            ev = r.get("E-value") or r.get("adj_p-value") or ""
            binl = r.get("bin_location") or ""
            if mid:
                centrimo_info[mid] = (ev, binl)
                if j == 0: centrimo_top = f"{mid} (E={ev}, bin={binl})"

# ---- 每个 motif 独立输出：命中明细 + 基因列表 ----
def safe(x): return re.sub(r'[^A-Za-z0-9._-]','_', x)
for label in motif_labels:
    hits = motif_hits[label]
    with open(f"{outdir}/motif_{safe(label)}_hits.tsv","w") as w:
        w.write("gene_id\tpeak_id\tfeature\tmotif_chr\tmotif_start\tmotif_end\tstrand\tp_value\tdist_to_TSS\n")
        for h in sorted(hits, key=lambda x:(str(x["gene"]), x["s"])):
            w.write(f'{h["gene"]}\t{h["pid"]}\t{h["feat"]}\t{h["chr"]}\t{h["s"]}\t{h["e"]}\t{h["strand"]}\t{h["p"]:.2e}\t{h["dist"]}\n')
    genes = sorted(set(h["gene"] for h in hits if h["gene"] and h["feat"] in target_features))
    open(f"{outdir}/motif_{safe(label)}_genes.txt","w").write("\n".join(genes)+("\n" if genes else ""))

# ---- 每 peak 的 motif 计数 ----
peak_cnt = collections.defaultdict(collections.Counter)
for label, hits in motif_hits.items():
    for h in hits: peak_cnt[h["pid"]][label] += 1

# ---- 候选表（每 peak 一行：peak 信息 + 各 motif 计数 + 判定）----
def score_of(pk):
    try: return float(pk["score"])
    except: return 0.0
rows_out = []; gene_best = {}
for pk in peaks:
    cnt = peak_cnt.get(pk["pid"], {})
    ncols = [str(cnt.get(l,0)) for l in motif_labels]
    any_m = 1 if sum(cnt.values()) > 0 else 0
    passf = pk["feat"] in target_features
    cand = 1 if (passf and any_m) else 0
    sc = score_of(pk)
    rows_out.append((pk, ncols, any_m, passf, cand, sc))
    if pk["gene"] and cand and sc > gene_best.get(pk["gene"], -1):
        gene_best[pk["gene"]] = sc
rows_out.sort(key=lambda x:(-x[4], -x[3], -x[5]))
with open(f"{outdir}/candidate_targets.tsv","w") as w:
    hdr = ["gene_id","peak_id","peak_chr","peak_start","peak_end","feature","dist_TSS_peak","peak_score"] \
          + [f"n_{l}" for l in motif_labels] + ["any_motif","pass_feature","candidate"]
    w.write("\t".join(hdr)+"\n")
    for pk,ncols,any_m,passf,cand,sc in rows_out:
        w.write("\t".join([pk["gene"],pk["pid"],pk["chr"],str(pk["s"]),str(pk["e"]),pk["feat"],
                str(pk["dist"]),str(pk["score"])]+ncols+[str(any_m),str(int(passf)),str(cand)])+"\n")

# ---- 高置信 = 候选基因（合适位置+motif），按强度降序 ----
hi = sorted(gene_best, key=lambda g:-gene_best[g])
open(f"{outdir}/high_confidence_genes.txt","w").write("\n".join(hi)+("\n" if hi else ""))

# ---- motif 汇总表（整合 CentriMo）----
with open(f"{outdir}/motif_summary.tsv","w") as w:
    w.write("motif\tn_hits\tn_peaks\tn_genes_target\tcentrimo_Evalue\tcentrimo_bin\tcentrally_enriched\n")
    for label in motif_labels:
        hits = motif_hits[label]
        npk = len(set(h["pid"] for h in hits))
        ng = len(set(h["gene"] for h in hits if h["gene"] and h["feat"] in target_features))
        ev, binl = centrimo_info.get(label, ("",""))
        w.write(f"{label}\t{len(hits)}\t{npk}\t{ng}\t{ev}\t{binl}\t{'yes' if ev else 'no'}\n")

# ---- 小结 ----
featdist = collections.Counter(pk["feat"] for pk in peaks)
print(f"  特征分布: {dict(featdist)}")
print(f"  合适位置: {target_features}")
for label in motif_labels:
    ng = len(set(h["gene"] for h in motif_hits[label] if h["gene"] and h["feat"] in target_features))
    ev,_ = centrimo_info.get(label, ("",""))
    tag = f"CentriMo中心富集(E={ev})" if ev else "非中心富集"
    print(f"    {label}: {len(motif_hits[label])} 命中 → {ng} 目标基因  [{tag}]")
print(f"  CentriMo 最中心: {centrimo_top or '(无)'}")
print(f"  候选(合适位置+motif)基因: {len(hi)}")
PYEOF

CAND="${OUT}/candidate_targets.tsv"
HICONF="${OUT}/high_confidence_genes.txt"
[[ -s "${CAND}" ]] || die "候选表生成失败"

# ─────────────── 自动: 导出 agriGO/g:Profiler 基因列表 ───────────────
echo ""
echo "[GO 导出] 生成 agriGO/g:Profiler 输入"
GO_OUT="${OUT}/go"; mkdir -p "${GO_OUT}"
if [[ -s "${HICONF}" ]]; then
    sort -u "${HICONF}" > "${GO_OUT}/target_candidates.raw.txt"
    sed -E 's/\.[0-9]+$//' "${GO_OUT}/target_candidates.raw.txt" | sort -u > "${GO_OUT}/target_candidates.stripped.txt"
fi
if [[ -s "${ANNO_OUT}/background_universe.txt" ]]; then
    cp "${ANNO_OUT}/background_universe.txt" "${GO_OUT}/background.raw.txt"
    sed -E 's/\.[0-9]+$//' "${ANNO_OUT}/background_universe.txt" | sort -u > "${GO_OUT}/background.stripped.txt"
fi
cat > "${GO_OUT}/README_GO.txt" <<EOF
GO 富集（外部执行）— 候选靶基因（合适位置+motif）
前景: target_candidates.stripped.txt
背景: background.stripped.txt (务必用作自定义背景)
陆生植物→agriGO；莱茵衣藻等→g:Profiler
回来后: chromaflow go-plot --input <结果> --type ${PIPELINE}
EOF
echo "  → ${GO_OUT}/target_candidates.stripped.txt"

echo ""
echo "============================================"
echo " 整合完成 (${PIPELINE})"
echo "============================================"
echo "  候选表(逐peak):  ${CAND}"
echo "  motif 汇总:      ${OUT}/motif_summary.tsv"
echo "  各 motif 明细:   ${OUT}/motif_<名>_hits.tsv / _genes.txt"
echo "  高置信基因:      ${HICONF}"
echo "  GO 导出:         ${GO_OUT}/"
echo "  (波峰图请单独跑: chromaflow tracks --gene-list ${HICONF})"
echo "[DONE] $(date)"
