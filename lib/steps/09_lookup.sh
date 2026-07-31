#!/usr/bin/env bash
# ============================================================================
# 09_lookup.sh — 先验基因反查（独立分析：直接从 05 注释 + 04 FIMO 自己判）
# ============================================================================
# 环境: mamba activate ChromaFlow_downstream  (出图需 pyGenomeTracks)
# 用法: chromaflow lookup --type chipseq --gene-list phosphate_genes.txt [--plot]
#       chromaflow lookup --type chipseq --gene-list X --peakset consensus
# ----------------------------------------------------------------------------
# 设计初衷：即使 candidate(06) 整体没筛出预期基因，也能拿先验列表独立诊断。
# 因此 lookup 不依赖 candidate 表，自己复刻"位置 + FIMO motif"判据，
# 对每个查询基因给出【可解释的证据分解】：有无 peak / 位置 / signalValue /
# 有无 FIMO motif / 哪个 motif → 你据此可视化时能对上"为什么是候选"。
# 判据与 candidate 一致：ChIPseeker 归属 + 位置(TARGET_FEATURES) + FIMO motif。
# ============================================================================
set -uo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true
GENE_LIST=""; DO_PLOT=false; OUTDIR_ARG=""; MIN_SIGNAL=""; LEVEL=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gene-list)  GENE_LIST="$2"; shift 2 ;;
        --plot)       DO_PLOT=true; shift ;;
        --min-signal) MIN_SIGNAL="$2"; shift 2 ;;
        --level)      LEVEL="$2"; shift 2 ;;
        --outdir)     OUTDIR_ARG="$2"; shift 2 ;;
        *) shift ;;
    esac
done
[[ -s "${GENE_LIST}" ]] || die "需 --gene-list 指向先验基因列表"
MIN_SIGNAL="${MIN_SIGNAL:-0}"
LEVEL="${LEVEL:-all}"
FPVAL="${FIMO_PVAL:-1e-4}"
FEATURES="${TARGET_FEATURES:-Promoter}"

init_pipeline_dirs "${PIPELINE}"
open_log "09_lookup" "${PIPELINE}"
log_section "Step 09: 先验基因反查 (${PIPELINE})"

resolve_peakset "${PIPELINE}" "${CHROMAFLOW_PEAKSET:-idr}"
ANNO_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/annotation${PEAKSET_SUFFIX}"
MOTIF_DIR="${RESULTS_DIR}/downstream/${PIPELINE}/motif${PEAKSET_SUFFIX}"
OUT="${OUTDIR_ARG:-$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/lookup${PEAKSET_SUFFIX}")}"
mkdir -p "${OUT}"

ANNO_TSV="${ANNO_OUT}/peak_annotation_full.tsv"
[[ -s "${ANNO_TSV}" ]] || die "缺 05 注释 ${ANNO_TSV}。请先: chromaflow annotate --type ${PIPELINE} ${CHROMAFLOW_PEAKSET:+--peakset ${CHROMAFLOW_PEAKSET}}"

# 04 的 FIMO（独立读，不经 candidate）；缺失则退化为"仅位置"分析
MC_DIR=$(find "${MOTIF_DIR}" -type d -name "meme-chip_output" 2>/dev/null | head -1)
[[ -z "${MC_DIR}" ]] && MC_DIR="${MOTIF_DIR}"
mapfile -t FIMO_TSVS < <(find "${MC_DIR}" -path "*/fimo_out_*/fimo.tsv" 2>/dev/null | sort)

echo "  peakset:  ${CHROMAFLOW_PEAKSET:-idr}  → $(basename "${PEAKSET_FILE}")"
echo "  注释:     ${ANNO_TSV}"
echo "  FIMO:     ${#FIMO_TSVS[@]} 个 $([[ ${#FIMO_TSVS[@]} -eq 0 ]] && echo '(无 → 仅位置分析，先跑 motif 可补 motif 证据)')"
echo "  查询列表: ${GENE_LIST}"
echo "  判据:     ChIPseeker 归属 + 位置(${FEATURES}) + FIMO motif(p<${FPVAL})"
echo "  min-signal: ${MIN_SIGNAL}$([[ "${MIN_SIGNAL}" == "0" ]] && echo ' (不过滤)')   level: ${LEVEL}"
echo ""

python3 - "${ANNO_TSV}" "${GENE_LIST}" "${OUT}" "${FEATURES}" "${FPVAL}" "${MIN_SIGNAL}" "${LEVEL}" "${FIMO_TSVS[@]}" <<'PYEOF'
import sys, csv, re, collections, os
anno_f, query_f, outdir, feat_s, pmax_s, minsig_s, level = sys.argv[1:8]
fimo_files = sys.argv[8:]
target_features = [x.strip() for x in feat_s.split(",") if x.strip()]
pmax = float(pmax_s); min_signal = float(minsig_s)

def norm_feat(a):
    a=a or ""
    if a.startswith("Promoter"): return "Promoter"
    if "5' UTR" in a: return "5' UTR"
    if "3' UTR" in a: return "3' UTR"
    if a.startswith("1st Exon"): return "1st Exon"
    if a.startswith("Other Exon") or a.startswith("Exon"): return "Other Exon"
    if a.startswith("1st Intron"): return "1st Intron"
    if a.startswith("Other Intron") or a.startswith("Intron"): return "Other Intron"
    if a.startswith("Downstream"): return "Downstream"
    return "Distal Intergenic"
# 只剥离「版本后缀」(如 .1 / .v6.1)，绝不能按第一个点切分：
# Phytozome ID 形如 Cre02.g075050_4532，按点切会塌缩成 "Cre02"（整条染色体同键）。
def strip(x): return re.sub(r'\.(v[\d.]+|\d+)$', '', (x or "").strip())

# 查询基因
query=[]
for line in open(query_f):
    g=line.strip()
    if g and not g.startswith("#"): query.append(g.split()[0])
query=list(dict.fromkeys(query))

# 05 注释 → peaks；建 by_chr(映射FIMO) 与 gene→peaks
rows=list(csv.DictReader(open(anno_f), delimiter="\t"))
peaks=[]; by_chr=collections.defaultdict(list)
by_gene=collections.defaultdict(list); by_bare=collections.defaultdict(list)
bare_ids=collections.defaultdict(set)   # 归一化键 → 原始基因ID集合（用于检测跨基因碰撞）
for i,r in enumerate(rows):
    ch=r.get("seqnames") or ""
    try: s=int(float(r.get("start")))-1; e=int(float(r.get("end")))
    except: continue
    def num(x):
        try: return float(x)
        except: return 0.0
    rec={"chr":ch,"s":s,"e":e,"pid":r.get("peak_id") or f"peak_{i}",
         "gene":r.get("geneId") or "","feat":norm_feat(r.get("annotation","")),
         "dist":r.get("distanceToTSS",""),"sig":num(r.get("peak_score")),
         "idr":r.get("peak_globalIDR","")}
    peaks.append(rec); by_chr[ch].append(rec)
    if rec["gene"]:
        by_gene[rec["gene"]].append(rec); by_bare[strip(rec["gene"])].append(rec)
        bare_ids[strip(rec["gene"])].add(rec["gene"])
for ch in by_chr: by_chr[ch].sort(key=lambda x:x["s"])

# 全基因组基因ID集合（来自 05 的 genes.bed）：用于区分
#   "该基因存在但无 peak"  vs  "该ID在注释中根本不存在(多半是ID格式/物种不匹配)"
all_genes=set(); all_genes_bare=set()
genes_bed=os.path.join(os.path.dirname(anno_f), "genes.bed")
if os.path.exists(genes_bed):
    for line in open(genes_bed):
        c=line.rstrip("\n").split("\t")
        if len(c)>=4 and c[3]:
            all_genes.add(c[3]); all_genes_bare.add(strip(c[3]))

# 04 FIMO → 每 peak 含哪些 motif（基因组坐标映射到 peak，与 06 同法）
peak_motifs=collections.defaultdict(set)
mode="?"
for ff in fimo_files:
    rws=[l for l in open(ff) if l.strip() and not l.startswith("#")]
    if not rws: continue
    for r in csv.DictReader(rws, delimiter="\t"):
        seq=(r.get("sequence_name") or "").strip()
        if not seq: continue
        try: p=float(r.get("p-value") or 1)
        except: p=1
        if p>pmax: continue
        try: fs=int(float(r.get("start"))); fe=int(float(r.get("stop")))
        except: continue
        m=re.match(r'^(.+):(\d+)-(\d+)', seq)
        if m: chrom=m.group(1); win=int(m.group(2)); gs=win+fs-1; ge=win+fe; mode="A"
        else: chrom=seq; gs=fs-1; ge=fe; mode="B"
        label=(r.get("motif_alt_id") or r.get("motif_id") or "motif").strip()
        for pk in by_chr.get(chrom, []):
            if pk["s"]>ge: break
            if gs<pk["e"] and pk["s"]<ge: peak_motifs[pk["pid"]].add(label)

# 每查询基因：证据分解
def tier_of(in_target, has_motif, has_peak):
    if not has_peak: return "none"
    if in_target and has_motif: return "candidate"
    if in_target: return "in_target"
    if has_motif: return "motif_only"
    return "peak_only"
tier_rank={"candidate":0,"in_target":1,"motif_only":2,"peak_only":3,"none":4}

out_rows=[]; n_peak=n_target=n_cand=0
bound_genes=[]; cand_genes=[]
n_exact=n_fallback=n_ambig=0
for q in query:
    recs = by_gene.get(q)
    if recs: n_exact+=1
    else:
        # 版本号归一后再试（双向都可能带/不带版本号，故不比较 qn 与 q）
        qn = strip(q)
        cand_recs = by_bare.get(qn)
        if cand_recs and len(bare_ids.get(qn, set())) > 1:
            # 归一化后仍对应多个不同基因 → 有歧义，不猜测
            cand_recs = None; n_ambig += 1
        recs = cand_recs
        if recs: n_fallback+=1
    recs = [r for r in (recs or []) if r["sig"]>=min_signal]
    if not recs:
        out_rows.append([q,"no","","-","no","","","no","-","none"]); continue
    # 每个 peak 的证据
    evald=[]
    for r in recs:
        mset=peak_motifs.get(r["pid"], set())
        in_t=r["feat"] in target_features
        hm=len(mset)>0
        evald.append((r, in_t, hm, mset, tier_of(in_t,hm,True)))
    # 取最佳 peak：tier 优先，其次 signalValue
    best=min(evald, key=lambda x:(tier_rank[x[4]], -x[0]["sig"]))
    r,in_t,hm,mset,tier=best
    motifs=",".join(sorted(mset)) if mset else "-"
    out_rows.append([q,"yes",r["pid"],r["feat"],"yes" if in_t else "no",
                     f'{r["sig"]:g}',str(r["idr"]),"yes" if hm else "no",motifs,tier])
    n_peak+=1
    if in_t: n_target+=1
    if in_t and hm: n_cand+=1
    bound_genes.append(q)
    if tier=="candidate": cand_genes.append(q)

# 输出（--level 聚焦）
def keep(row):
    if level=="candidate": return row[9]=="candidate"
    if level=="in_target": return row[4]=="yes"
    return True
hdr=["gene_id","has_peak","peak_id","feature","in_target","signalValue","globalIDR","has_motif","motifs","tier"]
with open(f"{outdir}/query_binding.tsv","w") as w:
    w.write("\t".join(hdr)+"\n")
    for row in sorted(out_rows, key=lambda r:(tier_rank[r[9]], -(float(r[5]) if r[5] else 0))):
        if keep(row): w.write("\t".join(row)+"\n")
open(f"{outdir}/bound_genes.txt","w").write("\n".join(bound_genes)+("\n" if bound_genes else ""))
open(f"{outdir}/candidate_genes.txt","w").write("\n".join(cand_genes)+("\n" if cand_genes else ""))

nq=len(query)
print(f"  查询 {nq} 基因（min-signal={min_signal:g}，FIMO 坐标模式 {mode}）：")
# 区分「ID 不在基因组注释中」(格式问题) 与「基因存在但无 peak」(真实生物学结果)
not_in_genome=[]
if all_genes:
    for q in query:
        if q not in all_genes and strip(q) not in all_genes_bare: not_in_genome.append(q)
print(f"    ID 匹配: 精确 {n_exact} / 版本号归一后 {n_fallback}"
      + (f" / 归一后有歧义(未采用) {n_ambig}" if n_ambig else ""))
if all_genes:
    n_absent=len(not_in_genome)
    print(f"    ID 校验: {nq-n_absent}/{nq} 存在于基因组注释；{n_absent} 个不存在")
    if n_absent:
        print(f"    ⚠ 以下 ID 在注释中找不到（多半是ID格式或物种不匹配，非'无结合'）: "
              f"{', '.join(not_in_genome[:5])}{' ...' if n_absent>5 else ''}")
        open(f"{outdir}/id_not_found.txt","w").write("\n".join(not_in_genome)+"\n")
print(f"    ① 有 peak(ChIPseeker 归属):        {n_peak}")
print(f"    ② 位置合适(∈TARGET_FEATURES):      {n_target}")
print(f"    ③ 候选级(位置合适 + FIMO motif):   {n_cand}  ← 最可信")
# 证据分解速览：各 tier 计数（帮助解释"为什么/为什么不是候选"）
tc=collections.Counter(r[9] for r in out_rows)
print(f"  证据分解: 候选级={tc['candidate']}  仅位置(无motif)={tc['in_target']}  "
      f"仅motif(位置不合适)={tc['motif_only']}  仅有peak={tc['peak_only']}  无peak={tc['none']}")
if cand_genes:
    print("  候选级基因(位置+motif，前 12):")
    for row in sorted([r for r in out_rows if r[9]=="candidate"], key=lambda r:-(float(r[5]) if r[5] else 0))[:12]:
        print(f"    {row[0]:22} feat={row[3]:16} sig={row[5]:>9} motif={row[8]}")
else:
    print("  ⚠ 无候选级基因：查询基因的 peak 要么位置不合适、要么无 FIMO motif（见 query_binding.tsv 的 tier 列）")
PYEOF

echo ""
echo "  产出: ${OUT}/query_binding.tsv（证据分解：has_peak/feature/in_target/signalValue/has_motif/motifs/tier）"
echo "        bound_genes.txt, candidate_genes.txt"

# 可选出图（优先候选级）
if [[ "${DO_PLOT}" == "true" ]]; then
    PLOT_LIST="${OUT}/candidate_genes.txt"; [[ -s "${PLOT_LIST}" ]] || PLOT_LIST="${OUT}/bound_genes.txt"
    if [[ -s "${PLOT_LIST}" ]] && command -v pyGenomeTracks &>/dev/null; then
        echo ""; echo "[波峰图] $(basename "${PLOT_LIST}") 中的基因"
        BW_DIR="${RESULTS_DIR}/${PIPELINE}/bwa/merged_library/bigwig"; [[ -d "${BW_DIR}" ]] || BW_DIR="${RESULTS_DIR}/${PIPELINE}"
        bash "${CHROMAFLOW_HOME}/lib/steps/plot_tracks.sh" \
            --gene-list "${PLOT_LIST}" --genes-bed "${ANNO_OUT}/genes.bed" --gtf "${GENOME_GTF}" \
            --bigwig-dir "${BW_DIR}" --peaks "${PEAKSET_FILE}" \
            --outdir "${OUT}/tracks" --flank "${TRACK_FLANK}" ${TRACK_MAX:+--ymax "${TRACK_MAX}"} || echo "  绘图返回非零"
    fi
fi
echo "[DONE] $(date)"
