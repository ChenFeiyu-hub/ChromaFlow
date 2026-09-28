#!/usr/bin/env python3
"""Candidate integration using shared motif evidence filters."""
import sys, csv, re, collections, os
from pathlib import Path
from annotation_policy import norm_feat, target_features as parse_features, location_rule
from motif_selection import Evidence, number
anno_f, outdir, pmax_s, centrimo, genes_bed, feat_s, maxtss_s = sys.argv[1:8]
fimo_files = sys.argv[8:]
max_tss_dist = int(maxtss_s or 0)
pmax = float(pmax_s)
target_features = parse_features(feat_s)
if max_tss_dist < 0: sys.exit('[ERROR] TARGET_MAX_TSS_DIST must be nonnegative')

# ---- 基因 TSS（算 motif 到 TSS 距离）----
# 坐标一律走 int(float(...))：老版本 05 写出的 genes.bed 里，整千万级坐标可能
# 是 "8e+06" 这种科学计数（R 的 write.table 对 double 用 as.character 所致），
# 直接 int() 会 ValueError。新版 05 已修，这里保持向后兼容。
def coord(x):
    return int(float(x))

gene_tss = {}
bad_lines = 0
if os.path.exists(genes_bed):
    for line in open(genes_bed):
        c = line.rstrip("\n").split("\t")
        if len(c) < 6: continue
        try:
            strand = c[5]; tss = coord(c[1]) if strand == "+" else coord(c[2]) - 1
        except (ValueError, TypeError):
            bad_lines += 1; continue
        gene_tss[c[3]] = (c[0], tss, strand)
    if bad_lines:
        print(f"  [WARN] Skipped {bad_lines} genes.bed rows with invalid coordinates")
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
    try: dnum = int(float(r.get("distanceToTSS") or "nan"))
    except: dnum = None
    # feature_geneId: peak 物理落在哪个基因的结构里（05 新增；老表则现场解析）
    fg = r.get("feature_geneId") or ""
    # R's write.table uses literal NA for missing character values. The
    # authoritative consistency flag distinguishes this from a real ID "NA".
    consistent = (r.get("gene_consistent") or "").strip().upper()
    if fg == "NA" and consistent == "TRUE":
        fg = ""
    if not fg:
        m = re.match(r"^.*\([^/]*/([^,]+),.*$", r.get("annotation","") or "")
        fg = m.group(1) if m else ""
    rec = {"chr":ch,"s":s,"e":e,"pid":pid,"gene":r.get("geneId") or "",
           "feat":norm_feat(r.get("annotation","")),
           "fgene":fg,
           "gene_consistent":consistent,
           "dist":r.get("distanceToTSS",""),"dnum":dnum,
           "score":r.get("peak_score","")}
    peaks.append(rec); by_chr[ch].append(rec)
for ch in by_chr: by_chr[ch].sort(key=lambda x:x["s"])

# Shared filtering happens before peak ranking or gene selection.
evidence = Evidence.from_env(fimo_files, pmax)
motif_hits = collections.defaultdict(list)
for h, pk in evidence.mapped(by_chr):
    motif_hits[h['motif_label']].append({"chr":h['chrom'],"s":h['start'],"e":h['end'],
        "strand":h['strand'],"p":h['p_value'],"gene":pk['gene'],"pid":pk['pid'],"feat":pk['feat'],
        "dist":dist_tss(pk['gene'], (h['start']+h['end'])//2),
        "pass_feature":location_rule(pk['feat'], pk['dnum'], target_features, max_tss_dist)[0],
        "evidence":dict(gene_id=pk['gene'], peak_id=pk['pid'], **h)})
motif_labels = sorted({k[1] for k in evidence.selected_models})
min_signal = number(os.environ.get('CFSEL_MIN_SIGNAL', '0'), 'minimum signal')
centrimo_info = {key[1]:evidence.models[key] for key in evidence.selected_models}

# Remove only owned per-motif exports from a previous filter in this directory.
for pattern in ('motif_*_hits.tsv', 'motif_*_genes.txt'):
    for stale in Path(outdir).glob(pattern):
        stale.unlink()
# ---- 每个 motif 独立输出：命中明细 + 基因列表 ----
signal_by_peak = {pk['pid']: number(pk['score'] or 0, 'peak signal') for pk in peaks}
def safe(x): return re.sub(r'[^A-Za-z0-9._-]','_', x)
for label in motif_labels:
    hits = motif_hits[label]
    with open(f"{outdir}/motif_{safe(label)}_hits.tsv","w") as w:
        w.write("gene_id\tpeak_id\tfeature\tmotif_chr\tmotif_start\tmotif_end\tstrand\tp_value\tdist_to_TSS\n")
        for h in sorted(hits, key=lambda x:(str(x["gene"]), x["s"])):
            w.write(f'{h["gene"]}\t{h["pid"]}\t{h["feat"]}\t{h["chr"]}\t{h["s"]}\t{h["e"]}\t{h["strand"]}\t{h["p"]:.2e}\t{h["dist"]}\n')
    genes = sorted(set(h["gene"] for h in hits if h["gene"] and h["pass_feature"] and signal_by_peak[h["pid"]] >= min_signal))
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
rescued = []; reasons = {}
for pk in peaks:
    cnt = peak_cnt.get(pk["pid"], {})
    ncols = [str(cnt.get(l,0)) for l in motif_labels]
    any_m = 1 if sum(cnt.values()) > 0 else 0
    # 位置判定。默认按 annotation 的结构特征筛。
    # 但 annotation 描述的可能是【邻近基因】的结构（peak 落在 A 的内含子里，
    # 而最近的 TSS 属于 B），此时按特征筛会漏掉「B 的启动子上游 peak」。
    # max_tss_dist > 0 时启用距离补救：只要 peak 在 geneId 的 TSS 上游该距离内，
    # 无论它压在谁的基因体里都算合适位置。
    passf, reason = location_rule(pk["feat"], pk["dnum"], target_features, max_tss_dist)
    if reason.startswith('upstream<='): rescued.append(pk)
    sc = score_of(pk)
    cand = 1 if (passf and any_m and sc >= min_signal) else 0
    reasons[pk["pid"]] = reason
    rows_out.append((pk, ncols, any_m, passf, cand, sc))
    if pk["gene"] and cand and sc > gene_best.get(pk["gene"], -1):
        gene_best[pk["gene"]] = sc
rows_out.sort(key=lambda x:(-x[4], -x[3], -x[5]))
with open(f"{outdir}/candidate_targets.tsv","w") as w:
    hdr = ["gene_id","peak_id","peak_chr","peak_start","peak_end","feature","feature_gene",
           "dist_TSS_peak","peak_score"] \
          + [f"n_{l}" for l in motif_labels] + ["any_motif","pass_feature","pass_reason","candidate"]
    w.write("\t".join(hdr)+"\n")
    for pk,ncols,any_m,passf,cand,sc in rows_out:
        w.write("\t".join([pk["gene"],pk["pid"],pk["chr"],str(pk["s"]),str(pk["e"]),pk["feat"],
                pk["fgene"] or "-", str(pk["dist"]),str(pk["score"])]
                +ncols+[str(any_m),str(int(passf)),reasons[pk["pid"]] or "-",str(cand)])+"\n")

# ---- 高置信 = 候选基因（合适位置+motif），按强度降序 ----
hi = sorted(gene_best, key=lambda g:-gene_best[g])
open(f"{outdir}/high_confidence_genes.txt","w").write("\n".join(hi)+("\n" if hi else ""))

# Correct positional-enrichment labels; local mode is not central enrichment.
with open(f"{outdir}/motif_summary.tsv","w") as w:
    w.write("motif\tn_hits\tn_peaks\tn_genes_target\tcentrimo_Evalue\tcentrimo_bin\tcentrimo_mode\tpositional_enrichment\tcentrally_enriched\n")
    for label in motif_labels:
        hits = motif_hits[label]; m = centrimo_info[label]
        npk = len({h['pid'] for h in hits})
        ng = len({h['gene'] for h in hits if h['gene'] and h['pass_feature'] and signal_by_peak[h['pid']] >= min_signal})
        w.write('\t'.join(map(str,[label,len(hits),npk,ng,m['centrimo_Evalue'],m['centrimo_bin'],m['centrimo_mode'],m['positional_enrichment'],m['centrally_enriched']]))+'\n')
selected_peaks = {pk['pid'] for pk,_,_,_,cand,_ in rows_out if cand}
selected = [h['evidence'] for hits in motif_hits.values() for h in hits if h['pid'] in selected_peaks]
parameters = dict(features=feat_s, min_signal=min_signal, max_tss_dist=max_tss_dist,
    promoter_up=int(os.environ['CFSEL_UP']), promoter_down=int(os.environ['CFSEL_DOWN']), level='candidate')
evidence.save(outdir, hi, selected, anno_f, os.environ['CFSEL_PEAKS'], parameters)
open(f"{outdir}/candidate_genes.txt","w").write(''.join(g+'\n' for g in hi))

# ---- 小结 ----
featdist = collections.Counter(pk["feat"] for pk in peaks)
print(f"  Feature distribution: {dict(featdist)}")
print(f"  Accepted locations: {target_features}")
n_mis = sum(1 for pk in peaks if (pk["gene_consistent"] == "FALSE" or
            (pk["gene_consistent"] not in ("TRUE", "FALSE") and
             pk["fgene"] and pk["fgene"] != pk["gene"])))
if n_mis:
    print(f"  NOTE: {n_mis} peaks overlap a different gene from the nearest-TSS gene (see feature_gene)")
if max_tss_dist > 0:
    print(f"  Distance rule: include peaks <= {max_tss_dist} bp upstream of TSS; adds {len(rescued)} peaks")
elif n_mis:
    print(f"  Review these loci before using --max-tss-dist or changing PROMOTER_UP and rerunning annotation")
for label in motif_labels:
    ng = len(set(h["gene"] for h in motif_hits[label] if h["gene"] and h["pass_feature"] and signal_by_peak[h["pid"]] >= min_signal))
    m = centrimo_info[label]
    print(f"    {label}: {len(motif_hits[label])} hits -> {ng} target genes [CentriMo {m['centrimo_mode']}: {m['positional_enrichment']}; E={m['centrimo_Evalue'] or 'NA'}; report threshold={evidence.display_max}]")
print(f"  Candidate genes (accepted location + motif): {len(hi)}")
