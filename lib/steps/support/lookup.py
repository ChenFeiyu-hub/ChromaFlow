#!/usr/bin/env python3
"""Exact-ID lookup with explicit namespace confirmation and regional evidence."""
import sys, csv, collections, os, math
from pathlib import Path
from annotation_policy import norm_feat, target_features as parse_features, location_rule, check_policy
from promoter_audit import read_genes, read_query, audit
from fimo_tracks import read_sizes
from motif_selection import Evidence
from peak_metadata import lookup_idr_values
anno_f, query_f, outdir, feat_s, pmax_s, minsig_s, level = sys.argv[1:8]
up, down, maxtss, confirmed, peak_file, chrom_sizes = sys.argv[8:14]
up, down, maxtss = int(up), int(down), int(maxtss)
fimo_files = sys.argv[14:]
status = Path(outdir)/'lookup.status'
status.write_text('checking_ids\n')
check_policy(anno_f, up, down)
target_features = parse_features(feat_s)
pmax = float(pmax_s); min_signal = float(minsig_s)

if not 0 < pmax <= 1 or not math.isfinite(min_signal) or min_signal < 0 or maxtss < 0:
    sys.exit('[ERROR] Invalid p-value, minimum signal or upstream-distance threshold')
if level not in ('all','candidate','in_target'): sys.exit('[ERROR] Invalid lookup level')
query = read_query(query_f)

# 05 注释 → peaks；建 by_chr(映射FIMO) 与 gene→peaks
rows=list(csv.DictReader(open(anno_f), delimiter="\t"))
idr_values = lookup_idr_values(rows, os.environ.get('PEAKSET_KIND', 'unknown'), peak_file)
if os.environ.get('PEAKSET_KIND', 'unknown') != 'idr':
    print('[METADATA] globalIDR=NA: this input is not an explicitly selected IDR peak set.')
peaks=[]; by_chr=collections.defaultdict(list)
by_gene=collections.defaultdict(list)
for i,r in enumerate(rows):
    ch=r.get("seqnames") or ""
    try: s=int(float(r.get("start")))-1; e=int(float(r.get("end")))
    except: continue
    def num(x):
        try: value = float(x); return value if math.isfinite(value) else 0.0
        except (ValueError, TypeError): return 0.0
    rec={"chr":ch,"s":s,"e":e,"pid":r.get("peak_id") or f"peak_{i}",
         "gene":r.get("geneId") or "","feat":norm_feat(r.get("annotation","")),
         "dist":r.get("distanceToTSS",""),"sig":num(r.get("peak_score")),
         "idr":idr_values[i]}
    peaks.append(rec); by_chr[ch].append(rec)
    if rec["gene"]:
        by_gene[rec["gene"]].append(rec)
for ch in by_chr: by_chr[ch].sort(key=lambda x:x["s"])

genes_bed = os.path.join(os.path.dirname(anno_f), 'genes.bed')
all_genes = set(read_genes(genes_bed))
matched = [g for g in query if g in all_genes]
not_in_genome = [g for g in query if g not in all_genes]
print('[ID CHECK] Exact, case-sensitive IDs only; no suffix stripping, aliases or regex matching.')
print('  Query IDs:      '+', '.join(query[:6]))
print('  Annotation IDs: '+', '.join(sorted(all_genes)[:6]))
print(f'  Exact matches: {len(matched)}/{len(query)}; absent: {len(not_in_genome)}')
for g in matched[:6]: print(f'  MATCH: {g} == {g}')
if not_in_genome: print('  UNMATCHED: '+', '.join(not_in_genome[:6]))
with open(f'{outdir}/id_check.tsv','w') as w:
    w.write('query_gene_id\tannotation_gene_id\tid_status\n')
    for g in query: w.write(f"{g}\t{g if g in all_genes else ''}\t{'exact_match' if g in all_genes else 'not_found'}\n")
Path(outdir,'id_not_found.txt').write_text(''.join(g+'\n' for g in not_in_genome), encoding='utf-8')
if not matched: sys.exit('[ERROR] No exact ID matches. Review id_check.tsv; IDs were not rewritten.')
if confirmed == 'true':
    print('[CONFIRM] Exact-ID comparison accepted via --yes; unmatched IDs remain unmatched.')
else:
    try: answer = input('Confirm these ID namespaces and continue exact matching? [y/N]: ').strip().lower()
    except EOFError: answer = ''
    if answer not in ('y','yes'):
        status.write_text('cancelled\n')
        print('[CANCELLED] No binding results generated. Unattended runs require explicit --yes.')
        sys.exit(2)
    print('[CONFIRM] User accepted exact-ID comparison.')
status.write_text('running\n')
Path(outdir,'query_genes.txt').write_text(''.join(g+'\n' for g in query), encoding='utf-8')
regional = audit(genes_bed, peak_file, anno_f, outdir, up, down, query, read_sizes(chrom_sizes))

# Identical motif filters to candidates, applied to every peak before ranking.
evidence = Evidence.from_env(fimo_files, pmax)
peak_motifs=collections.defaultdict(set)
mapped_hits=[]
for h, pk in evidence.mapped(by_chr):
    peak_motifs[pk['pid']].add(h['motif_label'])
    mapped_hits.append((h,pk))
mode=','.join(sorted({h['coordinate_mode'] for h in evidence.hits})) or 'no_hits'

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
n_exact=0
for q in query:
    recs = by_gene.get(q) if q in all_genes else None
    if recs: n_exact += 1
    recs = [r for r in (recs or []) if r["sig"]>=min_signal]
    if not recs:
        out_rows.append([q,"no","","-","no","","","no","-","none"]); continue
    # 每个 peak 的证据
    evald=[]
    for r in recs:
        mset=peak_motifs.get(r["pid"], set())
        in_t, _ = location_rule(r["feat"], r["dist"], target_features, maxtss)
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

for row in out_rows:
    q = row[0]
    row.extend(['exact_match' if q in all_genes else 'not_found', *map(str, regional[q][6:])])

# Output: keep every query unless --level explicitly filters evidence.
def keep(row):
    if level=="candidate": return row[9]=="candidate"
    if level=="in_target": return row[4]=="yes"
    return True
hdr=["gene_id","has_peak","peak_id","feature","in_target","signalValue","globalIDR","has_motif","motifs","tier","id_status","n_promoter_peaks","n_self_assigned","n_other_assigned","n_unannotated"]
with open(f"{outdir}/query_binding.tsv","w") as w:
    w.write("\t".join(hdr)+"\n")
    for row in sorted(out_rows, key=lambda r:(tier_rank[r[9]], -(float(r[5]) if r[5] else 0))):
        if keep(row): w.write("\t".join(row)+"\n")
open(f"{outdir}/bound_genes.txt","w").write("\n".join(bound_genes)+("\n" if bound_genes else ""))
open(f"{outdir}/candidate_genes.txt","w").write("\n".join(cand_genes)+("\n" if cand_genes else ""))

selected_genes = [r[0] for r in out_rows if keep(r) and r[0] in all_genes]
selected_set = set(selected_genes)
selected_evidence = [dict(gene_id=pk['gene'],peak_id=pk['pid'],**h) for h,pk in mapped_hits
    if pk['gene'] in selected_set and pk['sig'] >= min_signal and
    (level == 'all' or location_rule(pk['feat'],pk['dist'],target_features,maxtss)[0])]
evidence.save(outdir, selected_genes, selected_evidence, anno_f, peak_file,
    dict(features=feat_s,min_signal=min_signal,max_tss_dist=maxtss,
         promoter_up=up,promoter_down=down,level=level))
nq=len(query)
print(f"  Querying {nq} genes (min-signal={min_signal:g}; FIMO coordinate mode={mode}):")
print(f'    ID validation: exact={len(matched)}/{nq}; absent={len(not_in_genome)}; no ID transformation')
print(f"    ① Has assigned peak (ChIPseeker):        {n_peak}")
print(f"    ② Accepted location (TARGET_FEATURES):      {n_target}")
print(f"    ③ Candidate (accepted location + FIMO motif):   {n_cand}  ← candidate rule satisfied")
# 证据分解速览：各 tier 计数（帮助解释"为什么/为什么不是候选"）
tc=collections.Counter(r[9] for r in out_rows)
print(f"  Evidence categories: candidate={tc['candidate']}  location-only (no motif)={tc['in_target']}  "
      f"motif-only (outside accepted locations)={tc['motif_only']}  peak-only={tc['peak_only']}  no peak={tc['none']}")
if cand_genes:
    print("  Candidate genes (location + motif; first 12):")
    for row in sorted([r for r in out_rows if r[9]=="candidate"], key=lambda r:-(float(r[5]) if r[5] else 0))[:12]:
        print(f"    {row[0]:22} feat={row[3]:16} sig={row[5]:>9} motif={row[8]}")
else:
    print("  ⚠ No candidate genes: no peak satisfies both location and motif criteria (see tier in query_binding.tsv)")
status.write_text('complete\n')
