#!/usr/bin/env python3
"""
normalize_gff.py — 不依赖 AGAT 的轻量 GFF3 规范化器

目的与 AGAT 的 agat_convert_sp_gxf2gxf.pl 一致（但只做 ChIPseeker/TxDb 真正需要的部分）：
  1. 统一 feature 命名（5'UTR / UTR5 / five_prime_utr → five_prime_UTR）
  2. 缺 exon 行时，从 CDS(+UTR) 推断并补全 exon
  3. 缺 UTR 行时，用 exon − CDS 推断 five_prime_UTR / three_prime_UTR
  4. 缺 gene 行时，按转录本范围合成 gene，并补上 Parent 关系
  5. 基因/转录本坐标收敛到子特征的实际跨度，按坐标排序输出

局限：不做序列校验、不修复真正错乱的父子关系、不处理跨基因重叠合并。
GFF 结构本身很乱时，仍建议装 AGAT 跑 --engine agat。

用法: normalize_gff.py <in.gff3|in.gff3.gz> <out.gff3>
"""

import gzip
import sys
from collections import OrderedDict, defaultdict

GENE_TYPES = {"gene", "ncrna_gene", "pseudogene", "protein_coding_gene"}
TX_TYPES = {"mrna", "transcript", "trna", "rrna", "ncrna", "mirna", "snorna", "snrna",
            "lnc_rna", "lncrna", "pre_mirna", "pseudogenic_transcript", "primary_transcript"}
UTR_ALIASES = {
    "5'utr": "five_prime_UTR", "5utr": "five_prime_UTR", "utr5": "five_prime_UTR",
    "five_prime_utr": "five_prime_UTR", "five_prime_untranslated_region": "five_prime_UTR",
    "3'utr": "three_prime_UTR", "3utr": "three_prime_UTR", "utr3": "three_prime_UTR",
    "three_prime_utr": "three_prime_UTR", "three_prime_untranslated_region": "three_prime_UTR",
}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def parse_attrs(s):
    attrs = OrderedDict()
    for field in s.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
        else:                       # GTF 风格容错: key "value"
            parts = field.split(None, 1)
            if len(parts) != 2:
                continue
            k, v = parts[0], parts[1].strip('"')
        attrs[k.strip()] = v.strip()
    return attrs


def fmt_attrs(attrs):
    return ";".join("%s=%s" % (k, v) for k, v in attrs.items() if v != "")


class Feat(object):
    __slots__ = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attrs")

    def __init__(self, cols):
        self.seqid, self.source, self.type = cols[0], cols[1], cols[2]
        self.start, self.end = int(cols[3]), int(cols[4])
        self.score, self.strand, self.phase = cols[5], cols[6], cols[7]
        self.attrs = parse_attrs(cols[8]) if len(cols) > 8 else OrderedDict()

    @property
    def id(self):
        return self.attrs.get("ID", "")

    @property
    def parents(self):
        p = self.attrs.get("Parent", "")
        return [x for x in p.split(",") if x] if p else []

    def line(self):
        return "\t".join([self.seqid, self.source, self.type, str(self.start), str(self.end),
                          self.score, self.strand, self.phase, fmt_attrs(self.attrs)])


def clone(f, ftype, start, end, fid, parent):
    g = Feat([f.seqid, f.source, ftype, str(start), str(end), ".", f.strand, "."])
    g.attrs = OrderedDict([("ID", fid), ("Parent", parent)])
    return g


def subtract(blocks, holes):
    """区间差集：blocks - holes（闭区间，1-based）"""
    out = []
    for bs, be in blocks:
        cur = [(bs, be)]
        for hs, he in holes:
            nxt = []
            for s, e in cur:
                if he < s or hs > e:
                    nxt.append((s, e))
                    continue
                if s < hs:
                    nxt.append((s, min(e, hs - 1)))
                if e > he:
                    nxt.append((max(s, he + 1), e))
            cur = nxt
        out.extend(cur)
    return sorted(out)


def merge(blocks):
    out = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def main():
    if len(sys.argv) < 3:
        sys.exit("用法: normalize_gff.py <in.gff3> <out.gff3>")
    src, dst = sys.argv[1], sys.argv[2]

    feats, headers, auto = [], [], 0
    with opener(src) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region") or line.startswith("##gff-version"):
                    headers.append(line.rstrip("\n"))
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            try:
                f = Feat(cols)
            except ValueError:
                continue
            low = f.type.lower()
            f.type = UTR_ALIASES.get(low, f.type)
            if low in ("utr", "untranslated_region"):
                f.type = "UTR"                      # 稍后按 CDS 位置定 5'/3'
            feats.append(f)

    by_id = {}
    for f in feats:
        if f.id:
            by_id.setdefault(f.id, f)
    children = defaultdict(list)
    for f in feats:
        for p in f.parents:
            children[p].append(f)

    genes = [f for f in feats if f.type.lower() in GENE_TYPES]
    txs = [f for f in feats if f.type.lower() in TX_TYPES]

    # ── 无 transcript 层：把带 CDS/exon 的 gene 直接当转录本处理 ──
    if not txs:
        for g in genes:
            if children.get(g.id):
                tx = clone(g, "mRNA", g.start, g.end, g.id + ".t1", g.id)
                feats.append(tx)
                txs.append(tx)
                by_id[tx.id] = tx
                for c in children[g.id]:
                    c.attrs["Parent"] = tx.id
                children[tx.id] = children.pop(g.id)

    stats = defaultdict(int)
    new_feats = []

    for tx in txs:
        kids = children.get(tx.id, [])
        exons = [k for k in kids if k.type.lower() == "exon"]
        cds = [k for k in kids if k.type.lower() == "cds"]
        utrs = [k for k in kids if k.type in ("five_prime_UTR", "three_prime_UTR", "UTR")]
        minus = tx.strand == "-"

        # 1) 补 exon
        if not exons:
            blocks = merge([(k.start, k.end) for k in cds + utrs]) or [(tx.start, tx.end)]
            for i, (s, e) in enumerate(blocks, 1):
                auto += 1
                ex = clone(tx, "exon", s, e, "%s.exon%d" % (tx.id, i), tx.id)
                new_feats.append(ex)
                exons.append(ex)
                stats["exon_added"] += 1

        # 2) UTR：把无向 UTR 定向；缺失则由 exon − CDS 推断
        if cds:
            cds_blocks = merge([(k.start, k.end) for k in cds])
            cds_min, cds_max = cds_blocks[0][0], cds_blocks[-1][1]
            for u in utrs:
                if u.type == "UTR":
                    is_five = (u.end < cds_min) if not minus else (u.start > cds_max)
                    u.type = "five_prime_UTR" if is_five else "three_prime_UTR"
                    stats["utr_typed"] += 1
            has_five = any(u.type == "five_prime_UTR" for u in utrs)
            has_three = any(u.type == "three_prime_UTR" for u in utrs)
            if not (has_five and has_three):
                rest = subtract(merge([(e.start, e.end) for e in exons]), cds_blocks)
                n5 = n3 = 0
                for s, e in rest:
                    is_five = (e < cds_min) if not minus else (s > cds_max)
                    if is_five and has_five:
                        continue
                    if (not is_five) and has_three:
                        continue
                    if is_five:
                        n5 += 1
                        fid, ftype = "%s.utr5_%d" % (tx.id, n5), "five_prime_UTR"
                    else:
                        n3 += 1
                        fid, ftype = "%s.utr3_%d" % (tx.id, n3), "three_prime_UTR"
                    new_feats.append(clone(tx, ftype, s, e, fid, tx.id))
                    stats["utr_added"] += 1

        # 3) 转录本坐标收敛
        span = [(k.start, k.end) for k in exons + cds + utrs]
        if span:
            tx.start, tx.end = min(s for s, _ in span), max(e for _, e in span)

    feats.extend(new_feats)
    children = defaultdict(list)
    for f in feats:
        for p in f.parents:
            children[p].append(f)

    # ── 4) 补 gene 层 ──
    gene_ids = set(g.id for g in genes)
    synth = []
    for tx in txs:
        pids = [p for p in tx.parents if p in gene_ids]
        if pids:
            continue
        gid = tx.attrs.get("gene_id") or tx.attrs.get("geneID") or (tx.id + ".gene")
        if gid in by_id and by_id[gid].type.lower() in GENE_TYPES:
            tx.attrs["Parent"] = gid
            continue
        g = Feat([tx.seqid, tx.source, "gene", str(tx.start), str(tx.end), ".", tx.strand, "."])
        g.attrs = OrderedDict([("ID", gid)])
        if "Name" in tx.attrs:
            g.attrs["Name"] = tx.attrs["Name"]
        tx.attrs["Parent"] = gid
        synth.append(g)
        by_id[gid] = g
        gene_ids.add(gid)
        stats["gene_added"] += 1
    genes.extend(synth)
    feats.extend(synth)

    # ── 5) 基因坐标收敛到转录本跨度 ──
    tx_by_gene = defaultdict(list)
    for tx in txs:
        for p in tx.parents:
            tx_by_gene[p].append(tx)
    for g in genes:
        kids = tx_by_gene.get(g.id, [])
        if kids:
            g.start = min(k.start for k in kids)
            g.end = max(k.end for k in kids)

    # ── 6) 排序输出：gene → transcript → 子特征 ──
    rank = {"gene": 0, "ncrna_gene": 0, "pseudogene": 0}
    for t in TX_TYPES:
        rank[t] = 1

    def key(f):
        return (f.seqid, f.start, rank.get(f.type.lower(), 2), f.end, f.type)

    with open(dst, "w") as out:
        out.write("##gff-version 3\n")
        out.write("#!normalized-by ChromaFlow normalize_gff.py (builtin, AGAT-free)\n")
        for h in headers:
            if not h.startswith("##gff-version"):
                out.write(h + "\n")
        for f in sorted(feats, key=key):
            out.write(f.line() + "\n")

    print("  exon 补全: %d | UTR 推断: %d | UTR 定向: %d | gene 合成: %d"
          % (stats["exon_added"], stats["utr_added"], stats["utr_typed"], stats["gene_added"]))


if __name__ == "__main__":
    main()
