#!/usr/bin/env python3
"""Size each region's annotation tracks without hiding transcript isoforms."""
import argparse
import math
import os
from pathlib import Path
import re
from promoter_audit import IntervalIndex, read_genes, read_query


def transcript_index(gtf):
    intervals = {}
    with open(gtf, encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            match = re.search(r'transcript_id "([^"]+)"', c[8])
            if not match:
                continue
            key = (c[0], match[1])
            start, end = int(c[3])-1, int(c[4])
            old = intervals.get(key, (c[0], start, end))
            intervals[key] = (c[0], min(start, old[1]), max(end, old[2]))
    return IntervalIndex(intervals.values())


def bed_index(path):
    rows = []
    if path and Path(path).is_file():
        for line in Path(path).read_text(encoding='utf-8').splitlines():
            c = line.split('\t')
            if len(c) >= 3:
                rows.append((c[0], int(c[1]), int(c[2])))
    return IntervalIndex(rows)


def height(variable, count, minimum):
    value = os.environ.get(variable, 'auto')
    if value == 'auto':
        return max(minimum, 0.5*count+0.4)
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(variable+' must be auto or a positive height in cm')
    return result


def prepare(ini, gtf, genes_bed, query, flank, motifs, promoters, outdir):
    tx = transcript_index(gtf)
    motif = bed_index(motifs)
    promoter = bed_index(promoters)
    genes = read_genes(genes_bed)
    template = Path(ini).read_text(encoding='utf-8')
    out = Path(outdir)/'.layout'
    out.mkdir(exist_ok=True)
    for gid in read_query(query):
        if gid not in genes:
            continue
        chrom, start, end, _, strand = genes[gid]
        start, end = max(0, start-flank), end+flank
        counts = [len(idx.overlaps(chrom, start, end)) for idx in (tx,motif,promoter)]
        sizes = {'genes':height('TRACK_GENE_HEIGHT',counts[0],2),
                 'motif_hits':height('TRACK_MOTIF_HEIGHT',counts[1],1.2),
                 'query_promoters':max(0.8,0.5*counts[2]+0.2)}
        text = template
        for section, size in sizes.items():
            pattern = r'(\['+section+r'\][\s\S]*?\nheight = )[^\n]+'
            text = re.sub(pattern, lambda m: m[1]+f'{size:g}', text)
        (out/(gid+'.ini')).write_text(text, encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('ini','gtf','genes_bed','query','flank','motifs','promoters','outdir'):
        parser.add_argument(key, type=int if key == 'flank' else str)
    a = parser.parse_args()
    try:
        prepare(**vars(a))
    except (ValueError,OSError) as exc:
        parser.exit(1,'[ERROR] '+str(exc)+'\n')
