"""Audit all selected peaks against complete gene-level promoter intervals."""
import argparse
from bisect import bisect_left, bisect_right
import csv
import json
from pathlib import Path

from annotation_policy import sha256


def read_genes(path):
    genes = {}
    with open(path, encoding='utf-8') as stream:
        for no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith('#'): continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 6: raise ValueError(f'genes.bed row {no}: expected BED6')
            chrom, start, end, gid, _, strand = c[:6]
            start, end = int(start), int(end)
            if start < 0 or end <= start or strand not in ('+','-') or not gid:
                raise ValueError(f'genes.bed row {no}: invalid coordinates, ID or strand')
            if gid in genes: raise ValueError('Ambiguous duplicate gene ID: '+gid)
            genes[gid] = (chrom, start, end, gid, strand)
    if not genes: raise ValueError('No gene coordinates found')
    return genes


def read_query(path):
    # Plain IDs or a table with an explicit ID header. Never split punctuation,
    # remove version suffixes, change case, or infer names from other columns.
    lines = [x.rstrip('\r\n') for x in Path(path).read_text(encoding='utf-8-sig').splitlines()
             if x.strip() and not x.lstrip().startswith('#')]
    if not lines: raise ValueError('Gene list is empty')
    headers = ('gene_id', 'gene', 'gene_name', 'id')
    delimiter = '\t' if '\t' in lines[0] else ',' if ',' in lines[0] else None
    first = lines[0].split(delimiter) if delimiter else [lines[0].strip()]
    found = [i for i, value in enumerate(first) if value.strip().lower() in headers]
    if found:
        if len(found) != 1: raise ValueError('Gene-list table has multiple possible ID columns')
        rows = csv.reader(lines[1:], delimiter=delimiter or '\t')
        values = [row[found[0]].strip() for row in rows if row and len(row) > found[0]]
    else:
        if delimiter: raise ValueError('Multi-column gene lists need an explicit gene_id header')
        values = [line.strip() for line in lines]
    if any(not x or any(c.isspace() for c in x) for x in values):
        raise ValueError('Gene IDs must be nonempty and contain no whitespace')
    return list(dict.fromkeys(values))


class IntervalIndex:
    """Index (chrom,start,end,...) records, including long/nested intervals."""
    def __init__(self, rows):
        self.data = {}
        for row in rows: self.data.setdefault(row[0], []).append(row)
        self.starts = {}; self.ends = {}
        for chrom, items in self.data.items():
            items.sort(key=lambda r:(r[1],r[2]))
            self.starts[chrom] = [r[1] for r in items]
            high = -1; maxima = []
            for row in items:
                high = max(high, row[2]); maxima.append(high)
            self.ends[chrom] = maxima

    def overlaps(self, chrom, start, end):
        rows = self.data.get(chrom, [])
        lo = bisect_right(self.ends.get(chrom, []), start)
        hi = bisect_left(self.starts.get(chrom, []), end)
        return [r for r in rows[lo:hi] if r[2] > start]


def promoter(gene, up, down, sizes):
    chrom, start, end, gid, strand = gene
    tss = start if strand == '+' else end-1
    left = tss-up if strand == '+' else tss-down
    right = tss+down+1 if strand == '+' else tss+up+1
    left = max(0, left)
    if sizes:
        if chrom not in sizes or end > sizes[chrom]:
            raise ValueError('Gene is incompatible with reference sizes: '+gid)
        right = min(right, sizes[chrom])
    return left, right, tss


def audit(genes_bed, peaks_bed, annotation, outdir, up, down, query=None, sizes=None):
    if min(up, down) < 0: raise ValueError('Promoter extents must be nonnegative')
    genes = read_genes(genes_bed); sizes = sizes or {}
    requested = query if query is not None else sorted(genes)
    assignments = {}
    if annotation:
        with open(annotation, encoding='utf-8') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                key = (row['seqnames'], int(row['start'])-1, int(row['end']), row['peak_id'])
                assignments.setdefault(key, []).append((row.get('geneId',''), row.get('annotation','')))
    peaks = []
    with open(peaks_bed, encoding='utf-8') as stream:
        for no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(('#','track ','browser ')): continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 3: raise ValueError(f'Peak BED row {no}: fewer than three columns')
            start, end = int(c[1]), int(c[2])
            if start < 0 or end <= start: raise ValueError(f'Peak BED row {no}: invalid interval')
            peaks.append((c[0], start, end, c[3] if len(c)>3 else f'peak_{len(peaks)+1}'))
    peak_index = IntervalIndex(peaks); gene_index = IntervalIndex(genes.values())
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    summary = {}; details = []
    with (out/'promoter_regions.bed').open('w', encoding='utf-8', newline='') as bed:
        for gid in requested:
            if gid not in genes:
                summary[gid] = [gid,'id_not_found','','','','',0,0,0,0]
                continue
            gene = genes[gid]; chrom, _, _, _, strand = gene
            start, end, tss = promoter(gene, up, down, sizes)
            bed.write(f'{chrom}\t{start}\t{end}\t{gid}\t0\t{strand}\n')
            found = peak_index.overlaps(chrom, start, end)
            n_self = n_other = n_unknown = 0
            for peak in found:
                pairs = assignments.get(peak, [])
                assigned = sorted(set(g for g, a in pairs if g and (g != 'NA' or g in genes)))
                relation = 'self' if gid in assigned else 'other' if assigned else 'unannotated'
                n_self += relation == 'self'; n_other += relation == 'other'; n_unknown += relation == 'unannotated'
                bodies = sorted({g[3] for g in gene_index.overlaps(chrom, peak[1], peak[2])})
                local_bodies = sorted({g[3] for g in gene_index.overlaps(chrom, max(start,peak[1]), min(end,peak[2]))})
                details.append([gid,chrom,start,end,strand,tss+1,peak[3],peak[1],peak[2],
                                min(end,peak[2])-max(start,peak[1]), relation,
                                json.dumps(assigned, ensure_ascii=False),
                                json.dumps(sorted({a for g,a in pairs}), ensure_ascii=False),
                                json.dumps(bodies, ensure_ascii=False), json.dumps(local_bodies, ensure_ascii=False)])
            summary[gid] = [gid,'matched',chrom,start,end,strand,len(found),n_self,n_other,n_unknown]
    for name, header, rows in [
        ('promoter_peak_overlaps.tsv', ['gene_id','chrom','promoter_start','promoter_end','strand','tss_1based',
          'peak_id','peak_start','peak_end','overlap_bp','assignment_relation','assigned_gene_ids',
          'assigned_annotations','overlapped_gene_body_ids','promoter_overlap_gene_body_ids'], details),
        ('promoter_peak_summary.tsv', ['gene_id','id_status','chrom','promoter_start','promoter_end','strand',
          'n_promoter_peaks','n_self_assigned','n_other_assigned','n_unannotated'], summary.values())]:
        with (out/name).open('w', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(header); writer.writerows(rows)
    metadata = dict(promoter_up=up, promoter_down=down, tss_level='gene',
                    coordinates='BED 0-based half-open; inclusive strand-relative [-UP,+DOWN], including TSS',
                    peak_scope='entire selected peakset; any positive promoter overlap; no candidate filter',
                    peaks_sha256=sha256(peaks_bed), genes_sha256=sha256(genes_bed),
                    annotation_sha256=sha256(annotation) if annotation else None, queried_genes=len(requested))
    (out/'promoter_audit.json').write_text(json.dumps(metadata, indent=2)+'\n', encoding='utf-8')
    print(f'[PROMOTER] Audited {len(requested)} genes against {len(peaks)} peaks; {len(details)} overlaps; -{up}/+{down} bp')
    print('[PROMOTER] Regional overlaps are supplementary evidence; ChIPseeker assignments are unchanged.')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('genes-bed','peaks','outdir'): parser.add_argument('--'+name, required=True)
    parser.add_argument('--annotation', default='')
    parser.add_argument('--up', required=True, type=int); parser.add_argument('--down', required=True, type=int)
    parser.add_argument('--gene-list'); parser.add_argument('--chrom-sizes')
    args = parser.parse_args()
    try:
        from fimo_tracks import read_sizes
        audit(args.genes_bed, args.peaks, args.annotation, args.outdir, args.up, args.down,
              read_query(args.gene_list) if args.gene_list else None, read_sizes(args.chrom_sizes))
    except (ValueError, OSError, KeyError) as exc: parser.exit(1, '[ERROR] '+str(exc)+'\n')


if __name__ == '__main__': main()
