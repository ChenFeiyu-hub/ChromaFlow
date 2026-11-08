#!/usr/bin/env python3
"""Export existing FIMO hits as BED6 for genome tracks; do not rescan sequences."""
import argparse
import csv
import math
from pathlib import Path
import re


def fimo_rows(path, sources=None):
    """Read a FIMO TSV, including the completed, footer-only zero-hit format.

    FIMO 5.5.7 writes its header in print_site_as_tsv(), on the first hit.
    An empty scan therefore contains only the three completion comments.
    A blank file, arbitrary comments, or malformed columns are not zero hits.
    """
    comments = []
    required = {'motif_id', 'sequence_name', 'start', 'stop', 'strand', 'p-value'}
    with open(path, encoding='utf-8-sig', newline='') as stream:
        def lines():
            for line in stream:
                if not line.strip():
                    continue
                if line.startswith('#'):
                    comments.append(line.strip())
                else:
                    yield line
        reader = csv.DictReader(lines(), delimiter='\t')
        fields = reader.fieldnames
        if fields is None:
            completed = (
                len(comments) == 3
                and re.fullmatch(r'# FIMO \(Find Individual Motif Occurrences\): Version '
                                 r'\d+\.\d+\.\d+.* compiled on .+ at .+', comments[0])
                and comments[1].startswith('# The format of this file is described at ')
                and comments[1].endswith('/doc/fimo-output-format.html#tsv_results.')
                and re.match(r'# (?:\S*/)?fimo\s+\S+\s+\S+', comments[2])
            )
            if not completed:
                raise ValueError('Missing required result columns: ' + str(path)
                                 + '; no table header or recognized FIMO completion footer')
            # A companion GFF, when available, must agree. Never hide a lost TSV
            # table when the same run reports hits in its other output format.
            gff = Path(path).with_suffix('.gff')
            if gff.exists():
                header_seen = False
                with gff.open(encoding='utf-8') as handle:
                    for line in handle:
                        line = line.strip()
                        header_seen |= line == '##gff-version 3'
                        if line and not line.startswith('#'):
                            raise ValueError('Footer-only FIMO TSV conflicts with its GFF: ' + str(path))
                if not header_seen:
                    raise ValueError('Footer-only FIMO TSV conflicts with its GFF: ' + str(path))
                if sources is not None:
                    sources.append(gff.resolve())
            print('[FIMO] No reported hits (completed footer-only result): ' + str(path))
            return []
        missing = required - set(fields)
        if missing:
            raise ValueError('Missing required result columns: ' + str(path)
                             + '; missing=' + ','.join(sorted(missing))
                             + '; found=' + ','.join(fields))
        rows = list(reader)
        for rowno, row in enumerate(rows, 2):
            if None in row or any(row.get(key) is None for key in required):
                raise ValueError(f'Malformed FIMO table row {rowno}: {path}')
        return rows


def read_sizes(path):
    sizes = {}
    if path and Path(path).is_file():
        for line in Path(path).read_text().splitlines():
            cols = line.split('\t')
            if len(cols) >= 2:
                sizes[cols[0]] = int(cols[1])
    return sizes


def coordinates(sequence, start, stop, sizes):
    """FIMO coordinates are 1-based inclusive; BED is 0-based half-open."""
    lo, hi = min(start, stop), max(start, stop)
    if lo < 1:
        raise ValueError('FIMO start/stop must be positive')
    if sequence in sizes:
        chrom, left, right, mode = sequence, lo-1, hi, 'genomic'
    else:
        # BED getfasta -name+ may prefix the genomic window with a feature name.
        window = sequence.rsplit('::', 1)[-1]
        if window.endswith('(-)'):
            raise ValueError('Reverse-oriented window headers require explicit coordinate conversion')
        match = re.fullmatch(r'(.+):(\d+)-(\d+)(?:\([+-]\))?', window)
        if match:
            chrom, offset, end = match[1], int(match[2]), int(match[3])
            if hi > end-offset:
                raise ValueError('Relative hit exceeds its sequence window')
            left, right, mode = offset+lo-1, offset+hi, 'window_relative'
        else:
            chrom, left, right, mode = sequence, lo-1, hi, 'genomic'
    if sizes and chrom not in sizes:
        raise ValueError('Sequence is absent from the reference: '+chrom)
    if sizes and right > sizes[chrom]:
        raise ValueError('Hit exceeds reference sequence length: '+chrom)
    return chrom, left, right, mode


def export(files, output, audit, pmax=1e-4, sizes=None):
    if not 0 < pmax <= 1:
        raise ValueError('FIMO p-value threshold must be in (0, 1]')
    sizes = sizes or {}
    hits = {}
    raw_count = 0
    for path in files:
        for rowno, row in enumerate(fimo_rows(path), 2):
            try:
                p = float(row['p-value'])
                if not math.isfinite(p) or not 0 <= p <= 1:
                    raise ValueError('Invalid p-value')
                if p > pmax:
                    continue
                chrom, start, end, mode = coordinates(row['sequence_name'].strip(), int(row['start']), int(row['stop']), sizes)
                strand = row['strand'].strip()
                if strand not in ('+','-','.'):
                    raise ValueError('Invalid strand')
            except (ValueError, TypeError) as exc:
                raise ValueError(f'{path}: row {rowno}: {exc}') from exc
            motif = row['motif_id'].strip()
            alias = (row.get('motif_alt_id') or '').strip()
            label = alias if alias and alias != '.' else motif
            raw_count += 1
            key = (chrom, start, end, strand, motif, label)
            if key not in hits or p < hits[key][0]:
                hits[key] = (p, mode, str(path))
    # Parse everything before publishing so a bad row cannot look like success.
    with open(output, 'w', encoding='utf-8', newline='') as bed, open(audit, 'w', encoding='utf-8', newline='') as report:
        report.write('chrom\tstart\tend\tstrand\tmotif_id\tmotif_label\tp_value\tcoordinate_mode\tsource\n')
        for (chrom,start,end,strand,motif,label),(p,mode,path) in sorted(hits.items()):
            score = min(1000, round(-10*math.log10(max(p,1e-100))))
            bed.write(f'{chrom}\t{start}\t{end}\t{label} (p={p:.2g})\t{score}\t{strand}\n')
            report.write(f'{chrom}\t{start}\t{end}\t{strand}\t{motif}\t{label}\t{p:g}\t{mode}\t{path}\n')
    print(f'[MOTIF] {raw_count} passing rows; {len(hits)} unique genomic hits; p<={pmax:g}')
    print('[MOTIF] Existing scanned windows only; absence elsewhere does not establish absence of a motif.')
    return hits


def scan_windows(fasta, output, sizes):
    """Recover the extracted input windows, not BED regions getfasta skipped."""
    windows = set()
    with open(fasta, encoding='utf-8') as stream:
        for line in stream:
            if not line.startswith('>'):
                continue
            name = line[1:].strip().split()[0].rsplit('::', 1)[-1]
            match = re.fullmatch(r'(.+):(\d+)-(\d+)(?:\(\+\))?', name)
            if not match:
                raise ValueError('Cannot recover a genomic window from FASTA header: '+name)
            chrom, start, end = match[1], int(match[2]), int(match[3])
            if end <= start or (sizes and (chrom not in sizes or end > sizes[chrom])):
                raise ValueError('Invalid FASTA window: '+name)
            windows.add((chrom, start, end))
    if not windows:
        raise ValueError('No FASTA input windows found')
    with open(output, 'w', encoding='utf-8', newline='') as stream:
        for chrom, start, end in sorted(windows):
            stream.write(f'{chrom}\t{start}\t{end}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--audit', required=True, type=Path)
    parser.add_argument('--pval', type=float, default=1e-4)
    parser.add_argument('--chrom-sizes')
    parser.add_argument('--scan-fasta')
    parser.add_argument('--windows-output')
    parser.add_argument('fimo', nargs='+', type=Path)
    args = parser.parse_args()
    try:
        sizes = read_sizes(args.chrom_sizes)
        if args.scan_fasta:
            if not args.windows_output:
                raise ValueError('--scan-fasta requires --windows-output')
            scan_windows(args.scan_fasta, args.windows_output, sizes)
        export(args.fimo, args.output, args.audit, args.pval, sizes)
    except (ValueError, OSError) as exc:
        parser.exit(1, '[ERROR] '+str(exc)+'\n')


if __name__ == '__main__':
    main()
