#!/usr/bin/env python3
"""Shared, exact-model motif filtering and reproducible tracks selections."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
from annotation_policy import sha256
from fimo_tracks import coordinates, scan_windows


def number(value, name, maximum=None):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'Missing or invalid {name}: {value}') from exc
    if not math.isfinite(result) or result < 0 or (maximum is not None and result > maximum):
        raise ValueError(f'Invalid {name}: {value}')
    return result


def table(path, required=()):
    with open(path, encoding='utf-8') as stream:
        reader = csv.DictReader((line for line in stream if line.strip() and not line.startswith('#')), delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError('Missing required result columns: '+str(path))
        return list(reader)


def html_data(path):
    """Read MEME's JSON payload, never execute JavaScript from a report."""
    text = Path(path).read_text(encoding='utf-8')
    match = re.search(r'\bvar\s+data\s*=\s*', text)
    if not match:
        raise ValueError('MEME-ChIP report has no supported JSON payload: '+str(path))
    return json.JSONDecoder().raw_decode(text[match.end():])[0]


def write_tsv(path, fields, rows):
    with open(path, 'w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


class Evidence:
    def __init__(self, files, root, pmax=1e-4, motifs='', qmax=None,
                 emax=None, centrimo_max=None, display_max=0.05):
        self.files = [Path(f).resolve() for f in files]
        self.root = Path(root).resolve()
        self.sources = list(self.files)
        self.filters = {'motifs': [s.strip() for s in motifs.split(',') if s.strip()],
                        'pval': number(pmax, 'FIMO p-value threshold', 1),
                        'qval': None if qmax in (None, '') else number(qmax, 'FIMO q-value threshold', 1),
                        'motif_evalue': None if emax in (None, '') else number(emax, 'discovery E-value threshold'),
                        'centrimo_evalue': None if centrimo_max in (None, '') else number(centrimo_max, 'CentriMo E-value threshold')}
        if self.filters['pval'] == 0:
            raise ValueError('FIMO p-value threshold must be positive')
        if not self.files and (self.filters['motifs'] or any(self.filters[k] is not None for k in ('qval','motif_evalue','centrimo_evalue'))):
            raise ValueError('Requested motif filters require available FIMO results; run motif first')
        self.display_max = number(display_max, 'CentriMo report threshold')
        self.models = {}
        self.mode = 'unknown'
        report = self.root/'meme-chip.html'
        if report.is_file():
            data = html_data(report)
            self.sources.append(report)
            cmd = data.get('cmd', [])
            if cmd:
                self.mode = 'local' if '-centrimo-local' in cmd else 'central'
            for m in data.get('motifs', []):
                raw, label = m['id'], m.get('alt') or m['id']
                key = (raw, label)
                if key in self.models:
                    raise ValueError('Ambiguous motif identity in MEME-ChIP report: '+label)
                pwm = m.get('pwm')
                scanned = bool(m.get('fimo_gff'))
                self.models[key] = {'motif_id': raw, 'motif_label': label, 'program': m.get('prog', ''),
                    'discovery_evalue': m.get('evalue', ''),
                    'pwm_sha256': hashlib.sha256(json.dumps(pwm, sort_keys=True).encode()).hexdigest() if pwm else '',
                    'scan_status': 'scanned' if scanned else 'not_scanned',
                    'fimo_file': str(self.root/Path(m['fimo_gff']).parent/'fimo.tsv') if scanned else ''}
        all_hits = []
        for path in self.files:
            if path.parent.parent != self.root:
                raise ValueError('FIMO files must belong to one MEME-ChIP run: '+str(path))
            for row in table(path, ('motif_id','sequence_name','start','stop','strand','p-value')):
                raw = (row.get('motif_id') or '').strip()
                alias = (row.get('motif_alt_id') or '').strip()
                label = alias if alias and alias != '.' else raw
                if not raw:
                    raise ValueError('FIMO motif_id missing: '+str(path))
                key = (raw, label)
                model = self.models.setdefault(key, {'motif_id': raw, 'motif_label': label,
                    'program': '', 'discovery_evalue': '', 'pwm_sha256': ''})
                model.update(scan_status='scanned', fimo_file=str(path))
                p = number(row.get('p-value'), 'FIMO p-value', 1)
                q = row.get('q-value')
                q = None if q in (None, '', 'NA', '.') else number(q, 'FIMO q-value', 1)
                chrom, start, end, mode = coordinates(row.get('sequence_name', '').strip(), int(row['start']), int(row['stop']), {})
                strand = row.get('strand', '')
                if strand not in ('+', '-', '.') or not chrom:
                    raise ValueError('Invalid FIMO sequence/strand: '+str(path))
                all_hits.append(dict(chrom=chrom, start=start, end=end, strand=strand,
                    motif_id=raw, motif_label=label, p_value=p, q_value=q,
                    coordinate_mode=mode, source=str(path)))
        # Explicit selectors are exact within this run, not FIMO directory numbers.
        requested = set()
        for label in self.filters['motifs']:
            matches = [k for k in self.models if label in k]
            if len(matches) != 1:
                raise ValueError(f'Motif selector {label!r} is unknown or ambiguous; see motif_catalog.tsv / MEME-ChIP report')
            key = matches[0]
            if self.models[key]['scan_status'] != 'scanned' or not Path(self.models[key]['fimo_file']).is_file():
                raise ValueError(f'{label} has no available FIMO scan (it may be a non-representative cluster member). Scan that exact model before selecting it.')
            requested.add(key)
        centrimo = self.root/'centrimo_out/centrimo.tsv'
        crows = table(centrimo) if centrimo.is_file() else []
        if centrimo.is_file():
            self.sources.append(centrimo)
        self.selected_models = set()
        for key, m in self.models.items():
            matches = [r for r in crows if r.get('motif_id') == key[0] and
                       (r.get('motif_alt_id') or key[1]) in (key[1], '.')]
            if len(matches) > 1:
                raise ValueError('Ambiguous CentriMo records for motif: '+key[1])
            c = matches[0] if matches else {}
            ev = c.get('E-value', '')
            value = number(ev, 'CentriMo E-value') if ev != '' else None
            m.update(centrimo_Evalue=ev, centrimo_bin=c.get('bin_location', ''), centrimo_mode=self.mode,
                positional_enrichment=('not_reported' if value is None else
                    'passes_threshold' if value <= self.display_max else 'above_threshold'),
                centrimo_report_threshold=self.display_max,
                centrally_enriched=('not_assessed' if self.mode != 'central' or value is None else
                                   'yes' if value <= self.display_max else 'no'))
            if m['scan_status'] == 'scanned' and Path(m['fimo_file']).resolve() not in self.files:
                m['scan_status'] = 'scan_file_missing'
            selected = m['scan_status'] == 'scanned' and (not requested or key in requested)
            if selected and self.filters['motif_evalue'] is not None:
                if m['discovery_evalue'] == '':
                    raise ValueError('Discovery E-value unavailable for '+key[1]+'; restore meme-chip.html or omit --motif-evalue')
                selected = number(m['discovery_evalue'], 'discovery E-value') <= self.filters['motif_evalue']
            if selected and self.filters['centrimo_evalue'] is not None:
                selected = value is not None and value <= self.filters['centrimo_evalue']
            m['selected'] = int(selected)
            if selected:
                self.selected_models.add(key)
        labels = [k[1] for k in self.selected_models]
        if len(set(labels)) != len(labels):
            raise ValueError('Multiple selected models share a motif label; select a unique raw motif_id')
        self.hits = []
        seen = {}
        for h in all_hits:
            if (h['motif_id'], h['motif_label']) not in self.selected_models:
                continue
            # A requested q filter must never silently accept unknown q-values.
            if self.filters['qval'] is not None and h['q_value'] is None:
                raise ValueError('FIMO q-value missing for selected motif '+h['motif_label']+'; restore a scan with q-values or omit --qval')
            if h['p_value'] > self.filters['pval'] or (self.filters['qval'] is not None and h['q_value'] > self.filters['qval']):
                continue
            key = tuple(h[k] for k in ('chrom', 'start', 'end', 'strand', 'motif_id', 'motif_label'))
            if key not in seen or h['p_value'] < seen[key]['p_value']:
                seen[key] = h
        self.hits = [seen[k] for k in sorted(seen)]
        print(f'[FILTER] {len(self.selected_models)} motif models; {len(self.hits)} unique passing hits; {json.dumps(self.filters)}')
        print('[MOTIF] Filters apply to existing FIMO output only. Unscanned models/regions are not negatives; relaxing a scan threshold requires rescanning.')

    @classmethod
    def from_env(cls, files, pmax):
        return cls(files, os.environ.get('CFSEL_ROOT', '.'), pmax,
            os.environ.get('CFSEL_MOTIFS', ''), os.environ.get('CFSEL_QVAL'),
            os.environ.get('CFSEL_EVALUE'), os.environ.get('CFSEL_CENTRIMO'))

    def mapped(self, by_chr):
        for h in self.hits:
            for pk in by_chr.get(h['chrom'], []):
                if pk['s'] >= h['end']:
                    break
                if h['start'] < pk['e']:
                    yield h, pk

    def save(self, outdir, genes, evidence, annotation, peak_file, parameters):
        out = Path(outdir).resolve()
        fields = ['motif_id','motif_label','program','discovery_evalue','pwm_sha256','scan_status','fimo_file',
                  'selected','centrimo_Evalue','centrimo_bin','centrimo_mode','centrimo_report_threshold',
                  'positional_enrichment','centrally_enriched']
        write_tsv(out/'motif_catalog.tsv', fields, self.models.values())
        (out/'selected_genes.txt').write_text(''.join(g+'\n' for g in genes), encoding='utf-8')
        hit_fields = ['chrom','start','end','strand','motif_id','motif_label','p_value','q_value','coordinate_mode','source']
        write_tsv(out/'selected_evidence.tsv', ['gene_id','peak_id']+hit_fields, evidence)
        write_tsv(out/'motif_hits.tsv', hit_fields, self.hits)
        with (out/'motif_hits.bed').open('w', encoding='utf-8') as stream:
            for h in self.hits:
                stream.write(f"{h['chrom']}\t{h['start']}\t{h['end']}\t{h['motif_label']} (p={h['p_value']:.2g})\t0\t{h['strand']}\n")
        files = {name: {'file': name, 'sha256': sha256(out/name)} for name in
                 ('selected_genes.txt','selected_evidence.tsv','motif_hits.tsv','motif_hits.bed','motif_catalog.tsv')}
        fasta = self.root.parent/'peak_sequences.fa'
        if fasta.is_file():
            scan_windows(fasta, out/'motif_input_windows.bed', {})
            self.sources.append(fasta)
            files['motif_input_windows.bed'] = {'file': 'motif_input_windows.bed', 'sha256': sha256(out/'motif_input_windows.bed')}
        for extra in [Path(annotation).parent/'genes.bed', Path(annotation).parent/'annotation_policy.json', Path(os.environ.get('CFSEL_QUERY', '/__absent__'))]:
            if extra.is_file():
                self.sources.append(extra.resolve())
        sources = list(dict.fromkeys(self.sources+[Path(annotation).resolve(), Path(peak_file).resolve()]))
        result = {'schema': 1, 'status': 'complete', 'version': '1.2.7',
            'pipeline': os.environ.get('CFSEL_PIPELINE', 'chipseq'), 'peakset': os.environ.get('CHROMAFLOW_PEAKSET', 'idr'),
            'peak_file': str(Path(peak_file).resolve()), 'peak_sha256': sha256(peak_file),
            'filters': self.filters, 'parameters': parameters, 'files': files,
            'sources': [{'path': str(p), 'sha256': sha256(p)} for p in sources],
            'scope': 'Computational candidates; motif_hits includes all passing context hits, selected_evidence records selected gene/peak associations.'}
        tmp = out/'selection.json.tmp'
        tmp.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        tmp.replace(out/'selection.json')
        print('[SELECTION] '+str(out/'selection.json'))


def validate(path, peakset, pipeline, peak_file, up, down):
    path = Path(path).resolve()
    if path.is_dir():
        path = path/'selection.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    if data.get('schema') != 1 or data.get('status') != 'complete':
        raise ValueError('Unsupported or incomplete selection')
    if data['peakset'] != peakset or data['pipeline'] != pipeline or data['peak_sha256'] != sha256(peak_file):
        raise ValueError('Selection belongs to a different pipeline or peak set; use its original --peakset/config')
    if [data['parameters'].get('promoter_up'), data['parameters'].get('promoter_down')] != [up, down]:
        raise ValueError('Selection promoter parameters differ from this config')
    for src in data['sources']:
        if not Path(src['path']).is_file() or sha256(src['path']) != src['sha256']:
            raise ValueError('Selection input changed or is missing: '+src['path'])
    for item in data['files'].values():
        target = (path.parent/item['file']).resolve()
        if target.parent != path.parent or not target.is_file() or sha256(target) != item['sha256']:
            raise ValueError('Selection output changed or is missing: '+str(target))
    return path.parent


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Validate a candidates/lookup selection for tracks')
    parser.add_argument('selection')
    parser.add_argument('peakset')
    parser.add_argument('pipeline')
    parser.add_argument('peaks')
    parser.add_argument('up', type=int)
    parser.add_argument('down', type=int)
    a = parser.parse_args()
    try:
        print(validate(a.selection, a.peakset, a.pipeline, a.peaks, a.up, a.down))
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, '[ERROR] '+str(exc)+'\n')
