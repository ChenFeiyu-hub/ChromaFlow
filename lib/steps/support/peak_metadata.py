"""Strict non-IDR metadata correction; never change peak/gene/motif evidence."""
import csv
import hashlib
import io
import json
import math
from pathlib import Path

POLICY = 'non-idr-metadata-v1'


def sha_bytes(value):
    return hashlib.sha256(value).hexdigest()


def table_bytes(raw):
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8')), delimiter='\t')
    fields, rows = list(reader.fieldnames or []), list(reader)
    if not fields or len(fields) != len(set(fields)) or any(None in r or None in r.values() for r in rows):
        raise ValueError('Malformed annotation metadata table')
    return fields, rows


def non_idr_replacement(raw, policy_raw, peaks):
    """Return the exact allowed replacement, verified against selected BED input."""
    fields, rows = table_bytes(raw)
    if 'peak_globalIDR' not in fields or all(r['peak_globalIDR'] == 'NA' for r in rows):
        return None
    required = {'seqnames', 'start', 'end', 'peak_id', 'V7', 'peak_globalIDR'}
    if not required <= set(fields):
        raise ValueError('Legacy IDR correction requires original narrowPeak columns')
    source = {}
    for line in Path(peaks).read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
            continue
        c = line.split('\t')
        if len(c) != 10 or c[3] in source or int(c[1]) < 0 or int(c[2]) <= int(c[1]):
            raise ValueError('Legacy IDR correction requires unambiguous narrowPeak input')
        source[c[3]] = c
    seen = set()
    for row in rows:
        c = source.get(row['peak_id'])
        if c is None or row['peak_id'] in seen:
            raise ValueError('Annotation peak is absent or duplicated in selected input')
        seen.add(row['peak_id'])
        if (row['seqnames'] != c[0] or int(row['start']) != int(c[1]) + 1 or int(row['end']) != int(c[2])
                or not math.isfinite(float(row['V7'])) or float(row['V7']) != float(c[6])
                or row['peak_globalIDR'] != row['V7']):
            raise ValueError('Annotation does not match the known narrowPeak signal-as-IDR defect')
        row['peak_globalIDR'] = 'NA'
    policy = json.loads(policy_raw)
    if policy.get('annotation_sha256') != sha_bytes(raw):
        raise ValueError('Annotation policy identity does not match the saved table')
    out = io.StringIO()
    writer = csv.DictWriter(out, fields, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    corrected = out.getvalue().encode('utf-8')
    policy['annotation_sha256'] = sha_bytes(corrected)
    return corrected, (json.dumps(policy, indent=2) + '\n').encode('utf-8'), len(rows)


def lookup_idr_values(rows, peak_kind, peak_file):
    """Only an explicitly selected, matching ChromaFlow IDR BED8 has IDR values."""
    if peak_kind != 'idr':
        return ['NA'] * len(rows)
    source = {}
    for line in Path(peak_file).read_text(encoding='utf-8').splitlines():
        if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
            continue
        c = line.split('\t')
        if (len(c) != 8 or c[3] in source or not math.isfinite(float(c[6])) or float(c[6]) < 0
                or not math.isfinite(float(c[7])) or float(c[7]) < 1 or float(c[7]) != int(float(c[7]))):
            raise ValueError('IDR lookup requires the ChromaFlow BED8 schema')
        source[c[3]] = c
    values = []
    for row in rows:
        c = source.get(row.get('peak_id'))
        if (c is None or row.get('seqnames') != c[0] or int(row['start']) != int(c[1]) + 1
                or int(row['end']) != int(c[2]) or float(row.get('peak_globalIDR', 'nan')) != float(c[6])):
            raise ValueError('IDR annotation does not match the selected IDR peak set; rerun annotate')
        values.append(row['peak_globalIDR'])
    return values
