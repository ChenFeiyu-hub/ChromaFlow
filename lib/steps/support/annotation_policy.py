"""Shared annotation categories and configuration checks; never rewrite IDs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

FEATURES = ('Promoter', "5' UTR", "3' UTR", '1st Exon', 'Other Exon',
            '1st Intron', 'Other Intron', 'Downstream', 'Distal Intergenic')
PRIORITY = ['Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstream', 'Intergenic']


def norm_feat(annotation):
    text = (annotation or '').strip()
    if text.startswith('Promoter'): return 'Promoter'
    if text.startswith("5' UTR"): return "5' UTR"
    if text.startswith("3' UTR"): return "3' UTR"
    for kind in ('Exon', 'Intron'):
        if text.startswith('1st '+kind): return '1st '+kind
        if text.startswith('Other '+kind): return 'Other '+kind
        if text.startswith(kind):
            first = re.search(r',\s*'+kind.lower()+r'\s+1\s+of\s+\d+\s*\)$', text, re.I)
            return ('1st ' if first else 'Other ')+kind
    if text.startswith('Downstream'): return 'Downstream'
    if text in ('Intergenic', 'Distal Intergenic'): return 'Distal Intergenic'
    return 'Unknown'


def target_features(text):
    result = [x.strip() for x in text.split(',') if x.strip()]
    invalid = set(result)-set(FEATURES)
    if not result or invalid:
        raise ValueError('Invalid TARGET_FEATURES: '+str(sorted(invalid))+'; allowed: '+', '.join(FEATURES))
    return result


def location_rule(feature, distance, accepted, max_tss_dist=0):
    if max_tss_dist < 0: raise ValueError('TARGET_MAX_TSS_DIST must be nonnegative')
    if feature in accepted: return True, feature
    try: distance = float(distance)
    except (TypeError, ValueError): return False, ''
    if max_tss_dist > 0 and math.isfinite(distance) and -max_tss_dist <= distance <= 0:
        return True, f'upstream<={max_tss_dist}bp'
    return False, ''


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''): digest.update(block)
    return digest.hexdigest()


def check_policy(annotation, up, down):
    if min(up, down) < 0: raise ValueError('Promoter extents must be nonnegative')
    meta = Path(annotation).with_name('annotation_policy.json')
    if not meta.exists():
        print('[WARN] Annotation has no recorded policy; rerun annotate with 1.2.6 to verify its promoter settings.')
        return
    data = json.loads(meta.read_text())
    if data['promoter_up'] != up or data['promoter_down'] != down:
        raise ValueError('Promoter config differs from recorded annotation; rerun annotate before downstream analysis')
    if data['annotation_sha256'] != sha256(annotation):
        raise ValueError('Annotation changed after policy validation; rerun annotate')


def prepared_annotation(annotation, up, down):
    """Apply the existing feature policy in memory, without rewriting a file."""
    if min(up, down) < 0: raise ValueError('Promoter extents must be nonnegative')
    annotation = Path(annotation)
    with annotation.open(encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = list(reader.fieldnames or [])
        if not {'annotation','distanceToTSS','geneId'}.issubset(fields):
            raise ValueError('Annotation is missing required ChIPseeker columns')
        rows = list(reader)
    counts = {}
    for row in rows:
        feature = norm_feat(row['annotation'])
        distance = float(row['distanceToTSS'])
        if not math.isfinite(distance): raise ValueError('Nonfinite distanceToTSS in annotation')
        expected = -up <= distance <= down
        if (feature == 'Promoter') != expected:
            raise ValueError('Promoter priority/window mismatch for peak '+row.get('peak_id', '?'))
        row['feature_class'] = feature
        counts[feature] = counts.get(feature, 0)+1
    if 'feature_class' not in fields: fields.append('feature_class')
    return fields, rows, counts


def validate_annotation(annotation, up, down):
    annotation = Path(annotation)
    fields, rows, counts = prepared_annotation(annotation, up, down)
    temporary = annotation.with_suffix(annotation.suffix+'.tmp')
    with temporary.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    temporary.replace(annotation)
    data = dict(schema=1, chromaflow_version='1.2.6', promoter_up=up, promoter_down=down,
                level='gene', priority=PRIORITY, annotation_sha256=sha256(annotation),
                feature_counts=counts, promoter_priority_check='passed')
    annotation.with_name('annotation_policy.json').write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
    print('[POLICY] Promoter priority/window verified for', len(rows), 'peaks; feature classes:', counts)


def compare_annotation(annotation, reconstructed, up, down):
    """Compare the finished pipeline table with the reconstructed R table.

    Normal annotate adds feature_class after R exits and may quote TSV cells.
    Reproduce that same policy in memory, then compare parsed cells exactly.
    Do not ignore IDs, coordinates, scores, feature classes or row order.
    Neither input is edited, and numeric tolerances/ID normalization are absent.
    """
    fields, rows, _ = prepared_annotation(reconstructed, up, down)
    with Path(annotation).open(encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        expected_fields, expected_rows = list(reader.fieldnames or []), list(reader)
    if expected_fields != fields:
        raise ValueError('Annotation column mismatch after feature_class derivation; '
                         f'saved={expected_fields!r}; reconstructed={fields!r}')
    if len(expected_rows) != len(rows):
        raise ValueError(f'Annotation row count mismatch: saved={len(expected_rows)}, reconstructed={len(rows)}')
    differences = []
    for i, (old, new) in enumerate(zip(expected_rows, rows), 1):
        if None in old or None in new or any(old.get(k) is None or new.get(k) is None for k in fields):
            raise ValueError(f'Malformed annotation TSV row {i}')
        for key in fields:
            if old[key] != new[key]:
                differences.append(f'row {i}, peak {old.get("peak_id", "?")!r}, column {key}: '
                                   f'saved={old[key]!r}; reconstructed={new[key]!r}')
                if len(differences) == 5:
                    break
        if len(differences) == 5:
            break
    if differences:
        raise ValueError('Annotation content mismatch (first differences):\n' + '\n'.join(differences))
    print(f'[VERIFIED] Reconstructed annotation matches all {len(rows)} rows and {len(fields)} columns '
          'after the same feature_class policy; saved files unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('validate','check','compare'))
    parser.add_argument('annotation')
    parser.add_argument('--up', type=int, required=True)
    parser.add_argument('--down', type=int, required=True)
    parser.add_argument('--reconstructed', help='Raw R annotation TSV, required for compare; read-only')
    args = parser.parse_args()
    try:
        if args.mode == 'compare':
            if not args.reconstructed:
                parser.error('compare requires --reconstructed')
            compare_annotation(args.annotation, args.reconstructed, args.up, args.down)
        else:
            (validate_annotation if args.mode == 'validate' else check_policy)(args.annotation, args.up, args.down)
    except (ValueError, OSError, KeyError) as exc: parser.exit(1, '[ERROR] '+str(exc)+'\n')


if __name__ == '__main__': main()
