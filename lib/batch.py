#!/usr/bin/env python3
"""Manifest-driven DAP-seq orchestration. No species or accession heuristics.

Scientific work is delegated to the existing ChromaFlow commands. Only explicit
sample/replicate identities may connect shared nf-core output to downstream work.
"""
from __future__ import annotations

import argparse
import collections
import copy
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time

HOME = Path(__file__).resolve().parents[1]
CLI = HOME / 'bin/chromaflow'
VERSION = '1.2.9'
STEPS = ('reference', 'nfcore', 'annotate', 'motif', 'candidates')
MACS_PROCESS = 'NFCORE_CHIPSEQ:CHIPSEQ:BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER:MACS3_CALLPEAK'
TASK_FIELDS = ('task_id', 'process', 'tag', 'status', 'exit', 'workdir')
REPORT_POLICY = 'annotation-reports-v2-reference-bounds'
HEADER = ('analysis_id sample_id role replicate layout fastq_1 fastq_2 control '
          'control_replicate target condition run_accession').split()
CONFIG_KEYS = '''PROJECT_DIR INPUT_DIR OUTPUT_DIR REFERENCE_DIR RAW_DATA_DIR
SPECIES_SHORT SPECIES_NAME GENOME_VERSION GENOME_FASTA_SOURCE GENOME_GFF3_SOURCE
GENOME_COMBINED GENOME_GFF3 GENOME_GTF GENOME_NUCLEAR CHROM_SIZES
ANNOTATION_INFO MEME_DB BLACKLIST_BED GFF_TYPE_MAP SKIP_ORGANELLE
NFCORE_PROFILE NFCORE_CHIPSEQ_VERSION NFCORE_OFFLINE NFCORE_OFFLINE_BUNDLE
NFCORE_CONDA_MIRROR NFCORE_CONDA_TIMEOUT SKIP_SPP MAX_CPUS MAX_MEMORY MAX_TIME
PROMOTER_UP PROMOTER_DOWN TARGET_DIST_KB MAKE_PLOTS FIMO_PVAL TARGET_FEATURES
TARGET_MAX_TSS_DIST MEME_SEARCHSIZE MOTIF_IDS FIMO_QVAL MOTIF_EVALUE_MAX
CENTRIMO_EVALUE_MAX KEGG_ORGANISM MACS_GSIZE READ_LENGTH
BATCH_ENV_MODE BATCH_CONDA_EXE BATCH_ENV_REF BATCH_ENV_RUN BATCH_ENV_ANNOTATE
BATCH_ENV_MOTIF BATCH_ENV_CANDIDATES BATCH_ALLOW_SE_WITH_PE_CONTROL
BATCH_FASTQ_SHA256 BATCH_MEME_THREADS'''.split()
DEFAULTS = dict(SPECIES_SHORT='Species', NFCORE_PROFILE='singularity',
                NFCORE_CHIPSEQ_VERSION='2.1.0', NFCORE_OFFLINE='false',
                MAX_CPUS='8', MAX_MEMORY='32.GB', MAX_TIME='72.h',
                PROMOTER_UP='1000', PROMOTER_DOWN='200', TARGET_DIST_KB='3',
                MAKE_PLOTS='true', FIMO_PVAL='1e-4',
                TARGET_FEATURES="Promoter,5' UTR,1st Exon", KEGG_ORGANISM='',
                BATCH_ENV_MODE='conda', BATCH_ENV_REF='ChromaFlow_ref',
                BATCH_ENV_ANNOTATE='ChromaFlow_R', BATCH_ENV_MOTIF='ChromaFlow_motif',
                BATCH_ENV_CANDIDATES='ChromaFlow_downstream',
                BATCH_ALLOW_SE_WITH_PE_CONTROL='false', BATCH_FASTQ_SHA256='false')
ID = re.compile(r'^[A-Za-z][A-Za-z0-9_]*$')


class BatchError(RuntimeError):
    pass


def fail(message):
    raise BatchError(message)


def now():
    return dt.datetime.now().astimezone().isoformat(timespec='seconds')


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=True).encode()).hexdigest()


def resume_code_ids(code):
    """Preserve stage signatures only across checked, compatible hotfixes.

    The installer records both entire runtime trees, not success stamps. Any
    later code change disables this compatibility mapping. Inputs, parameters,
    environments and output hashes still undergo the normal resume checks.
    """
    current = digest(code)
    fallback = dict.fromkeys(STEPS, current)
    release = HOME / 'release/compatibility.json'
    if release.is_file():
        try:
            data = json.loads(release.read_text(encoding='utf-8'))
            stages = data['step_code_sha256']
            if (data['schema'] != 1 or data['version'] != VERSION
                    or data['runtime'] != code or set(stages) != set(STEPS)
                    or any(not isinstance(v, str) or not re.fullmatch(r'[0-9a-f]{64}', v)
                           for v in stages.values())):
                raise ValueError('release runtime or stage identities do not match')
            print(f'[RELEASE] {VERSION}: verified compatible stage identities; normal input/environment/output checks remain active', flush=True)
            return dict(stages)
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            print(f'[WARN] Release compatibility is not applicable ({exc}); using current code identity', flush=True)
            return fallback
    patch_root = HOME / '.chromaflow-hotfixes'
    fimo = patch_root / 'fimo-zero-hits-v1/manifest.json'
    read_length = patch_root / 'read-length-v1/manifest.json'
    peak_isolation = patch_root / 'peak-isolation-v1/manifest.json'
    resume_reports = patch_root / 'resume-reports-v1/manifest.json'
    report_guard = patch_root / 'report-table-guard-v1/manifest.json'
    tss_metadata = patch_root / 'tss-idr-metadata-v1/manifest.json'

    def verified_before(manifest, expected, patch_id, allowed, compatible, added=frozenset()):
        data = json.loads(manifest.read_text(encoding='utf-8'))
        before, after = data['before_code'], data['after_code']
        changed = {name for name in set(before) | set(after) if before.get(name) != after.get(name)}
        if (data['schema'] != 1 or data['patch'] != patch_id
                or after != expected or not changed or not changed <= allowed
                or set(before) - set(after) or set(after) - set(before) != set(added)
                or data['compatible_steps'] != list(compatible)
                or any(not re.fullmatch(r'[0-9a-f]{64}', value) for value in before.values())):
            raise ValueError(patch_id + ': runtime identity mismatch')
        return before

    try:
        previous = code
        if tss_metadata.is_file():
            if not report_guard.is_file():
                raise ValueError('tss-idr-metadata-v1: missing predecessor manifest')
            previous = verified_before(tss_metadata, previous, 'tss-idr-metadata-v1',
                {'lib/batch.py', 'lib/common.sh', 'lib/steps/05_peak_annotation.R',
                 'lib/steps/05_peak_annotation.sh', 'lib/steps/tss_windows.R',
                 'lib/steps/support/peak_metadata.py', 'lib/steps/support/lookup.py'}, STEPS,
                added={'lib/steps/tss_windows.R', 'lib/steps/support/peak_metadata.py'})
        if report_guard.is_file():
            if not resume_reports.is_file():
                raise ValueError('report-table-guard-v1: missing predecessor manifest')
            previous = verified_before(report_guard, previous, 'report-table-guard-v1',
                {'lib/batch.py', 'lib/steps/05_peak_annotation.R', 'lib/steps/support/annotation_policy.py'}, STEPS)
        if resume_reports.is_file():
            if not peak_isolation.is_file():
                raise ValueError('resume-reports-v1: missing predecessor manifest')
            previous = verified_before(resume_reports, previous, 'resume-reports-v1',
                {'lib/batch.py', 'lib/help.sh', 'lib/steps/05_peak_annotation.R',
                 'lib/steps/annotation_plots.R', 'lib/steps/support/motif_selection.py'},
                STEPS, added={'lib/steps/annotation_plots.R'})
        if peak_isolation.is_file():
            if not read_length.is_file():
                raise ValueError('peak-isolation-v1: missing predecessor manifest')
            previous = verified_before(peak_isolation, previous, 'peak-isolation-v1', {'lib/batch.py'}, STEPS)
        if read_length.is_file():
            if not fimo.is_file():
                raise ValueError('read-length-v1: missing predecessor manifest')
            previous = verified_before(read_length, previous, 'read-length-v1',
                {'lib/batch.py', 'lib/steps/02_run_nfcore.sh', 'lib/steps/02_run_nfcore_conda.sh'}, STEPS)
        if not fimo.is_file():
            return fallback
        before = verified_before(fimo, previous, 'fimo-zero-hits-v1',
            {'lib/batch.py', 'lib/steps/support/fimo_tracks.py',
             'lib/steps/support/motif_selection.py'}, STEPS[:-1])
        stages = dict.fromkeys(STEPS, digest(previous))
        for step in STEPS[:-1]:
            stages[step] = digest(before)
        if tss_metadata.is_file():
            print('[HOTFIX] tss-idr-metadata-v1: verified predecessor chain; reference-bound reports and checked metadata correction', flush=True)
        elif report_guard.is_file():
            print('[HOTFIX] report-table-guard-v1: verified complete predecessor chain; scientific stages retain normal input/output checks', flush=True)
        elif resume_reports.is_file():
            print('[HOTFIX] resume-reports-v1: verified complete predecessor chain; scientific stages retain normal input/output checks', flush=True)
        elif peak_isolation.is_file():
            print('[HOTFIX] fimo-zero-hits-v1 + read-length-v1 + peak-isolation-v1: verified resume compatibility; all stages retain normal input/output checks', flush=True)
        elif read_length.is_file():
            print('[HOTFIX] fimo-zero-hits-v1 + read-length-v1: verified resume compatibility; all stages retain normal input/output checks', flush=True)
        else:
            print('[HOTFIX] fimo-zero-hits-v1: verified upstream resume compatibility; candidates use updated code', flush=True)
        return stages
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        print(f'[WARN] Hotfix resume metadata is not applicable ({exc}); verifying with current code identity', flush=True)
        return fallback


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.tmp-{os.getpid()}')
    temp.write_bytes(data.encode('utf-8'))
    temp.replace(path)


def save(path, data):
    write(path, json.dumps(data, indent=2, ensure_ascii=True) + '\n')


def read_json(path, default=None):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        return default


def tsv(path, header, rows):
    out = io.StringIO()
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerow(header)
    writer.writerows(rows)
    write(path, out.getvalue())


def clean_env():
    # Never inherit a different project's selectors or Nextflow runtime.
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('CHROMAFLOW_', 'CF_BATCH_', 'CFSEL_', 'NXF_', 'NFCORE_'))}
    env.update(PYTHONUTF8='1', PYTHONIOENCODING='utf-8', LC_ALL='C.UTF-8',
               NXF_SYNTAX_PARSER='v1')
    return env


def load_config(path):
    """Source the user's trusted Bash config, exporting only documented fields.

    Do not run common.sh's directory discovery: a batch must not adopt another
    project's manifest or derived read length. Relative inputs use PROJECT_DIR.
    """
    path = Path(path).resolve()
    if not path.is_file():
        fail(f'Config not found: {path}')
    script = '''set -e
set +u
source "$1" >&2
: "${PROJECT_DIR:=$PWD}"
: "${INPUT_DIR:=$PROJECT_DIR/inputs}"
: "${OUTPUT_DIR:=$PROJECT_DIR}"
: "${REFERENCE_DIR:=$INPUT_DIR/reference}"
: "${RAW_DATA_DIR:=$PROJECT_DIR/raw_data}"
source "$1" >&2
shift
for key in "$@"; do printf '%s\\0%s\\0' "$key" "${!key-}"; done
for item in "${ORGANELLE_SOURCES[@]}"; do printf 'ORGANELLE_SOURCES[]\\0%s\\0' "$item"; done
'''
    result = subprocess.run(['bash', '-c', script, 'batch-config', str(path), *CONFIG_KEYS],
                            env=clean_env(), cwd=path.parent, capture_output=True)
    if result.returncode:
        fail('Cannot load config: ' + result.stderr.decode(errors='replace'))
    parts = result.stdout.decode('utf-8').split('\0')
    if parts[-1] != '' or len(parts) % 2 != 1:
        fail('Config produced invalid output; use Bash variable assignments only')
    cfg = dict(DEFAULTS)
    cfg['ORGANELLE_SOURCES'] = []
    for key, value in zip(parts[0:-1:2], parts[1:-1:2]):
        if key == 'ORGANELLE_SOURCES[]':
            cfg['ORGANELLE_SOURCES'].append(value)
        elif value or key not in DEFAULTS:
            cfg[key] = value
    for key in ('PROJECT_DIR', 'INPUT_DIR', 'OUTPUT_DIR', 'REFERENCE_DIR', 'RAW_DATA_DIR'):
        p = Path(cfg[key]).expanduser()
        cfg[key] = str((path.parent / p).resolve())
    for key in ('GENOME_FASTA_SOURCE', 'GENOME_GFF3_SOURCE', 'GENOME_COMBINED',
                'GENOME_GFF3', 'GENOME_GTF', 'GENOME_NUCLEAR', 'CHROM_SIZES',
                'ANNOTATION_INFO', 'MEME_DB', 'BLACKLIST_BED', 'GFF_TYPE_MAP',
                'NFCORE_OFFLINE_BUNDLE'):
        if cfg.get(key):
            cfg[key] = str((Path(cfg['PROJECT_DIR']) / Path(cfg[key]).expanduser()).resolve())
    if cfg['NFCORE_CHIPSEQ_VERSION'] != '2.1.0':
        fail('Batch output mapping is validated for nf-core/chipseq 2.1.0 only')
    if cfg['NFCORE_PROFILE'] not in ('singularity', 'apptainer', 'docker', 'podman', 'conda', 'mamba'):
        fail('Batch requires one explicit engine profile (no institutional/combined profiles)')
    if cfg['BATCH_ENV_MODE'] not in ('conda', 'current'):
        fail('BATCH_ENV_MODE must be conda or current')
    for key in ('NFCORE_OFFLINE', 'BATCH_ALLOW_SE_WITH_PE_CONTROL', 'BATCH_FASTQ_SHA256'):
        if cfg[key] not in ('true', 'false'):
            fail(f'{key} must be true or false')
    if not re.fullmatch(r'[1-9][0-9]*', cfg['MAX_CPUS']):
        fail('MAX_CPUS must be a positive integer')
    if not re.fullmatch(r'[1-9][0-9]*(?:\.[0-9]+)?\.(?:MB|GB|TB)', cfg['MAX_MEMORY']):
        fail('MAX_MEMORY must use a positive value such as 64.GB')
    cfg['BATCH_MEME_THREADS'] = cfg.get('BATCH_MEME_THREADS') or str(min(8, int(cfg['MAX_CPUS'])))
    if not cfg['BATCH_MEME_THREADS'].isdigit() or not 1 <= int(cfg['BATCH_MEME_THREADS']) <= int(cfg['MAX_CPUS']):
        fail('BATCH_MEME_THREADS must be between 1 and MAX_CPUS')
    cfg['BATCH_ENV_RUN'] = cfg.get('BATCH_ENV_RUN') or (
        'ChromaFlow_runtime' if cfg['NFCORE_OFFLINE'] == 'true' else 'nfcore_env')
    for key in ('PROMOTER_UP', 'PROMOTER_DOWN', 'TARGET_MAX_TSS_DIST', 'MEME_SEARCHSIZE'):
        if cfg.get(key) and not re.fullmatch('[0-9]+', cfg[key]):
            fail(f'{key} must be a nonnegative integer')
    for key in ('FIMO_PVAL', 'FIMO_QVAL', 'MOTIF_EVALUE_MAX', 'CENTRIMO_EVALUE_MAX', 'TARGET_DIST_KB'):
        if not cfg.get(key):
            continue
        try:
            number = float(cfg[key])
        except ValueError:
            fail(f'{key} must be numeric')
        if not math.isfinite(number) or number < 0 or (key == 'FIMO_PVAL' and number == 0) or (
                key in ('FIMO_PVAL', 'FIMO_QVAL') and number > 1):
            fail(f'{key} is out of range')
    if cfg['MAKE_PLOTS'] not in ('true', 'false'):
        fail('MAKE_PLOTS must be true or false')
    sys.path.insert(0, str(HOME/'lib/steps/support'))
    from annotation_policy import target_features
    target_features(cfg['TARGET_FEATURES'])
    if cfg.get('READ_LENGTH'):
        fail('Remove global READ_LENGTH for batch; use optional per-row read_length overrides')
    if cfg.get('MACS_GSIZE') and not re.fullmatch(r'[1-9][0-9]*', cfg['MACS_GSIZE']):
        fail('Batch MACS_GSIZE overrides must be positive integers, not species aliases')
    if cfg['NFCORE_OFFLINE'] == 'true':
        if cfg['NFCORE_PROFILE'] not in ('singularity', 'apptainer'):
            fail('Offline batch requires singularity or apptainer')
        if not cfg.get('NFCORE_OFFLINE_BUNDLE'):
            fail('Set NFCORE_OFFLINE_BUNDLE to the installed chipseq-2.1.0 bundle')
    return cfg


def snapshot(path, cfg):
    lines = ['# Generated by ChromaFlow batch. Edit the shared project config instead.']
    for key, value in sorted(cfg.items()):
        if key == 'ORGANELLE_SOURCES':
            lines.append(key + '=(' + ' '.join(shlex.quote(v) for v in value) + ')')
        else:
            lines.append(f'{key}={shlex.quote(str(value))}')
    write(path, '\n'.join(lines) + '\n')
    return Path(path)


def fastq_info(path, full_hash=False):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        fail(f'FASTQ missing or empty: {path}')
    with path.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    lengths, identifiers = set(), []
    try:
        with (gzip.open(path, 'rt') if compressed else path.open()) as f:
            for _ in range(2000):
                name = f.readline()
                if not name:
                    break
                seq, plus, qual = f.readline().rstrip(), f.readline(), f.readline().rstrip()
                if not name.startswith('@') or not plus.startswith('+') or not seq or len(seq) != len(qual):
                    fail(f'Malformed FASTQ record in first 2000 reads: {path}')
                lengths.add(len(seq))
                identifiers.append(re.sub(r'/[12]$', '', name.split()[0]))
    except (OSError, UnicodeError, EOFError) as e:
        fail(f'Cannot read FASTQ {path}: {e}')
    if not lengths:
        fail(f'No FASTQ records: {path}')
    st = path.stat()
    result = dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns,
                  sampled_reads=len(identifiers), min_length=min(lengths), max_length=max(lengths))
    if full_hash:
        result['sha256'] = sha(path)
    return result, identifiers


def load_samples(path, cfg, only=None):
    path = Path(path).resolve()
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            fail('Sample table has missing or duplicate column names')
        missing = set(HEADER) - set(reader.fieldnames)
        if missing:
            fail('Sample table is missing columns: ' + ', '.join(sorted(missing)))
        rows = list(reader)
    if not rows:
        fail('Sample table is empty')
    samples, analyses, seen_files = {}, {}, set()
    for n, row in enumerate(rows, 2):
        if None in row or any(v is None for v in row.values()):
            fail(f'Malformed TSV row {n}')
        if any(v != v.strip() for v in row.values()):
            fail(f'Leading/trailing whitespace in TSV row {n}; IDs are never normalized')
        sid, aid = row['sample_id'], row['analysis_id']
        if not ID.fullmatch(sid) or '_REP' in sid:
            fail(f'Invalid sample_id {sid!r}: use letters, digits, underscores; reserve _REP')
        if row['role'] not in ('target', 'control') or row['layout'] not in ('SE', 'PE'):
            fail(f'Row {n}: role must be target/control and layout SE/PE')
        if not re.fullmatch('[1-9][0-9]*', row['replicate']):
            fail(f'Row {n}: replicate must be a positive biological-replicate number')
        if row['role'] == 'target':
            if aid != sid or not row['target'] or not row['condition']:
                fail(f'Row {n}: target rows require analysis_id=sample_id, target and condition')
            signature = (row['target'], row['condition'], row['layout'])
            if aid in analyses and analyses[aid] != signature:
                fail(f'{aid}: target/condition/layout must agree across replicates')
            analyses[aid] = signature
            if not row['control'] or not re.fullmatch('[1-9][0-9]*', row['control_replicate']):
                fail(f'{aid}: an explicit control and control_replicate are required')
        elif aid or row['control'] or row['control_replicate']:
            fail(f'Control row {sid}: leave analysis_id/control/control_replicate empty')
        signature = (row['role'], row['layout'])
        if sid in samples and samples[sid] != signature:
            fail(f'{sid}: conflicting role or layout')
        samples[sid] = signature
        if row['layout'] == 'PE' and not row['fastq_2']:
            fail(f'{sid}: PE requires both mates; no automatic conversion to SE')
        if row['layout'] == 'SE' and row['fastq_2']:
            fail(f'{sid}: SE requires an empty fastq_2')
        for key in ('fastq_1', 'fastq_2'):
            if not row[key]:
                continue
            p = (Path(cfg['RAW_DATA_DIR']) / row[key]).resolve()
            if any(c.isspace() or c in ',\"' for c in str(p)) or not str(p).endswith(('.fastq.gz', '.fq.gz')):
                fail(f'nf-core 2.1.0 requires .fastq.gz/.fq.gz paths without whitespace, commas or quotes: {p}')
            if str(p) in seen_files:
                fail(f'FASTQ is assigned more than once: {p}; list shared controls once')
            seen_files.add(str(p))
            row[key] = str(p)
    for sid in samples:
        reps = sorted({int(r['replicate']) for r in rows if r['sample_id'] == sid})
        if reps != list(range(1, len(reps) + 1)):
            fail(f'{sid}: biological replicate numbers must be consecutive from 1')
        for rep in reps:
            assignments = {(r['control'], r['control_replicate']) for r in rows
                           if r['sample_id'] == sid and int(r['replicate']) == rep}
            if len(assignments) != 1:
                fail(f'{sid} replicate {rep}: technical runs must share one control assignment')
    if not analyses:
        fail('No target analyses in sample table')
    # Exact selection only; never use a substring/regular expression for IDs.
    chosen = set(only.split(',')) if only else set(analyses)
    if chosen - set(analyses):
        fail('Unknown --only analysis IDs: ' + ', '.join(sorted(chosen - set(analyses))))
    target_rows = [r for r in rows if r['analysis_id'] in chosen]
    refs = {(r['control'], r['control_replicate']) for r in target_rows}
    controls = [r for r in rows if r['role'] == 'control' and r['sample_id'] in {x[0] for x in refs}]
    control_keys = {(r['sample_id'], r['replicate']) for r in controls}
    if refs - control_keys:
        fail(f'Control sample/replicate not found: {sorted(refs - control_keys)}')
    for r in target_rows:
        layout = samples[r['control']][1]
        if r['layout'] == 'PE' and layout != 'PE':
            fail(f"{r['analysis_id']}: PE treatment requires a PE control for MACS3 BAMPE")
        if r['layout'] == 'SE' and layout == 'PE' and cfg['BATCH_ALLOW_SE_WITH_PE_CONTROL'] != 'true':
            fail(f"{r['analysis_id']}: SE treatment/PE control requires explicit BATCH_ALLOW_SE_WITH_PE_CONTROL=true")
    info, lengths = {}, collections.defaultdict(set)
    for row in target_rows + controls:
        one, ids1 = fastq_info(row['fastq_1'], cfg['BATCH_FASTQ_SHA256'] == 'true')
        info[row['fastq_1']] = one
        detected = one['max_length']
        if row['fastq_2']:
            two, ids2 = fastq_info(row['fastq_2'], cfg['BATCH_FASTQ_SHA256'] == 'true')
            info[row['fastq_2']] = two
            if ids1 != ids2:
                fail(f"{row['sample_id']}: mate names/order differ in first 2000 reads")
            detected = max(detected, two['max_length'])
        override = row.get('read_length', '')
        if override and not re.fullmatch('[1-9][0-9]*', override):
            fail('Optional read_length must be a positive integer')
        row['effective_read_length'] = int(override or detected)
        if row['role'] == 'target':
            lengths[row['analysis_id']].add(row['effective_read_length'])
    for aid, values in lengths.items():
        if len(values) != 1:
            fail(f'{aid}: different read lengths across runs; supply reviewed per-row read_length overrides')
    return target_rows + controls, info, {a: next(iter(v)) for a, v in lengths.items()}, sorted(analyses)


def source_inputs(cfg):
    sources = {}
    for name, fallback in (('GENOME_FASTA_SOURCE', 'GENOME_COMBINED'), ('GENOME_GFF3_SOURCE', 'GENOME_GFF3')):
        value = cfg.get(name) or cfg.get(fallback)
        if not value or not Path(value).is_file() or Path(value).stat().st_size == 0:
            fail(f'Set {name} to an existing local reference file')
        sources[name] = dict(path=value, sha256=sha(value))
    for name in ('GFF_TYPE_MAP', 'BLACKLIST_BED', 'ANNOTATION_INFO', 'MEME_DB'):
        if cfg.get(name):
            if not Path(cfg[name]).is_file():
                fail(f'{name} does not exist: {cfg[name]}')
            sources[name] = dict(path=cfg[name], sha256=sha(cfg[name]))
    for value in cfg['ORGANELLE_SOURCES']:
        label, sep, filename = value.partition(':')
        if not sep or not filename or not Path(filename).is_file():
            fail('Batch ORGANELLE_SOURCES must use name:/absolute/local/file; no network downloads')
        sources['organelle:' + label] = dict(path=filename, sha256=sha(filename))
    return sources


def plan_project(args):
    cfg = load_config(args.config)
    rows, fastqs, lengths, all_ids = load_samples(args.samples, cfg, args.only)
    sources = source_inputs(cfg)
    groups, analyses = {}, {}
    for aid, k in sorted(lengths.items()):
        targets = [r for r in rows if r['analysis_id'] == aid]
        # Split by treatment layout and k; each group has one MACS effective size.
        gid = targets[0]['layout'].lower() + '_k' + str(k)
        groups.setdefault(gid, dict(id=gid, k=k, analyses=[], rows=[]))['analyses'].append(aid)
        analyses[aid] = dict(id=aid, group=gid, target=targets[0]['target'],
                            condition=targets[0]['condition'],
                            biological_replicates=sorted({int(r['replicate']) for r in targets}),
                            rows=targets)
    for g in groups.values():
        targets = [r for r in rows if r['analysis_id'] in g['analyses']]
        controls = {r['control'] for r in targets}
        # Include all replicates of a selected control sample so nf-core's
        # consecutive-replicate invariant remains true even for --only.
        g['rows'] = targets + [r for r in rows if r['role'] == 'control' and r['sample_id'] in controls]
    code = {p.relative_to(HOME).as_posix(): sha(p) for d in ('lib', 'bin')
            for p in sorted((HOME/d).rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
    return dict(schema=1, version=VERSION, created=now(), config=cfg,
                config_path=str(Path(args.config).resolve()),
                samples_path=str(Path(args.samples).resolve()), samples_sha256=sha(args.samples),
                all_analysis_ids=all_ids, analyses=analyses, groups=groups, sources=sources,
                fastqs=fastqs, code_sha256=digest(code), step_code_sha256=resume_code_ids(code))


def env_prefix(cfg, step):
    if cfg['BATCH_ENV_MODE'] == 'current':
        return []
    conda = cfg.get('BATCH_CONDA_EXE') or os.environ.get('CONDA_EXE') or shutil.which('conda')
    if not conda or not Path(conda).is_file():
        fail('Conda executable not found. Set BATCH_CONDA_EXE or initialize Conda first')
    name = cfg['BATCH_ENV_' + dict(reference='REF', nfcore='RUN', annotate='ANNOTATE',
                                 motif='MOTIF', candidates='CANDIDATES')[step]]
    return [str(Path(conda).resolve()), 'run', '--no-capture-output',
            '-p' if '/' in name else '-n', name]


def command(cfg, step, subcommand, config, *args):
    return env_prefix(cfg, step) + ['bash', str(CLI), subcommand, '-c', str(config), *args]


def stable_environment(identity):
    """Exclude one known, explicitly ignored MPI loader diagnostic, not versions.

    Keep the raw probe in the environment log. Unknown warnings, executable
    paths, versions and nonzero probe exits remain significant.
    """
    if not isinstance(identity, dict) or not isinstance(identity.get('tools'), str):
        return identity
    result = dict(identity)
    ignored = re.compile(r'^\[[^\]\r\n:]+:\d+\] mca_base_component_repository_open: '
                         r'unable to open mca_btl_openib: librdmacm\.so\.1: '
                         r'cannot open shared object file: No such file or directory \(ignored\)(?:\r?\n|$)', re.M)
    result['tools'] = ignored.sub('', identity['tools'])
    return result


def preflight(cfg, steps=STEPS):
    required = dict(reference=['python3', 'samtools', 'gffread', 'bedtools', 'unique-kmers.py'],
                    nfcore=['python3', 'java'], annotate=['python3', 'Rscript'],
                    motif=['python3', 'meme-chip', 'fasta-get-markov', 'bedtools'], candidates=['python3'])
    if cfg.get('MACS_GSIZE'):
        required['reference'].remove('unique-kmers.py')
    if cfg['NFCORE_OFFLINE'] != 'true':
        required['nfcore'].append('nextflow')
    engine = cfg['NFCORE_PROFILE']
    reports = {}
    for step in steps:
        probe = '''set -eu
printf 'CONDA_PREFIX=%s\\n' "${CONDA_PREFIX:-current}"
for tool in "$@"; do command -v "$tool" || exit 1; done
python3 --version
'''
        if step == 'nfcore':
            probe += 'java -version\n'
            if engine in ('apptainer', 'singularity'):
                probe += '(command -v apptainer && apptainer --version) || (command -v singularity && singularity --version)\n'
            elif engine in ('docker', 'podman'):
                probe += engine + ' --version\n'
            if cfg['NFCORE_OFFLINE'] == 'true':
                probe += "printf '[BUNDLE_RUNTIME_CHECK]\\n'\n"
                probe += 'bash ' + shlex.quote(str(CLI)) + ' offline check --bundle ' + shlex.quote(cfg['NFCORE_OFFLINE_BUNDLE']) + ' --runtime-check\n'
            else:
                probe += 'nextflow -version\n'
        elif step == 'reference':
            probe += 'samtools --version | head -1\ngffread --version\nbedtools --version\n'
            if not cfg.get('MACS_GSIZE'):
                probe += 'unique-kmers.py --version\n'
        elif step == 'annotate':
            probe += '''Rscript -e 'for (p in c("ChIPseeker","GenomicFeatures","rtracklayer","ggplot2")) { if (!requireNamespace(p, quietly=TRUE)) stop(paste("Missing R package",p)); cat(p, as.character(packageVersion(p)), "\\n") }'\n'''
        elif step == 'motif':
            probe += 'meme -version\n'
        cmd = env_prefix(cfg, step) + ['bash', '-c', probe, 'batch-check', *required[step]]
        print(f'[CHECK] {step}: {shlex.join(cmd[:len(env_prefix(cfg, step))]) or "current environment"}', flush=True)
        p = subprocess.run(cmd, env=clean_env(), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(p.stdout, end='', flush=True)
        if p.returncode:
            fail(f'Environment check failed for {step}; check BATCH_ENV_* and setup --mode offline')
        # Runtime-check logs contain temporary directories and run names. Keep
        # them for diagnostics, but do not turn them into changing cache keys.
        identity = dict(environment=env_prefix(cfg, step), tools=p.stdout.split('[BUNDLE_RUNTIME_CHECK]')[0])
        if step == 'motif':
            identity = stable_environment(identity)
        if step == 'nfcore' and cfg['NFCORE_OFFLINE'] == 'true':
            identity['bundle_manifest_sha256'] = sha(Path(cfg['NFCORE_OFFLINE_BUNDLE'])/'chromaflow-bundle.json')
        reports[step] = dict(identity=identity, log=p.stdout)
    return reports


def nf_sheet(path, rows):
    tsv_rows = []
    for r in rows:
        tsv_rows.append([r['sample_id'], r['fastq_1'], r['fastq_2'], r['replicate'],
                         r['analysis_id'], r['control'], r['control_replicate']])
    out = io.StringIO()
    w = csv.writer(out, lineterminator='\n')
    w.writerow(['sample', 'fastq_1', 'fastq_2', 'replicate', 'antibody', 'control', 'control_replicate'])
    w.writerows(tsv_rows)
    write(path, out.getvalue())


def artifact_files(paths):
    result = {}
    for path in paths:
        p = Path(path)
        if not p.exists():
            fail(f'Required output missing: {p}')
        files = sorted(q for q in p.rglob('*') if q.is_file()) if p.is_dir() else [p]
        if not files:
            fail(f'Output directory contains no files: {p}')
        for f in files:
            result[str(f)] = sha(f)
    return result


def valid_record(record, signature):
    if not record or record.get('status') != 'complete' or signature not in (
            record.get('signature'), record.get('resume_signature')):
        return False
    try:
        return bool(record['outputs']) and all(Path(p).is_file() and sha(p) == h
                                                for p, h in record['outputs'].items())
    except (OSError, KeyError):
        return False


def nfcore_executor_config(cfg, trace):
    """Keep a failed peak task out of output channels while other samples finish.

    Native exit zero is required before interpreting the trace. Only this exact
    process gets an override; ChromaFlow reports ignored failures and exits one.
    The scientific MACS command and all other process strategies are unchanged.
    """
    return (f"executor {{ name = 'local'; cpus = {cfg['MAX_CPUS']}; memory = {cfg['MAX_MEMORY']}; queueSize = {cfg['MAX_CPUS']} }}\n"
            f"process {{ withName: '{MACS_PROCESS}' {{ errorStrategy = 'ignore' }} }}\n"
            "workflow.failOnIgnore = false\n"
            "trace {\n    enabled = true\n    overwrite = false\n"
            f"    file = {json.dumps(str(trace))}\n"
            f"    fields = '{','.join(TASK_FIELDS)}'\n}}\n")


def peak_task_records(trace, expected):
    """Require fresh, unambiguous task outcomes for every requested replicate."""
    trace = Path(trace)
    if not trace.is_file():
        fail(f'Nextflow task trace is missing: {trace}')
    found = {}
    with trace.open(encoding='utf-8') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not set(TASK_FIELDS) <= set(fields):
            fail(f'Invalid Nextflow task trace header: {trace}')
        for row in reader:
            if None in row or any(row.get(k) is None for k in TASK_FIELDS):
                fail(f'Malformed Nextflow task trace row: {trace}')
            if row['process'] != MACS_PROCESS:
                continue
            tag = row['tag']
            if tag not in expected or tag in found:
                fail(f'Unexpected or duplicate MACS3 task in trace: {tag}; see {trace}')
            if row['status'] not in ('COMPLETED', 'CACHED', 'FAILED'):
                fail(f'Unfinished MACS3 task: {tag} ({row["status"]}); see {trace}')
            if row['status'] != 'FAILED' and row['exit'] != '0':
                fail(f'Inconsistent MACS3 task exit status: {tag}; see {trace}')
            found[tag] = row
    missing = set(expected) - set(found)
    if missing:
        fail(f'MACS3 tasks absent from current trace: {", ".join(sorted(missing))}; see {trace}')
    return found


def peak_task_failure(row, work_root, trace):
    """Classify a model error only from its own task stderr, without changing it."""
    result = dict(status='peak_calling_failed', task=row['tag'], task_exit=row['exit'],
                  reason=f'{row["tag"]}: MACS3 failed (exit {row["exit"]}); see {trace}')
    work = Path(row['workdir'])
    if not work.is_absolute() or not work.resolve().is_relative_to(Path(work_root).resolve()):
        return result
    error = work/'.command.err'
    if not error.is_file() or not error.resolve().is_relative_to(Path(work_root).resolve()):
        return result
    try:
        with error.open('rb') as stream:
            stream.seek(max(0, error.stat().st_size - 128 * 1024))
            text = stream.read().decode('utf-8', errors='replace')
    except OSError:
        return result
    result['diagnostic_log'] = str(error)
    result['reason'] = f'{row["tag"]}: MACS3 failed (exit {row["exit"]}); see {error}'
    match = re.search(r'MACS3 needs at least (\d+) paired peaks[^\n]*can only find (\d+)', text)
    if match and row['exit'] == '1':
        minimum, observed = map(int, match.groups())
        result.update(status='model_failed', minimum_model_pairs=minimum, observed_model_pairs=observed,
            reason=f'{row["tag"]}: MACS3 model has {observed} paired peaks; requires {minimum}; see {error}')
    return result


class Runner:
    def __init__(self, plan, resume, until):
        self.plan, self.cfg, self.resume, self.until = plan, plan['config'], resume, until
        self.root = Path(self.cfg['OUTPUT_DIR'])
        self.state = self.root / '.chromaflow/batch'
        self.results = self.root / 'results'
        self.logs = self.root / 'logs'
        self.state.mkdir(parents=True, exist_ok=True)
        self.status_path = self.state / 'status.json'
        self.records = read_json(self.status_path, {})
        self.envs = {}
        self.run_id = dt.datetime.now().strftime('%Y%m%d_%H%M%S') + f'_{os.getpid()}_{time.time_ns()}'
        self.errors = []
        self.completed_this_run = set()
        self.signature_aliases = {}
        self.legacy_motif_envs = None

    def signature(self, step, inputs):
        environment = self.envs.get(step, '')
        if isinstance(environment, dict):
            environment = environment['identity']
        code = self.plan.get('step_code_sha256', {}).get(step, self.plan['code_sha256'])
        if step == 'motif':
            environment = stable_environment(environment)
        def key(env):
            return digest(dict(step=step, inputs=inputs, code=code, environment=env))
        result = key(environment)
        if step == 'motif' and self.resume:
            if self.legacy_motif_envs is None:
                self.legacy_motif_envs = []
                for path in sorted((self.state/'environments').glob('*.json')):
                    try:
                        old = read_json(path).get('motif', {})
                        if isinstance(old, dict) and 'identity' in old:
                            self.legacy_motif_envs.append(old['identity'])
                    except (OSError, ValueError, AttributeError):
                        print(f'[WARN] Cannot read historical environment record: {path}', flush=True)
            self.signature_aliases[result] = {key(old) for old in self.legacy_motif_envs
                if old != environment and stable_environment(old) == environment}
        return result

    def execute(self, cmd, log, cwd):
        log.parent.mkdir(parents=True, exist_ok=True)
        cwd.mkdir(parents=True, exist_ok=True)
        print(f'[LOG] {log}', flush=True)
        env = clean_env()
        # Serial stages share an explicit CPU budget; prevent BLAS oversubscription.
        env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        with log.open('w', encoding='utf-8') as f:
            def emit(line):
                print(line, end='', flush=True)
                f.write(line); f.flush()
            emit(f'[START] {now()}\n[VERSION] {VERSION}\n[COMMAND] {shlex.join(cmd)}\n')
            p = subprocess.Popen(cmd, cwd=cwd, env=env, text=True, errors='replace', start_new_session=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            try:
                for line in p.stdout:
                    emit(line)
                rc = p.wait()
            except BaseException:
                os.killpg(p.pid, signal.SIGTERM)
                try:
                    p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL); p.wait()
                raise
            finally:
                p.stdout.close()
            emit(f'[END] {now()} exit_code={rc}\n')
        if rc:
            fail(f'Command exited {rc}; see {log}')

    def stage(self, key, signature, action, allow_partial=False):
        old = self.records.get(key)
        if (self.resume and old and old.get('signature') in self.signature_aliases.get(signature, set())
                and valid_record(old, old['signature'])):
            old['resume_signature'] = signature
            old['resume_identity_policy'] = 'ignored-openib-diagnostic-v1'
            save(self.status_path, self.records)
            print(f'[RESUME] {key}: migrated ignored MPI diagnostic; inputs, tools and outputs verified', flush=True)
        if (self.resume or key in self.completed_this_run) and valid_record(old, signature):
            print(f'[RESUME] {key}: verified outputs', flush=True)
            return old
        reason = ('resume not requested' if not self.resume else
                  'no previous record' if not old else
                  'previous status is ' + str(old.get('status')) if old.get('status') != 'complete' else
                  'inputs, parameters, environment or code changed' if signature not in (old.get('signature'), old.get('resume_signature')) else
                  'recorded output is missing or changed')
        print(f'[RUN] {key}: {reason}', flush=True)
        self.records[key] = dict(status='running', signature=signature, started=now())
        save(self.status_path, self.records)
        try:
            paths, details = action()
            status = 'partial' if allow_partial and details.get('failures') else 'complete'
            record = dict(status=status, signature=signature, finished=now(),
                          outputs=artifact_files(paths), details=details)
            self.records[key] = record
            self.completed_this_run.add(key)
            save(self.status_path, self.records)
            return record
        except (BatchError, OSError, ValueError) as e:
            self.records[key] = dict(status='failed', signature=signature, finished=now(), error=str(e))
            save(self.status_path, self.records)
            raise

    def archive(self, path):
        """Move only generated output under this project; never touch raw/reference sources."""
        path = Path(path)
        if not path.exists():
            return
        allowed = (self.results / 'analyses', self.results / 'nfcore')
        if path.is_symlink() or not any(path.resolve().is_relative_to(p.resolve()) and path != p for p in allowed):
            fail(f'Refusing to archive outside owned batch output: {path}')
        dest = self.state / 'archive' / self.run_id / path.relative_to(self.results)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            fail(f'Archive destination already exists: {dest}')
        path.rename(dest)
        print(f'[ARCHIVE] {dest}', flush=True)

    def shared_reference(self):
        keys = ('GENOME_FASTA_SOURCE', 'GENOME_GFF3_SOURCE', 'GFF_TYPE_MAP')
        identity = {k: v for k, v in self.plan['sources'].items() if k in keys or k.startswith('organelle:')}
        sig = self.signature('reference', [identity, self.cfg['SPECIES_SHORT'], self.cfg.get('SKIP_ORGANELLE')])
        ref = self.root / 'reference/batch' / sig[:16]
        rt = self.state / 'reference' / sig[:16]
        cfg = dict(self.cfg, OUTPUT_DIR=str(rt), REFERENCE_DIR=str(ref),
                   GENOME_FASTA_SOURCE=str(ref/'nuclear.fa'),
                   GENOME_GFF3_SOURCE=self.plan['sources']['GENOME_GFF3_SOURCE']['path'],
                   GENOME_NUCLEAR=str(ref/'nuclear.fa'),
                   CF_BATCH_LOG_DIR=str(self.logs/'batch/reference'), MACS_GSIZE='', READ_LENGTH='')
        for k in ('GENOME_COMBINED', 'GENOME_GFF3', 'GENOME_GTF', 'CHROM_SIZES'):
            cfg[k] = ''
        cfile = snapshot(rt / 'reference.config', cfg)
        def action():
            # Indexes must be written beside an owned FASTA, never beside a
            # potentially read-only/shared source outside this project.
            ref.mkdir(parents=True, exist_ok=True)
            original = Path(self.plan['sources']['GENOME_FASTA_SOURCE']['path'])
            with original.open('rb') as f:
                compressed = f.read(2) == b'\x1f\x8b'
            staged = ref/'nuclear.fa.partial'
            with (gzip.open(original, 'rb') if compressed else original.open('rb')) as source:
                with staged.open('wb') as target:
                    shutil.copyfileobj(source, target, 4 * 1024 * 1024)
            staged.replace(ref/'nuclear.fa')
            self.execute(command(cfg, 'reference', 'normalize-gff', cfile, '--engine', 'builtin', '--force'),
                         self.logs/f'batch/{self.run_id}_normalize.log', rt)
            self.execute(command(cfg, 'reference', 'reference', cfile, '--no-measure', '--strict', '--force'),
                         self.logs/f'batch/{self.run_id}_reference.log', rt)
            manifest = rt / '.chromaflow/reference.manifest'
            values = dict(line.split('=', 1) for line in manifest.read_text().splitlines() if '=' in line)
            out = {k: values[k] for k in ('GENOME_COMBINED', 'GENOME_GFF3', 'GENOME_GTF', 'CHROM_SIZES')}
            self.check_sequence_names(out)
            return [*out.values(), out['GENOME_COMBINED'] + '.fai'], out
        rec = self.stage('reference', sig, action)
        return rec['details'], sig

    @staticmethod
    def check_sequence_names(ref):
        names = {s.split('\t')[0] for s in Path(ref['CHROM_SIZES']).read_text().splitlines()}
        with Path(ref['GENOME_GFF3']).open() as f:
            gff_names = {s.split('\t')[0] for s in f if s.strip() and not s.startswith('#')}
        if gff_names - names:
            fail('GFF3 seqids absent from FASTA: ' + ', '.join(sorted(gff_names-names)[:20]))

    def gsize(self, ref, refsig, k):
        sig = self.signature('reference', [refsig, k, self.cfg.get('MACS_GSIZE')])
        out = self.state / 'gsize' / (sig + '.txt')
        def action():
            if self.cfg.get('MACS_GSIZE'):
                value, method = self.cfg['MACS_GSIZE'], 'explicit_config'
                write(out, value + ' ' + method + '\n')
            else:
                # Fail on a missing/failed estimator: never silently substitute non-N size.
                script = '''set -euo pipefail
source "$1/lib/common.sh"
command -v unique-kmers.py >/dev/null
result=$(measure_effective_genome_size "$2" "$3")
[[ "$result" == *'unique-kmers'* ]] || { echo '[ERROR] Unique-kmer estimation failed'; exit 1; }
printf '%s\\n' "$result" > "$4"
'''
                out.parent.mkdir(parents=True, exist_ok=True)
                cmd = env_prefix(self.cfg, 'reference') + ['bash', '-c', script, 'batch-gsize',
                      str(HOME), ref['GENOME_COMBINED'], str(k), str(out)]
                self.execute(cmd, self.logs/f'batch/{self.run_id}_gsize_k{k}.log', self.state)
                value, method = out.read_text().strip().split(' ', 1)
            if not value.isdigit() or int(value) <= 0:
                fail('Invalid effective genome size')
            return [out], dict(value=value, method=method, k=k, reference_signature=refsig)
        return self.stage(f'gsize/k{k}', sig, action)['details']

    def group_run(self, group, ref, refsig, index):
        gid = group['id']
        rt = self.state / 'runs' / gid
        dest = self.results / 'nfcore' / gid
        executor = rt / 'executor.config'
        trace = self.logs/'batch'/gid/(self.run_id+'_tasks.tsv')
        write(executor, nfcore_executor_config(self.cfg, trace))
        size = self.gsize(ref, refsig, group['k'])
        cfg = dict(self.cfg, **ref, OUTPUT_DIR=str(rt),
                   MACS_GSIZE=size['value'], READ_LENGTH=str(group['k']),
                   CF_BATCH_NFCORE_OUT=str(dest), CF_BATCH_EXECUTOR_CONFIG=str(executor),
                   CF_BATCH_WORK_DIR=str(self.root/'work/batch'/gid),
                   CF_BATCH_LOG_DIR=str(self.logs/'batch'/gid), CF_BATCH_TASK_TRACE=str(trace))
        cfile = snapshot(rt / 'nfcore.config', cfg)
        nf_sheet(rt / 'samplesheets/samplesheet_dapseq.csv', group['rows'])
        # Input signatures include reference content, explicit control assignments,
        # observed FASTQ identities and effective parameters, not a success stamp.
        selected = {r[k]: self.plan['fastqs'][r[k]] for r in group['rows']
                    for k in ('fastq_1', 'fastq_2') if r[k]}
        params = {k: cfg.get(k, '') for k in ('NFCORE_PROFILE', 'NFCORE_OFFLINE',
                  'NFCORE_OFFLINE_BUNDLE', 'MAX_CPUS', 'MAX_MEMORY', 'MAX_TIME', 'SKIP_SPP',
                  'NFCORE_CONDA_MIRROR', 'NFCORE_CONDA_TIMEOUT', 'BLACKLIST_BED')}
        sig = self.signature('nfcore', [group['rows'], selected, refsig, size, params,
                                      self.plan['sources'].get('BLACKLIST_BED')])
        def action():
            self.archive(dest)
            trace.parent.mkdir(parents=True, exist_ok=True)
            if trace.exists():
                fail(f'Refusing to reuse a task trace from another execution: {trace}')
            extra = ['--save_reference', '--aligner', 'bwa', '--skip_deseq2_qc']
            if index:
                extra += ['--bwa_index', str(index)]
            self.execute(command(cfg, 'nfcore', 'run', cfile, '--type', 'dapseq', '--yes', '--', *extra),
                         self.logs/f'batch/{self.run_id}_{gid}_nfcore.log', rt)
            expected = {f'{aid}_REP{rep}' for aid in group['analyses']
                        for rep in self.plan['analyses'][aid]['biological_replicates']}
            tasks = peak_task_records(trace, expected)
            paths = [trace]
            mapping, failures = {}, {}
            for aid in group['analyses']:
                a = self.plan['analyses'][aid]
                failed = [peak_task_failure(tasks[f'{aid}_REP{rep}'], cfg['CF_BATCH_WORK_DIR'], trace)
                          for rep in a['biological_replicates'] if tasks[f'{aid}_REP{rep}']['status'] == 'FAILED']
                if failed:
                    failures[aid] = dict(status='model_failed' if all(t['status'] == 'model_failed' for t in failed)
                                         else 'peak_calling_failed', candidates='NA', peak_count='NA',
                                         reason='; '.join(t['reason'] for t in failed), tasks=failed, trace=str(trace))
                    continue
                base = dest / 'bwa/merged_library/macs3/narrow_peak'
                peaks = [base/f'{aid}_REP{rep}_peaks.narrowPeak' for rep in a['biological_replicates']]
                consensus = base/'consensus'/aid/f'{aid}.consensus_peaks.bed'
                missing = [p for p in peaks + ([consensus] if len(peaks) > 1 else []) if not p.is_file()]
                if missing:
                    failures[aid] = dict(status='nfcore_incomplete', candidates='NA', peak_count='NA',
                        reason='Successful peak task but expected output is missing: '+', '.join(map(str, missing)), trace=str(trace))
                    continue
                # A singleton commonly has no native consensus in nf-core 2.1.0.
                kind = 'replicate_consensus' if len(peaks) > 1 else 'single_replicate_no_idr'
                source = consensus if len(peaks) > 1 else peaks[0]
                mapping[aid] = dict(replicates=list(map(str, peaks)), peaks=str(source), peakset_kind=kind,
                                    nfcore_group=gid, biological_replicates=len(peaks), idr=False)
                paths.extend(peaks + ([source] if source not in peaks else []))
            manifest = rt/'outputs.json'
            save(manifest, mapping)
            audit = rt/'peak_tasks.json'
            save(audit, dict(trace=str(trace), tasks=tasks, failures=failures))
            paths.extend([manifest, audit])
            indexes = sorted((dest/'genome/index').rglob('*.bwt'))
            index_dir = ''
            if len(indexes) == 1 and all(indexes[0].with_suffix(s).is_file() for s in ('.amb', '.ann', '.pac', '.sa')):
                index_dir = str(indexes[0].parent)
                paths.extend(indexes[0].with_suffix(s) for s in ('.amb', '.ann', '.bwt', '.pac', '.sa'))
            if failures:
                print(f'[PARTIAL] {gid}: {len(mapping)} analyses have accepted peak outputs; {len(failures)} failed/incomplete; batch will return nonzero', flush=True)
            return paths, dict(mapping=mapping, failures=failures, bwa_index=index_dir, reference_signature=refsig, trace=str(trace))
        record = self.stage('nfcore/' + gid, sig, action, allow_partial=True)
        return record['details'], sig

    def downstream(self, aid, mapping, ref, nf_sig):
        root = self.results/'analyses'/aid
        rt = self.state/'analyses'/aid
        peaks = rt/'selected_peaks.bed'
        files = [Path(p) for p in mapping['replicates']]
        sig = digest([nf_sig, mapping, {str(p): sha(p) for p in files}])
        # Write a bounded per-analysis manifest even for no-peak/too-small cases.
        root.mkdir(parents=True, exist_ok=True)
        save(root/'inputs.json', dict(analysis=self.plan['analyses'][aid], **mapping))
        source = Path(mapping['peaks'])
        if mapping['peakset_kind'] == 'single_replicate_no_idr':
            lines = []
            for line in source.read_text().splitlines():
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                cols = line.split('\t')
                if len(cols) != 10 or int(cols[2]) <= int(cols[1]):
                    fail(f'Malformed narrowPeak output: {source}')
                cols[4] = cols[6]  # Use MACS signalValue while preserving the summit column.
                lines.append('\t'.join(cols))
            write(peaks, '\n'.join(lines) + ('\n' if lines else ''))
        else:
            combined = rt/'replicate_signal.narrowPeak'
            write(combined, ''.join(p.read_text() for p in files))
            # This utility has a CLI-style main; invoke using the orchestrator's Python.
            if not source.stat().st_size:
                write(peaks, '')
            else:
                p = subprocess.run([sys.executable, str(HOME/'lib/steps/support/consensus_signal.py'),
                                    str(source), str(combined), str(peaks)], capture_output=True, text=True)
                if p.returncode:
                    fail('Consensus signal extraction failed: ' + p.stderr)
        n = sum(1 for line in peaks.read_text().splitlines() if line.strip())
        if n == 0:
            return dict(status='no_peaks', peak_count=0, candidates='NA', reason='No MACS peaks; downstream not run')
        cfg = dict(self.cfg, **ref, OUTPUT_DIR=str(rt), CF_BATCH_ANALYSIS_DIR=str(root),
                   CF_BATCH_PEAK_FILE=str(peaks), CF_BATCH_LOG_DIR=str(self.logs/'analyses'/aid),
                   DEFAULT_PEAKSET='consensus', MEME_THREADS=self.cfg['BATCH_MEME_THREADS'])
        cfile = snapshot(rt/'analysis.config', cfg)
        previous = sig
        for step in ('annotate', 'motif', 'candidates'):
            if STEPS.index(step) > STEPS.index(self.until):
                return dict(status='stopped_after_' + self.until, peak_count=n, candidates='NA')
            if step == 'motif' and n < 10:
                return dict(status='insufficient_peaks', peak_count=n, candidates='NA',
                            reason='Fewer than 10 peaks; motif discovery skipped by ChromaFlow; annotation is available')
            # Keep filter changes from needlessly recomputing peak calling or motifs.
            keys = dict(annotate=['PROMOTER_UP', 'PROMOTER_DOWN', 'TARGET_DIST_KB', 'MAKE_PLOTS', 'KEGG_ORGANISM'],
                        motif=['MEME_SEARCHSIZE', 'BATCH_MEME_THREADS'],
                        candidates=['TARGET_FEATURES', 'TARGET_MAX_TSS_DIST', 'FIMO_PVAL', 'FIMO_QVAL',
                                    'MOTIF_IDS', 'MOTIF_EVALUE_MAX', 'CENTRIMO_EVALUE_MAX'])
            step_sig = self.signature(step, [previous, {k: cfg.get(k, '') for k in keys[step]},
                                             self.plan['sources'].get('ANNOTATION_INFO') if step == 'annotate' else
                                             self.plan['sources'].get('MEME_DB') if step == 'motif' else None])
            dest = root / ('annotation' if step == 'annotate' else step)
            def action(step=step, dest=dest):
                self.archive(dest)
                self.execute(command(cfg, step, step, cfile, '--type', 'dapseq', '--peakset', 'consensus'),
                             self.logs/f'analyses/{aid}/{self.run_id}_{step}.log', rt)
                required = dict(annotate=[dest/'peak_annotation_full.tsv', dest/'genes.bed'],
                                motif=[dest/'meme-chip_output/meme-chip.html'],
                                candidates=[dest/'candidate_targets.tsv', dest/'selection.json', dest/'selected_genes.txt'])
                for p in required[step]:
                    if not p.is_file():
                        fail(f'{step} reported success but output is missing: {p}')
                if step == 'motif' and not list((dest/'meme-chip_output').glob('fimo_out_*/fimo.tsv')):
                    fail('MEME-ChIP produced no FIMO tables; candidates cannot be evaluated')
                return [dest], dict(directory=str(dest), peak_count=n)
            rec = self.stage(f'{aid}/{step}', step_sig, action)
            # A verified legacy environment alias retains the original dependency
            # identity, so unchanged motif outputs do not invalidate candidates.
            previous = self.annotation_dependency(aid, rec) if step == 'annotate' else digest([rec['signature'], rec['outputs']])
        count = len((root/'candidates/selected_genes.txt').read_text().splitlines())
        return dict(status='complete', peak_count=n, candidates=count, config=str(cfile))

    def summary(self, outcomes):
        out = self.results/'summary'
        table, genes, motifs = [], [], []
        for aid, a in self.plan['analyses'].items():
            r = outcomes.get(aid, dict(status='not_run', candidates='NA', peak_count='NA'))
            save(self.results/'analyses'/aid/'analysis_status.json', dict(updated=now(), **r))
            table.append([aid, a['target'], a['condition'], a['group'], len(a['biological_replicates']),
                          r['status'], r.get('peak_count', 'NA'), r.get('candidates', 'NA'), r.get('reason', '')])
            if r['status'] != 'complete':
                continue
            d = self.results/'analyses'/aid/'candidates'
            genes.extend([aid, a['target'], gene] for gene in (d/'selected_genes.txt').read_text().splitlines())
            with (d/'motif_summary.tsv').open() as f:
                for row in csv.DictReader(f, delimiter='\t'):
                    motifs.append([aid, aid + '::' + row['motif'], row['motif'], row['n_hits'],
                                   row['n_genes_target'], row['centrimo_Evalue'], row['centrimo_mode']])
        tsv(out/'batch_status.tsv', ['analysis_id', 'target', 'condition', 'nfcore_group',
            'biological_replicates', 'status', 'peak_count', 'candidate_genes', 'reason'], table)
        tsv(out/'candidate_genes.tsv', ['analysis_id', 'target', 'gene_id'], genes)
        tsv(out/'motifs.tsv', ['analysis_id', 'scoped_motif_id', 'motif_id', 'n_hits',
                              'n_genes_target', 'centrimo_Evalue', 'centrimo_mode'], motifs)
        save(out/'batch_summary.json', dict(version=VERSION, updated=now(),
             scope='Current invocation selection only; unrun/failed results are NA, not negative binding evidence',
             analyses=outcomes, errors=self.errors, plan=str(self.state/'plan.json')))
        print(f'[REPORT] {out / "batch_status.tsv"}', flush=True)

    def refresh_selection_version(self, aid):
        """Correct only the known 1.2.8 reporter typo, after full output checks."""
        record = self.records.get(aid + '/candidates')
        if not record or record.get('status') != 'complete':
            return
        path = self.results/'analyses'/aid/'candidates/selection.json'
        if not valid_record(record, record['signature']):
            fail(f'{aid}: candidate outputs changed; metadata was not edited')
        original = path.read_bytes()
        data = json.loads(original)
        # Reused results retain the version that produced them. A report redraw
        # must not relabel valid 1.2.8 candidate evidence as newly computed.
        if data.get('version') in (VERSION, '1.2.8'):
            return
        if (self.plan.get('version') != '1.2.8'
                or data.get('version') != '1.2.7' or data.get('status') != 'complete'):
            fail(f'{aid}: unexpected selection version; metadata was not edited')
        backup = self.state/f'report_metadata_backup/{self.run_id}/{aid}-selection.json'
        write(backup, original.decode('utf-8'))
        data['version'] = '1.2.8'
        prior_record = copy.deepcopy(record)
        try:
            save(path, data)
            record['outputs'][str(path)] = sha(path)
            record['metadata_correction'] = dict(field='selection.json:version', previous='1.2.7',
                current='1.2.8', backup=str(backup), reason='Known reporter typo in the recorded 1.2.8 run')
            save(self.status_path, self.records)
        except BaseException:
            path.write_bytes(original)
            self.records[aid + '/candidates'] = prior_record
            raise
        print(f'[METADATA] {aid}: selection version corrected to 1.2.8; scientific evidence and filters preserved', flush=True)

    def annotation_dependency(self, aid, record):
        """Preserve a prior dependency key only for a verified metadata-only edit."""
        correction = record.get('non_idr_metadata_correction')
        if not correction:
            return digest([record['signature'], record['outputs']])
        from steps.support.peak_metadata import POLICY, non_idr_replacement
        if correction.get('policy') != POLICY or record['outputs'] != correction['after_outputs']:
            fail(f'{aid}: annotation metadata correction receipt does not match outputs')
        annotation = self.results/'analyses'/aid/'annotation/peak_annotation_full.tsv'
        policy = annotation.with_name('annotation_policy.json')
        backup = Path(correction['backup'])
        raw, policy_raw = (backup/annotation.name).read_bytes(), (backup/policy.name).read_bytes()
        before = correction['before_outputs']
        if sha(backup/annotation.name) != before.get(str(annotation)) or sha(backup/policy.name) != before.get(str(policy)):
            fail(f'{aid}: annotation metadata backup changed')
        changes = {p for p in set(before) | set(record['outputs']) if before.get(p) != record['outputs'].get(p)}
        if changes != {str(annotation), str(policy)}:
            fail(f'{aid}: metadata receipt contains an unrelated output change')
        self.non_idr_input(aid)
        result = non_idr_replacement(raw, policy_raw, self.state/'analyses'/aid/'selected_peaks.bed')
        if result is None or annotation.read_bytes() != result[0] or policy.read_bytes() != result[1]:
            fail(f'{aid}: metadata correction differs from the permitted field-only transformation')
        return digest([record['signature'], before])

    def non_idr_input(self, aid):
        data = read_json(self.results/'analyses'/aid/'inputs.json', {})
        if data.get('idr') is not False or data.get('peakset_kind') not in ('single_replicate_no_idr', 'replicate_consensus'):
            fail(f'{aid}: explicit non-IDR provenance is required for metadata correction')

    def repair_annotation_metadata(self, aid):
        """Repair only the verified legacy narrowPeak V7 label, with a backup."""
        from steps.support.peak_metadata import POLICY, non_idr_replacement
        record = self.records[aid + '/annotate']
        if not valid_record(record, record['signature']):
            fail(f'{aid}: saved annotation outputs changed; metadata correction refused')
        if record.get('non_idr_metadata_correction'):
            self.annotation_dependency(aid, record)
            return
        annotation = self.results/'analyses'/aid/'annotation/peak_annotation_full.tsv'
        policy = annotation.with_name('annotation_policy.json')
        self.non_idr_input(aid)
        raw, policy_raw = annotation.read_bytes(), policy.read_bytes()
        result = non_idr_replacement(raw, policy_raw, self.state/'analyses'/aid/'selected_peaks.bed')
        if result is None:
            return
        # Selection manifests also pin annotation sources for tracks. Refresh
        # only these two source identities, after validating every original
        # input/output; candidate rows, genes and thresholds remain untouched.
        candidate = self.records.get(aid + '/candidates')
        selection = self.results/'analyses'/aid/'candidates/selection.json'
        selection_raw = selection_new = None
        prior_candidate = copy.deepcopy(candidate)
        if candidate and candidate.get('status') == 'complete':
            if not valid_record(candidate, candidate['signature']):
                fail(f'{aid}: candidate outputs changed; annotation metadata correction refused')
            selection_raw = selection.read_bytes()
            data = json.loads(selection_raw)
            corrected_sources = {str(annotation): hashlib.sha256(result[0]).hexdigest(),
                                 str(policy): hashlib.sha256(result[1]).hexdigest()}
            found = []
            for item in data['sources']:
                if sha(item['path']) != item['sha256']:
                    fail(f'{aid}: selection source changed; metadata correction refused: ' + item['path'])
                if item['path'] in corrected_sources:
                    found.append(item['path'])
                    item['sha256'] = corrected_sources[item['path']]
            if sorted(found) != sorted(corrected_sources):
                fail(f'{aid}: selection does not unambiguously record both annotation sources')
            for item in data['files'].values():
                target = (selection.parent/item['file']).resolve()
                if target.parent != selection.parent.resolve() or sha(target) != item['sha256']:
                    fail(f'{aid}: selection output changed; metadata correction refused')
            selection_new = (json.dumps(data, indent=2) + '\n').encode('utf-8')
        backup = self.state/f'annotation_metadata_backup/{self.run_id}/{aid}'
        write(backup/annotation.name, raw.decode('utf-8'))
        write(backup/policy.name, policy_raw.decode('utf-8'))
        prior = copy.deepcopy(record)
        save(backup/'record.json', prior)
        if selection_raw is not None:
            write(backup/'selection.json', selection_raw.decode('utf-8'))
            save(backup/'candidates_record.json', prior_candidate)
        try:
            write(annotation, result[0].decode('utf-8'))
            write(policy, result[1].decode('utf-8'))
            record['outputs'][str(annotation)] = sha(annotation)
            record['outputs'][str(policy)] = sha(policy)
            record['non_idr_metadata_correction'] = dict(policy=POLICY, backup=str(backup),
                before_outputs=prior['outputs'], after_outputs=copy.deepcopy(record['outputs']),
                rows=result[2], field='peak_globalIDR', value='NA', corrected=now())
            if selection_new is not None:
                write(selection, selection_new.decode('utf-8'))
                candidate['outputs'][str(selection)] = sha(selection)
                candidate['non_idr_metadata_sources'] = dict(policy=POLICY, backup=str(backup/'selection.json'),
                    scope='Only the two verified annotation source identities were refreshed; selection evidence unchanged')
            self.annotation_dependency(aid, record)
            save(self.status_path, self.records)
        except BaseException:
            annotation.write_bytes(raw); policy.write_bytes(policy_raw)
            self.records[aid + '/annotate'] = prior
            if selection_raw is not None:
                selection.write_bytes(selection_raw)
                self.records[aid + '/candidates'] = prior_candidate
            raise
        print(f'[METADATA] {aid}: {result[2]} non-IDR rows corrected to globalIDR=NA; backup={backup}', flush=True)

    def reports(self, only=None):
        """Render separate reports from a saved run without changing its results."""
        selected = list(self.plan['analyses'])
        if only:
            wanted = only.split(',')
            if len(set(wanted)) != len(wanted) or set(wanted) - set(selected):
                fail('Unknown or duplicate --only analysis IDs')
            selected = wanted
        reference = self.records.get('reference')
        if not reference or not valid_record(reference, reference['signature']):
            fail('Saved reference outputs are missing or changed; report refresh refused')
        for item in self.plan['sources'].values():
            if not Path(item['path']).is_file() or sha(item['path']) != item['sha256']:
                fail('Saved reference/annotation/database source changed: ' + item['path'])
        envs = preflight(self.cfg, ('annotate',))
        report_records_path = self.state/'reports.json'
        report_records = read_json(report_records_path, {})
        errors = []
        for aid in selected:
            record = self.records.get(aid + '/annotate')
            if not record or record.get('status') != 'complete':
                print(f'[SKIP] {aid}: no completed annotation; no report was invented', flush=True)
                continue
            try:
                if not valid_record(record, record['signature']):
                    fail('Saved annotation outputs are missing or changed')
                self.repair_annotation_metadata(aid)
                analysis = self.results/'analyses'/aid
                rt = self.state/'analyses'/aid
                peaks = rt/'selected_peaks.bed'
                annotation = analysis/'annotation/peak_annotation_full.tsv'
                gff = reference['details']['GENOME_GFF3']
                signature = digest([REPORT_POLICY, record['outputs'], reference['outputs'], sha(peaks),
                    envs, self.cfg, {name: sha(HOME/'lib/steps'/name) for name in
                    ('05_peak_annotation.R', 'annotation_plots.R', 'plot_tss.R', 'tss_windows.R')}])
                dest = analysis/'reports/annotation'
                if valid_record(report_records.get(aid), signature):
                    print(f'[RESUME] {aid}/reports: verified outputs', flush=True)
                else:
                    stage_root = rt/'report-staging'
                    stage_root.mkdir(parents=True, exist_ok=True)
                    with tempfile.TemporaryDirectory(prefix='annotation-', dir=stage_root) as temporary:
                        staged = Path(temporary)
                        cmd = env_prefix(self.cfg, 'annotate') + ['Rscript', str(HOME/'lib/steps/05_peak_annotation.R'),
                            gff, str(peaks), str(staged), self.cfg['PROMOTER_UP'], self.cfg['PROMOTER_DOWN'],
                            self.cfg.get('ANNOTATION_INFO', ''), self.cfg['TARGET_DIST_KB'], 'dapseq', '',
                            self.cfg['MAKE_PLOTS'], str(annotation), reference['details']['CHROM_SIZES'], 'non_idr']
                        self.execute(cmd, self.logs/f'analyses/{aid}/{self.run_id}_reports.log', rt)
                        if not (staged/'plot_status.tsv').is_file():
                            fail('Report renderer returned without plot_status.tsv')
                        save(staged/'report.json', dict(version=VERSION, policy=REPORT_POLICY, created=now(),
                            annotation=str(annotation), annotation_sha256=sha(annotation),
                            chrom_sizes=reference['details']['CHROM_SIZES'], peak_kind='non_idr',
                            scope='Diagnostic plots; checked non-IDR metadata correction; peak/gene assignments and selection filters preserved'))
                        # Publish only after reconstruction and all requested plots succeed.
                        self.archive(dest)
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        staged.rename(dest)
                    report_records[aid] = dict(status='complete', signature=signature,
                                              outputs=artifact_files([dest]), finished=now())
                    save(report_records_path, report_records)
                    print(f'[REPORT] {dest}', flush=True)
                self.refresh_selection_version(aid)
            except (BatchError, OSError, ValueError, KeyError) as exc:
                errors.append(f'{aid}: {exc}')
                print('[ERROR] ' + errors[-1], flush=True)
        print(f'[REPORTS] {len(selected)} analyses considered; {len(errors)} errors; batch_status.tsv preserved', flush=True)
        return 1 if errors else 0

    def run(self):
        outcomes = {}
        try:
            self.envs = preflight(self.cfg, STEPS[:STEPS.index(self.until)+1])
            save(self.state/f'environments/{self.run_id}.json', self.envs)
            save(self.state/'plan.json', self.plan)
            save(self.state/f'plans/{self.run_id}.json', self.plan)
            ref, refsig = self.shared_reference()
            if self.until == 'reference':
                outcomes = {a: dict(status='stopped_after_reference', candidates='NA') for a in self.plan['analyses']}
                return 0
            index = None
            for gid, group in self.plan['groups'].items():
                try:
                    details, sig = self.group_run(group, ref, refsig, index)
                    if details['bwa_index']:
                        index = details['bwa_index']
                except (BatchError, OSError, ValueError) as e:
                    self.errors.append(f'{gid}: {e}')
                    print('[ERROR] ' + self.errors[-1], flush=True)
                    outcomes.update({a: dict(status='nfcore_failed', candidates='NA', reason=str(e)) for a in group['analyses']})
                    continue
                for aid in group['analyses']:
                    if aid in details.get('failures', {}):
                        outcomes[aid] = details['failures'][aid]
                        self.errors.append(f'{aid}: {outcomes[aid]["reason"]}')
                        print('[ERROR] ' + self.errors[-1], flush=True)
                        self.summary(outcomes)
                        continue
                    if self.until == 'nfcore':
                        outcomes[aid] = dict(status='stopped_after_nfcore', candidates='NA')
                        continue
                    try:
                        outcomes[aid] = self.downstream(aid, details['mapping'][aid], ref, sig)
                    except (BatchError, OSError, ValueError) as e:
                        self.errors.append(f'{aid}: {e}')
                        print('[ERROR] ' + self.errors[-1], flush=True)
                        outcomes[aid] = dict(status='downstream_failed', candidates='NA', reason=str(e))
                    self.summary(outcomes)
        except (BatchError, OSError, ValueError) as e:
            self.errors.append(str(e))
            outcomes.update({a: dict(status='shared_step_failed', candidates='NA', reason=str(e))
                             for a in self.plan['analyses'] if a not in outcomes})
            raise
        finally:
            self.summary(outcomes)
            counts = collections.Counter(r['status'] for r in outcomes.values())
            print('[SUMMARY] ' + '; '.join(f'{status}={count}' for status, count in sorted(counts.items())), flush=True)
        return 1 if self.errors else 0


def init_project(args):
    project = Path(args.directory).resolve()
    for name in ('chromaflow.config', 'samples.tsv'):
        if (project/name).exists():
            fail(f'Refusing to overwrite {project/name}')
    project.mkdir(parents=True, exist_ok=True)
    for name in ('raw_data', 'reference', 'logs', 'results'):
        (project/name).mkdir(exist_ok=True)
    config = (HOME/'templates/batch.config').read_text()
    config = config.replace('PROJECT_DIR="/data/projects/MyDAPbatch"', 'PROJECT_DIR=' + shlex.quote(str(project)))
    write(project/'chromaflow.config', config)
    write(project/'samples.tsv', (HOME/'templates/samples_batch.tsv').read_text())
    print(f'[INIT] {project}\nEdit chromaflow.config and samples.tsv, then run chromaflow batch plan.')


class Tee:
    def __init__(self, terminal, logfile):
        self.terminal, self.logfile = terminal, logfile

    def write(self, text):
        self.terminal.write(text)
        self.logfile.write(text)
        self.logfile.flush()
        return len(text)

    def flush(self):
        self.terminal.flush()
        self.logfile.flush()


def main(argv=None):
    parser = argparse.ArgumentParser(description='ChromaFlow 1.2.8 batch DAP-seq: shared nf-core, isolated downstream through candidates.',
        epilog='One project config + one TSV. Read docs/BATCH.md before preparing samples. No lookup/tracks/IDR are run automatically.')
    commands = parser.add_subparsers(dest='action', required=True)
    init = commands.add_parser('init', help='Create shared config and example sample table; never overwrite')
    init.add_argument('directory')
    for name in ('plan', 'run', 'status', 'reports'):
        p = commands.add_parser(name)
        p.add_argument('-c', '--config', default=os.environ.get('CHROMAFLOW_CONFIG', 'chromaflow.config'))
        if name == 'status':
            continue
        if name == 'reports':
            p.add_argument('--only', help='Comma-separated exact analysis IDs from the saved run')
            continue
        p.add_argument('--samples', default='samples.tsv')
        p.add_argument('--only', help='Comma-separated exact analysis IDs; required controls are included')
        if name == 'plan':
            p.add_argument('--check-env', action='store_true', help='Also check each configured execution environment')
        else:
            p.add_argument('--resume', action='store_true', help='Reuse stages only after input signatures and output hashes agree')
            p.add_argument('--until', choices=STEPS, default='candidates', help='Stop after a named stage (default: candidates)')
    args = parser.parse_args(argv)
    if args.action == 'init':
        init_project(args)
        return 0
    if args.action == 'status':
        cfg = load_config(args.config)
        path = Path(cfg['OUTPUT_DIR'])/'results/summary/batch_status.tsv'
        if not path.is_file():
            fail('No batch run report exists yet; run batch plan or batch run first')
        print(path.read_text(), end='')
        print(f'[STATE] {Path(cfg["OUTPUT_DIR"])/".chromaflow/batch/status.json"}')
        return 0
    cfg = load_config(args.config)
    logfile = Path(cfg['OUTPUT_DIR'])/'logs/batch'/(dt.datetime.now().strftime('%Y%m%d_%H%M%S') + f'_{os.getpid()}_{args.action}.log')
    logfile.parent.mkdir(parents=True, exist_ok=True)
    terminal_out, terminal_err = sys.stdout, sys.stderr
    with logfile.open('w', encoding='utf-8') as stream:
        sys.stdout, sys.stderr = Tee(terminal_out, stream), Tee(terminal_err, stream)
        try:
            print(f'[START] {now()}\n[VERSION] {VERSION}\n[LOG] {logfile}\n[COMMAND] {shlex.join([str(CLI), "batch", *(argv or sys.argv[1:])])}', flush=True)
            rc = dispatch(args)
            print(f'[END] {now()} exit_code={rc}', flush=True)
            return rc
        except BaseException as e:
            print(f'[ERROR] {e}\n[END] {now()} exit_code={130 if isinstance(e, KeyboardInterrupt) else 1}', file=sys.stderr, flush=True)
            raise
        finally:
            sys.stdout, sys.stderr = terminal_out, terminal_err


def dispatch(args):
    if args.action == 'reports':
        cfg = load_config(args.config)
        state = Path(cfg['OUTPUT_DIR'])/'.chromaflow/batch'
        plan = read_json(state/'plan.json')
        if not plan or cfg != plan.get('config'):
            fail('Reports require the unchanged config from the saved batch run')
        import fcntl
        with (state/'run.lock').open('a+') as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                fail('Another batch run is active for this OUTPUT_DIR')
            return Runner(plan, True, 'candidates').reports(args.only)
    plan = plan_project(args)
    root = Path(plan['config']['OUTPUT_DIR'])
    state = root/'.chromaflow/batch'
    print(f'[PLAN] {len(plan["analyses"])} analyses; {len(plan["groups"])} compatible nf-core groups; serial execution')
    for gid, g in plan['groups'].items():
        ncontrol = len({(r['sample_id'], r['replicate']) for r in g['rows'] if r['role'] == 'control'})
        print(f'[GROUP] {gid}: {len(g["analyses"])} targets; {ncontrol} control libraries; k={g["k"]}')
    if args.action == 'plan':
        save(state/'proposed_plan.json', plan)
        for gid, g in plan['groups'].items():
            nf_sheet(state/f'proposed/{gid}.csv', g['rows'])
        if args.check_env:
            save(state/'environment_check.json', preflight(plan['config']))
        print(f'[PLAN] {state/"proposed_plan.json"}\n[OK] No analysis was started. FASTQ checks sampled up to 2000 records per file.')
        return 0
    state.mkdir(parents=True, exist_ok=True)
    # fcntl automatically releases the lock after interruption. Do not delete a
    # persistent lock inode: that could allow two concurrent project writers.
    import fcntl
    with (state/'run.lock').open('a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            fail('Another batch run is active for this OUTPUT_DIR')
        return Runner(plan, args.resume, args.until).run()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (BatchError, OSError, ValueError, KeyError) as e:
        print('[ERROR] ' + str(e), file=sys.stderr, flush=True)
        sys.exit(1)
    except KeyboardInterrupt:
        print('[INTERRUPTED] Resume using chromaflow batch run --resume', file=sys.stderr)
        sys.exit(130)
