#!/usr/bin/env python3
"""Release regression, upgrade/resume and report acceptance; stdlib only.

No runtime patch, result cleanup, parameter changes or new content hashing.
ChromaFlow itself still performs its normal cache validation.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter
from contextlib import contextmanager
import csv
import datetime as dt
import io
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tarfile
import time

DEFAULT_PROJECT = '.'
DEFAULT_CF = str(Path(__file__).resolve().parents[1]/'bin/chromaflow')
KIT_VERSION = 'release-1.2.9-acceptance-v1'
RELEASE = '1.2.9'
PLAN_FIELDS = ('config', 'config_path', 'samples_path', 'samples_sha256',
               'all_analysis_ids', 'analyses', 'groups', 'sources', 'fastqs', 'step_code_sha256')
SCIENCE_FIELDS = ('status', 'peak_count', 'candidates')
MODEL_FIELDS = ('status', 'task', 'task_exit', 'minimum_model_pairs', 'observed_model_pairs')
PLOTS = ('anno_pie.png', 'dist_to_tss.png', 'anno_upset.png', 'tss_avgprof.png')
NATIVE_CORE = ('TRIMGALORE', 'BWA_MEM', 'PICARD_MARKDUPLICATES')
TABLES = {'dist_to_tss.png': 'dist_to_tss.tsv', 'anno_upset.png': 'annotation_intersections.tsv',
          'tss_avgprof.png': 'tss_avgprof.tsv'}


def now():
    return dt.datetime.now().astimezone().isoformat(timespec='seconds')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')


def tsv(path, columns, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def inside(path, root):
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except ValueError:
        return False


def file_stat(path):
    try:
        st = Path(path).stat()
        return dict(exists=Path(path).is_file(), size=st.st_size, mtime_ns=st.st_mtime_ns)
    except FileNotFoundError:
        return dict(exists=False)


@contextmanager
def lock(path):
    import fcntl
    with Path(path).open('a+') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another batch/test is active. Wait for it to finish: ' + str(path))
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def snapshot(project, dest):
    state = project/'.chromaflow/batch'
    dest.mkdir(parents=True, exist_ok=True)
    for src, name in ((state/'status.json', 'status.json'), (state/'plan.json', 'plan.json'),
                      (state/'reports.json', 'reports.json'),
                      (project/'results/summary/batch_summary.json', 'batch_summary.json'),
                      (project/'results/summary/batch_status.tsv', 'batch_status.tsv')):
        if src.is_file():
            shutil.copy2(src, dest/name)
    records = read(dest/'status.json')
    # Stat existing downstream artifacts, not FASTQs/BAMs/reference files.
    metadata = {key: {p: file_stat(p) for p in rec.get('outputs', {})}
                for key, rec in records.items()
                if key.rsplit('/', 1)[-1] in ('annotate', 'motif', 'candidates')}
    save(dest/'output_stats.json', metadata)
    audits, traces = {}, {}
    for gid in read(dest/'plan.json')['groups']:
        path = state/'runs'/gid/'peak_tasks.json'
        if path.is_file():
            audits[gid] = read(path)
        trace = records.get('nfcore/' + gid, {}).get('details', {}).get('trace') or audits.get(gid, {}).get('trace')
        if trace and inside(trace, project/'logs') and Path(trace).is_file():
            with Path(trace).open(encoding='utf-8', newline='') as stream:
                reader = csv.DictReader(stream, delimiter='\t')
                if {'process', 'tag', 'status', 'exit'} <= set(reader.fieldnames or []):
                    traces[gid] = dict(path=trace, rows=list(reader))
    save(dest/'peak_tasks.json', audits)
    save(dest/'native_traces.json', traces)


def baseline_check(plan, summary, records):
    selected = set(plan['analyses'])
    if selected != set(plan['all_analysis_ids']) or selected != set(summary['analyses']):
        raise ValueError('Saved run is not the complete sample selection. Do not use --only or a smoke TSV for this test.')
    if not any(r['status'] == 'complete' for r in summary['analyses'].values()):
        raise ValueError('No completed analyses are available for resume acceptance.')
    for aid, outcome in summary['analyses'].items():
        status = outcome['status']
        if status not in ('complete', 'insufficient_peaks', 'no_peaks', 'model_failed'):
            raise ValueError('Resolve/review this unexpected baseline status first: ' + aid + ': ' + status)
        stages = ('annotate', 'motif', 'candidates') if status == 'complete' else ('annotate',) if status == 'insufficient_peaks' else ()
        for step in stages:
            rec = records.get(aid + '/' + step, {})
            if rec.get('status') != 'complete' or not rec.get('outputs'):
                raise ValueError('Missing complete stage record: ' + aid + '/' + step)
        if status == 'model_failed' and not outcome.get('tasks'):
            raise ValueError('Missing model-failure evidence: ' + aid)


def stream_command(command, cwd, logfile):
    print('[COMMAND] ' + shlex.join(map(str, command)), flush=True)
    with logfile.open('w', encoding='utf-8') as out:
        out.write('[COMMAND] ' + shlex.join(map(str, command)) + '\n')
        out.flush()
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(('CHROMAFLOW_', 'CF_BATCH_', 'CFSEL_', 'NXF_', 'NFCORE_'))}
        env.update(PYTHONUNBUFFERED='1', PYTHONUTF8='1', LC_ALL='C.UTF-8')
        with subprocess.Popen(list(map(str, command)), cwd=cwd, env=env, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, errors='replace') as proc:
            try:
                for line in proc.stdout:
                    print(line, end='', flush=True)
                    out.write(line)
                    out.flush()
                code = proc.wait()
            except KeyboardInterrupt:
                # The foreground process group also receives Ctrl-C. Wait for
                # ChromaFlow to stop its child tasks and flush its own status.
                print('[INTERRUPTED] Waiting for ChromaFlow to stop.', flush=True)
                proc.wait()
                raise
        out.write(f'[CAPTURE_EXIT] {code}\n')
    return code


def stage_events(text):
    events = {}
    for kind, key, reason in re.findall(r'^\[(RESUME|RUN)\] ([^:\s]+): (.*)$', text, re.M):
        events.setdefault(key, []).append(dict(event=kind, reason=reason))
    return events


def audit_resume(before, after, log_text, batch_exit):
    old, new = read(before/'status.json'), read(after/'status.json')
    bplan, aplan = read(before/'plan.json'), read(after/'plan.json')
    bsum, asum = read(before/'batch_summary.json')['analyses'], read(after/'batch_summary.json')['analyses']
    bstats, astats = read(before/'output_stats.json'), read(after/'output_stats.json')
    events, issues, rows = stage_events(log_text), [], []
    changed_plan = [key for key in PLAN_FIELDS if bplan.get(key) != aplan.get(key)]
    if changed_plan:
        issues.append('Run plan changed: ' + ', '.join(changed_plan))
    for key, rec in sorted(old.items()):
        # Partial native groups are deliberately retried. Complete scientific
        # stages and complete native groups must reuse their records.
        if rec.get('status') != 'complete':
            continue
        if '/' in key and not key.startswith(('gsize/', 'nfcore/')) and key.split('/')[0] not in bplan['analyses']:
            continue
        current = new.get(key, {})
        ev = events.get(key, [])
        problems = []
        if not any(e['event'] == 'RESUME' for e in ev):
            problems.append('no current RESUME event')
        if any(e['event'] == 'RUN' for e in ev):
            problems.append('stage executed again')
        for field in ('status', 'signature', 'finished', 'outputs'):
            if rec.get(field) != current.get(field):
                problems.append('record changed: ' + field)
        if key in bstats:
            if bstats[key] != astats.get(key):
                problems.append('output size/mtime/file presence changed')
            if any(not v['exists'] for v in bstats[key].values()):
                problems.append('baseline output missing')
        verdict = 'REVIEW' if problems else 'PASS'
        rows.append(dict(stage=key, result=verdict, reason='; '.join(problems),
                         events=' | '.join(e['event'] + ': ' + e['reason'] for e in ev)))
        issues.extend(key + ': ' + p for p in problems)
    for aid in sorted(set(bsum) | set(asum)):
        b, a = bsum.get(aid, {}), asum.get(aid, {})
        if any(b.get(k) != a.get(k) for k in SCIENCE_FIELDS):
            issues.append('Analysis outcome/count changed: ' + aid)
        if b.get('status') == 'model_failed':
            def failures(value):
                return sorted(tuple(str(t.get(k, '')) for k in MODEL_FIELDS) for t in value.get('tasks', []))
            if failures(b) != failures(a):
                issues.append('Model-failure evidence changed: ' + aid)
    expected_rc = 1 if any(r['status'] == 'model_failed' for r in bsum.values()) else 0
    if batch_exit != expected_rc:
        issues.append(f'Batch exit changed: expected {expected_rc}, observed {batch_exit}')
    # Native cache checks use this invocation's task trace, never an old audit.
    prior_audits, new_audits = read(before/'peak_tasks.json'), read(after/'peak_tasks.json')
    native, core, process_counts = [], [], []
    prior_traces, new_traces = read(before/'native_traces.json'), read(after/'native_traces.json')
    for gid in bplan['groups']:
        key = 'nfcore/' + gid
        if old.get(key, {}).get('status') == 'complete':
            native.append(dict(group=gid, task='(whole group)', result='RESUMED_GROUP', status='', reason='Checked above'))
            continue
        audit = new_audits.get(gid, {})
        if (not any(e['event'] == 'RUN' for e in events.get(key, [])) or not audit.get('tasks')
                or not audit.get('trace') or audit.get('trace') == prior_audits.get(gid, {}).get('trace')):
            issues.append('Missing fresh native task evidence: ' + gid)
            continue
        # Require every previously successful peak call to be CACHED, including
        # zero-peak successes; do not count the known failed task as a cache miss.
        old_tasks = prior_audits.get(gid, {}).get('tasks', {})
        for tag, task in old_tasks.items():
            current = audit['tasks'].get(tag, {})
            expected = 'CACHED' if task['status'] in ('CACHED', 'COMPLETED') else task['status']
            good = current.get('status') == expected and str(current.get('exit')) == str(task.get('exit'))
            native.append(dict(group=gid, task=tag, result='PASS' if good else 'REVIEW',
                               status=current.get('status', 'MISSING'), reason='Expected ' + expected))
            if not good:
                issues.append('Native peak task was not reused/changed: ' + tag)
        if not old_tasks:
            issues.append('Missing baseline native task evidence: ' + gid)
        prior, current = prior_traces.get(gid, {}), new_traces.get(gid, {})
        if not prior.get('rows') or not current.get('rows') or current['path'] != audit['trace']:
            issues.append('Missing complete native trace for cache review: ' + gid)
            continue
        counts = Counter((row['process'], row['status']) for row in current['rows'])
        process_counts.extend(dict(group=gid, process=proc, status=status, count=count)
                              for (proc, status), count in sorted(counts.items()))
        for previous in prior['rows']:
            if previous['process'].split(':')[-1] not in NATIVE_CORE or previous['status'] not in ('CACHED', 'COMPLETED'):
                continue
            matches = [row for row in current['rows']
                       if (row['process'], row['tag']) == (previous['process'], previous['tag'])]
            good = len(matches) == 1 and matches[0]['status'] == 'CACHED' and matches[0]['exit'] == '0'
            core.append(dict(group=gid, process=previous['process'], tag=previous['tag'],
                             result='PASS' if good else 'REVIEW', status=','.join(row['status'] for row in matches) or 'MISSING'))
            if not good:
                issues.append('Native trimming/alignment/dedup task was not reused: ' + previous['process'] + ' (' + previous['tag'] + ')')
        if not any(row['group'] == gid for row in core):
            issues.append('No baseline trimming/alignment/dedup tasks available: ' + gid)
    return dict(result='PASS' if not issues else 'REVIEW', batch_exit=batch_exit, expected_batch_exit=expected_rc,
                issues=issues, stages=rows, native_peak_tasks=native, native_core_tasks=core, native_process_counts=process_counts,
                before_counts=dict(Counter(r['status'] for r in bsum.values())),
                after_counts=dict(Counter(r['status'] for r in asum.values())),
                scope='Resume evidence only; does not turn existing biological/model failures into successes')


def png_dimensions(path):
    with path.open('rb') as stream:
        head = stream.read(24)
        stream.seek(-12, os.SEEK_END)
        tail = stream.read(12)
    if head[:8] != b'\x89PNG\r\n\x1a\n' or head[12:16] != b'IHDR' or tail != b'\0\0\0\0IEND\xaeB`\x82':
        raise ValueError('Invalid PNG header or end marker')
    width, height = struct.unpack('>II', head[16:24])
    if not width or not height:
        raise ValueError('Empty image dimensions')
    return width, height


def audit_tss_metadata(project, aid, folder, records):
    """Independent interval-count and denominator check; no file rewriting."""
    annotation = project/'results/analyses'/aid/'annotation/peak_annotation_full.tsv'
    with annotation.open(encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if any(row.get('peak_globalIDR') != 'NA' for row in rows):
        raise ValueError('Non-IDR annotation still has a populated/missing globalIDR field')
    info = read(folder/'report.json')
    if info.get('policy') != 'annotation-reports-v2-reference-bounds' or info.get('peak_kind') != 'non_idr':
        raise ValueError('Report has not been refreshed with reference-bound TSS policy')
    sizes_path = Path(records['reference']['details']['CHROM_SIZES'])
    sizes = dict(line.split('\t') for line in sizes_path.read_text().splitlines() if line.strip())
    sizes = {name: int(length) for name, length in sizes.items()}
    peak_starts, peak_ends = {}, {}
    for line in (project/'.chromaflow/batch/analyses'/aid/'selected_peaks.bed').read_text().splitlines():
        if not line.strip() or line.startswith(('#', 'track ', 'browser ')): continue
        c = line.split('\t')
        peak_starts.setdefault(c[0], []).append(int(c[1])+1)
        peak_ends.setdefault(c[0], []).append(int(c[2]))
    for group in (peak_starts, peak_ends):
        for values in group.values(): values.sort()
    with (folder/'tss_windows.tsv').open(encoding='utf-8') as stream:
        windows = list(csv.DictReader(stream, delimiter='\t'))
    seen, counts, recovered = set(), Counter(), 0
    for row in windows:
        ch, start, end = row['seqnames'], int(row['start']), int(row['end'])
        key = (ch, start, end, row['strand'])
        if key in seen or end-start != 6000 or row['strand'] not in ('+', '-'):
            raise ValueError('Duplicate or invalid TSS window')
        seen.add(key)
        if int(row['chrom_length']) != sizes[ch] or int(row['tss']) != start+3000:
            raise ValueError('TSS coordinates/length differ from reference bounds')
        overlaps = bisect_right(peak_starts.get(ch, []), end) - bisect_left(peak_ends.get(ch, []), start)
        expected = 'outside_reference' if start < 1 or end > sizes[ch] else ('retained' if overlaps else 'no_peak_overlap')
        if row['status'] != expected or int(row['overlapping_peaks']) != overlaps:
            raise ValueError('TSS retention or overlap count differs from selected input peaks')
        legacy = expected == 'retained' and end > peak_ends.get(ch, [0])[-1]
        if row['legacy_end_bound_exclusion'] != ('TRUE' if legacy else 'FALSE'):
            raise ValueError('Incorrect legacy boundary exclusion flag')
        recovered += legacy
        counts[expected] += 1
    with (folder/'tss_window_summary.tsv').open(encoding='utf-8') as stream:
        summary = {r['metric']: int(r['windows']) for r in csv.DictReader(stream, delimiter='\t')}
    expected_summary = {s: counts[s] for s in ('retained', 'no_peak_overlap', 'outside_reference')}
    expected_summary['legacy_end_bound_exclusion'] = recovered
    if summary != expected_summary: raise ValueError('TSS window summary differs from the audit')
    profile = folder/'tss_avgprof.tsv'
    if counts['retained']:
        with profile.open(encoding='utf-8') as stream:
            curve = list(csv.DictReader(stream, delimiter='\t'))
        if ([int(r['distance_bp']) for r in curve] != list(range(-3000, 3001)) or
            any(int(r['windows']) != counts['retained'] or not math.isfinite(float(r['mean_peak_overlap']))
                or float(r['mean_peak_overlap']) < 0 for r in curve)):
            raise ValueError('TSS numerical profile has an invalid denominator, grid or value')
        if not any(float(r['mean_peak_overlap']) > 0 for r in curve):
            raise ValueError('Retained TSS windows unexpectedly produce a zero curve')
    elif profile.exists():
        raise ValueError('A zero-window report retained a stale profile')


def audit_reports(project, plan, records):
    rows, issues = [], []
    plots_enabled = str(plan['config'].get('MAKE_PLOTS', 'true')).lower() == 'true'
    for aid in plan['analyses']:
        rec = records.get(aid + '/annotate', {})
        if rec.get('status') != 'complete':
            rows.append(dict(analysis_id=aid, plot='', status='not_applicable', reason='No complete annotation'))
            continue
        folder = project/'results/analyses'/aid/'reports/annotation'
        try:
            info = read(folder/'report.json')
            annotation = str(project/'results/analyses'/aid/'annotation/peak_annotation_full.tsv')
            if info.get('annotation') != annotation or not info.get('annotation_sha256') or info['annotation_sha256'] != rec['outputs'].get(annotation):
                raise ValueError('Report provenance differs from the recorded annotation')
            with (folder/'plot_status.tsv').open(encoding='utf-8', newline='') as stream:
                table = csv.DictReader(stream, delimiter='\t')
                if not {'plot', 'status', 'reason'} <= set(table.fieldnames or []):
                    raise ValueError('Malformed plot_status.tsv')
                entries = list(table)
            names = [entry['plot'] for entry in entries]
            if len(names) != len(set(names)) or (plots_enabled and set(names) != set(PLOTS)):
                raise ValueError('Missing, duplicate or unexpected annotation plot entries')
            if not plots_enabled:
                issues.append(aid + ': MAKE_PLOTS=false; visual acceptance is not available')
            else:
                audit_tss_metadata(project, aid, folder, records)
            for entry in entries:
                name, status = entry['plot'], entry['status']
                note = dict(analysis_id=aid, plot=name, status=status, reason=entry['reason'], width='', height='')
                if name not in PLOTS or status not in ('written', 'skipped'):
                    issues.append(aid + ': unexpected/failed plot: ' + name)
                elif status == 'written':
                    note['width'], note['height'] = png_dimensions(folder/name)
                    companion = TABLES.get(name)
                    if companion and (not (folder/companion).is_file() or (folder/companion).stat().st_size == 0):
                        issues.append(aid + ': missing numerical table: ' + companion)
                elif not entry['reason'] or (folder/name).exists():
                    issues.append(aid + ': unexplained skip or stale image: ' + name)
                rows.append(note)
        except (OSError, ValueError, KeyError, struct.error) as exc:
            issues.append(aid + ': ' + str(exc))
            rows.append(dict(analysis_id=aid, plot='', status='review', reason=str(exc)))
    return dict(result='FILES_OK_VISUAL_REVIEW_PENDING' if not issues else 'REVIEW', issues=issues, plots=rows,
                counts=dict(Counter(row['status'] for row in rows)),
                scope='File/status/provenance checks only; labels, scales and scientific interpretation need visual review')


def validate_release(cf):
    root = cf.parent.parent
    sys.path.insert(0, str(root/'lib'))
    import batch
    data = read(root/'release/compatibility.json')
    code = {p.relative_to(root).as_posix(): batch.sha(p) for folder in ('lib', 'bin')
            for p in sorted((root/folder).rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}
    if batch.VERSION != RELEASE or data.get('version') != RELEASE or data.get('runtime') != code:
        raise ValueError('Release installation differs from its recorded runtime; no analysis started')
    if batch.resume_code_ids(code) != data['step_code_sha256']:
        raise ValueError('Release stage compatibility check failed')
    return batch


def regression_checks(cf, plan, directory):
    batch = validate_release(cf)
    root, cfg = cf.parent.parent, plan['config']
    checks = [
        ('python-regression', batch.env_prefix(cfg, 'reference') +
         ['python3', str(root/'tests/run_suite.py'), '--json', str(directory/'python-regression.json')]),
        ('r-tss-windows', batch.env_prefix(cfg, 'annotate') +
         ['Rscript', str(root/'tests/test_tss_windows.R'), str(root)]),
        ('r-annotation-plots', batch.env_prefix(cfg, 'annotate') +
         ['Rscript', str(root/'tests/test_annotation_plots.R'), str(root), str(directory/'r-plots')]),
        ('r-report-guard', batch.env_prefix(cfg, 'annotate') +
         ['Rscript', str(root/'tests/test_report_annotation_guard.R'), str(root), 'python3', str(directory/'r-guard')])]
    results = []
    for name, cmd in checks:
        rc = stream_command(cmd, directory, directory/(name+'.log'))
        results.append(dict(test=name, exit=rc))
    save(directory/'regressions.json', results)
    if any(r['exit'] for r in results):
        raise ValueError('Release regressions failed or skipped; analysis resume was not started')


def refresh_reports(cf, project, config, directory):
    command = [cf, 'batch', 'reports', '-c', config]
    first = stream_command(command, project, directory/'reports-first.console.log')
    if first:
        raise ValueError('Report refresh failed; inspect reports-first.console.log')
    state = project/'.chromaflow/batch/reports.json'
    baseline = read(state)
    stats = {p: file_stat(p) for rec in baseline.values() for p in rec.get('outputs', {})}
    second = stream_command(command, project, directory/'reports-repeat.console.log')
    current = {p: file_stat(p) for rec in read(state).values() for p in rec.get('outputs', {})}
    if second or read(state) != baseline or current != stats:
        raise ValueError('Repeated report refresh did not preserve verified report files')
    result = dict(first_exit=first, repeat_exit=second, status='PASS',
                  report_records_and_file_stats_preserved=True)
    save(directory/'report-idempotence.json', result)
    return result


def run(args):
    project, cf = Path(args.project).resolve(), Path(args.cf).resolve()
    state = project/'.chromaflow/batch'
    if not (state/'plan.json').is_file():
        raise ValueError('No saved batch plan: ' + str(state))
    with lock(state/'resume-acceptance.lock'):
        stamp = dt.datetime.now().strftime('%Y%m%d-%H%M%S') + '-' + str(os.getpid())
        directory = project/'diagnostics'/('release-acceptance-' + stamp)
        directory.mkdir(parents=True, exist_ok=False)
        meta = dict(kit=KIT_VERSION, project=str(project), cf=str(cf), directory=str(directory),
                    started=now(), state='preflight', started_ns=time.time_ns())
        save(directory/'run.json', meta)
        save(project/'diagnostics/latest-release-acceptance.json', dict(directory=str(directory)))
        print('[TEST_DIRECTORY] ' + str(directory), flush=True)
        try:
            with lock(state/'run.lock'):
                snapshot(project, directory/'before')
                plan = read(directory/'before/plan.json')
                baseline_check(plan, read(directory/'before/batch_summary.json'), read(directory/'before/status.json'))
                if Path(plan['config']['OUTPUT_DIR']).resolve() != project:
                    raise ValueError('Saved OUTPUT_DIR differs from --project')
                validate_release(cf)
                config, samples = Path(plan['config_path']), Path(plan['samples_path'])
                shutil.copy2(config, directory/'chromaflow.config')
                shutil.copy2(samples, directory/'samples.tsv')
            base = [cf, 'batch']
            meta['plan_exit'] = stream_command(base + ['plan', '-c', config, '--samples', samples, '--check-env'],
                                               project, directory/'plan.console.log')
            if meta['plan_exit']:
                raise ValueError('Environment/input plan check failed; batch run was not started.')
            proposed = read(state/'proposed_plan.json')
            differences = [k for k in PLAN_FIELDS if plan.get(k) != proposed.get(k)]
            save(directory/'plan-comparison.json', dict(changed_fields=differences))
            shutil.copy2(state/'proposed_plan.json', directory/'proposed_plan.json')
            if differences:
                raise ValueError('Unchanged-input test refused; plan fields differ: ' + ', '.join(differences))
            regression_checks(cf, plan, directory)
            meta['state'] = 'running'; save(directory/'run.json', meta)
            batch_exit = stream_command(base + ['run', '-c', config, '--samples', samples, '--resume', '--until', 'candidates'],
                                        project, directory/'resume.console.log')
            meta['batch_exit'] = batch_exit
            print(f'[BATCH_EXIT] {batch_exit}', flush=True)
            meta['reports'] = refresh_reports(cf, project, config, directory)
            with lock(state/'run.lock'):
                snapshot(project, directory/'after')
                result = audit_resume(directory/'before', directory/'after',
                                      (directory/'resume.console.log').read_text(encoding='utf-8'), batch_exit)
                reports = audit_reports(project, read(state/'plan.json'), read(state/'status.json'))
            save(directory/'resume_acceptance.json', result)
            tsv(directory/'stage_reuse.tsv', ['stage', 'result', 'reason', 'events'], result['stages'])
            tsv(directory/'native_peak_reuse.tsv', ['group', 'task', 'result', 'status', 'reason'], result['native_peak_tasks'])
            tsv(directory/'native_core_reuse.tsv', ['group', 'process', 'tag', 'result', 'status'], result['native_core_tasks'])
            tsv(directory/'native_process_counts.tsv', ['group', 'process', 'status', 'count'], result['native_process_counts'])
            save(directory/'report_acceptance.json', reports)
            tsv(directory/'report_inventory.tsv', ['analysis_id', 'plot', 'status', 'reason', 'width', 'height'], reports['plots'])
            meta.update(state='finished', resume_result=result['result'], report_files_result=reports['result'])
            print('[RESUME_ACCEPTANCE] ' + result['result'], flush=True)
            print('[ANALYSIS_COUNTS] ' + '; '.join(f'{k}={v}' for k, v in sorted(result['after_counts'].items())), flush=True)
            print('[REPORT_FILES] ' + reports['result'], flush=True)
            for issue in result['issues'] + reports['issues']:
                print('[REVIEW] ' + issue, flush=True)
            code = 0 if result['result'] == 'PASS' and not reports['issues'] else 1
        except KeyboardInterrupt:
            meta.update(state='interrupted', error='Interrupted by user'); code = 130
        except (OSError, ValueError, KeyError) as exc:
            meta.update(state='incomplete', error=str(exc)); code = 2
            print('[ERROR] ' + str(exc), flush=True)
        meta.update(finished=now(), test_exit=code)
        save(directory/'run.json', meta)
        print(f'[ACCEPTANCE_EXIT] {code}\n[FEEDBACK] Run the pack command even if this test needs review.', flush=True)
        return code


def collect_files(project, test_dir, cf):
    """Explicit allowlist; never recurse through work, raw_data or references."""
    files, excluded = {}, []
    def add(path, arcname=None, root=None):
        root = root or project
        path = Path(path)
        if not path.is_file():
            return
        if not inside(path, root):
            excluded.append(dict(path=str(path), reason='symlink resolves outside the allowed root'))
            return
        name = arcname or ('project/' + path.relative_to(project).as_posix())
        files[name] = path
    for path in sorted(test_dir.rglob('*')):
        if path.suffix in ('.json', '.tsv', '.log', '.config', '.txt', '.png'):
            add(path, 'test/' + path.relative_to(test_dir).as_posix())
    state = project/'.chromaflow/batch'
    for relative in ('status.json', 'plan.json', 'reports.json', 'proposed_plan.json', 'environment_check.json'):
        add(state/relative)
    for path in (state/'environments').glob('*.json'):
        add(path)
    for path in (state/'annotation_metadata_backup').rglob('*'):
        if path.suffix in ('.tsv', '.json'): add(path)
    reference = read(state/'status.json').get('reference', {}).get('details', {})
    if reference.get('CHROM_SIZES'):
        add(Path(reference['CHROM_SIZES']))
    for path in (state/'runs').glob('*/*'):
        if path.name in ('outputs.json', 'peak_tasks.json', 'nfcore.config', 'executor.config'):
            add(path)
    for path in (project/'results/summary').glob('*'):
        if path.suffix in ('.tsv', '.json'):
            add(path)
    # All batch/downstream logs include prior report failures and successful
    # report reconstruction, even when the latest reports command was cached.
    for folder in ('logs/batch', 'logs/analyses'):
        for path in (project/folder).rglob('*'):
            if path.suffix in ('.log', '.tsv', '.json'):
                add(path)
    plan = read(state/'plan.json')
    for aid in plan['analyses']:
        root = project/'results/analyses'/aid
        for relative in ('analysis_status.json', 'inputs.json', 'annotation/peak_annotation_full.tsv',
                         'annotation/genomic_distribution.tsv', 'annotation/genes.bed',
                         'annotation/annotation_policy.json', 'candidates/selection.json'):
            add(root/relative)
        for path in (root/'reports/annotation').glob('*'):
            if path.suffix in ('.png', '.tsv', '.json'):
                add(path)
        for name in ('analysis.config', 'selected_peaks.bed'):
            add(state/'analyses'/aid/name)
    # Only command diagnostics for failed MACS tasks; never their BAMs.
    for key, record in read(state/'status.json').items():
        if not key.startswith('nfcore/'):
            continue
        for failure in record.get('details', {}).get('failures', {}).values():
            for task in failure.get('tasks', []):
                diagnostic = task.get('diagnostic_log')
                if diagnostic and inside(diagnostic, project/'work/batch'):
                    for name in ('.command.err', '.command.out', '.command.sh', '.exitcode'):
                        add(Path(diagnostic).parent/name)
    installation = cf.parent.parent
    add(installation/'release/compatibility.json', 'installation/compatibility.json', installation)
    add(Path(__file__), 'kit/acceptance.py', Path(__file__).parent)
    add(installation/'docs/RELEASE_ACCEPTANCE.md', 'kit/RELEASE_ACCEPTANCE.md', installation)
    return files, excluded


def pack(args):
    project, cf = Path(args.project).resolve(), Path(args.cf).resolve()
    state = project/'.chromaflow/batch'
    if args.test_dir:
        directory = Path(args.test_dir).resolve()
    else:
        directory = Path(read(project/'diagnostics/latest-release-acceptance.json')['directory']).resolve()
    if not inside(directory, project/'diagnostics') or not (directory/'run.json').is_file():
        raise ValueError('Select an acceptance test directory under this project/diagnostics.')
    meta = read(directory/'run.json')
    if meta.get('project') != str(project):
        raise ValueError('The selected test belongs to a different project.')
    cf = Path(meta['cf']).resolve()
    with lock(state/'resume-acceptance.lock'), lock(state/'run.lock'):
        if meta.get('state') == 'running':
            print('[WARN] Interrupted/incomplete test: collecting available evidence.', flush=True)
        reports = audit_reports(project, read(state/'plan.json'), read(state/'status.json'))
        save(directory/'report_acceptance_at_pack.json', reports)
        tsv(directory/'report_inventory_at_pack.tsv', ['analysis_id', 'plot', 'status', 'reason', 'width', 'height'], reports['plots'])
        files, excluded = collect_files(project, directory, cf)
        archive_name = 'ChromaFlow-' + RELEASE + '-acceptance-' + dt.datetime.now().strftime('%Y%m%d-%H%M%S') + '-' + str(os.getpid())
        archive = Path(args.outdir).expanduser().resolve()/(archive_name + '.tar.gz')
        archive.parent.mkdir(parents=True, exist_ok=True)
        manifest = dict(kit=KIT_VERSION, created=now(), project=str(project), test=str(directory),
                        report_files_result=reports['result'], excluded=excluded,
                        files=[dict(member=name, size=path.stat().st_size) for name, path in sorted(files.items())],
                        note='No FASTQs, BAMs, reference FASTA, containers or new transfer checksum files.')
        payload = (json.dumps(manifest, indent=2) + '\n').encode('utf-8')
        partial = archive.with_name(archive.name + '.partial')
        if archive.exists():
            raise ValueError('Refusing to overwrite an existing archive: ' + str(archive))
        with tarfile.open(partial, 'x:gz', dereference=True) as tar:
            member = tarfile.TarInfo(archive_name + '/MANIFEST.json')
            member.size = len(payload); member.mtime = int(time.time()); member.mode = 0o644
            tar.addfile(member, io.BytesIO(payload))
            for name, path in sorted(files.items()):
                tar.add(path, arcname=archive_name + '/' + name, recursive=False)
        partial.rename(archive)
    print(f'[ARCHIVE] {archive}\n[SIZE_MIB] {archive.stat().st_size / 1048576:.2f}\n[FILES] {len(files)}')
    print('[REPORT_FILES] ' + reports['result'])
    if excluded:
        print(f'[WARN] {len(excluded)} external symlink targets excluded; see MANIFEST.json')
    print('[PACK_EXIT] 0 (collection succeeded; not a scientific acceptance result)')
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='action', required=True)
    for name in ('run', 'pack'):
        sub = commands.add_parser(name)
        sub.add_argument('--project', default=DEFAULT_PROJECT)
        sub.add_argument('--cf', default=DEFAULT_CF)
        if name == 'pack':
            sub.add_argument('--test-dir', help='Defaults to the latest test directory recorded by this helper')
            sub.add_argument('--outdir', default='.')
    args = parser.parse_args()
    try:
        return run(args) if args.action == 'run' else pack(args)
    except (OSError, ValueError, KeyError) as exc:
        print('[ERROR] ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
