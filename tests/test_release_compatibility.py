"""Release upgrades retain reviewed caches, but never altered evidence."""
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'lib'))
import batch as b


def runtime_tree():
    return {p.relative_to(ROOT).as_posix(): b.sha(p) for folder in ('lib', 'bin')
            for p in sorted((ROOT/folder).rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}


class ReleaseCompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.home = Path(self.tmp.name)
        self.file = self.home/'release/compatibility.json'
        self.data = b.read_json(ROOT/'release/compatibility.json')
        self.code = runtime_tree()

    def test_packaged_release_runtime_is_exact(self):
        self.assertEqual(self.data['version'], b.VERSION)
        self.assertEqual(self.data['runtime'], self.code)
        self.assertEqual(b.resume_code_ids(self.code), self.data['step_code_sha256'])

    def test_fresh_install_needs_no_previous_hotfix_directories(self):
        b.save(self.file, self.data)
        with patch.object(b, 'HOME', self.home):
            self.assertEqual(b.resume_code_ids(self.code), self.data['step_code_sha256'])
        self.assertFalse((self.home/'.chromaflow-hotfixes').exists())

    def test_later_runtime_edits_disable_release_reuse(self):
        b.save(self.file, self.data)
        with patch.object(b, 'HOME', self.home):
            for name in ('lib/steps/support/candidates.py', 'lib/steps/05_peak_annotation.R', 'bin/chromaflow'):
                changed = dict(self.code, **{name: '0'*64})
                self.assertEqual(set(b.resume_code_ids(changed).values()), {b.digest(changed)})

    def test_bad_manifest_never_uses_partial_stage_aliases(self):
        cases = [[], {}, dict(self.data, version='unreviewed'),
                 dict(self.data, step_code_sha256={'reference': '0'*64}),
                 dict(self.data, step_code_sha256=dict.fromkeys(b.STEPS, 'not-a-digest'))]
        with patch.object(b, 'HOME', self.home):
            for value in cases:
                b.save(self.file, value)
                self.assertEqual(set(b.resume_code_ids(self.code).values()), {b.digest(self.code)})
            self.file.write_text('invalid json')
            self.assertEqual(set(b.resume_code_ids(self.code).values()), {b.digest(self.code)})

    def test_prior_result_version_is_preserved_without_relabeled_evidence(self):
        plan = dict(version=b.VERSION, config={'OUTPUT_DIR': str(self.home)}, code_sha256='fixture')
        runner = b.Runner(plan, True, 'candidates')
        selection = runner.results/'analyses/A/candidates/selection.json'
        b.save(selection, {'version': '1.2.8', 'status': 'complete'})
        runner.records['A/candidates'] = dict(status='complete', signature='fixture', outputs={str(selection): b.sha(selection)})
        before = selection.read_bytes(), copy.deepcopy(runner.records), selection.stat().st_mtime_ns
        runner.refresh_selection_version('A')
        self.assertEqual(before, (selection.read_bytes(), runner.records, selection.stat().st_mtime_ns))
        selection.write_text('{"version":"1.2.8","status":"changed"}')
        with self.assertRaisesRegex(b.BatchError, 'outputs changed'):
            runner.refresh_selection_version('A')


@unittest.skipUnless(os.name == 'posix', 'Linux/WSL Bash required')
class ReleaseResumeTests(unittest.TestCase):
    def test_upgrade_reuses_completed_stages_but_honors_new_filters(self):
        from test_batch import BatchTests, FixtureRunner
        fixture = BatchTests(methodName='test_shared_control_and_independent_ids')
        fixture.setUp(); self.addCleanup(fixture.doCleanups)
        current = fixture.plan()
        old = copy.deepcopy(current)
        old.update(version='1.2.8', code_sha256='previous-full-runtime')
        first = FixtureRunner(old, False, 'candidates')
        self.assertEqual(first.run(), 0)
        original = copy.deepcopy(first.records)
        upgraded = FixtureRunner(current, True, 'candidates')
        self.assertEqual(upgraded.run(), 0)
        self.assertEqual(upgraded.calls, [])
        self.assertEqual(upgraded.records, original)
        changed = copy.deepcopy(current)
        changed['config']['FIMO_PVAL'] = '1e-6'
        rerun = FixtureRunner(changed, True, 'candidates')
        self.assertEqual(rerun.run(), 0)
        self.assertEqual(rerun.calls, [('candidates', 'TF_A'), ('candidates', 'TF_B')])


if __name__ == '__main__':
    unittest.main()
