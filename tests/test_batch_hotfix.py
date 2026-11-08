"""Only verified, compatible hotfix chains may retain stage cache identities."""
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib'))
import batch as b


class HotfixCodeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.home = Path(self.tmp.name)
        self.marker = self.home/'.chromaflow-hotfixes/fimo-zero-hits-v1/manifest.json'
        self.marker.parent.mkdir(parents=True)
        self.before = {name: 'a'*64 for name in ('bin/chromaflow', 'lib/batch.py',
            'lib/steps/02_run_nfcore.sh', 'lib/steps/02_run_nfcore_conda.sh',
            'lib/steps/05_peak_annotation.R', 'lib/common.sh', 'lib/steps/05_peak_annotation.sh',
            'lib/steps/support/lookup.py',
            'lib/steps/support/annotation_policy.py',
            'lib/steps/support/fimo_tracks.py', 'lib/steps/support/motif_selection.py')}
        self.after = {name: ('b' if name in ('lib/batch.py', 'lib/steps/support/fimo_tracks.py',
                        'lib/steps/support/motif_selection.py') else 'a')*64 for name in self.before}
        self.data = dict(schema=1, patch='fimo-zero-hits-v1', before_code=self.before,
                         after_code=self.after, compatible_steps=list(b.STEPS[:-1]))
        self.home_patch = patch.object(b, 'HOME', self.home)
        self.home_patch.start(); self.addCleanup(self.home_patch.stop)

    def save(self):
        self.marker.write_text(json.dumps(self.data))

    def test_verified_patch_preserves_only_unaffected_stages(self):
        self.save()
        result = b.resume_code_ids(self.after)
        self.assertEqual(result['candidates'], b.digest(self.after))
        self.assertTrue(all(result[s] == b.digest(self.before) for s in b.STEPS[:-1]))

    def test_later_code_changes_disable_compatibility(self):
        self.save()
        changed = dict(self.after, **{'bin/chromaflow': 'c'*64})
        self.assertEqual(set(b.resume_code_ids(changed).values()), {b.digest(changed)})

    def test_unreviewed_changes_cannot_be_aliased(self):
        self.after['bin/chromaflow'] = 'c'*64
        self.save()
        self.assertEqual(set(b.resume_code_ids(self.after).values()), {b.digest(self.after)})

    def test_missing_or_malformed_manifest_fails_closed_for_reuse(self):
        for text in (None, '{}', 'invalid json', '[]'):
            if text is not None:
                self.marker.write_text(text)
            self.assertEqual(set(b.resume_code_ids(self.after).values()), {b.digest(self.after)})

    def read_length_patch(self):
        self.save()
        code = dict(self.after)
        for name in ('lib/batch.py', 'lib/steps/02_run_nfcore.sh', 'lib/steps/02_run_nfcore_conda.sh'):
            code[name] = 'c'*64
        marker = self.home/'.chromaflow-hotfixes/read-length-v1/manifest.json'
        marker.parent.mkdir(parents=True)
        data = dict(schema=1, patch='read-length-v1', before_code=self.after,
                    after_code=code, compatible_steps=list(b.STEPS))
        marker.write_text(json.dumps(data))
        return code, marker, data

    def test_read_length_patch_preserves_verified_fimo_stage_identities(self):
        code, _, _ = self.read_length_patch()
        result = b.resume_code_ids(code)
        self.assertEqual(result['candidates'], b.digest(self.after))
        self.assertTrue(all(result[s] == b.digest(self.before) for s in b.STEPS[:-1]))

    def test_missing_or_mismatched_predecessor_invalidates_entire_chain(self):
        code, _, _ = self.read_length_patch()
        self.marker.unlink()
        self.assertEqual(set(b.resume_code_ids(code).values()), {b.digest(code)})
        self.data['after_code'] = dict(self.after, **{'lib/batch.py': 'd'*64})
        self.save()
        self.assertEqual(set(b.resume_code_ids(code).values()), {b.digest(code)})

    def test_unreviewed_or_later_edits_cannot_reuse_chained_cache(self):
        code, marker, data = self.read_length_patch()
        later = dict(code, **{'bin/chromaflow': 'd'*64})
        self.assertEqual(set(b.resume_code_ids(later).values()), {b.digest(later)})
        data['after_code'] = later
        marker.write_text(json.dumps(data))
        self.assertEqual(set(b.resume_code_ids(later).values()), {b.digest(later)})

    def test_invalid_read_length_manifest_cannot_partially_reuse_cache(self):
        code, marker, data = self.read_length_patch()
        data['compatible_steps'] = ['reference']
        marker.write_text(json.dumps(data))
        self.assertEqual(set(b.resume_code_ids(code).values()), {b.digest(code)})
        marker.write_text('[]')
        self.assertEqual(set(b.resume_code_ids(code).values()), {b.digest(code)})

    def test_peak_isolation_patch_preserves_the_predecessor_stage_identities(self):
        previous, _, _ = self.read_length_patch()
        expected = b.resume_code_ids(previous)
        current = dict(previous, **{'lib/batch.py': 'd'*64})
        marker = self.home/'.chromaflow-hotfixes/peak-isolation-v1/manifest.json'
        marker.parent.mkdir(parents=True)
        data = dict(schema=1, patch='peak-isolation-v1', before_code=previous,
                    after_code=current, compatible_steps=list(b.STEPS))
        marker.write_text(json.dumps(data))
        self.assertEqual(b.resume_code_ids(current), expected)
        unrelated = dict(current, **{'bin/chromaflow': 'e'*64})
        data['after_code'] = unrelated
        marker.write_text(json.dumps(data))
        self.assertEqual(set(b.resume_code_ids(unrelated).values()), {b.digest(unrelated)})

    def test_peak_isolation_manifest_requires_its_read_length_predecessor(self):
        previous, predecessor, _ = self.read_length_patch()
        current = dict(previous, **{'lib/batch.py': 'd'*64})
        marker = self.home/'.chromaflow-hotfixes/peak-isolation-v1/manifest.json'
        marker.parent.mkdir(parents=True)
        marker.write_text(json.dumps(dict(schema=1, patch='peak-isolation-v1', before_code=previous,
            after_code=current, compatible_steps=list(b.STEPS))))
        predecessor.unlink()
        self.assertEqual(set(b.resume_code_ids(current).values()), {b.digest(current)})

    def test_report_patch_allows_only_its_reviewed_new_helper(self):
        previous, _, _ = self.read_length_patch()
        peak = dict(previous, **{'lib/batch.py': 'd'*64})
        marker = self.home/'.chromaflow-hotfixes/peak-isolation-v1/manifest.json'
        marker.parent.mkdir(parents=True)
        marker.write_text(json.dumps(dict(schema=1, patch='peak-isolation-v1', before_code=previous,
            after_code=peak, compatible_steps=list(b.STEPS))))
        expected = b.resume_code_ids(peak)
        current = dict(peak)
        for name in ('lib/batch.py', 'lib/steps/05_peak_annotation.R',
                     'lib/steps/annotation_plots.R', 'lib/steps/support/motif_selection.py'):
            current[name] = 'e'*64
        report = self.home/'.chromaflow-hotfixes/resume-reports-v1/manifest.json'
        report.parent.mkdir(parents=True)
        data = dict(schema=1, patch='resume-reports-v1', before_code=peak,
                    after_code=current, compatible_steps=list(b.STEPS))
        report.write_text(json.dumps(data))
        self.assertEqual(b.resume_code_ids(current), expected)
        guarded = dict(current)
        for name in ('lib/batch.py', 'lib/steps/05_peak_annotation.R', 'lib/steps/support/annotation_policy.py'):
            guarded[name] = 'f'*64
        guard = self.home/'.chromaflow-hotfixes/report-table-guard-v1/manifest.json'
        guard.parent.mkdir(parents=True)
        guard.write_text(json.dumps(dict(schema=1, patch='report-table-guard-v1', before_code=current,
            after_code=guarded, compatible_steps=list(b.STEPS))))
        self.assertEqual(b.resume_code_ids(guarded), expected)
        repaired = dict(guarded)
        for name in ('lib/batch.py', 'lib/common.sh', 'lib/steps/05_peak_annotation.R',
                     'lib/steps/05_peak_annotation.sh', 'lib/steps/support/lookup.py',
                     'lib/steps/tss_windows.R', 'lib/steps/support/peak_metadata.py'):
            repaired[name] = '1'*64
        repair = self.home/'.chromaflow-hotfixes/tss-idr-metadata-v1/manifest.json'
        repair.parent.mkdir(parents=True)
        repair.write_text(json.dumps(dict(schema=1, patch='tss-idr-metadata-v1', before_code=guarded,
            after_code=repaired, compatible_steps=list(b.STEPS))))
        self.assertEqual(b.resume_code_ids(repaired), expected)
        unrelated_repair = dict(repaired, **{'lib/steps/support/motif_selection.py': '2'*64})
        self.assertEqual(set(b.resume_code_ids(unrelated_repair).values()), {b.digest(unrelated_repair)})
        repair.unlink()
        unreviewed_guard = dict(guarded, **{'bin/chromaflow': 'f'*64})
        self.assertEqual(set(b.resume_code_ids(unreviewed_guard).values()), {b.digest(unreviewed_guard)})
        guard.unlink()
        unreviewed = dict(current, **{'lib/extra.py': 'f'*64})
        data['after_code'] = unreviewed
        report.write_text(json.dumps(data))
        self.assertEqual(set(b.resume_code_ids(unreviewed).values()), {b.digest(unreviewed)})
        data['after_code'] = current
        report.write_text(json.dumps(data))
        marker.unlink()
        self.assertEqual(set(b.resume_code_ids(current).values()), {b.digest(current)})


@unittest.skipUnless(os.name == 'posix', 'Linux/WSL Bash required')
class HotfixResumeTests(unittest.TestCase):
    def test_upgrade_reuses_upstream_and_reruns_candidate_stages(self):
        from test_batch import BatchTests, FixtureRunner
        fixture = BatchTests(methodName='test_shared_control_and_independent_ids')
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        old = fixture.plan()
        old['code_sha256'] = 'old-reviewed-runtime'
        old['step_code_sha256'] = dict.fromkeys(b.STEPS, old['code_sha256'])
        failed = FixtureRunner(old, False, 'candidates')
        failed.fail_analysis = 'TF_A'
        self.assertEqual(failed.run(), 1)

        patched = copy.deepcopy(old)
        patched['code_sha256'] = 'patched-reader-runtime'
        patched['step_code_sha256']['candidates'] = patched['code_sha256']
        resumed = FixtureRunner(patched, True, 'candidates')
        self.assertEqual(resumed.run(), 0, resumed.errors)
        self.assertEqual(resumed.calls, [('candidates', 'TF_A'), ('candidates', 'TF_B')])
        again = FixtureRunner(patched, True, 'candidates')
        self.assertEqual(again.run(), 0)
        self.assertEqual(again.calls, [])
        # A changed upstream output still invalidates reuse despite the alias.
        (fixture.root/'results/analyses/TF_A/annotation/genes.bed').write_text('changed\n')
        changed = FixtureRunner(patched, True, 'candidates')
        self.assertEqual(changed.run(), 0, changed.errors)
        self.assertIn(('annotate', 'TF_A'), changed.calls)


if __name__ == '__main__':
    unittest.main()
