"""Cache migration must distinguish unstable diagnostics from scientific changes."""
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

WARNING = ('[host:123] mca_base_component_repository_open: unable to open mca_btl_openib: '
           'librdmacm.so.1: cannot open shared object file: No such file or directory (ignored)\n')


class EnvironmentTests(unittest.TestCase):
    def setUp(self):
        self.identity = dict(environment=['conda', 'run', '-n', 'motif'],
                             tools='CONDA_PREFIX=/env\n/env/bin/meme-chip\nPython 3.12.13\n'+WARNING+'5.5.7\n')

    def test_only_known_ignored_diagnostic_is_excluded(self):
        expected = dict(self.identity, tools=self.identity['tools'].replace(WARNING, ''))
        self.assertEqual(b.stable_environment(self.identity), expected)
        for raw in (self.identity['tools'].replace('123', '987654'), self.identity['tools'].replace('host:', 'another-host:')):
            self.assertEqual(b.stable_environment(dict(self.identity, tools=raw)), expected)
        for old, new in [('5.5.7', '5.5.8'), ('/env/bin', '/other/bin'), ('librdmacm.so.1', 'libmpi.so.1'),
                         ('(ignored)', '(fatal)')]:
            changed = dict(self.identity, tools=self.identity['tools'].replace(old, new))
            self.assertNotEqual(b.stable_environment(changed), expected)
        self.assertIn(WARNING, self.identity['tools'])

    def test_saved_signature_can_migrate_but_changed_outputs_or_inputs_cannot(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = dict(config=dict(OUTPUT_DIR=str(root)), code_sha256='runtime')
            runner = b.Runner(plan, True, 'candidates')
            old = self.identity
            new = dict(old, tools=old['tools'].replace('123', '456'))
            runner.envs = dict(motif=dict(identity=new))
            b.save(runner.state/'environments/old.json', dict(motif=dict(identity=old)))
            output = root/'motif.txt'; output.write_text('real result')
            inputs = ['same inputs', 'same parameters']
            oldsig = b.digest(dict(step='motif', inputs=inputs, code='runtime', environment=old))
            runner.records['A/motif'] = dict(status='complete', signature=oldsig, outputs={str(output): b.sha(output)})
            canonical = runner.signature('motif', inputs)
            with patch.object(b, 'artifact_files', side_effect=AssertionError('must not rerun')):
                record = runner.stage('A/motif', canonical, lambda: self.fail('reran motif'))
            self.assertEqual(record['signature'], oldsig)
            self.assertEqual(record['resume_signature'], canonical)
            self.assertTrue(b.valid_record(record, canonical))
            self.assertFalse(b.valid_record(record, runner.signature('motif', inputs + ['changed'])))
            for environment in (dict(new, tools=new['tools'].replace('5.5.7', '6.0.0')),
                                dict(new, environment=['another environment'])):
                runner.envs['motif']['identity'] = environment
                self.assertFalse(b.valid_record(record, runner.signature('motif', inputs)))
            output.write_text('changed result')
            self.assertFalse(b.valid_record(record, canonical))


@unittest.skipUnless(os.name == 'posix', 'Linux/WSL Bash required')
class ResumeIntegrationTests(unittest.TestCase):
    def fixture(self):
        from test_batch import BatchTests
        f = BatchTests(methodName='test_shared_control_and_independent_ids')
        f.setUp(); self.addCleanup(f.doCleanups)
        return f

    def test_pid_migration_keeps_both_motif_and_candidates_and_repeat_resume(self):
        from test_batch import FixtureRunner
        f = self.fixture()
        plan = f.plan()
        def env(pid):
            data = {s: 'fixture environment' for s in b.STEPS}
            data['motif'] = dict(identity=dict(tools='meme 5.5.7\n'+WARNING.replace('123', str(pid))), log='raw')
            return data
        class LegacyRunner(FixtureRunner):
            def signature(self, step, inputs):
                if step != 'motif':
                    return super().signature(step, inputs)
                return b.digest(dict(step=step, inputs=inputs, code=self.plan['step_code_sha256'][step],
                                     environment=self.envs[step]['identity']))
        with patch.object(b, 'preflight', return_value=env(101)):
            before = LegacyRunner(plan, False, 'candidates')
            self.assertEqual(before.run(), 0)
        original = copy.deepcopy(before.records)
        for pid in (202, 303):
            with patch.object(b, 'preflight', return_value=env(pid)):
                after = FixtureRunner(plan, True, 'candidates')
                self.assertEqual(after.run(), 0)
            self.assertEqual(after.calls, [])
            for aid in plan['analyses']:
                for step in ('annotate', 'motif', 'candidates'):
                    self.assertEqual(after.records[f'{aid}/{step}']['outputs'], original[f'{aid}/{step}']['outputs'])
        with patch.object(b, 'preflight', return_value=env(404)):
            changed = copy.deepcopy(plan)
            changed['config']['MEME_SEARCHSIZE'] = '100'
            rerun = FixtureRunner(changed, True, 'candidates')
            self.assertEqual(rerun.run(), 0)
            self.assertEqual(rerun.calls, [('motif', 'TF_A'), ('candidates', 'TF_A'), ('motif', 'TF_B'), ('candidates', 'TF_B')])

    def test_reports_are_separate_idempotent_and_do_not_invalidate_resume(self):
        from test_batch import FixtureRunner
        f = self.fixture(); plan = f.plan()
        plan['version'] = '1.2.8'  # Historical reporter typo belongs to this run.
        before = FixtureRunner(plan, False, 'candidates')
        self.assertEqual(before.run(), 0)
        original = copy.deepcopy(before.records)
        # Reproduce the old reporter typo, keeping the saved output hash honest.
        selection = f.root/'results/analyses/TF_A/candidates/selection.json'
        data = b.read_json(selection); data['version'] = '1.2.7'; b.save(selection, data)
        before.records['TF_A/candidates']['outputs'][str(selection)] = b.sha(selection)
        b.save(before.status_path, before.records)
        statuses = (f.root/'results/summary/batch_status.tsv').read_bytes()
        renders = []
        def render(cmd, log, cwd):
            self.assertEqual(cmd[0], 'Rscript')
            staged = Path(cmd[4]); renders.append(staged)
            self.assertTrue(Path(cmd[-3]).is_file())
            self.assertTrue(Path(cmd[-2]).is_file())
            self.assertEqual(cmd[-1], 'non_idr')
            b.tsv(staged/'plot_status.tsv', ['plot', 'status', 'reason'], [['tss_avgprof.png', 'skipped', 'No retained windows']])
            b.write(staged/'anno_upset.png', 'fixture plot')
        runner = b.Runner(plan, True, 'candidates')
        with patch.object(runner, 'execute', side_effect=render):
            self.assertEqual(runner.reports(), 0)
            self.assertEqual(runner.reports(), 0)
        self.assertEqual(len(renders), 2)
        self.assertEqual(b.read_json(selection)['version'], '1.2.8')
        corrected = b.read_json(selection); corrected['version'] = '1.2.7'
        self.assertEqual(corrected, data)
        self.assertEqual((f.root/'results/summary/batch_status.tsv').read_bytes(), statuses)
        resumed = FixtureRunner(plan, True, 'candidates')
        self.assertEqual(resumed.run(), 0)
        self.assertEqual(resumed.calls, [])
        for aid in plan['analyses']:
            for step in ('annotate', 'motif'):
                self.assertEqual(resumed.records[f'{aid}/{step}'], original[f'{aid}/{step}'])
        # A modified upstream artifact must block report refresh and metadata edits.
        (f.root/'results/analyses/TF_A/annotation/genes.bed').write_text('modified')
        with patch.object(runner, 'execute', side_effect=AssertionError('must refuse')):
            self.assertEqual(runner.reports('TF_A'), 1)

    def test_failed_report_render_keeps_scientific_results_and_previous_report(self):
        from test_batch import FixtureRunner
        f = self.fixture(); plan = f.plan()
        runner = FixtureRunner(plan, False, 'candidates')
        self.assertEqual(runner.run(), 0)
        before = copy.deepcopy(runner.records)
        report = f.root/'results/analyses/TF_A/reports/annotation/old.txt'
        b.write(report, 'previous diagnostic')
        with patch.object(runner, 'execute', side_effect=b.BatchError('Reconstructed annotation differs')):
            self.assertEqual(runner.reports('TF_A'), 1)
        self.assertEqual(report.read_text(), 'previous diagnostic')
        self.assertEqual(runner.records, before)


if __name__ == '__main__':
    unittest.main()
