"""Acceptance gates must reject missing evidence, rerendering and failed tests."""
import importlib.util
import argparse
import os
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('release_acceptance', ROOT/'tools/acceptance.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class AcceptanceWorkflowTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name).resolve()
        self.project = self.root/'project'
        self.directory = self.project/'diagnostics/test'
        self.directory.mkdir(parents=True)
        self.cf = ROOT/'bin/chromaflow'

    def report_fixture(self):
        path = self.project/'results/analyses/A/reports/annotation/image.png'
        path.parent.mkdir(parents=True)
        path.write_bytes(b'fixture')
        state = self.project/'.chromaflow/batch/reports.json'
        a.save(state, {'A': {'status': 'complete', 'outputs': {str(path): 'existing-identity'}}})
        return state, path

    def test_packaged_runtime_gate_accepts_the_complete_release(self):
        self.assertEqual(a.validate_release(self.cf).VERSION, a.RELEASE)

    def test_report_refresh_and_repeat_preserve_records_and_file_stats(self):
        self.report_fixture()
        with patch.object(a, 'stream_command', return_value=0) as run:
            result = a.refresh_reports(self.cf, self.project, self.project/'config', self.directory)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(run.call_count, 2)
        self.assertTrue(all(call.args[0][1:3] == ['batch', 'reports'] for call in run.call_args_list))

    def test_failed_report_command_is_not_reported_as_success(self):
        self.report_fixture()
        with patch.object(a, 'stream_command', return_value=1) as run:
            with self.assertRaisesRegex(ValueError, 'refresh failed'):
                a.refresh_reports(self.cf, self.project, self.project/'config', self.directory)
        self.assertEqual(run.call_count, 1)
        self.assertFalse((self.directory/'report-idempotence.json').exists())

    def test_unexpected_rewrite_on_repeat_is_rejected(self):
        _, path = self.report_fixture()
        calls = []
        def run(*args):
            calls.append(args)
            if len(calls) == 2:
                path.write_bytes(b'changed report content')
            return 0
        with patch.object(a, 'stream_command', side_effect=run):
            with self.assertRaisesRegex(ValueError, 'did not preserve'):
                a.refresh_reports(self.cf, self.project, self.project/'config', self.directory)

    def test_regression_failure_is_recorded_and_rejected(self):
        plan = dict(config={'BATCH_ENV_MODE': 'current'})
        with patch.object(a, 'stream_command', side_effect=[1, 0, 0, 0]) as run:
            with self.assertRaisesRegex(ValueError, 'regressions failed or skipped'):
                a.regression_checks(self.cf, plan, self.directory)
        self.assertEqual(run.call_count, 4)
        self.assertEqual(a.read(self.directory/'regressions.json')[0]['exit'], 1)

    def saved_project(self):
        state, _ = self.report_fixture()
        state = state.parent
        config = self.project/'chromaflow.config'; config.write_text('fixture config\n')
        samples = self.project/'samples.tsv'; samples.write_text('fixture samples\n')
        plan = dict(version='1.2.8', config={'OUTPUT_DIR': str(self.project), 'MAKE_PLOTS': 'true'},
                    config_path=str(config), samples_path=str(samples),
                    analyses={'A': {}}, all_analysis_ids=['A'], groups={'se_k100': {}},
                    sources={}, fastqs={}, step_code_sha256={}, samples_sha256='fixture')
        records = {}
        for key in ('reference', 'nfcore/se_k100', 'A/annotate', 'A/motif', 'A/candidates'):
            path = self.project/(key.replace('/', '-')+'.txt'); path.write_text('unchanged result\n')
            records[key] = dict(status='complete', signature='same-'+key, finished='previous-run',
                                outputs={str(path): 'existing-output-identity'})
        a.save(state/'status.json', records)
        a.save(state/'plan.json', plan)
        a.save(self.project/'results/summary/batch_summary.json',
               {'analyses': {'A': {'status': 'complete', 'peak_count': 20, 'candidates': 2}}})
        return state, plan, records

    @unittest.skipUnless(os.name == 'posix', 'Linux project locks required')
    def test_run_and_pack_keep_scientific_state_and_separate_acceptance_results(self):
        state, plan, records = self.saved_project()
        calls = []
        def producer(command, cwd, logfile):
            calls.append(command[2])
            content = ''
            if command[2] == 'plan':
                a.save(state/'proposed_plan.json', dict(plan, version='1.2.9'))
            elif command[2] == 'run':
                content = ''.join('[RESUME] '+key+': verified outputs\n' for key in records)
            logfile.write_text(content)
            return 0
        reports = dict(result='FILES_OK_VISUAL_REVIEW_PENDING', plots=[], issues=[])
        args = argparse.Namespace(project=str(self.project), cf=str(self.cf))
        with patch.object(a, 'stream_command', side_effect=producer), \
             patch.object(a, 'regression_checks') as regressions, \
             patch.object(a, 'audit_reports', return_value=reports):
            self.assertEqual(a.run(args), 0)
            self.assertEqual(calls, ['plan', 'run', 'reports', 'reports'])
            regressions.assert_called_once()
            latest = a.read(self.project/'diagnostics/latest-release-acceptance.json')
            directory = Path(latest['directory'])
            meta = a.read(directory/'run.json')
            self.assertEqual((meta['test_exit'], meta['batch_exit'], meta['resume_result']), (0, 0, 'PASS'))
            args.test_dir = None; args.outdir = str(self.root/'feedback')
            self.assertEqual(a.pack(args), 0)
        self.assertEqual(a.read(state/'status.json'), records)
        archive = next((self.root/'feedback').glob('*.tar.gz'))
        with tarfile.open(archive) as stream:
            names = stream.getnames()
        self.assertTrue(any(n.endswith('/test/run.json') for n in names))
        self.assertTrue(any(n.endswith('/test/report-idempotence.json') for n in names))
        self.assertFalse(any(n.endswith(('.bam', '.fastq.gz', '.sha256')) for n in names))

    @unittest.skipUnless(os.name == 'posix', 'Linux project locks required')
    def test_changed_plan_stops_before_regressions_or_analysis(self):
        state, plan, _ = self.saved_project()
        def producer(command, cwd, logfile):
            a.save(state/'proposed_plan.json', dict(plan, samples_sha256='changed-input'))
            logfile.write_text('changed plan\n')
            return 0
        args = argparse.Namespace(project=str(self.project), cf=str(self.cf))
        with patch.object(a, 'stream_command', side_effect=producer) as commands, \
             patch.object(a, 'regression_checks') as regressions:
            self.assertEqual(a.run(args), 2)
            self.assertEqual(commands.call_count, 1)
            regressions.assert_not_called()


if __name__ == '__main__':
    unittest.main()
