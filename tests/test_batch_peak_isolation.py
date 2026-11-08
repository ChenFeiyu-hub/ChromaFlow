"""A failed peak task must never become a successful or zero-peak analysis."""
import csv
import json
import os
from pathlib import Path
import tempfile
import unittest

import test_batch as fixtures

b = fixtures.b
FixtureRunner = fixtures.FixtureRunner
config_values = fixtures.config_values


class TaskTraceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.trace = self.root/'tasks.tsv'
        self.rows = [['1', b.MACS_PROCESS, 'TF_A_REP1', 'CACHED', '0', str(self.root/'task')]]

    def save(self):
        b.tsv(self.trace, b.TASK_FIELDS, self.rows)

    def test_cached_success_and_unrelated_processes_are_distinct(self):
        self.rows.append(['2', 'OTHER', 'TF_A_REP1', 'COMPLETED', '0', '-'])
        self.save()
        self.assertEqual(set(b.peak_task_records(self.trace, {'TF_A_REP1'})), {'TF_A_REP1'})

    def test_missing_header_missing_task_and_partial_task_are_errors(self):
        with self.assertRaisesRegex(b.BatchError, 'trace is missing'):
            b.peak_task_records(self.trace, {'TF_A_REP1'})
        self.trace.write_text('status\nCOMPLETED\n')
        with self.assertRaisesRegex(b.BatchError, 'header'):
            b.peak_task_records(self.trace, {'TF_A_REP1'})
        self.save()
        with self.assertRaisesRegex(b.BatchError, 'absent from current trace'):
            b.peak_task_records(self.trace, {'TF_A_REP1', 'TF_B_REP1'})
        self.rows[0][3] = 'ABORTED'; self.save()
        with self.assertRaisesRegex(b.BatchError, 'Unfinished'):
            b.peak_task_records(self.trace, {'TF_A_REP1'})

    def test_duplicate_unexpected_and_contradictory_records_are_errors(self):
        self.rows *= 2; self.save()
        with self.assertRaisesRegex(b.BatchError, 'duplicate'):
            b.peak_task_records(self.trace, {'TF_A_REP1'})
        self.rows = self.rows[:1]; self.save()
        with self.assertRaisesRegex(b.BatchError, 'Unexpected'):
            b.peak_task_records(self.trace, {'TF_B_REP1'})
        self.rows[0][4] = '1'; self.save()
        with self.assertRaisesRegex(b.BatchError, 'Inconsistent'):
            b.peak_task_records(self.trace, {'TF_A_REP1'})

    def test_model_failure_requires_actual_error_in_own_task_directory(self):
        work = self.root/'work'; task = work/'aa/bb'
        task.mkdir(parents=True)
        error = task/'.command.err'
        error.write_text('MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can only find 40!\n')
        row = dict(tag='TF_A_REP1', workdir=str(task), exit='1')
        failure = b.peak_task_failure(row, work, self.trace)
        self.assertEqual(failure['status'], 'model_failed')
        self.assertEqual(failure['observed_model_pairs'], 40)
        self.assertEqual(failure['minimum_model_pairs'], 100)
        self.assertEqual(b.peak_task_failure(row, self.root/'another-project', self.trace)['status'], 'peak_calling_failed')
        error.write_text('BAM file could not be read\n')
        self.assertEqual(b.peak_task_failure(row, work, self.trace)['status'], 'peak_calling_failed')

    def test_executor_only_changes_exact_peak_process_and_has_current_trace(self):
        config = b.nfcore_executor_config(dict(MAX_CPUS='8', MAX_MEMORY='32.GB'), self.trace)
        self.assertIn("withName: '"+b.MACS_PROCESS+"'", config)
        self.assertEqual(config.count('errorStrategy'), 1)
        self.assertNotIn('nomodel', config)
        self.assertNotIn('extsize', config)
        self.assertIn('workflow.failOnIgnore = false', config)
        self.assertIn('overwrite = false', config)


class PartialRunner(FixtureRunner):
    fail_tags = {'TF_A_REP1'}
    model_error = True
    missing_output = False
    global_error = False
    incomplete_trace = False

    def execute(self, cmd, log, cwd):
        super().execute(cmd, log, cwd)
        if cmd[cmd.index(str(b.CLI))+1] != 'run':
            return
        cfg = config_values(cmd[cmd.index('-c')+1])
        trace = Path(cfg['CF_BATCH_TASK_TRACE'])
        with trace.open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        for row in rows:
            if row['tag'] not in self.fail_tags:
                continue
            peak = Path(cfg['CF_BATCH_NFCORE_OUT'])/'bwa/merged_library/macs3/narrow_peak'/(row['tag']+'_peaks.narrowPeak')
            peak.unlink()
            if self.missing_output:
                continue
            row.update(status='FAILED', exit='1')
            error = Path(row['workdir'])/'.command.err'
            b.write(error, 'MACS3 needs at least 100 paired peaks at + and - strand to build the model, but can only find 40!\n'
                    if self.model_error else 'Invalid BAM fixture\n')
        if self.incomplete_trace:
            rows = rows[:1]
        b.tsv(trace, b.TASK_FIELDS, [[r[k] for k in b.TASK_FIELDS] for r in rows])
        if self.global_error:
            b.fail('Unrelated nf-core process failed; native exit was nonzero')


@unittest.skipUnless(os.name == 'posix', 'Linux/WSL Bash required')
class IsolationTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.BatchTests(methodName='test_shared_control_and_independent_ids')
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root

    def run_partial(self, resume=False, until='candidates', **options):
        runner = PartialRunner(self.fixture.plan(), resume, until)
        for key, value in options.items():
            setattr(runner, key, value)
        rc = runner.run()
        report = json.loads((self.root/'results/summary/batch_summary.json').read_text())
        return runner, rc, report['analyses']

    def test_model_failure_does_not_block_same_group_and_resume_reuses_success(self):
        first, rc, rows = self.run_partial()
        self.assertEqual(rc, 1)
        self.assertEqual(rows['TF_A']['status'], 'model_failed')
        self.assertEqual(rows['TF_A']['candidates'], 'NA')
        self.assertEqual(rows['TF_B']['status'], 'complete')
        self.assertEqual(first.records['nfcore/se_k50']['status'], 'partial')
        self.assertNotIn(('annotate', 'TF_A'), first.calls)
        genes = (self.root/'results/summary/candidate_genes.tsv').read_text()
        self.assertNotIn('TF_A', genes)
        self.assertIn('TF_B', genes)
        again, rc, _ = self.run_partial(resume=True)
        self.assertEqual(rc, 1)
        self.assertEqual(again.calls, [('run', 'reference')])
        fixed, rc, rows = self.run_partial(resume=True, fail_tags=set())
        self.assertEqual(rc, 0, fixed.errors)
        self.assertEqual(rows['TF_A']['status'], 'complete')
        self.assertNotIn(('candidates', 'TF_B'), fixed.calls)
        reused, rc, _ = self.run_partial(resume=True, fail_tags=set())
        self.assertEqual(rc, 0)
        self.assertEqual(reused.calls, [])

    def test_all_failed_tasks_are_reported_without_fake_zero_peak_results(self):
        _, rc, rows = self.run_partial(fail_tags={'TF_A_REP1', 'TF_B_REP1'})
        self.assertEqual(rc, 1)
        self.assertTrue(all(r['status'] == 'model_failed' and r['peak_count'] == 'NA' for r in rows.values()))

    def test_unrelated_failure_and_missing_trace_do_not_accept_partial_group(self):
        for option in ('global_error', 'incomplete_trace'):
            with self.subTest(option=option):
                runner, rc, rows = self.run_partial(**{option: True})
                self.assertEqual(rc, 1)
                self.assertTrue(all(r['status'] == 'nfcore_failed' for r in rows.values()))
                self.assertFalse(any(s in ('annotate', 'motif', 'candidates') for s, _ in runner.calls))

    def test_other_peak_errors_and_missing_success_outputs_remain_failures(self):
        _, rc, rows = self.run_partial(model_error=False)
        self.assertEqual((rc, rows['TF_A']['status']), (1, 'peak_calling_failed'))
        _, rc, rows = self.run_partial(missing_output=True)
        self.assertEqual((rc, rows['TF_A']['status']), (1, 'nfcore_incomplete'))
        self.assertEqual(rows['TF_B']['status'], 'complete')

    def test_until_nfcore_still_reports_failed_samples_and_nonzero_exit(self):
        runner, rc, rows = self.run_partial(until='nfcore')
        self.assertEqual(rc, 1)
        self.assertEqual(rows['TF_A']['status'], 'model_failed')
        self.assertEqual(rows['TF_B']['status'], 'stopped_after_nfcore')
        self.assertFalse(any(s == 'annotate' for s, _ in runner.calls))

    def test_one_failed_replicate_invalidates_the_whole_analysis(self):
        extra = dict(self.fixture.rows[0], replicate='2', fastq_1='TF_A_rep2.fastq.gz')
        self.fixture.fastq(extra['fastq_1'])
        self.fixture.rows.append(extra); self.fixture.save_rows()
        runner, rc, rows = self.run_partial(fail_tags={'TF_A_REP2'})
        self.assertEqual((rc, rows['TF_A']['status']), (1, 'model_failed'))
        self.assertNotIn(('annotate', 'TF_A'), runner.calls)
        self.assertEqual(rows['TF_B']['status'], 'complete')


if __name__ == '__main__':
    unittest.main()
