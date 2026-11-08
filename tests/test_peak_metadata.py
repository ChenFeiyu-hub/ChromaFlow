"""Signal/IDR semantics and narrowly verified cache-preserving metadata repair."""
import copy
import csv
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'lib'))
from steps.support.peak_metadata import non_idr_replacement, lookup_idr_values, sha_bytes, table_bytes
import batch as b


class MetadataTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.peaks = self.root/'peaks.bed'
        self.peaks.write_text('Chr2\t100\t120\tp\t3.21\t.\t3.21\t7\t6\t10\n')
        self.raw = b'seqnames\tstart\tend\tpeak_id\tgeneId\tV7\tpeak_globalIDR\nChr2\t101\t120\tp\tG_123\t3.21\t3.21\n'
        self.policy = json.dumps(dict(annotation_sha256=sha_bytes(self.raw), promoter_up=1500)).encode()

    def test_only_idr_label_is_corrected_and_other_cells_are_exact(self):
        result = non_idr_replacement(self.raw, self.policy, self.peaks)
        fields, rows = table_bytes(self.raw)
        new_fields, new_rows = table_bytes(result[0])
        self.assertEqual(fields, new_fields)
        rows[0]['peak_globalIDR'] = 'NA'
        self.assertEqual(new_rows, rows)
        self.assertEqual(json.loads(result[1])['annotation_sha256'], sha_bytes(result[0]))
        self.assertIsNone(non_idr_replacement(result[0], result[1], self.peaks))
        for raw in (self.raw.replace(b'\t101\t', b'\t102\t'),
                    self.raw.replace(b'3.21\n', b'4.21\n')):
            with self.assertRaises(ValueError): non_idr_replacement(raw, self.policy, self.peaks)

    def test_lookup_requires_explicit_idr_and_matching_bed8(self):
        rows = table_bytes(self.raw)[1]
        for kind in ('non_idr', 'unknown', ''):
            self.assertEqual(lookup_idr_values(rows, kind, self.peaks), ['NA'])
        with self.assertRaises(ValueError): lookup_idr_values(rows, 'idr', self.peaks)
        self.peaks.write_text('Chr2\t100\t120\tp\t9\t+\t3.21\t2\n')
        self.assertEqual(lookup_idr_values(rows, 'idr', self.peaks), ['3.21'])
        rows[0]['start'] = '102'
        with self.assertRaises(ValueError): lookup_idr_values(rows, 'idr', self.peaks)

    def runner(self):
        runner = b.Runner(dict(config=dict(OUTPUT_DIR=str(self.root))), True, 'candidates')
        directory = runner.results/'analyses/A/annotation'
        b.write(directory/'peak_annotation_full.tsv', self.raw.decode())
        b.write(directory/'annotation_policy.json', self.policy.decode())
        b.write(directory/'genes.bed', 'unchanged genes')
        b.save(directory.parent/'inputs.json', dict(idr=False, peakset_kind='single_replicate_no_idr'))
        b.write(runner.state/'analyses/A/selected_peaks.bed', self.peaks.read_text())
        runner.records['A/annotate'] = dict(status='complete', signature='old', outputs=b.artifact_files([directory]))
        return runner

    def test_repair_preserves_dependency_and_rejects_modified_receipt_or_backups(self):
        runner = self.runner()
        old = copy.deepcopy(runner.records['A/annotate'])
        dependency = b.digest([old['signature'], old['outputs']])
        runner.repair_annotation_metadata('A')
        rec = runner.records['A/annotate']
        self.assertTrue(b.valid_record(rec, 'old'))
        self.assertEqual(runner.annotation_dependency('A', rec), dependency)
        runner.repair_annotation_metadata('A')
        backup = Path(rec['non_idr_metadata_correction']['backup'])
        (backup/'peak_annotation_full.tsv').write_bytes(b'changed')
        with self.assertRaises(b.BatchError): runner.annotation_dependency('A', rec)

    def test_wrong_provenance_and_failed_state_commit_cannot_corrupt_files(self):
        runner = self.runner(); old = copy.deepcopy(runner.records)
        p = runner.results/'analyses/A/inputs.json'
        b.save(p, dict(idr=True, peakset_kind='single_replicate_no_idr'))
        with self.assertRaises(b.BatchError): runner.repair_annotation_metadata('A')
        b.save(p, dict(idr=False, peakset_kind='single_replicate_no_idr'))
        save = b.save
        def fail_state(path, data):
            if path == runner.status_path: raise OSError('simulated state write failure')
            return save(path, data)
        with patch.object(b, 'save', side_effect=fail_state):
            with self.assertRaises(OSError): runner.repair_annotation_metadata('A')
        self.assertEqual(runner.records, old)
        self.assertTrue(b.valid_record(old['A/annotate'], 'old'))


@unittest.skipUnless(os.name == 'posix', 'Linux fixture CLI required')
class MetadataResumeTests(unittest.TestCase):
    def test_report_repair_preserves_motif_and_candidates_resume(self):
        from test_batch import BatchTests, FixtureRunner
        fixture = BatchTests(methodName='test_shared_control_and_independent_ids')
        fixture.setUp(); self.addCleanup(fixture.doCleanups)
        class Legacy(FixtureRunner):
            def execute(self, cmd, log, cwd):
                super().execute(cmd, log, cwd)
                if self.calls[-1][0] == 'annotate':
                    aid = self.calls[-1][1]
                    table = self.results/'analyses'/aid/'annotation/peak_annotation_full.tsv'
                    fields, rows = table_bytes(table.read_bytes())
                    fields += ['V7', 'peak_globalIDR']
                    for row in rows: row['V7'] = row['peak_globalIDR'] = row['peak_score']
                    b.tsv(table, fields, [[row[f] for f in fields] for row in rows])
                    policy = table.with_name('annotation_policy.json')
                    data = b.read_json(policy); data['annotation_sha256'] = b.sha(table); b.save(policy, data)
        plan = fixture.plan()
        before = Legacy(plan, False, 'candidates')
        self.assertEqual(before.run(), 0)
        old = copy.deepcopy(before.records)
        for aid in plan['analyses']: before.repair_annotation_metadata(aid)
        for _ in range(2):
            after = FixtureRunner(plan, True, 'candidates')
            self.assertEqual(after.run(), 0)
            self.assertEqual(after.calls, [])
            for aid in plan['analyses']:
                for step in ('motif', 'candidates'):
                    record = after.records[f'{aid}/{step}']
                    if step == 'motif':
                        self.assertEqual(record, old[f'{aid}/{step}'])
                    else:
                        self.assertEqual(record['signature'], old[f'{aid}/{step}']['signature'])
                        changed = {name for name, value in record['outputs'].items()
                                   if old[f'{aid}/{step}']['outputs'][name] != value}
                        selection = after.results/'analyses'/aid/'candidates/selection.json'
                        self.assertEqual(changed, {str(selection)})
                        from steps.support import motif_selection
                        peaks = after.state/'analyses'/aid/'selected_peaks.bed'
                        self.assertEqual(motif_selection.validate(str(selection), 'consensus', 'dapseq', str(peaks),
                            int(plan['config']['PROMOTER_UP']), int(plan['config']['PROMOTER_DOWN'])), selection.parent)
        # Scientific changes still invalidate cache normally.
        (fixture.root/'results/analyses/TF_A/annotation/genes.bed').write_text('changed')
        changed = FixtureRunner(plan, True, 'candidates')
        self.assertEqual(changed.run(), 0)
        self.assertIn(('annotate', 'TF_A'), changed.calls)


if __name__ == '__main__': unittest.main()
