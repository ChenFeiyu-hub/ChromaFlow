"""Independent server acceptance rejects false bounds and misleading metadata."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

BASE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('acceptance', BASE/'tools/acceptance.py')
a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)


class AcceptanceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.project = Path(self.tmp.name)
        self.folder = self.project/'results/analyses/A/reports/annotation'; self.folder.mkdir(parents=True)
        self.anno = self.project/'results/analyses/A/annotation/peak_annotation_full.tsv'
        self.anno.parent.mkdir(parents=True)
        self.anno.write_text('peak_id\tpeak_globalIDR\np\tNA\n')
        sizes = self.project/'chrom.sizes'; sizes.write_text('Chr1\t10000\n')
        self.records = {'reference': {'details': {'CHROM_SIZES': str(sizes)}}}
        p = self.project/'.chromaflow/batch/analyses/A/selected_peaks.bed'; p.parent.mkdir(parents=True)
        p.write_text('Chr1\t3999\t4001\tp\t2\t.\t2\t3\t4\t1\n')
        a.save(self.folder/'report.json', dict(policy='annotation-reports-v2-reference-bounds', peak_kind='non_idr'))
        self.windows = [dict(seqnames='Chr1', start=1000, end=7000, strand='-', window_id=1, tss=4000,
            chrom_length=10000, overlapping_peaks=1, status='retained', legacy_end_bound_exclusion='TRUE')]
        self.write_windows()
        a.tsv(self.folder/'tss_window_summary.tsv', ['metric', 'windows'], [
            {'metric': k, 'windows': v} for k, v in dict(retained=1, no_peak_overlap=0, outside_reference=0,
                                                      legacy_end_bound_exclusion=1).items()])
        a.tsv(self.folder/'tss_avgprof.tsv', ['distance_bp', 'mean_peak_overlap', 'windows'],
            [dict(distance_bp=i, mean_peak_overlap=1 if -1 <= i <= 0 else 0, windows=1) for i in range(-3000, 3001)])

    def write_windows(self):
        a.tsv(self.folder/'tss_windows.tsv', list(self.windows[0]), self.windows)

    def check(self):
        a.audit_tss_metadata(self.project, 'A', self.folder, self.records)

    def test_valid_recovered_window_and_na_metadata(self): self.check()

    def test_false_last_peak_boundary_is_rejected(self):
        self.windows[0]['status'] = 'outside_reference'; self.write_windows()
        with self.assertRaisesRegex(ValueError, 'retention'): self.check()

    def test_mislabeled_idr_is_rejected(self):
        self.anno.write_text('peak_id\tpeak_globalIDR\np\t2\n')
        with self.assertRaisesRegex(ValueError, 'globalIDR'): self.check()

    def test_invented_counts_and_duplicate_windows_are_rejected(self):
        self.windows[0]['overlapping_peaks'] = 2; self.write_windows()
        with self.assertRaisesRegex(ValueError, 'overlap count'): self.check()
        self.windows[0]['overlapping_peaks'] = 1
        self.windows *= 2; self.write_windows()
        with self.assertRaisesRegex(ValueError, 'Duplicate'): self.check()

    def test_wrong_denominator_and_old_report_are_rejected(self):
        path = self.folder/'tss_avgprof.tsv'
        path.write_text(path.read_text().replace('\t1\n', '\t2\n'))
        with self.assertRaisesRegex(ValueError, 'denominator'): self.check()
        a.save(self.folder/'report.json', dict(policy='annotation-reports-v1'))
        with self.assertRaisesRegex(ValueError, 'refreshed'): self.check()


if __name__ == '__main__': unittest.main()
