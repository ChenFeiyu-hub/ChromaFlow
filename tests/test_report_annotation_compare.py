"""Finished annotate tables include the Python policy column after R exits."""
import copy
import csv
import io
from pathlib import Path
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib/steps/support'))
import annotation_policy as p


class ReportAnnotationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.saved = self.root/'saved.tsv'
        self.reconstructed = self.root/'raw.tsv'
        self.fields = ['seqnames','start','end','peak_id','geneId','annotation','distanceToTSS','peak_score','gene_consistent']
        self.rows = [dict(zip(self.fields, row)) for row in [
            ['Chr1','100','200','P1','gene_1','Promoter (<=1kb)','-10','2.123456','TRUE'],
            ['Chr2','500','700','P2','gene "two"','Intron (T2/gene_two, intron 1 of 3)','3000','4.5678','FALSE']]]
        self.write(self.reconstructed, self.fields, self.rows)
        self.saved.write_bytes(self.reconstructed.read_bytes())
        with redirect_stdout(io.StringIO()):
            p.validate_annotation(self.saved, 1500, 0)

    def write(self, path, fields, rows):
        with path.open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)

    def test_same_annotation_after_real_policy_step_is_accepted_read_only(self):
        self.assertNotEqual(self.saved.read_bytes(), self.reconstructed.read_bytes())
        before = {f.name: f.read_bytes() for f in self.root.iterdir()}
        p.compare_annotation(self.saved, self.reconstructed, 1500, 0)
        self.assertEqual({f.name: f.read_bytes() for f in self.root.iterdir()}, before)

    def test_every_scientific_column_is_still_checked(self):
        for key, value in [('seqnames','ChrX'),('start','101'),('end','201'),('peak_id','changed'),
                           ('geneId','changed_gene'),('distanceToTSS','-11'),('peak_score','2.123457'),
                           ('gene_consistent','FALSE'),('annotation','Promoter (1-2kb)')]:
            with self.subTest(column=key):
                rows = copy.deepcopy(self.rows); rows[0][key] = value
                self.write(self.reconstructed, self.fields, rows)
                with self.assertRaisesRegex(ValueError, 'column '+key):
                    p.compare_annotation(self.saved, self.reconstructed, 1500, 0)

    def test_saved_feature_class_is_not_ignored(self):
        with self.saved.open() as stream:
            reader = csv.DictReader(stream, delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
        rows[0]['feature_class'] = 'Distal Intergenic'
        self.write(self.saved, fields, rows)
        with self.assertRaisesRegex(ValueError, 'column feature_class'):
            p.compare_annotation(self.saved, self.reconstructed, 1500, 0)

    def test_row_count_order_and_extra_columns_are_errors(self):
        self.write(self.reconstructed, self.fields, self.rows[:1])
        with self.assertRaisesRegex(ValueError, 'row count'):
            p.compare_annotation(self.saved, self.reconstructed, 1500, 0)
        self.write(self.reconstructed, self.fields, self.rows[::-1])
        with self.assertRaisesRegex(ValueError, 'content mismatch'):
            p.compare_annotation(self.saved, self.reconstructed, 1500, 0)
        self.write(self.reconstructed, self.fields+['unexpected'], self.rows)
        with self.assertRaisesRegex(ValueError, 'column mismatch'):
            p.compare_annotation(self.saved, self.reconstructed, 1500, 0)

    def test_changed_promoter_policy_is_not_accepted(self):
        with self.assertRaisesRegex(ValueError, 'Promoter priority/window mismatch'):
            p.compare_annotation(self.saved, self.reconstructed, 0, 0)


if __name__ == '__main__':
    unittest.main()
