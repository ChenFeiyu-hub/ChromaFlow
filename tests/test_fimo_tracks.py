"""Coordinate and evidence-scope checks for motif track export."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

PATH = Path(__file__).resolve().parents[1]/'lib/steps/support/fimo_tracks.py'
spec = importlib.util.spec_from_file_location('fimo_tracks', PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FimoTracksTests(unittest.TestCase):
    def test_coordinate_conventions_and_boundaries(self):
        sizes = {'chromosome_01': 1000, 'plastome': 500}
        self.assertEqual(module.coordinates('plastome', 1, 8, sizes), ('plastome', 0, 8, 'genomic'))
        self.assertEqual(module.coordinates('chromosome_01', 1000, 996, sizes), ('chromosome_01', 995, 1000, 'genomic'))
        self.assertEqual(module.coordinates('peak::chromosome_01:100-300', 1, 8, sizes),
                         ('chromosome_01', 100, 108, 'window_relative'))
        self.assertEqual(module.coordinates('chromosome_01:100-300(+)', 193, 200, sizes),
                         ('chromosome_01', 292, 300, 'window_relative'))
        for args in [('chromosome_01', 0, 8), ('other', 1, 8),
                     ('chromosome_01', 999, 1001), ('chromosome_01:100-300', 199, 205),
                     ('chromosome_01:100-300(-)', 1, 8)]:
            with self.assertRaises(ValueError):
                module.coordinates(*args, sizes)

    def test_filter_deduplicate_and_do_not_publish_malformed_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); f = p/'fimo.tsv'; bed = p/'hits.bed'; audit = p/'hits.tsv'
            header = 'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\n'
            f.write_text(header+'M1\tDof\tchr1\t101\t108\t-\t0.0001\n'
                         +'M1\tDof\tchr1\t101\t108\t-\t0.00001\n'
                         +'M2\t.\tchr1\t111\t115\t+\t0.5\n')
            hits = module.export([f], bed, audit, sizes={'chr1': 1000})
            self.assertEqual(len(hits), 1)
            self.assertEqual(bed.read_text(), 'chr1\t100\t108\tDof (p=1e-05)\t50\t-\n')
            before = bed.read_bytes()
            f.write_text(header+'M1\tDof\tchr1\t1\t8\t+\tNaN\n')
            with self.assertRaises(ValueError): module.export([f], bed, audit)
            self.assertEqual(bed.read_bytes(), before)
            f.write_text(header)
            module.export([f], bed, audit)
            self.assertEqual(bed.read_text(), '')

    def test_input_windows_come_from_extracted_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); f = p/'input.fa'; bed = p/'windows.bed'
            f.write_text('>p2::chr1:200-400\nACGT\n>p1::chr1:100-300\nAAAA\n')
            module.scan_windows(f, bed, {'chr1': 1000})
            self.assertEqual(bed.read_text(), 'chr1\t100\t300\nchr1\t200\t400\n')
            f.write_text('>unknown-header\nACGT\n')
            with self.assertRaises(ValueError): module.scan_windows(f, bed, {})


if __name__ == '__main__':
    unittest.main()
