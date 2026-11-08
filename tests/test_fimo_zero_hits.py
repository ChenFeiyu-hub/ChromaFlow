"""Regression for MEME/FIMO 5.5.7 footer-only zero-hit output."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / 'lib/steps/support'
sys.path.insert(0, str(SUPPORT))
from fimo_tracks import fimo_rows, export
from motif_selection import Evidence, table

HEADER = 'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n'
FOOTER = ('\n# FIMO (Find Individual Motif Occurrences): Version 5.5.7 compiled on Nov 12 2025 at 00:25:07\n'
          '# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.\n'
          '# fimo --verbosity 1 --oc fimo_out_5 --motif 3-ATAAAG streme.xml peak_sequences.fa\n')


class FimoZeroHitsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.empty = self.root / 'fimo_out_5/fimo.tsv'
        self.empty.parent.mkdir()
        self.empty.write_text(FOOTER, encoding='utf-8')
        self.empty.with_suffix('.gff').write_text('##gff-version 3\n', encoding='utf-8')

    def test_reproduces_old_parser_failure_and_accepts_completed_empty_result(self):
        with self.assertRaisesRegex(ValueError, 'Missing required result columns'):
            table(self.empty, ('motif_id', 'sequence_name'))
        sources = []
        self.assertEqual(fimo_rows(self.empty, sources), [])
        self.assertEqual(sources, [self.empty.with_suffix('.gff').resolve()])

    def test_malformed_and_incomplete_outputs_still_fail(self):
        for bad in ('', '\n\n', '# no hits\n', FOOTER.rsplit('# fimo', 1)[0],
                    FOOTER + '# motif_id\twrong\n', 'bad\theader\n',
                    HEADER + 'M\tAlias\tchr1\t10\n', HEADER + 'a\t' * 10 + 'extra\n'):
            with self.subTest(text=bad):
                self.empty.write_text(bad, encoding='utf-8')
                with self.assertRaises(ValueError):
                    fimo_rows(self.empty)

    def test_companion_hits_cannot_be_mistaken_for_zero_hits(self):
        self.empty.with_suffix('.gff').write_text('##gff-version 3\nchr1\tfimo\tmotif\t1\t6\t.\t+\t.\tID=M\n')
        with self.assertRaisesRegex(ValueError, 'conflicts with its GFF'):
            fimo_rows(self.empty)

    def test_header_only_and_populated_formats_still_work(self):
        self.empty.write_text(HEADER + FOOTER)
        self.assertEqual(fimo_rows(self.empty), [])
        self.empty.write_text('\ufeff' + HEADER + 'M\tAlias\tchr1\t10\t15\t+\t8\t1e-5\t.01\tAAAAAA\n' + FOOTER,
                              encoding='utf-8')
        self.assertEqual(len(fimo_rows(self.empty)), 1)

    def prepare_models(self, positive=True):
        files = [self.empty]
        motifs = [dict(id='3-ATAAAG', alt='STREME-3', prog='streme', evalue='1e-6',
                       fimo_gff='fimo_out_5/fimo.gff')]
        if positive:
            full = self.root/'fimo_out_1/fimo.tsv'
            full.parent.mkdir()
            full.write_text(HEADER + 'ACGTAC\tMEME-1\tChr1\t110\t115\t+\t8\t1e-5\t.01\tACGTAC\n')
            files.insert(0, full)
            motifs.append(dict(id='ACGTAC', alt='MEME-1', prog='meme', evalue='1e-8',
                               fimo_gff='fimo_out_1/fimo.gff'))
        (self.root/'meme-chip.html').write_text('<script>var data = ' + json.dumps(dict(cmd=['meme-chip'], motifs=motifs)) + ';</script>')
        return files

    def test_zero_hit_model_remains_catalogued_and_other_hits_are_preserved(self):
        e = Evidence(self.prepare_models(), self.root)
        self.assertEqual(len(e.hits), 1)
        self.assertEqual(e.hits[0]['motif_label'], 'MEME-1')
        self.assertIn(('3-ATAAAG', 'STREME-3'), e.selected_models)
        self.assertEqual(e.models[('3-ATAAAG', 'STREME-3')]['scan_status'], 'scanned')
        empty = Evidence([self.empty], self.root, motifs='STREME-3', qmax=.05)
        self.assertEqual(empty.hits, [])

    def test_tracks_export_accepts_empty_and_preserves_positive_hits(self):
        files = self.prepare_models()
        bed, audit = self.root/'hits.bed', self.root/'audit.tsv'
        self.assertEqual(len(export(files, bed, audit)), 1)
        self.assertIn('MEME-1', bed.read_text())
        self.assertEqual(export([self.empty], bed, audit), {})
        self.assertEqual(bed.read_text(), '')

    def candidate_run(self, positive):
        files = self.prepare_models(positive)
        annotation = self.root/'annotation.tsv'
        annotation.write_text('seqnames\tstart\tend\tpeak_id\tgeneId\tannotation\tdistanceToTSS\tpeak_score\n'
                              'Chr1\t101\t200\tpeak1\tGENE_1\tPromoter\t-100\t5\n')
        peaks = self.root/'peaks.bed'
        peaks.write_text('Chr1\t100\t200\tpeak1\t5\t+\n')
        genes = self.root/'genes.bed'
        genes.write_text('Chr1\t200\t1000\tGENE_1\t0\t+\n')
        out = self.root/'candidates'; out.mkdir()
        env = {k: v for k, v in os.environ.items() if not k.startswith('CFSEL_')}
        env.update(CFSEL_ROOT=str(self.root), CFSEL_PEAKS=str(peaks), CFSEL_UP='1500',
                   CFSEL_DOWN='0', CFSEL_PIPELINE='dapseq', CHROMAFLOW_PEAKSET='consensus', PYTHONUTF8='1')
        result = subprocess.run([sys.executable, str(SUPPORT/'candidates.py'), str(annotation), str(out),
                                 '1e-4', str(self.root/'centrimo.tsv'), str(genes),
                                 "Promoter,5' UTR,1st Exon", '0', *map(str, files)],
                                env=env, text=True, encoding='utf-8', capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(json.loads((out/'selection.json').read_text())['status'], 'complete')
        with (out/'motif_summary.tsv').open() as handle:
            rows = {row['motif']: row for row in csv.DictReader(handle, delimiter='\t')}
        self.assertEqual(rows['STREME-3']['n_hits'], '0')
        self.assertEqual((out/'selected_genes.txt').read_text(), 'GENE_1\n' if positive else '')
        return out

    def test_candidates_integrate_mixed_hit_and_zero_hit_models(self):
        self.candidate_run(True)

    def test_candidates_all_zero_hits_complete_with_empty_gene_list(self):
        self.candidate_run(False)


if __name__ == '__main__':
    unittest.main()
