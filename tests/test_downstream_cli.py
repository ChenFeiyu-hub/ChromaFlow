"""Exercise logging, failure propagation and peak-set isolation via the real CLI."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.name == 'posix' and shutil.which('bash'), 'Linux Bash required')
class DownstreamTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='chromaflow downstream ')
        self.addCleanup(self.tmp.cleanup)
        self.project = Path(self.tmp.name)
        self.tools = self.project / 'tools'
        self.tools.mkdir()
        self.config = self.project / 'project.config'
        self.config.write_text('PROJECT_DIR=' + shlex.quote(str(self.project)) + '\n'
                               'OUTPUT_DIR="${PROJECT_DIR}"\nSPECIES_SHORT=Test\n')
        self.reference = self.project / 'inputs/reference'
        self.reference.mkdir(parents=True)
        (self.reference / 'Test_genes.gff3').write_text('chr1\ttest\texon\t1\t9\t.\t+\t.\tID=e\n')
        self.down = self.project / 'results/downstream/chipseq'
        self.peaks = self.project / 'results/chipseq/bwa/merged_library/macs3/narrow_peak'
        self.peaks.mkdir(parents=True)
        (self.peaks / 'consensus').mkdir()
        (self.peaks / 'consensus/Test.consensus_peaks.bed').write_text('chr1\t500\t600\tcons\t0\t+\n')
        self.stub('Rscript', '''echo annotation-stdout; echo annotation-stderr >&2
if [[ "${TEST_EXIT:-0}" != 0 ]]; then exit "$TEST_EXIT"; fi
python3 - "$@" <<'PY'
import sys
from pathlib import Path
out=Path(sys.argv[4]); out.mkdir(parents=True,exist_ok=True)
bed=Path(sys.argv[3]).read_text().splitlines()[0].split('\\t')
(out/'peak_annotation_full.tsv').write_text('seqnames\\tstart\\tend\\tpeak_id\\tgeneId\\tannotation\\tdistanceToTSS\\tpeak_score\\n'
    +f'{bed[0]}\\t{int(bed[1])+1}\\t{bed[2]}\\t{bed[3]}\\tGENE1\\tPromoter\\t0\\t5\\n')
(out/'genes.bed').write_text('chr1\\t100\\t500\\tGENE1\\t0\\t+\\n')
PY''')

    def stub(self, name, body):
        path = self.tools / name
        path.write_text('#!/usr/bin/env bash\n' + body + '\n')
        path.chmod(0o755)

    def cli(self, *args, exit_code=0, stdin_text=''):
        env = {k: v for k, v in os.environ.items() if not k.startswith(('CHROMAFLOW_', 'NFCORE_', 'NXF_'))}
        env.update(PATH=str(self.tools) + ':/usr/bin:/bin', TEST_EXIT=str(exit_code))
        return subprocess.run(['bash', str(ROOT / 'bin/chromaflow'), *args, '-c', str(self.config)],
                              env=env, cwd=self.project, text=True, input=stdin_text, capture_output=True, timeout=30)

    def idr_fixture(self):
        path = self.down / 'idr'
        path.mkdir(parents=True, exist_ok=True)
        (path / 'idr_peaks_threshold_0.01.bed').write_text('chr1\t100\t200\tidr\t8\t+\t3\t1\n')

    def logs(self):
        return sorted((self.project / 'logs/downstream/chipseq').glob('*.log'))

    def test_success_logged_without_empty_feature_directories(self):
        self.idr_fixture()
        result = self.cli('annotate')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        log = self.logs()[0].read_text()
        for text in ('annotation-stdout', 'annotation-stderr', '[VERSION] 1.2.7', '[COMMAND]', 'exit_code=0'):
            self.assertIn(text, log)
        self.assertFalse((self.down / 'signal').exists())
        self.assertFalse((self.down / 'motif').exists())
        self.assertEqual(len(list((self.project / 'logs').glob('*.log'))), 0)
        result = self.cli('annotate')
        self.assertEqual(result.returncode, 0)
        self.assertEqual(len(self.logs()), 2)

    def test_tool_failure_is_logged_and_returned(self):
        self.idr_fixture()
        result = self.cli('annotate', exit_code=17)
        self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
        log = self.logs()[0].read_text()
        self.assertIn('annotation-stderr', log)
        self.assertIn('exit_code=17', log)
        self.assertNotIn('Annotation completed', log)

    def test_early_input_failure_is_logged(self):
        result = self.cli('annotate')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Requested IDR peak set is missing', self.logs()[0].read_text())
        self.assertNotIn('annotation-stdout', result.stdout)

    def test_failed_idr_cannot_be_consumed(self):
        self.idr_fixture()
        (self.down / 'idr/idr.status').write_text('failed\n')
        result = self.cli('annotate')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('latest IDR run did not complete', result.stdout)

    def test_failed_idr_cannot_export_go_from_old_annotation(self):
        self.idr_fixture()
        (self.down / 'idr/idr.status').write_text('failed\n')
        annotation = self.down/'annotation'; annotation.mkdir()
        (annotation/'promoter_target_genes.txt').write_text('OLD_RESULT\n')
        result = self.cli('go-prep')
        self.assertNotEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('latest IDR run did not complete', result.stdout)
        self.assertFalse((self.down/'go/target_promoter.raw.txt').exists())

    def test_candidates_respect_r_missing_gene_and_consistency_flag(self):
        self.idr_fixture()
        annotation = self.down/'annotation'; annotation.mkdir()
        (annotation/'peak_annotation_full.tsv').write_text(
            'seqnames\tstart\tend\tpeak_id\tgeneId\tannotation\tdistanceToTSS\tpeak_score\tfeature_geneId\tgene_consistent\n'
            'chr1\t101\t200\tp1\tG1\tPromoter\t-100\t5\tNA\tTRUE\n'
            'chr1\t301\t400\tp2\tG2\tExon\t-50\t10\tOTHER\tFALSE\n'
            'chr1\t501\t600\tp3\tG3\tExon\t-50\t10\tNA\tFALSE\n')
        motif = self.down/'motif/meme-chip_output/fimo_out_1'; motif.mkdir(parents=True)
        (motif/'fimo.tsv').write_text('motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\n'
                                    'M1\tDof\tchr1\t120\t128\t+\t0.00001\n')
        r = self.cli('candidates', '--features', 'Promoter')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        self.assertIn('NOTE: 2 peaks overlap a different gene', r.stdout)
        self.assertEqual((self.down/'candidates/high_confidence_genes.txt').read_text(), 'G1\n')

    def test_motif_tracks_include_all_requested_genes_and_record_missing(self):
        with self.config.open('a') as stream: stream.write('PROMOTER_UP=100\nPROMOTER_DOWN=0\n')
        annotation = self.down/'annotation_consensus'; annotation.mkdir(parents=True)
        (annotation/'genes.bed').write_text('chr1\t100\t500\tGENE1\t0\t+\tGeneOne\n'
                                           'chr1\t700\t1000\tGENE2\t0\t-\tGeneTwo\n')
        (self.reference/'Test_genes.gtf').write_text('fixture\n')
        bw = self.project/'results/chipseq/bwa/merged_library/bigwig'; bw.mkdir(parents=True)
        (bw/'test.bigWig').write_text('fixture')
        self.stub('pyGenomeTracks', 'echo render "$@"')
        motif = self.down/'motif_consensus/meme-chip_output/fimo_out_1'; motif.mkdir(parents=True)
        (motif/'fimo.tsv').write_text('motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\n'
                                    'M1\tDof\tchr1\t520\t528\t+\t0.00001\n')
        (self.down/'motif_consensus/peak_sequences.fa').write_text('>p1::chr1:450-650\nACGT\n')
        query = self.project/'genes.txt'; query.write_text('GENE1\nGENE2\nUNKNOWN\n')
        r = self.cli('tracks', '--peakset', 'consensus', '--gene-list', str(query), '--motifs', '--flank', '150')
        self.assertEqual(r.returncode, 0, r.stdout+r.stderr)
        out = self.down/'tracks_consensus'
        ini = (out/'tracks.ini').read_text()
        self.assertIn('[motif_hits]', ini)
        self.assertIn('[motif_scanned_windows]', ini)
        self.assertTrue((out/'motif_hits.bed').read_text().startswith('chr1\t519\t528\tDof'))
        report = (out/'plot_status.tsv').read_text()
        self.assertIn('GENE1\tplotted\tchr1:0-650', report)
        self.assertIn('GENE2\tplotted\tchr1:550-1150', report)
        self.assertIn('UNKNOWN\tmissing_coordinates', report)
        self.assertFalse((self.down/'tracks').exists())

    def test_consensus_annotation_uses_separate_directory_and_log(self):
        result = self.cli('annotate', '--peakset', 'consensus')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((self.down / 'annotation_consensus').is_dir())
        self.assertFalse((self.down / 'annotation').exists())
        self.assertIn('annotate_consensus', self.logs()[0].name)

    def test_signal_and_go_preserve_peakset_selection(self):
        self.idr_fixture()
        self.stub('computeMatrix', 'echo matrix-stdout; echo matrix-stderr >&2; exit "${TEST_EXIT:-0}"')
        self.stub('plotHeatmap', 'exit 0')
        self.stub('plotProfile', 'exit 0')
        bw = self.project / 'results/chipseq/bwa/merged_library/bigwig'
        bw.mkdir(parents=True)
        (bw / 'test.bigWig').write_text('fixture')
        for peakset, suffix, start in [('idr', '', '100'), ('consensus', '_consensus', '500')]:
            result = self.cli('signal', '--peakset', peakset)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((self.down / ('signal'+suffix) / 'peaks.bed').read_text().split('\t')[1], start)
            annotation = self.down / ('annotation'+suffix)
            annotation.mkdir(exist_ok=True)
            (annotation / 'promoter_target_genes.txt').write_text(peakset + '\n')
            # No background is valid: it must warn, not call an undefined function.
            result = self.cli('go-prep', '--peakset', peakset)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((self.down / ('go'+suffix) / 'target_promoter.raw.txt').read_text(), peakset+'\n')
        result = self.cli('signal', exit_code=23)
        self.assertEqual(result.returncode, 23, result.stdout + result.stderr)
        self.assertIn('matrix-stderr', result.stdout)

    def test_tracks_renderer_failure_and_success_are_logged(self):
        self.idr_fixture()
        annotation=self.down/'annotation'
        annotation.mkdir()
        (annotation/'genes.bed').write_text('chr1\t100\t500\tGENE1\t0\t+\tGeneOne\n')
        (self.reference/'Test_genes.gtf').write_text('fixture\n')
        bw=self.project/'results/chipseq/bwa/merged_library/bigwig'
        bw.mkdir(parents=True)
        (bw/'test.bigWig').write_text('fixture')
        self.stub('pyGenomeTracks','echo renderer-stdout; echo renderer-stderr >&2; exit "${TEST_EXIT:-0}"')
        for code in (0,25):
            r=self.cli('tracks','--genes','GENE1',exit_code=code)
            if code:
                self.assertNotEqual(r.returncode,0,r.stdout+r.stderr)
                self.assertNotIn('[DONE]',r.stdout)
            else:
                self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            self.assertIn('renderer-stderr',r.stdout)

    def test_spacing_selects_consensus_motifs_and_logs_bedtools_error(self):
        for i in (1,2):
            directory=self.down/f'motif_consensus/meme-chip_output/fimo_out_{i}'
            directory.mkdir(parents=True)
            (directory/'fimo.tsv').write_text('motif_id\tmotif_alt\tsequence\tstart\tstop\tstrand\tscore\tpvalue\tqvalue\tsequence\n'
                                             f'M{i}\tM{i}\tchr1\t520\t525\t+\t5\t0.0001\t0.01\tACGTA\n')
        self.stub('bedtools','echo spacing-bedtools-error >&2; exit 26')
        r=self.cli('spacing','--peakset','consensus')
        self.assertEqual(r.returncode,26,r.stdout+r.stderr)
        self.assertIn('spacing-bedtools-error',r.stdout)
        self.assertTrue((self.down/'motif_consensus/motif_spacing').is_dir())
        self.assertFalse((self.down/'motif').exists())

    def test_lookup_without_motifs_retains_location_only_mode(self):
        self.idr_fixture()
        annotation=self.down/'annotation';annotation.mkdir()
        (annotation/'peak_annotation_full.tsv').write_text('seqnames\tstart\tend\tpeak_id\tgeneId\tannotation\tdistanceToTSS\tpeak_score\n'
                                                          'chr1\t101\t200\tidr\tGENE1\tPromoter\t-100\t5\n')
        (annotation/'genes.bed').write_text('chr1\t100\t500\tGENE1\t0\t+\tGENE1\n')
        gene_list=self.project/'query.txt';gene_list.write_text('GENE1\n')
        r=self.cli('lookup','--gene-list',str(gene_list),'--yes')
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        self.assertIn('location-only',r.stdout)
        self.assertIn('in_target',(self.down/'lookup/query_binding.tsv').read_text())

    def lookup_fixture(self):
        self.idr_fixture()
        annotation=self.down/'annotation'; annotation.mkdir()
        (annotation/'peak_annotation_full.tsv').write_text('seqnames\tstart\tend\tpeak_id\tgeneId\tannotation\tdistanceToTSS\tpeak_score\n'
            'chr1\t101\t200\tidr\tCre03.g146207_4532\tPromoter\t-100\t5\n')
        (annotation/'genes.bed').write_text('chr1\t200\t500\tCre03.g146207_4532\t0\t+\n'
                                          'chr1\t600\t900\tG.v6.1\t0\t+\n')
        query=self.project/'query.txt'; query.write_text('Cre03.g146207_4532\n4532\nCre03.g146207\nG\nG.v6.1\n')
        return query

    def test_lookup_y_and_n_confirmation_and_no_automatic_id_changes(self):
        query=self.lookup_fixture()
        r=self.cli('lookup','--gene-list',str(query),stdin_text='n\n')
        self.assertEqual(r.returncode,2,r.stdout+r.stderr)
        self.assertIn('[CANCELLED]',r.stdout)
        self.assertFalse((self.down/'lookup/query_binding.tsv').exists())
        self.assertEqual((self.down/'lookup/lookup.status').read_text(),'cancelled\n')
        r=self.cli('lookup','--gene-list',str(query),stdin_text='y\n')
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        self.assertIn('MATCH: Cre03.g146207_4532 == Cre03.g146207_4532',r.stdout)
        import csv
        rows={r['gene_id']:r for r in csv.DictReader((self.down/'lookup/query_binding.tsv').open(),delimiter='\t')}
        for gid in ('4532','Cre03.g146207','G'):
            self.assertEqual(rows[gid]['id_status'],'not_found')
            self.assertEqual(rows[gid]['has_peak'],'no')
        self.assertEqual(rows['G.v6.1']['id_status'],'exact_match')
        self.assertEqual(rows['G.v6.1']['has_peak'],'no')

    def test_lookup_eof_cancels_and_yes_does_not_allow_zero_matches(self):
        query=self.lookup_fixture()
        r=self.cli('lookup','--gene-list',str(query))
        self.assertEqual(r.returncode,2,r.stdout+r.stderr)
        query.write_text('4532\nG\n')
        r=self.cli('lookup','--gene-list',str(query),'--yes')
        self.assertNotEqual(r.returncode,0)
        self.assertIn('No exact ID matches',r.stdout)
        self.assertFalse((self.down/'lookup/query_binding.tsv').exists())

    def test_candidates_and_lookup_share_first_exon_selection(self):
        self.idr_fixture()
        with self.config.open('a') as w: w.write('PROMOTER_DOWN=0\nTARGET_FEATURES="Promoter,5\' UTR,1st Exon"\n')
        annotation=self.down/'annotation'; annotation.mkdir()
        (annotation/'peak_annotation_full.tsv').write_text('seqnames\tstart\tend\tpeak_id\tgeneId\tannotation\tdistanceToTSS\tpeak_score\n'
            'chr1\t101\t200\tidr\tG1\tExon (T1/G1, exon 1 of 3)\t10\t5\n'
            'chr1\t501\t600\tp2\tG2\tExon (T2/G2, exon 2 of 3)\t10\t8\n')
        (annotation/'genes.bed').write_text('chr1\t90\t400\tG1\t0\t+\nchr1\t490\t800\tG2\t0\t+\n')
        (self.down/'idr/idr_peaks_threshold_0.01.bed').write_text('chr1\t100\t200\tidr\t5\t+\nchr1\t500\t600\tp2\t8\t+\n')
        motif=self.down/'motif/meme-chip_output/fimo_out_1'; motif.mkdir(parents=True)
        (motif/'fimo.tsv').write_text('motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\n'
            'M1\tM1\tchr1\t110\t118\t+\t0.00001\nM1\tM1\tchr1\t510\t518\t+\t0.00001\n')
        query=self.project/'query'; query.write_text('G1\nG2\n')
        for features, expected in [(None,'G1'), ('Other Exon','G2')]:
            if features:
                with self.config.open('a') as w: w.write('TARGET_FEATURES="Other Exon"\n')
            for cmd,args in [('candidates',()),('lookup',('--gene-list',str(query),'--yes'))]:
                r=self.cli(cmd,*args)
                self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            self.assertEqual((self.down/'candidates/high_confidence_genes.txt').read_text(),expected+'\n')
            self.assertEqual((self.down/'lookup/candidate_genes.txt').read_text(),expected+'\n')

    def test_tracks_show_other_assigned_peaks_and_no_name_aliases(self):
        query=self.lookup_fixture()
        (self.reference/'Test_genes.gtf').write_text('fixture\n')
        bw=self.project/'results/chipseq/bwa/merged_library/bigwig'; bw.mkdir(parents=True)
        (bw/'sample.bw').write_text('fixture')
        self.stub('pyGenomeTracks','echo rendered "$@"')
        query.write_text('G.v6.1\nG\n')
        r=self.cli('tracks','--gene-list',str(query),'--flank','10')
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        import csv
        rows=list(csv.DictReader((self.down/'tracks/promoter_peak_overlaps.tsv').open(),delimiter='\t'))
        self.assertTrue(any(r['gene_id']=='G.v6.1' and r['assignment_relation']=='other' for r in rows))
        self.assertIn('G\tmissing_coordinates',(self.down/'tracks/plot_status.tsv').read_text())
        self.assertIn('idr (sig=8)',(self.down/'tracks/peaks.bed').read_text())


    def test_shared_filters_select_non_best_signal_peak_and_drive_tracks(self):
        query=self.lookup_fixture()
        query.write_text('Cre03.g146207_4532\nG.v6.1\n')
        anno=self.down/'annotation/peak_annotation_full.tsv'
        with anno.open('a') as w:
            w.write('chr1\t301\t350\tp2\tCre03.g146207_4532\tPromoter\t-50\t2\n')
        motif=self.down/'motif/meme-chip_output/fimo_out_1';motif.mkdir(parents=True)
        (motif/'fimo.tsv').write_text('motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\tq-value\n'
            'M1\tMEME-1\tchr1\t120\t128\t+\t0.00001\t0.01\n'
            'M4\tMEME-4\tchr1\t310\t318\t+\t0.00001\t0.02\n')
        (self.reference/'Test_genes.gtf').write_text('fixture\n')
        bw=self.project/'results/chipseq/bwa/merged_library/bigwig';bw.mkdir(parents=True)
        (bw/'sample.bw').write_text('fixture')
        self.stub('pyGenomeTracks','echo rendered "$@"')
        r=self.cli('candidates','--motif','MEME-4','--qval','0.05')
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        r=self.cli('lookup','--gene-list',str(query),'--yes','--level','candidate','--motif','MEME-4','--qval','0.05','--plot')
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        out=self.down/'lookup'
        self.assertEqual((out/'selected_genes.txt').read_text(),'Cre03.g146207_4532\n')
        self.assertEqual((out/'selected_genes.txt').read_text(),(self.down/'candidates/selected_genes.txt').read_text())
        self.assertIn('\tp2\t',(out/'query_binding.tsv').read_text())
        self.assertNotIn('MEME-1',(out/'motif_hits.bed').read_text())
        self.assertIn('MEME-4',(out/'motif_hits.bed').read_text())
        self.assertNotIn('G.v6.1',(out/'tracks/plot_status.tsv').read_text())
        self.assertIn('idr (sig=8)',(out/'tracks/peaks.bed').read_text())
        with anno.open('a') as w: w.write('\n')
        r=self.cli('tracks','--selection',str(out))
        self.assertNotEqual(r.returncode,0,r.stdout+r.stderr)
        self.assertIn('Selection input changed',r.stdout)

    def test_project_peakset_default_and_explicit_override(self):
        self.idr_fixture()
        with self.config.open('a') as w: w.write('DEFAULT_PEAKSET=consensus\n')
        r=self.cli('annotate')
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        self.assertTrue((self.down/'annotation_consensus/peak_annotation_full.tsv').exists())
        r=self.cli('annotate','--peakset','idr')
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        self.assertTrue((self.down/'annotation/peak_annotation_full.tsv').exists())


if __name__ == '__main__':
    unittest.main()
