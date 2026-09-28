import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'lib/steps/support'))
from motif_selection import Evidence, validate
from track_layout import prepare


class MotifSelectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fimo = self.root/'fimo_out_3/fimo.tsv'
        self.fimo.parent.mkdir()
        self.header = 'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\tq-value\n'
        self.fimo.write_text(self.header+'ACGT\tMEME-4\tchr1:100-200\t10\t15\t+\t0.000079\t0.2\n')
        self.report = {'cmd':['meme-chip','-centrimo-local'], 'motifs':[
            {'id':'ACGT','alt':'MEME-4','prog':'meme','evalue':'1.8e-9','pwm':[[0.25]*4], 'fimo_gff':'fimo_out_3/fimo.gff'},
            {'id':'AAAA','alt':'MEME-5','prog':'meme','evalue':'1e-7','pwm':[[1,0,0,0]]}]}
        self.save_report()
        (self.root/'centrimo_out').mkdir()
        (self.root/'centrimo_out/centrimo.tsv').write_text('motif_id\tmotif_alt_id\tE-value\tbin_location\nACGT\tMEME-4\t1.8\t0\n')

    def save_report(self):
        (self.root/'meme-chip.html').write_text('<script>var data = '+json.dumps(self.report)+';</script>')

    def evidence(self, **kwargs):
        return Evidence([self.fimo], self.root, **kwargs)

    def test_model_selector_is_not_fimo_directory_number(self):
        e = self.evidence(motifs='MEME-4')
        self.assertEqual(len(e.hits),1)
        self.assertEqual((e.hits[0]['start'],e.hits[0]['end']),(109,115))
        self.assertEqual(e.hits[0]['q_value'],0.2)
        self.assertEqual(len(e.models[('ACGT','MEME-4')]['pwm_sha256']),64)
        with self.assertRaisesRegex(ValueError,'no available FIMO scan'):
            self.evidence(motifs='MEME-5')
        with self.assertRaisesRegex(ValueError,'unknown or ambiguous'):
            self.evidence(motifs='MEME-3')

    def test_p_q_model_e_and_centrimo_are_distinct_filters(self):
        self.assertEqual(len(self.evidence(pmax=1e-5).hits),0)
        self.assertEqual(len(self.evidence(qmax=.05).hits),0)
        self.assertEqual(len(self.evidence(emax=1e-8).hits),1)
        self.assertEqual(len(self.evidence(emax=1e-10).hits),0)
        self.assertEqual(len(self.evidence(centrimo_max=.05).hits),0)
        self.assertEqual(len(self.evidence(centrimo_max=2).hits),1)

    def test_unknown_q_and_nan_p_fail_without_false_success(self):
        self.fimo.write_text(self.header+'ACGT\tMEME-4\tchr1\t10\t15\t+\t1e-5\t\n')
        with self.assertRaisesRegex(ValueError,'q-value missing'):
            self.evidence(qmax=.05)
        self.fimo.write_text(self.header+'ACGT\tMEME-4\tchr1\t10\t15\t+\tnan\t.01\n')
        with self.assertRaisesRegex(ValueError,'Invalid FIMO p-value'):
            self.evidence()

    def test_requested_filters_without_scans_and_bad_header_fail(self):
        with self.assertRaisesRegex(ValueError,'require available FIMO'):
            Evidence([], self.root, qmax=.05)
        self.fimo.write_text('bad\theader\n')
        with self.assertRaisesRegex(ValueError,'Missing required result columns'):
            self.evidence()

    def test_centrimo_local_threshold_and_missing_are_not_central_claims(self):
        m=self.evidence().models[('ACGT','MEME-4')]
        self.assertEqual(m['positional_enrichment'],'above_threshold')
        self.assertEqual(m['centrally_enriched'],'not_assessed')
        (self.root/'centrimo_out/centrimo.tsv').write_text('motif_id\tmotif_alt_id\tE-value\tbin_location\nACGT\tMEME-4\t1e-9\t0\n')
        m=self.evidence().models[('ACGT','MEME-4')]
        self.assertEqual(m['positional_enrichment'],'passes_threshold')
        self.assertEqual(m['centrally_enriched'],'not_assessed')
        self.report['cmd']=['meme-chip']; self.save_report()
        self.assertEqual(self.evidence().models[('ACGT','MEME-4')]['centrally_enriched'],'yes')
        (self.root/'centrimo_out/centrimo.tsv').unlink()
        self.assertEqual(self.evidence().models[('ACGT','MEME-4')]['positional_enrichment'],'not_reported')

    def test_selection_rejects_changed_inputs_outputs_and_other_branch(self):
        out=self.root/'out'; out.mkdir()
        anno=self.root/'annotation.tsv'; anno.write_text('fixture\n')
        peaks=self.root/'peaks.bed'; peaks.write_text('chr1\t100\t200\tp\n')
        e=self.evidence()
        with patch.dict(os.environ, {'CHROMAFLOW_PEAKSET':'consensus','CFSEL_PIPELINE':'chipseq'}):
            e.save(out,['G'],[],anno,peaks,dict(promoter_up=3000,promoter_down=0))
        args=(out,'consensus','chipseq',peaks,3000,0)
        self.assertEqual(validate(*args),out)
        with self.assertRaisesRegex(ValueError,'different pipeline or peak set'):
            validate(out,'idr','chipseq',peaks,3000,0)
        with self.assertRaisesRegex(ValueError,'promoter parameters'):
            validate(out,'consensus','chipseq',peaks,1500,0)
        (out/'selected_genes.txt').write_text('OTHER\n')
        with self.assertRaisesRegex(ValueError,'output changed'):
            validate(*args)
        (out/'selected_genes.txt').write_text('G\n')
        self.fimo.write_text(self.fimo.read_text()+'\n')
        with self.assertRaisesRegex(ValueError,'input changed'):
            validate(*args)

    def test_adaptive_height_keeps_all_neighbor_transcripts(self):
        ini=self.root/'tracks.ini'; ini.write_text('[genes]\nheight = 2\nmerge_transcripts = false\n[motif_hits]\nheight = 1.2\n')
        gtf=self.root/'genes.gtf'; gtf.write_text(''.join(f'chr1\ts\texon\t101\t200\t.\t+\t.\tgene_id "G"; transcript_id "T{i}";\n' for i in range(16)))
        genes=self.root/'genes.bed'; genes.write_text('chr1\t100\t200\tG\t0\t+\n')
        query=self.root/'query';query.write_text('G\n')
        with patch.dict(os.environ, {'TRACK_GENE_HEIGHT':'auto','TRACK_MOTIF_HEIGHT':'auto'}):
            prepare(ini,gtf,genes,query,3000,'','',self.root)
        result=(self.root/'.layout/G.ini').read_text()
        self.assertIn('height = 8.4',result)
        self.assertIn('merge_transcripts = false',result)
