import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

SUPPORT = Path(__file__).resolve().parents[1]/'lib/steps/support'
sys.path.insert(0, str(SUPPORT))
from annotation_policy import norm_feat, target_features, location_rule, validate_annotation, check_policy
from promoter_audit import audit, promoter, IntervalIndex, read_genes, read_query


class AnnotationPolicyTests(unittest.TestCase):
    def test_exon_and_intron_rank_not_gene_name(self):
        cases = {
            'Exon (TX1/GENE.v6.1, exon 1 of 12)': '1st Exon',
            'Exon (TX1/exon1_gene, exon 10 of 12)': 'Other Exon',
            'Intron (TXminus/G-, intron 1 of 3)': '1st Intron',
            'Intron (TXminus/G-, intron 2 of 3)': 'Other Intron',
            '1st Exon': '1st Exon', 'Other Exon': 'Other Exon',
            'Promoter (1-2kb)': 'Promoter', "5' UTR": "5' UTR", "3' UTR": "3' UTR",
            'Distal Intergenic': 'Distal Intergenic', 'unexpected': 'Unknown'}
        for raw, expected in cases.items(): self.assertEqual(norm_feat(raw), expected)
        features = target_features("Promoter,5' UTR,1st Exon")
        self.assertTrue(location_rule(norm_feat(next(iter(cases))), -2500, features)[0])
        self.assertFalse(location_rule('Other Exon', -500, features)[0])
        self.assertTrue(location_rule('Other Exon', -2000, features, 2000)[0])
        self.assertFalse(location_rule('Other Exon', -2001, features, 2000)[0])
        with self.assertRaises(ValueError): target_features('Promotor')

    def test_annotation_priority_and_configuration_change_guard(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'peak_annotation_full.tsv'
            header = 'geneId\tpeak_id\tannotation\tdistanceToTSS\n'
            p.write_text(header+'G\tp\tPromoter (1-2kb)\t-2000\nG\tq\tExon (T/G, exon 1 of 2)\t100\n')
            validate_annotation(p, 2000, 0)
            with p.open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual(rows[1]['feature_class'], '1st Exon')
            check_policy(p, 2000, 0)
            with self.assertRaisesRegex(ValueError,'differs'): check_policy(p, 1500, 0)
            p.write_text(p.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'changed'): check_policy(p, 2000, 0)
            p.write_text(header+'G\tp\tIntron (T/G, intron 1 of 2)\t-1500\n')
            with self.assertRaisesRegex(ValueError,'priority/window'): validate_annotation(p, 2000, 0)

    def test_strand_boundaries_and_sequence_ends(self):
        self.assertEqual(promoter(('c',1000,1200,'P','+'),200,0,{}), (800,1001,1000))
        self.assertEqual(promoter(('c',2000,2200,'M','-'),200,0,{}), (2199,2400,2199))
        self.assertEqual(promoter(('c',2,10,'P','+'),200,3,{'c':100}), (0,6,2))
        self.assertEqual(promoter(('c',80,99,'M','-'),20,3,{'c':100}), (95,100,98))
        with self.assertRaises(ValueError): promoter(('c',80,101,'M','-'),20,3,{'c':100})

    def test_overlap_index_keeps_long_nested_intervals(self):
        idx = IntervalIndex([('c',0,10000,'long'),('c',40,50,'nested'),('c',900,950,'short'),('d',0,10000,'other')])
        self.assertEqual({r[3] for r in idx.overlaps('c',940,960)}, {'long','short'})
        self.assertEqual({r[3] for r in idx.overlaps('c',950,960)}, {'long'})

    def test_promoter_audit_keeps_other_assigned_unannotated_and_all_strands(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); genes = root/'genes.bed'; peaks = root/'peaks.bed'; anno = root/'anno.tsv'
            genes.write_text('c\t1000\t1200\tP\t0\t+\n'
                             'c\t850\t950\tB\t0\t+\n'
                             'c\t700\t750\tOutsidePromoter\t0\t+\n'
                             'c\t2000\t2200\tM\t0\t-\n')
            peaks.write_text('c\t799\t800\toutside\n'
                             'c\t740\t850\tspans_window_edge\n'
                             'c\t800\t801\tleft\n'
                             'c\t900\t910\tother_gene_peak\n'
                             'c\t1000\t1001\ttss\n'
                             'c\t1001\t1002\tafter\n'
                             'c\t2399\t2400\tminus_boundary\n'
                             'c\t2400\t2401\tminus_outside\n')
            anno.write_text('seqnames\tstart\tend\tpeak_id\tgeneId\tannotation\n'
                            'c\t901\t910\tother_gene_peak\tB\tExon (T/B, exon 1 of 2)\n'
                            'c\t1001\t1001\ttss\tP\tPromoter\n')
            summary = audit(genes,peaks,anno,root/'out',200,0,['P','M','Missing'],{'c':5000})
            self.assertEqual(summary['P'][6:], [4,1,1,2])
            self.assertEqual(summary['M'][6:], [1,0,0,1])
            self.assertEqual(summary['Missing'][1], 'id_not_found')
            with (root/'out/promoter_peak_overlaps.tsv').open() as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            other = next(r for r in rows if r['peak_id']=='other_gene_peak')
            self.assertEqual(other['assignment_relation'],'other')
            self.assertEqual(json.loads(other['assigned_gene_ids']),['B'])
            self.assertEqual(json.loads(other['overlapped_gene_body_ids']),['B'])
            self.assertEqual(json.loads(other['promoter_overlap_gene_body_ids']),['B'])
            edge = next(r for r in rows if r['peak_id']=='spans_window_edge')
            self.assertEqual(json.loads(edge['overlapped_gene_body_ids']),['OutsidePromoter'])
            self.assertEqual(json.loads(edge['promoter_overlap_gene_body_ids']),[])
            self.assertEqual(anno.read_text().count('other_gene_peak'),1)

    def test_query_ids_remain_literal_and_duplicate_coordinates_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'query'
            p.write_text('Cre03.g146207_4532\nCre03.g146207_4532.v6.1\n4532\nCre03.g146207_4532\n')
            self.assertEqual(read_query(p), ['Cre03.g146207_4532','Cre03.g146207_4532.v6.1','4532'])
            p.write_text('gene_id\tscore\nA.1\t2\nB_4532\t3\n')
            self.assertEqual(read_query(p), ['A.1','B_4532'])
            p.write_text('c\t1\t5\tG\t0\t+\nc\t10\t15\tG\t0\t+\n')
            with self.assertRaisesRegex(ValueError,'duplicate'): read_genes(p)


if __name__ == '__main__': unittest.main()
