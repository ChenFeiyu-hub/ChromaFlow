"""Batch contracts using miniature FASTQs and controlled workflow producers.

Candidate integration runs the real Bash CLI and Python evidence implementation;
expensive reference, nf-core, R and MEME producers use explicit fixtures.
"""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'lib'))
import batch as b


def config_values(path):
    values = {}
    for s in Path(path).read_text().splitlines():
        if not s or s.startswith('#'):
            continue
        k, _, v = s.partition('=')
        values[k] = shlex.split(v)[0] if v != "''" else ''
    return values


class FixtureRunner(b.Runner):
    def __init__(self, *args):
        super().__init__(*args)
        self.calls = []
        self.fail_analysis = None
        self.npeaks = 12

    def execute(self, cmd, log, cwd):
        step = cmd[cmd.index(str(b.CLI)) + 1]
        cfg = config_values(cmd[cmd.index('-c') + 1])
        aid = Path(cfg.get('CF_BATCH_ANALYSIS_DIR', 'reference')).name
        self.calls.append((step, aid))
        if step == 'candidates':
            if aid == self.fail_analysis:
                b.fail('Intentional candidate failure')
            return super().execute(cmd, log, cwd)
        if step == 'normalize-gff':
            ref = Path(cfg['REFERENCE_DIR']); ref.mkdir(parents=True, exist_ok=True)
            out = ref/'genes.normalized.gff3'
            out.write_text('chr1\ttest\tgene\t1\t10000\t.\t+\t.\tID=G\n')
            b.write(Path(cfg['OUTPUT_DIR'])/'.chromaflow/reference.manifest', 'GENOME_GFF3='+str(out)+'\n')
        elif step == 'reference':
            ref = Path(cfg['REFERENCE_DIR'])
            paths = dict(GENOME_COMBINED=ref/'genome.fa', GENOME_GFF3=ref/'genes.normalized.gff3',
                         GENOME_GTF=ref/'genes.gtf', CHROM_SIZES=ref/'chrom.sizes')
            paths['GENOME_COMBINED'].write_text('>chr1\n'+'ACGT'*2500+'\n')
            Path(str(paths['GENOME_COMBINED'])+'.fai').write_text('chr1\t10000\t6\t10000\t10001\n')
            paths['CHROM_SIZES'].write_text('chr1\t10000\n')
            paths['GENOME_GTF'].write_text('fixture GTF\n')
            b.write(Path(cfg['OUTPUT_DIR'])/'.chromaflow/reference.manifest',
                    ''.join(f'{k}={p}\n' for k, p in paths.items()))
        elif step == 'run':
            out = Path(cfg['CF_BATCH_NFCORE_OUT'])
            base = out/'bwa/merged_library/macs3/narrow_peak'
            base.mkdir(parents=True, exist_ok=True)
            with (Path(cfg['OUTPUT_DIR'])/'samplesheets/samplesheet_dapseq.csv').open() as f:
                rows = list(csv.DictReader(f))
            for row in rows:
                if not row['antibody']:
                    continue
                sample, rep = row['sample'], row['replicate']
                offset = 100 if sample == 'TF_A' else 5000
                signal = 3 if sample == 'TF_A' else 100
                np = base/f'{sample}_REP{rep}_peaks.narrowPeak'
                np.write_text(''.join(f'chr1\t{offset+i*100}\t{offset+i*100+50}\t{sample}_{i}\t999\t.\t{signal}\t9\t8\t25\n'
                                      for i in range(self.npeaks)))
            for sample in {r['sample'] for r in rows if r['antibody']}:
                if len({r['replicate'] for r in rows if r['sample'] == sample}) > 1:
                    source = base/f'{sample}_REP1_peaks.narrowPeak'
                    b.write(base/'consensus'/sample/f'{sample}.consensus_peaks.bed',
                            ''.join('\t'.join(line.split('\t')[:4]+['0', '+'])+'\n' for line in source.read_text().splitlines()))
            # An unrelated consensus must never be discovered by a broad glob.
            b.write(base/'consensus/Other/Other.consensus_peaks.bed', 'chr1\t1\t2\twrong\t0\t+\n')
            for suffix in ('.amb', '.ann', '.bwt', '.pac', '.sa'):
                b.write(out/'genome/index/bwa'/('genome.fa'+suffix), 'bwa fixture\n')
            if cfg.get('CF_BATCH_TASK_TRACE'):
                tags = sorted({f'{r["sample"]}_REP{r["replicate"]}' for r in rows if r['antibody']})
                b.tsv(cfg['CF_BATCH_TASK_TRACE'], b.TASK_FIELDS,
                      [[i, b.MACS_PROCESS, tag, 'COMPLETED', '0', str(Path(cfg['CF_BATCH_WORK_DIR'])/tag)]
                       for i, tag in enumerate(tags, 1)])
        elif step == 'annotate':
            dest = Path(cfg['CF_BATCH_ANALYSIS_DIR'])/'annotation'
            dest.mkdir(parents=True, exist_ok=True)
            peaks = Path(cfg['CF_BATCH_PEAK_FILE']).read_text().splitlines()
            lines = ['seqnames\tstart\tend\tpeak_id\tgeneId\tannotation\tdistanceToTSS\tpeak_score']
            for p in peaks:
                c = p.split('\t')
                lines.append(f'{c[0]}\t{int(c[1])+1}\t{c[2]}\t{c[3]}\tGENE_{aid}\tPromoter\t0\t{c[4]}')
            b.write(dest/'peak_annotation_full.tsv', '\n'.join(lines)+'\n')
            b.write(dest/'genes.bed', f'chr1\t1\t9000\tGENE_{aid}\t0\t+\n')
            subprocess.run([sys.executable, str(ROOT/'lib/steps/support/annotation_policy.py'), 'validate',
                            str(dest/'peak_annotation_full.tsv'), '--up', cfg['PROMOTER_UP'], '--down', cfg['PROMOTER_DOWN']], check=True)
        elif step == 'motif':
            dest = Path(cfg['CF_BATCH_ANALYSIS_DIR'])/'motif/meme-chip_output'
            peak = Path(cfg['CF_BATCH_PEAK_FILE']).read_text().splitlines()[0].split('\t')
            start, end = peak[1:3]
            b.write(dest/'fimo_out_1/fimo.tsv',
                    'motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tp-value\tq-value\n'
                    +f'ACGT\tMEME-1\tchr1:{start}-{end}\t10\t15\t+\t1e-5\t.01\n')
            data = dict(cmd=['meme-chip', '-centrimo-local'], motifs=[dict(id='ACGT', alt='MEME-1',
                        prog='meme', evalue='1e-6', pwm=[[.25]*4], fimo_gff='fimo_out_1/fimo.gff')])
            b.write(dest/'meme-chip.html', '<script>var data = '+json.dumps(data)+';</script>')
        else:
            raise AssertionError(step)


@unittest.skipUnless(os.name == 'posix', 'Linux/WSL Bash required')
class BatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='chromaflow-batch-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.raw = self.root/'raw'; self.raw.mkdir()
        self.fa = self.root/'source.fa'; self.fa.write_text('>chr1\nACGT\n')
        self.gff = self.root/'source.gff3'; self.gff.write_text('##gff-version 3\n')
        self.cfg = self.root/'shared config.config'
        self.cfg.write_text(f'PROJECT_DIR={shlex.quote(str(self.root))}\nOUTPUT_DIR="${{PROJECT_DIR}}"\n'
                            'RAW_DATA_DIR="${PROJECT_DIR}/raw"\nSPECIES_SHORT=Test\n'
                            'GENOME_FASTA_SOURCE="${PROJECT_DIR}/source.fa"\nGENOME_GFF3_SOURCE="${PROJECT_DIR}/source.gff3"\n'
                            'BATCH_ENV_MODE=current\nMACS_GSIZE=10000\n')
        self.rows = []
        for sid, role in [('TF_A', 'target'), ('TF_B', 'target'), ('Input', 'control')]:
            fn = sid+'.fastq.gz'; self.fastq(fn)
            self.rows.append(dict(zip(b.HEADER, [sid if role == 'target' else '', sid, role, '1', 'SE', fn, '',
                 'Input' if role == 'target' else '', '1' if role == 'target' else '',
                 sid if role == 'target' else '', 'col', sid+'_RUN'])))
        self.samples = self.root/'samples.tsv'
        self.save_rows()
        self.args = argparse.Namespace(config=str(self.cfg), samples=str(self.samples), only=None)
        self.env = patch.object(b, 'preflight', return_value={s: 'fixture environment' for s in b.STEPS})
        self.env.start(); self.addCleanup(self.env.stop)

    def fastq(self, name, k=50, record='read1'):
        with gzip.open(self.raw/name, 'wt') as f:
            f.write(f'@{record}\n'+('A'*k)+'\n+\n'+('I'*k)+'\n')

    def save_rows(self):
        b.tsv(self.samples, b.HEADER, [[r.get(k, '') for k in b.HEADER] for r in self.rows])

    def plan(self):
        return b.plan_project(self.args)

    def run_fixture(self, plan=None, resume=False, fail=None, npeaks=12):
        r = FixtureRunner(plan or self.plan(), resume, 'candidates')
        r.fail_analysis, r.npeaks = fail, npeaks
        code = r.run()
        return r, code

    def test_shared_control_and_independent_ids(self):
        p = self.plan()
        self.assertEqual(list(p['groups']), ['se_k50'])
        self.assertEqual(len(p['groups']['se_k50']['rows']), 3)
        self.assertEqual(set(p['analyses']), {'TF_A', 'TF_B'})
        self.args.only = 'TF_B'
        self.assertEqual(len(self.plan()['groups']['se_k50']['rows']), 2)
        self.args.only = 'TF'
        with self.assertRaisesRegex(b.BatchError, 'Unknown --only'):
            self.plan()

    def test_environment_diagnostic_paths_do_not_invalidate_runtime_identity(self):
        runner = FixtureRunner(self.plan(), True, 'candidates')
        identity = dict(tools='java17 /same/path', bundle_manifest_sha256='fixed')
        runner.envs = dict(nfcore=dict(identity=identity, log='/tmp/run-111'))
        first = runner.signature('nfcore', ['fixed inputs'])
        runner.envs['nfcore']['log'] = '/tmp/run-222'
        self.assertEqual(first, runner.signature('nfcore', ['fixed inputs']))
        runner.envs['nfcore']['identity']['bundle_manifest_sha256'] = 'changed'
        self.assertNotEqual(first, runner.signature('nfcore', ['fixed inputs']))

    def test_shared_estimate_is_reused_within_one_invocation_without_resume(self):
        runner = FixtureRunner(self.plan(), False, 'candidates')
        calls = []
        out = self.root/'estimate.txt'
        def action():
            calls.append(1); out.write_text('1000\n')
            return [out], dict(value=1000)
        runner.stage('gsize/k50', 'same-reference-and-k', action)
        runner.stage('gsize/k50', 'same-reference-and-k', action)
        self.assertEqual(len(calls), 1)

    def test_orphan_pe_missing_control_and_no_id_normalization(self):
        self.rows[0]['layout'] = 'PE'; self.save_rows()
        with self.assertRaisesRegex(b.BatchError, 'PE requires both'):
            self.plan()
        self.rows[0]['layout'] = 'SE'; self.rows[0]['control_replicate'] = '2'; self.save_rows()
        with self.assertRaisesRegex(b.BatchError, 'Control sample/replicate not found'):
            self.plan()
        self.rows[0]['control_replicate'] = '1'; self.rows[0]['sample_id'] = ' TF_A'; self.save_rows()
        with self.assertRaisesRegex(b.BatchError, 'never normalized'):
            self.plan()

    def test_duplicate_files_mismatched_mates_and_mixed_control_policy(self):
        self.rows[1]['fastq_1'] = self.rows[0]['fastq_1']; self.save_rows()
        with self.assertRaisesRegex(b.BatchError, 'assigned more than once'):
            self.plan()
        self.rows[1]['fastq_1'] = 'TF_B.fastq.gz'
        self.rows[2]['layout'] = 'PE'; self.rows[2]['fastq_2'] = 'Input_2.fastq.gz'; self.fastq('Input_2.fastq.gz', record='wrong'); self.save_rows()
        with self.assertRaisesRegex(b.BatchError, 'BATCH_ALLOW_SE'):
            self.plan()
        with self.cfg.open('a') as f: f.write('BATCH_ALLOW_SE_WITH_PE_CONTROL=true\n')
        with self.assertRaisesRegex(b.BatchError, 'mate names/order differ'):
            self.plan()

    def test_grouping_by_k_not_control_or_longest_project_read(self):
        self.fastq('TF_B.fastq.gz', k=100)
        p = self.plan()
        self.assertEqual(set(p['groups']), {'se_k50', 'se_k100'})
        self.assertEqual(p['groups']['se_k50']['analyses'], ['TF_A'])

    def test_nonstandard_length_is_measured_and_retained_in_batch_metadata(self):
        self.cfg.write_text(self.cfg.read_text().replace('MACS_GSIZE=10000\n', ''))
        self.fastq('TF_A.fastq.gz', k=130)
        self.fastq('TF_B.fastq.gz', k=100)
        self.fastq('Input.fastq.gz', k=100)
        tools = self.root/'tools'; tools.mkdir()
        calls = self.root/'estimator-calls.txt'
        estimator = tools/'unique-kmers.py'
        estimator.write_text('#!/usr/bin/env python3\nimport os, sys\nfrom pathlib import Path\n'
            'k = int(sys.argv[sys.argv.index("-k") + 1])\n'
            'with open(os.environ["TEST_KMER_CALLS"], "a") as stream: stream.write(str(k) + "\\n")\n'
            'assert Path(sys.argv[-1]).is_file()\n'
            'print(f"Total estimated number of unique {k}-mers: {k * 1000}")\n')
        estimator.chmod(0o755)

        class EstimatorRunner(FixtureRunner):
            def execute(self, cmd, log, cwd):
                if 'batch-gsize' in cmd:
                    return b.Runner.execute(self, cmd, log, cwd)
                return super().execute(cmd, log, cwd)

        with patch.dict(os.environ, PATH=str(tools)+os.pathsep+os.environ['PATH'], TEST_KMER_CALLS=str(calls)):
            plan = self.plan()
            self.assertEqual(set(plan['groups']), {'se_k100', 'se_k130'})
            runner = EstimatorRunner(plan, False, 'nfcore')
            self.assertEqual(runner.run(), 0, runner.errors)
            for k in (100, 130):
                info = runner.records[f'gsize/k{k}']['details']
                self.assertEqual((info['k'], info['method'], info['value']),
                                 (k, f'unique-kmers(k={k})', str(k*1000)))
                cfg = config_values(self.root/f'.chromaflow/batch/runs/se_k{k}/nfcore.config')
                self.assertEqual(cfg['READ_LENGTH'], str(k))
                self.assertEqual(cfg['MACS_GSIZE'], str(k*1000))
            again = EstimatorRunner(self.plan(), True, 'nfcore')
            self.assertEqual(again.run(), 0, again.errors)
            self.assertEqual(again.calls, [])
        self.assertEqual(sorted(calls.read_text().splitlines()), ['100', '130'])

    def test_changed_fastq_rebuilds_group_and_new_source_does_not_reuse_reference(self):
        self.run_fixture()
        self.fastq('TF_A.fastq.gz', record='new_read_name')
        runner, code = self.run_fixture(resume=True)
        self.assertEqual(code, 0, runner.errors)
        self.assertIn(('run', 'reference'), runner.calls)
        self.assertNotIn(('reference', 'reference'), runner.calls)
        self.fa.write_text('>chr1\nACGTACGT\n')
        runner, code = self.run_fixture(resume=True)
        self.assertEqual(code, 0, runner.errors)
        self.assertIn(('reference', 'reference'), runner.calls)

    def test_singletons_no_native_consensus_no_cross_contamination_and_real_candidates(self):
        r, code = self.run_fixture()
        self.assertEqual(code, 0, r.errors)
        for aid, expected_signal in [('TF_A', '3'), ('TF_B', '100')]:
            root = self.root/'results/analyses'/aid
            self.assertEqual((root/'candidates/selected_genes.txt').read_text(), 'GENE_'+aid+'\n')
            peaks = (self.root/'.chromaflow/batch/analyses'/aid/'selected_peaks.bed').read_text()
            self.assertEqual(peaks.splitlines()[0].split('\t')[4], expected_signal)
            self.assertNotIn('wrong', peaks)
            self.assertEqual(json.loads((root/'inputs.json').read_text())['peakset_kind'], 'single_replicate_no_idr')
        motifs = (self.root/'results/summary/motifs.tsv').read_text()
        self.assertIn('TF_A::MEME-1', motifs); self.assertIn('TF_B::MEME-1', motifs)
        self.assertFalse((self.root/'results/downstream').exists())
        resumed, code = self.run_fixture(resume=True)
        self.assertEqual(code, 0, resumed.errors)
        self.assertEqual(resumed.calls, [])

    def test_failure_does_not_stop_other_analysis_and_resume_only_failed_step(self):
        r, code = self.run_fixture(fail='TF_A')
        self.assertEqual(code, 1)
        summary = json.loads((self.root/'results/summary/batch_summary.json').read_text())
        self.assertEqual(summary['analyses']['TF_A']['candidates'], 'NA')
        self.assertEqual(summary['analyses']['TF_B']['candidates'], 1)
        resumed, code = self.run_fixture(resume=True)
        self.assertEqual(code, 0, resumed.errors)
        self.assertEqual(resumed.calls, [('candidates', 'TF_A')])

    def test_changed_filter_and_tampered_output_invalidate_reuse(self):
        self.run_fixture()
        with self.cfg.open('a') as f: f.write('FIMO_PVAL=1e-6\n')
        r, code = self.run_fixture(resume=True)
        self.assertEqual(code, 0, r.errors)
        self.assertEqual(r.calls, [('candidates', 'TF_A'), ('candidates', 'TF_B')])
        self.assertEqual((self.root/'results/analyses/TF_A/candidates/selected_genes.txt').read_text(), '')
        output = self.root/'results/analyses/TF_A/candidates/selected_genes.txt'
        output.write_text('CONTAMINATED\n')
        r, code = self.run_fixture(resume=True)
        self.assertEqual(code, 0, r.errors)
        self.assertEqual(r.calls, [('candidates', 'TF_A')])
        self.assertEqual(output.read_text(), '')

    def test_zero_and_low_peak_counts_are_na_not_zero_candidate_evidence(self):
        r, code = self.run_fixture(npeaks=3)
        self.assertEqual(code, 0, r.errors)
        summary = json.loads((self.root/'results/summary/batch_summary.json').read_text())
        self.assertEqual(summary['analyses']['TF_A']['status'], 'insufficient_peaks')
        self.assertEqual(summary['analyses']['TF_A']['candidates'], 'NA')
        self.assertFalse(any(c[0] in ('motif', 'candidates') for c in r.calls))

    def test_multireplicate_uses_exact_consensus_and_shared_signal(self):
        extra = dict(self.rows[0], replicate='2', fastq_1='TF_A_rep2.fastq.gz')
        self.fastq(extra['fastq_1']); self.rows.append(extra); self.save_rows()
        r, code = self.run_fixture()
        self.assertEqual(code, 0, r.errors)
        info = json.loads((self.root/'results/analyses/TF_A/inputs.json').read_text())
        self.assertEqual(info['peakset_kind'], 'replicate_consensus')
        self.assertEqual(info['biological_replicates'], 2)
        self.assertEqual((self.root/'results/analyses/TF_A/candidates/selected_genes.txt').read_text(), 'GENE_TF_A\n')

    def test_init_does_not_overwrite_and_cli_dispatch(self):
        dest = self.root/'new project'
        p = subprocess.run(['bash', str(b.CLI), 'batch', 'init', str(dest)], text=True, capture_output=True)
        self.assertEqual(p.returncode, 0, p.stdout+p.stderr)
        self.assertTrue((dest/'samples.tsv').is_file())
        p = subprocess.run(['bash', str(b.CLI), 'batch', 'init', str(dest)], text=True, capture_output=True)
        self.assertNotEqual(p.returncode, 0)
        p = subprocess.run(['bash', str(b.CLI), 'help', 'batch'], text=True, capture_output=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn('candidates', p.stdout)


if __name__ == '__main__':
    unittest.main()
