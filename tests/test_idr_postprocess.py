"""Real bedtools processing with deterministic IDR output (not an IDR fit test).

Set CHROMAFLOW_TEST_BEDTOOLS_IMAGE to an existing local bedtools SIF if the tool
is not installed. This test never downloads containers.
"""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
IMAGE=os.environ.get('CHROMAFLOW_TEST_BEDTOOLS_IMAGE', '')
BEDTOOLS=shutil.which('bedtools')
APPTAINER=shutil.which('apptainer')


@unittest.skipUnless(os.name=='posix' and (BEDTOOLS or (APPTAINER and Path(IMAGE).is_file())), 'Real bedtools or a local bedtools SIF required')
class IdrPostprocessTests(unittest.TestCase):
    def run_fixture(self, broad=False, groups=('TF','TF'), failure=0, empty=False, intervals=None):
        tmp=tempfile.TemporaryDirectory(prefix='chromaflow-idr-')
        self.addCleanup(tmp.cleanup)
        p=Path(tmp.name)
        tools=p/'tools';tools.mkdir()
        peak_type='broadPeak' if broad else 'narrowPeak'
        pipeline='atacseq' if broad else 'chipseq'
        peaks=p/f'results/{pipeline}/bwa/merged_library/macs3'
        (peaks/'narrow_peak/consensus').mkdir(parents=True)
        directory=peaks/('broad_peak' if broad else 'narrow_peak');directory.mkdir(exist_ok=True)
        intervals = intervals or [('chr1',100,200)]
        for i,group in enumerate(groups,1):
            (directory/f'{group}_REP{i}_peaks.{peak_type}').write_text(
                ''.join(f'{chrom}\t{start}\t{end}\tpeak{i}_{j}\t100\t.\t{i*5}\t900\t800'+('' if broad else '\t50')+'\n'
                        for j,(chrom,start,end) in enumerate(intervals)))
        stub=tools/'idr'
        stub.write_text('''#!/usr/bin/env python3
import os,sys
from pathlib import Path
print('idr-tool-stdout'); print('idr-tool-stderr',file=sys.stderr)
if int(os.environ['TEST_FAIL']): sys.exit(int(os.environ['TEST_FAIL']))
args=sys.argv
out=Path(args[args.index('--output-file')+1])
rows=[]
for line in Path(args[args.index('--samples')+1]).read_text().splitlines():
    row=line.split('\\t')[:3]+['merged','1000','.','-1','900','-1']
    if args[args.index('--input-file-type')+1]=='narrowPeak': row+=['50']
    row+=['0.1','3'] # Deliberately distinct local and global IDR.
    rows.append('\\t'.join(row)+'\\n')
out.write_text('' if os.environ['TEST_EMPTY']=='1' else ''.join(rows))
''')
        stub.chmod(0o755)
        wrapper=tools/'bedtools'
        command=shlex.quote(BEDTOOLS) if BEDTOOLS else shlex.join([APPTAINER,'exec','--cleanenv',IMAGE,'bedtools'])
        wrapper.write_text('#!/usr/bin/env bash\nexec '+command+' "$@"\n');wrapper.chmod(0o755)
        cfg=p/'config'
        cfg.write_text('PROJECT_DIR='+shlex.quote(str(p))+'\nOUTPUT_DIR="${PROJECT_DIR}"\nSPECIES_SHORT=Test\n')
        env={k:v for k,v in os.environ.items() if not k.startswith(('CHROMAFLOW_','NFCORE_','NXF_'))}
        env.update(PATH=str(tools)+':/usr/bin:/bin',TEST_FAIL=str(failure),TEST_EMPTY='1' if empty else '0')
        result=subprocess.run(['bash',str(ROOT/'bin/chromaflow'),'idr','--type',pipeline,'-c',str(cfg)],
                              env=env,text=True,capture_output=True,timeout=60)
        return result,p/f'results/downstream/{pipeline}/idr',p

    def test_narrow_and_broad_global_idr_and_macs_signal(self):
        for broad in (False,True):
            r,out,_=self.run_fixture(broad=broad)
            self.assertEqual(r.returncode,0,r.stdout+r.stderr)
            row=(out/'idr_peaks_threshold_0.01.bed').read_text().strip().split('\t')
            self.assertEqual(len(row),8)
            self.assertEqual(float(row[4]),10) # MACS col 7, not p-value 900.
            self.assertEqual(float(row[6]),3) # Correct global IDR column for both formats.
            self.assertEqual(int(row[7]),1)
            self.assertEqual((out/'idr.status').read_text(),'complete\n')

    def test_fit_failure_and_empty_output_are_not_success(self):
        for args in ({'failure':21},{'empty':True}):
            r,out,p=self.run_fixture(**args)
            self.assertNotEqual(r.returncode,0,r.stdout+r.stderr)
            self.assertEqual((out/'idr.status').read_text(),'failed\n')
            log=next((p/'logs/downstream/chipseq').glob('*.log')).read_text()
            self.assertIn('idr-tool-stderr',log)
            self.assertIn('exit_code=',log)

    def test_order_after_bedtools_annotate_across_bins_and_contigs(self):
        intervals=[('chloroplast',100,150), ('chromosome_01',100,200),
                   ('chromosome_01',16300,16480), ('chromosome_01',20000,20100),
                   ('chromosome_04',743255,743773), ('contig_33',800,900)]
        r,out,_=self.run_fixture(intervals=intervals)
        self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        rows=[line.split('\t') for line in (out/'idr_peaks_threshold_0.01.bed').read_text().splitlines()]
        self.assertEqual([(r[0],int(r[1]),int(r[2])) for r in rows],sorted(intervals))
        self.assertTrue(all(float(r[4])==10 and int(r[7])==1 for r in rows))

    def test_distinct_conditions_are_not_paired(self):
        r,out,_=self.run_fixture(groups=('Control','Treatment'))
        self.assertNotEqual(r.returncode,0)
        self.assertIn('one biological condition',r.stdout)
        self.assertNotIn('idr-tool-stdout',r.stdout)
        self.assertEqual((out/'idr.status').read_text(),'failed\n')
