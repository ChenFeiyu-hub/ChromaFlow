"""Run the actual Bash entry point with stub executables; no analysis/downloads."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

from test_offline_bundle import make_bundle, b


ROOT = Path(__file__).resolve().parents[1]
BASH = shutil.which('bash')


def shell_path(path):
    path = Path(path).resolve().as_posix()
    if os.name == 'nt':
        return '/' + path[0].lower() + path[2:]
    return path


@unittest.skipUnless(BASH, 'Bash is required for CLI integration tests')
class RunTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='.test-run-', dir=ROOT)
        self.project = Path(self.tmp.name).resolve() / 'project with spaces'
        self.project.mkdir()
        # Pre-create fixture directories in Python: the Windows sandbox's MSYS
        # mkdir cannot traverse protected user-profile ancestors to create them.
        for directory in ('inputs/reference', 'inputs/raw', 'outputs/logs', 'outputs/work',
                          'outputs/.singularity_cache', 'outputs/.chromaflow'):
            (self.project / directory).mkdir(parents=True, exist_ok=True)
        tools = self.project / 'stubs'
        tools.mkdir()
        self.tools = tools
        nextflow = tools / 'nextflow'
        nextflow.write_text('''#!/usr/bin/env bash
if [[ "$1" == "-version" ]]; then echo 'nextflow version 25.10.4'; exit 0; fi
printf '%s\\0' "$@" > "$TEST_RECORD"
if [[ "$1" == "-log" ]]; then touch "$2"; fi
exit "${TEST_EXIT:-0}"
''', encoding='utf-8')
        nextflow.chmod(0o755)
        apptainer = tools / 'apptainer'
        apptainer.write_text('#!/usr/bin/env bash\necho "apptainer version 1.3.0"\n', encoding='utf-8')
        apptainer.chmod(0o755)
        for name in ('genome.fa', 'genes.gtf'):
            (self.project / name).write_text('fixture\n', encoding='utf-8')
        sheets = self.project / 'outputs/samplesheets'
        sheets.mkdir(parents=True)
        (sheets / 'samplesheet_dapseq.csv').write_text('sample,fastq_1,fastq_2\n', encoding='utf-8')
        project = shell_path(self.project)
        (self.project / 'chromaflow.config').write_text(
            'PROJECT_DIR=' + shlex.quote(project) + '\nSPECIES_SHORT=Os\n'
            'GENOME_COMBINED=' + shlex.quote(project + '/genome.fa') + '\n'
            'GENOME_GTF=' + shlex.quote(project + '/genes.gtf') + '\n'
            'MACS_GSIZE=100\nREAD_LENGTH=150\nNFCORE_PROFILE=singularity\n', encoding='utf-8')
        self.record = self.project / 'args.bin'

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, extra=(), exit_code=0):
        env = os.environ.copy()
        for key in list(env):
            if key.startswith(('CHROMAFLOW_', 'NFCORE_', 'NXF_')):
                env.pop(key)
        env.update(TEST_RECORD=shell_path(self.record), TEST_EXIT=str(exit_code))
        args = ['bash', shell_path(ROOT / 'bin/chromaflow'), 'run', '--type', 'dapseq',
                '--profile', 'apptainer', '--yes'] + list(extra)
        command = 'export PATH=' + shlex.quote(shell_path(self.tools)) + ':/usr/bin:/bin:"$PATH"\n' + shlex.join(args)
        return subprocess.run([BASH, '-c', command], env=env, cwd=self.project,
                              capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=30)

    def recorded_args(self):
        return self.record.read_bytes().decode('utf-8').rstrip('\0').split('\0')

    def test_profile_override_and_explicit_debug_log(self):
        result = self.cli(['--skip-spp'])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        args = self.recorded_args()
        self.assertEqual(args[args.index('-profile') + 1], 'apptainer')
        self.assertEqual(args[args.index('-r') + 1], '2.1.0')
        self.assertIn('nf-core/chipseq', args)
        self.assertIn('--skip_spp', args)
        self.assertIn('-resume', args)
        self.assertTrue(list((self.project / 'outputs/logs').glob('nextflow_dapseq_*.log')))

    def test_failure_exit_status_and_log_location(self):
        result = self.cli(exit_code=17)
        self.assertEqual(result.returncode, 17, result.stdout + result.stderr)
        self.assertIn('17', result.stdout + result.stderr)
        self.assertIn('nextflow_dapseq_', result.stdout + result.stderr)

    def test_dry_run_never_executes_analysis(self):
        result = self.cli(['--dry-run'])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse(self.record.exists())
        self.assertIn('[DRY-RUN]', result.stdout)

    def test_offline_without_bundle_fails_before_nextflow(self):
        result = self.cli(['--offline'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--bundle', result.stdout + result.stderr)
        self.assertFalse(self.record.exists())

    def test_offline_uses_local_runtime_without_revision_lookup(self):
        bundle = self.project / 'software bundle'
        manifest = make_bundle(bundle)
        shutil.copy2(self.tools / 'nextflow', bundle / 'bin/nextflow')
        manifest['files'] = b.inventory(bundle, ['payload', 'bin', 'nextflow/framework', 'nextflow/plugins'])
        b.write_json(bundle / b.MANIFEST, manifest)
        python = self.tools / 'python3'
        python.write_text('#!/usr/bin/env bash\nexec ' + shlex.quote(Path(sys.executable).as_posix()) + ' "$@"\n', encoding='utf-8')
        python.chmod(0o755)
        result = self.cli(['--offline', '--bundle', shell_path(bundle)])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        args = self.recorded_args()
        self.assertNotIn('-r', args)
        self.assertNotIn('nf-core/chipseq', args)
        self.assertIn('-C', args)
        self.assertIn((bundle / 'payload/2_1_0').as_posix(), args)
        envfile = self.project / 'outputs/.chromaflow/offline-runtime.env'
        self.assertIn('NXF_OFFLINE=true', envfile.read_text(encoding='utf-8'))


if __name__ == '__main__':
    unittest.main()
