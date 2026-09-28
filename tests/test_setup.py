"""Verify installation plans and propagation of package-solver failures."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(os.name == 'posix' and shutil.which('bash'), 'Linux Bash required')
class SetupTests(unittest.TestCase):
    def plan(self, *args):
        return subprocess.run(['bash', str(ROOT/'bin/chromaflow'), 'setup', '--dry-run', *args],
                              capture_output=True, text=True, timeout=15)

    def test_modes_have_distinct_runtime_dependencies(self):
        for mode, expected, absent in [('offline', 'ChromaFlow_runtime', ['nf-core=', 'nextflow=']),
                                      ('prepare', 'nf-core=3.3.2', ['nextflow=']),
                                      ('conda', 'nextflow=25.10.4', ['nf-core='])]:
            r = self.plan('--mode', mode, '--components', 'runtime')
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn(expected, r.stdout)
            self.assertIn('openjdk=17', r.stdout)
            for token in absent + ['ChromaFlow_R', 'meme=', 'deeptools=', 'idr=']:
                self.assertNotIn(token, r.stdout)

    def test_downstream_only_and_optional_components(self):
        r = self.plan('--components', 'downstream', '--skip-r', '--with-agat')
        self.assertEqual(r.returncode, 0, r.stderr)
        for name in ['ChromaFlow_ref', 'ChromaFlow_idr', 'ChromaFlow_motif', 'ChromaFlow_downstream', 'ChromaFlow_agat']:
            self.assertIn(name, r.stdout)
        for name in ['create -y -n nfcore_env', 'create -y -n ChromaFlow_runtime', 'install_chipseeker.sh', 'homer', 'csvtk']:
            self.assertNotIn(name, r.stdout)

    def test_help_and_invalid_arguments(self):
        r = subprocess.run(['bash', str(ROOT/'bin/chromaflow'), 'setup', '--help'], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0)
        self.assertIn('--components', r.stdout)
        self.assertNotEqual(self.plan('--mode', 'typo').returncode, 0)

    def test_solver_failure_is_not_reported_as_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            solver = Path(tmp)/'mamba'
            solver.write_text('#!/usr/bin/env bash\nif [[ "$1" == env ]]; then exit 0; fi\necho solver-failed >&2; exit 19\n')
            solver.chmod(0o755)
            env = dict(os.environ, PATH=tmp+':/usr/bin:/bin')
            env.pop('CHROMAFLOW_DRY_RUN', None)
            r = subprocess.run(['bash', str(ROOT/'setup_environments.sh'), '--components', 'runtime', '--solver', 'mamba'],
                               env=env, capture_output=True, text=True)
            self.assertEqual(r.returncode, 19, r.stdout+r.stderr)
            self.assertNotIn('[DONE]', r.stdout)
