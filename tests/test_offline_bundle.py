import argparse
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / 'lib/offline_bundle.py'
spec = importlib.util.spec_from_file_location('bundle', MODULE)
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)


def put(root, name, content='fixture\n'):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding='utf-8')
    return path


def payload(root):
    put(root, '2_1_0/main.nf')
    put(root, '2_1_0/nextflow.config', '''params {
    custom_config_base = "${projectDir}/../configs"
}
plugins {
    id 'nf-validation@1.1.3'
}
manifest {
    version = '2.1.0'
}
''')
    put(root, 'configs/nfcore_custom.config')
    put(root, 'singularity-images/tool.img', 'container\n')


def make_bundle(root):
    payload(root / 'payload')
    put(root, 'bin/nextflow', '#!/usr/bin/env bash\nexit 0\n').chmod(0o755)
    put(root, 'nextflow/framework/25.10.4/nextflow-25.10.4-one.jar')
    put(root, 'nextflow/plugins/nf-validation-1.1.3/lib/plugin.jar')
    data = dict(format=1, complete=True, pipeline='chipseq', revision='2.1.0',
                nextflow_version='25.10.4', mode='singularity', architecture=b.platform.machine(),
                system='Linux', workflow='payload/2_1_0', configs='payload/configs',
                plugins=['nf-validation@1.1.3'], images='payload/singularity-images')
    data['files'] = b.inventory(root, ['payload', 'bin', 'nextflow/framework', 'nextflow/plugins'])
    b.write_json(root / b.MANIFEST, data)
    return data


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name).resolve()
        self.root = self.base / 'bundle with spaces'
        self.manifest = make_bundle(self.root)

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_plugin_fails_before_launch(self):
        (self.root / 'nextflow/plugins/nf-validation-1.1.3/lib/plugin.jar').unlink()
        with self.assertRaisesRegex(ValueError, 'Missing bundle file'):
            b.validate(self.root)

    def test_plugin_ready_recognizes_nextflow_installed_directory(self):
        # Nextflow 25.10.4 installs id@version as an id-version directory.
        self.assertTrue(b.plugin_ready(self.root, 'nf-validation@1.1.3'))
        self.assertFalse(b.plugin_ready(self.root, 'nf-validation@1.1.4'))

    def test_plugin_ready_rejects_wrong_directory_and_empty_jar(self):
        root = self.base / 'incomplete'
        put(root, 'nextflow/plugins/nf-validation@1.1.3/lib/plugin.jar')
        self.assertFalse(b.plugin_ready(root, 'nf-validation@1.1.3'))
        put(root, 'nextflow/plugins/nf-validation-1.1.3/lib/plugin.jar', '')
        self.assertFalse(b.plugin_ready(root, 'nf-validation@1.1.3'))

    def test_resume_after_plugin_check_failure_reuses_downloads(self):
        (self.root / b.MANIFEST).unlink()
        b.write_json(self.root / '.prepare.json', {
            key: self.manifest[key] for key in (
                'pipeline', 'revision', 'nextflow_version', 'mode', 'architecture', 'system')
        })
        commands = []

        def existing_runtime(command, **kwargs):
            commands.append(command)
            self.assertEqual(command, [str(self.root / 'bin/nextflow'), '-version'])
            return 'nextflow version 25.10.4\n'

        args = argparse.Namespace(bundle=self.root, pipeline='chipseq', revision='2.1.0',
                                  nextflow_version='25.10.4', mode='singularity', parallel_downloads=4)
        with patch.object(b.platform, 'system', return_value='Linux'), \
                patch.object(b.shutil, 'which', return_value='/bin/tool'), \
                patch.object(b, 'download', side_effect=AssertionError('Unexpected download')), \
                patch.object(b, 'run', side_effect=existing_runtime):
            b.prepare(args)
        self.assertEqual(len(commands), 1)
        b.validate(self.root, verify=True)

    def test_version_and_pipeline_mismatch(self):
        for kw in ({'revision': '2.0.0'}, {'pipeline': 'atacseq'}):
            with self.assertRaises(ValueError):
                b.validate(self.root, **kw)

    def test_atac_prepare_pins_upstream_unversioned_plugin(self):
        (self.root / b.MANIFEST).unlink()
        workflow = self.root / 'payload/2_1_2'
        (self.root / 'payload/2_1_0').rename(workflow)
        config = workflow / 'nextflow.config'
        config.write_text(config.read_text().replace("version = '2.1.0'", "version = '2.1.2'")
                          .replace("id 'nf-validation@1.1.3'", "id 'nf-validation' // upstream declaration"))
        state = {key: self.manifest[key] for key in (
            'pipeline', 'revision', 'nextflow_version', 'mode', 'architecture', 'system')}
        state.update(pipeline='atacseq', revision='2.1.2')
        b.write_json(self.root / '.prepare.json', state)
        args = argparse.Namespace(bundle=self.root, pipeline='atacseq', revision='2.1.2',
                                  nextflow_version='25.10.4', mode='singularity', parallel_downloads=4)
        with patch.object(b.platform, 'system', return_value='Linux'), \
                patch.object(b.shutil, 'which', return_value='/bin/tool'), \
                patch.object(b, 'run', return_value='nextflow version 25.10.4\n'):
            b.prepare(args)
        manifest = b.validate(self.root, verify=True, pipeline='atacseq', revision='2.1.2')
        self.assertEqual(manifest['plugins'], ['nf-validation@1.1.3'])
        self.assertIn("id 'nf-validation@1.1.3' // upstream declaration", config.read_text())

    def test_release_plugin_pin_preserves_explicit_versions(self):
        text = "plugins {\n    id 'nf-validation@1.0.0'\n}\n"
        self.assertEqual(b.pin_release_plugins(text, 'atacseq', '2.1.2'), text)

    def test_unknown_release_unversioned_plugin_still_rejected(self):
        text = "plugins {\n    id 'nf-validation'\n}\n"
        text = b.pin_release_plugins(text, 'atacseq', '9.9.9')
        workflow = self.base / 'unknown'
        put(workflow, 'nextflow.config', text)
        with self.assertRaisesRegex(ValueError, 'pinned plugin'):
            b.plugins_in(workflow)

    def test_same_size_corruption_detected_by_sha256(self):
        path = self.root / 'payload/singularity-images/tool.img'
        original = path.read_bytes()
        path.write_bytes(b'X' * len(original))
        with self.assertRaisesRegex(ValueError, 'SHA256'):
            b.validate(self.root, verify=True)

    def test_manifest_cannot_escape_bundle(self):
        self.manifest['files']['../outside'] = {'size': 1, 'sha256': 'x'}
        b.write_json(self.root / b.MANIFEST, self.manifest)
        with self.assertRaisesRegex(ValueError, 'Invalid path'):
            b.validate(self.root)

    def test_pack_relocates_without_caches_or_failed_attempts(self):
        put(self.root, '.download-cache/unused.img')
        put(self.root, '.download-attempts/partial/private.txt')
        archive = self.base / 'bundle.tar'
        b.pack(argparse.Namespace(bundle=self.root, archive=archive))
        with tarfile.open(archive) as tar:
            self.assertFalse(any('.download-cache' in x or 'private' in x for x in tar.getnames()))
            tar.extractall(self.base / 'moved', filter='data')
        moved = self.base / 'moved' / self.root.name
        self.root.resolve().relative_to(self.base.resolve())
        shutil.rmtree(self.root)
        b.validate(moved, verify=True)
        b.activate(argparse.Namespace(bundle=moved, pipeline='chipseq', revision='2.1.0',
                                      profile='apptainer', state=self.base / 'state',
                                      output=self.base / 'activate.env'))
        config = (self.base / 'state/nfcore_offline.config').read_text(encoding='utf-8')
        self.assertIn(moved.as_posix(), config)
        env = (self.base / 'activate.env').read_text(encoding='utf-8')
        self.assertIn('NXF_OFFLINE=true', env)
        self.assertIn('NXF_PLUGINS_DEFAULT=false', env)

    def test_refuse_overwriting_archive_and_archive_inside_bundle(self):
        archive = put(self.base, 'existing.tar')
        for dest in (archive, self.root / 'recursive.tar'):
            with self.assertRaises(ValueError):
                b.pack(argparse.Namespace(bundle=self.root, archive=dest))

    def test_conda_and_remote_profiles_are_not_certified_offline(self):
        for profile in ('conda', 'singularity,test', 'wave'):
            with self.assertRaisesRegex(ValueError, 'Complete offline execution'):
                b.activate(argparse.Namespace(bundle=self.root, pipeline='chipseq', revision='2.1.0',
                                              profile=profile, state=self.base / 'state',
                                              output=self.base / 'env'))

    def test_runtime_env_does_not_reuse_unrelated_paths(self):
        with patch.dict(os.environ, {'NXF_BIN': '/wrong.jar', 'NXF_PLUGINS_DIR': '/wrong',
                                     'NXF_PLUGINS_REGISTRY_URL': 'https://wrong.example'}):
            env = b.bundle_env(self.root, True, '25.10.4')
        self.assertNotIn('NXF_BIN', env)
        self.assertNotIn('NXF_PLUGINS_REGISTRY_URL', env)
        self.assertEqual(env['NXF_PLUGINS_DIR'], str(self.root / 'nextflow/plugins'))

    def test_prepare_download_flags_and_fixed_plugins(self):
        root = self.base / 'fresh'
        commands = []

        def fake_run(command, **kwargs):
            commands.append(command)
            self.assertEqual(kwargs['env']['NXF_VER'], '25.10.4')
            if command[0] == 'nf-core':
                self.assertEqual(command[command.index('--download-configuration') + 1], 'yes')
                self.assertNotIn('--force', command)
                payload(Path(command[command.index('--outdir') + 1]))
            elif 'install' in command:
                self.assertEqual(command[-1], 'nf-validation@1.1.3')
                put(root, 'nextflow/plugins/nf-validation-1.1.3/lib/plugin.jar')
            elif '-version' in command:
                put(root, 'nextflow/framework/25.10.4/nextflow-25.10.4-one.jar')
                return 'nextflow version 25.10.4\n'

        args = argparse.Namespace(bundle=root, pipeline='chipseq', revision='2.1.0',
                                  nextflow_version='25.10.4', mode='singularity', parallel_downloads=4)
        with patch.object(b.platform, 'system', return_value='Linux'), \
                patch.object(b.shutil, 'which', return_value='/bin/tool'), \
                patch.object(b, 'download', side_effect=lambda url, dest: put(dest.parent, dest.name)), \
                patch.object(b, 'run', side_effect=fake_run), \
                patch.object(b.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, 'help')):
            b.prepare(args)
            b.prepare(args)  # idempotent: no second download
        self.assertEqual(sum(c[0] == 'nf-core' for c in commands), 1)
        b.validate(root, verify=True)

    def test_failed_prepare_never_marks_complete(self):
        root = self.base / 'failed'
        args = argparse.Namespace(bundle=root, pipeline='chipseq', revision='2.1.0',
                                  nextflow_version='25.10.4', mode='singularity', parallel_downloads=4)
        with patch.object(b.platform, 'system', return_value='Linux'), \
                patch.object(b.shutil, 'which', return_value='/bin/tool'), \
                patch.object(b, 'download', side_effect=OSError('network unavailable')):
            with self.assertRaises(OSError):
                b.prepare(args)
        self.assertFalse((root / b.MANIFEST).exists())
        self.assertTrue((root / '.prepare.json').exists())


if __name__ == '__main__':
    unittest.main()
