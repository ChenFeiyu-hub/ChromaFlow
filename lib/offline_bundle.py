#!/usr/bin/env python3
"""Prepare and validate relocatable nf-core bundles. Python standard library only."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request

FORMAT = 1
MANIFEST = 'chromaflow-bundle.json'
# ATAC 2.1.2 declares nf-validation without a version. Pin the offline copy only.
RELEASE_PLUGIN_PINS = {('atacseq', '2.1.2'): {'nf-validation': '1.1.3'}}


def fail(message):
    raise ValueError(message)


def write_json(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    os.replace(tmp, path)


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def contained(root, relative):
    p = root / relative
    if Path(relative).is_absolute() or '..' in Path(relative).parts:
        fail('Invalid path in offline bundle: ' + str(relative))
    try:
        p.resolve().relative_to(root.resolve())
    except ValueError:
        fail('Bundle path or symlink points outside the bundle: ' + str(relative))
    return p


def run(command, env=None, cwd=None, capture=False):
    print('[CMD] ' + shlex.join([str(x) for x in command]), flush=True)
    result = subprocess.run(command, env=env, cwd=cwd, check=True,
                            stdout=subprocess.PIPE if capture else None,
                            text=True)
    return result.stdout if capture else None


def bundle_env(root, offline, version, state_dir=None):
    env = os.environ.copy()
    # Do not inherit an unrelated user's runtime, plugin or container cache.
    for key in ('NXF_BIN', 'NXF_JAR', 'NXF_LAUNCHER', 'NXF_PACK', 'NXF_CLASSPATH',
                'NXF_PLUGINS_TEST_REPOSITORY', 'NXF_PLUGINS_REGISTRY_URL', 'NXF_ASSETS'):
        env.pop(key, None)
    env.update(NXF_HOME=str(state_dir or root / 'nextflow'),
               NXF_DIST=str(root / 'nextflow/framework'),
               NXF_PLUGINS_DIR=str(root / 'nextflow/plugins'),
               NXF_VER=version, NXF_OFFLINE=str(offline).lower(),
               NXF_DISABLE_CHECK_LATEST='true', NXF_PLUGINS_DEFAULT='false',
               NXF_SYNTAX_PARSER='v1')
    return env


def download(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_suffix('.part')
    print('[DOWNLOAD] ' + url, flush=True)
    request = urllib.request.Request(url, headers={'User-Agent': 'ChromaFlow-offline'})
    with urllib.request.urlopen(request, timeout=120) as response, part.open('wb') as out:
        shutil.copyfileobj(response, out)
    os.replace(part, destination)


def pin_release_plugins(text, pipeline, revision):
    pins = RELEASE_PLUGIN_PINS.get((pipeline, revision), {})
    if not pins:
        return text

    def pin_declaration(match):
        prefix, quote, name = match.groups()
        if name not in pins:
            return match.group(0)
        print('[PIN] ' + name + '@' + pins[name], flush=True)
        return prefix + quote + name + '@' + pins[name] + quote

    def pin_block(match):
        return re.sub(r'''(?m)^([ \t]*id[ \t]+)(['"])([\w.-]+)\2''',
                      pin_declaration, match.group(0))

    return re.sub(r'(?m)^\s*plugins\s*\{[^}]+\}', pin_block, text)


def plugins_in(workflow):
    text = (workflow / 'nextflow.config').read_text(encoding='utf-8')
    blocks = re.findall(r'(?m)^\s*plugins\s*\{([^}]+)\}', text)
    specs = sorted(set(re.findall(r"\bid\s+['\"]([^'\"]+)['\"]", '\n'.join(blocks))))
    if not specs or any(not re.fullmatch(r'[\w.-]+@[\w.-]+', x) for x in specs):
        fail('Cannot identify pinned plugin declarations; check workflow/nextflow.config.')
    return specs


def plugin_directory(root, spec):
    # Plugin declarations use id@version; Nextflow installs them as id-version.
    name, version = spec.split('@', 1)
    return root / 'nextflow/plugins' / (name + '-' + version)


def plugin_ready(root, spec):
    folder = plugin_directory(root, spec)
    return folder.is_dir() and any(p.stat().st_size > 0 for p in folder.rglob('*.jar'))


def safe_links(root, subdir):
    # nf-core's image aliases must survive copying the bundle to another path.
    for p in subdir.rglob('*'):
        if p.is_symlink():
            target = p.resolve(strict=True)
            try:
                target.relative_to(root)
            except ValueError:
                fail('Download contains an external symlink; refusing to create a relocatable bundle: ' + str(p))
            if os.path.isabs(os.readlink(p)):
                p.unlink()
                p.symlink_to(os.path.relpath(target, p.parent), target_is_directory=target.is_dir())


def inventory(root, roots):
    files = {}
    for relative in roots:
        folder = root / relative
        paths = [folder] if folder.is_file() else sorted(folder.rglob('*'))
        for p in paths:
            if p.is_symlink() and p.is_dir():
                fail('Directory symlinks are not supported: ' + str(p))
            if p.is_file():
                name = p.relative_to(root).as_posix()
                contained(root, name)
                if p.is_symlink():
                    files[name] = {'link': os.readlink(p)}
                else:
                    files[name] = {'size': p.stat().st_size, 'sha256': digest(p)}
    return files


def validate(root, verify=False, pipeline=None, revision=None):
    root = root.resolve()
    if not (root / MANIFEST).is_file():
        fail('Offline bundle is incomplete; missing ' + str(root / MANIFEST))
    manifest = read_json(root / MANIFEST)
    if manifest.get('format') != FORMAT or manifest.get('complete') is not True:
        fail('Unsupported bundle format or incomplete preparation.')
    if pipeline and manifest['pipeline'] != pipeline:
        fail('Bundle pipeline does not match --type.')
    if revision and manifest['revision'] != revision:
        fail('Bundle revision ' + manifest['revision'] + ' differs from configured revision ' + revision + '.')
    required = [manifest['workflow'] + '/main.nf', manifest['workflow'] + '/nextflow.config',
                manifest['configs'] + '/nfcore_custom.config', 'bin/nextflow',
                'nextflow/framework/' + manifest['nextflow_version'] + '/nextflow-' +
                manifest['nextflow_version'] + '-one.jar']
    entries = manifest.get('files', {})
    if not entries or any(name not in entries for name in required):
        fail('Bundle manifest is missing workflow, configuration, or Nextflow runtime files.')
    if any(contained(root, name).stat().st_size == 0 for name in required):
        fail('A required bundle file is empty.')
    for name, item in entries.items():
        p = contained(root, name)
        if not p.is_file():
            fail('Missing bundle file or broken symlink: ' + name)
        if 'link' in item:
            if not p.is_symlink() or os.readlink(p) != item['link'] or os.path.isabs(item['link']):
                fail('Bundle symlink changed or is not relocatable: ' + name)
        elif p.stat().st_size != item['size']:
            fail('Bundle file size mismatch: ' + name)
        elif verify and digest(p) != item['sha256']:
            fail('Bundle SHA256 mismatch: ' + name)
    workflow = contained(root, manifest['workflow'])
    if plugins_in(workflow) != manifest['plugins']:
        fail('Pipeline plugin declarations differ from the bundle manifest.')
    for spec in manifest['plugins']:
        if not plugin_ready(root, spec):
            fail('Missing plugin files: ' + spec)
    if manifest['mode'] == 'singularity':
        image_dir = contained(root, manifest['images'])
        if not any(p.suffix in ('.img', '.sif') and p.stat().st_size for p in image_dir.iterdir()):
            fail('Bundle contains no container images.')
    return manifest


def write_runtime_config(root, manifest, path, engine):
    def quote(value):
        value = value.as_posix() if isinstance(value, Path) else str(value)
        return "'" + value.replace('\\', '\\\\').replace("'", "\\'") + "'"
    lines = ['// Generated by ChromaFlow offline',
             'params.custom_config_base = ' + quote(root / manifest['configs']),
             'wave.enabled = false', 'tower.enabled = false', 'fusion.enabled = false']
    if engine in ('singularity', 'apptainer'):
        images = quote(root / manifest['images'])
        lines += ['singularity.cacheDir = ' + images, 'apptainer.cacheDir = ' + images,
                  'singularity.ociAutoPull = false', 'apptainer.ociAutoPull = false']
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def prepare(args):
    if platform.system() != 'Linux':
        fail('Run prepare on a connected Linux machine; offline execution also requires Linux.')
    if not shutil.which('java') or not shutil.which('nf-core'):
        fail('Preparation requires Java and the nf-core CLI; activate nfcore_env first.')
    if not re.fullmatch(r'\d+\.\d+\.\d+(?:-edge)?', args.nextflow_version):
        fail('--nextflow-version must be pinned, for example 25.10.4.')
    if not re.fullmatch(r'\d+\.\d+\.\d+', args.revision):
        fail('--revision must be a fixed release, for example 2.1.0.')
    root = args.bundle.resolve()
    spec = {'pipeline': args.pipeline, 'revision': args.revision,
            'nextflow_version': args.nextflow_version, 'mode': args.mode,
            'architecture': platform.machine(), 'system': platform.system()}
    state = root / '.prepare.json'
    if root.exists() and any(root.iterdir()) and not state.is_file():
        fail('Destination is nonempty and is not a ChromaFlow bundle; use a new directory.')
    root.mkdir(parents=True, exist_ok=True)
    if state.exists():
        if read_json(state) != spec:
            fail('This directory belongs to another revision or mode; use a new directory for upgrades.')
    else:
        write_json(state, spec)
    if (root / MANIFEST).exists():
        validate(root, verify=True)
        print('[OK] Offline bundle is complete; reusing it.')
        return
    launcher = root / 'bin/nextflow'
    if not launcher.exists():
        download('https://github.com/nextflow-io/nextflow/releases/download/v' +
                 args.nextflow_version + '/nextflow', launcher)
    launcher.chmod(0o755)
    env = bundle_env(root, False, args.nextflow_version)
    env['PATH'] = str(root / 'bin') + os.pathsep + env.get('PATH', '')
    cache = root / '.download-cache'
    cache.mkdir(exist_ok=True)
    env['NXF_SINGULARITY_CACHEDIR'] = str(cache)
    env['NXF_APPTAINER_CACHEDIR'] = str(cache)
    version_text = run([str(launcher), '-version'], env=env, cwd=root, capture=True)
    if not re.search(r'\bversion\s+' + re.escape(args.nextflow_version) + r'(?:\s|$)', version_text):
        fail('Installed Nextflow version differs from the requested version: ' + version_text.strip())
    print(version_text.strip(), flush=True)
    payload = root / 'payload'
    if not payload.exists():
        # Failed attempts are retained; retries reuse the separate image cache.
        attempts = root / '.download-attempts'
        attempts.mkdir(exist_ok=True)
        target = attempts / ('attempt-' + str(time.time_ns()))
        help_result = subprocess.run(['nf-core', 'pipelines', 'download', '--help'],
                                     stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, env=env, text=True)
        modern = help_result.returncode == 0
        command = ['nf-core', 'pipelines', 'download'] if modern else ['nf-core', 'download']
        command += ['nf-core/' + args.pipeline, '--revision', args.revision,
                    '--outdir', str(target), '--compress', 'none',
                    '--container-system', 'singularity' if args.mode == 'singularity' else 'none']
        # nf-core 3.x pipelines download uses an explicit yes/no value;
        # the deprecated 2.x command used a boolean switch.
        command += ['--download-configuration', 'yes'] if modern else ['--download-configuration']
        if args.mode == 'singularity':
            command += ['--container-cache-utilisation', 'copy', '--parallel-downloads', str(args.parallel_downloads)]
        run(command, env=env, cwd=root)
        if not target.is_dir():
            fail('nf-core did not create the expected download directory.')
        # Normalize links before the downloaded directory changes location.
        safe_links(target.resolve(), target)
        target.rename(payload)
    workflows = [p.parent for p in payload.rglob('main.nf') if (p.parent / 'nextflow.config').is_file()
                 and p.parent.parent == payload]
    if len(workflows) != 1:
        fail('Cannot uniquely identify the downloaded pipeline; inspect the payload directory.')
    workflow = workflows[0]
    text = (workflow / 'nextflow.config').read_text(encoding='utf-8')
    if not re.search(r"(?m)^\s*version\s*=\s*['\"]" + re.escape(args.revision) + r"['\"]", text):
        fail('Downloaded pipeline revision differs from the requested revision.')
    configs = payload / 'configs'
    if not (configs / 'nfcore_custom.config').is_file():
        fail('Download is missing nf-core/configs.')
    # Ensure even older downloaders cannot leave the remote shared-config URL active.
    old = 'https://raw.githubusercontent.com/nf-core/configs/${params.custom_config_version}'
    text = text.replace(old, '${projectDir}/../configs')
    if old in text or not re.search(r'custom_config_base\s*=\s*[\"\'][^\n]*projectDir[^\n]*configs', text):
        fail('Cannot confirm that shared configurations point inside the bundle.')
    text = pin_release_plugins(text, args.pipeline, args.revision)
    (workflow / 'nextflow.config').write_text(text, encoding='utf-8')
    plugins = plugins_in(workflow)
    for plugin in plugins:
        if not plugin_ready(root, plugin):
            run([str(launcher), 'plugin', 'install', plugin], env=env, cwd=root)
        if not plugin_ready(root, plugin):
            fail('Plugin installation produced no nonempty JAR; expected directory: ' + str(plugin_directory(root, plugin)))
    images = None
    if args.mode == 'singularity':
        candidates = [p for p in payload.iterdir() if p.is_dir() and p.name.startswith('singularity-')
                      and any(x.suffix in ('.img', '.sif') for x in p.iterdir())]
        if len(candidates) != 1:
            fail('Container download is incomplete or image directory cannot be identified.')
        images = candidates[0].relative_to(root).as_posix()
    safe_links(root, payload)
    manifest = dict(spec, format=FORMAT, complete=True,
                    workflow=workflow.relative_to(root).as_posix(), configs='payload/configs',
                    plugins=plugins, images=images,
                    created=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    manifest['files'] = inventory(root, ['payload', 'bin', 'nextflow/framework', 'nextflow/plugins'])
    # Publish the completion marker only after all downloads and inventory succeed.
    write_json(root / MANIFEST, manifest)
    try:
        validate(root, verify=True)
    except (ValueError, OSError, KeyError):
        os.replace(root / MANIFEST, root / '.failed-manifest.json')
        raise
    print('[OK] Prepared: ' + str(root))
    if args.mode == 'code-only':
        print('[NOTE] code-only includes source, runtime and plugins, but no analysis tools; it is not a complete offline bundle.')
    print('After transfer, run offline check --verify --runtime-check, then validate a small local dataset.')


def check(args):
    root = args.bundle.resolve()
    m = validate(root, args.verify)
    print('[OK] ' + m['pipeline'] + ' ' + m['revision'] + '; Nextflow ' + m['nextflow_version'])
    print('[OK] Plugins: ' + ', '.join(m['plugins']))
    print('[OK] File manifest' + (' and SHA256' if args.verify else ', file sizes and symlinks') + ' passed')
    print('[INFO] Bundle mode: ' + m['mode'] + '; Build architecture: ' + m['architecture'])
    if m['mode'] == 'code-only':
        print('[NOTE] This bundle lacks analysis tools and cannot be used with run --offline.')
    if args.runtime_check:
        if platform.system() != 'Linux':
            fail('--runtime-check requires Linux, Java and the selected container engine.')
        if platform.machine() != m['architecture']:
            fail('Runtime and preparation machines have different CPU architectures.')
        engine = 'apptainer' if shutil.which('apptainer') else 'singularity'
        if m['mode'] == 'singularity' and not shutil.which(engine):
            fail('Apptainer/Singularity was not found.')
        state = args.state_dir.resolve() if args.state_dir else Path(tempfile.mkdtemp(prefix='chromaflow-offline-check-'))
        state.mkdir(parents=True, exist_ok=True)
        print('[LOG] Offline startup check log: ' + str(state / 'nextflow.log'))
        config = state / 'offline.config'
        write_runtime_config(root, m, config, engine if m['mode'] == 'singularity' else 'conda')
        env = bundle_env(root, True, m['nextflow_version'], state / 'home')
        run([str(root / 'bin/nextflow'), '-log', str(state / 'nextflow.log'),
             '-C', str(root / m['workflow'] / 'nextflow.config') + ',' + str(config),
             'run', str(root / m['workflow']), '--help'], env=env, cwd=state)
        print('[OK] Offline startup and plugin loading passed; complete data analysis has not yet been validated.')


def pack(args):
    root = args.bundle.resolve()
    m = validate(root, verify=True)
    archive = args.archive.resolve()
    try:
        archive.relative_to(root)
    except ValueError:
        pass
    else:
        fail('Archive must be outside the bundle directory.')
    if archive.exists():
        fail('Archive already exists; refusing to overwrite: ' + str(archive))
    archive.parent.mkdir(parents=True, exist_ok=True)
    part = archive.with_name(archive.name + '.partial')
    if part.exists():
        fail('An incomplete archive already exists; use another output filename: ' + str(part))
    # Images are compressed internally. Plain tar is faster and preserves links/modes.
    with tarfile.open(part, 'w') as out:
        for name in [MANIFEST] + sorted(m['files']):
            out.add(root / name, arcname=root.name + '/' + name, recursive=False)
    os.replace(part, archive)
    checksum = digest(archive)
    archive.with_name(archive.name + '.sha256').write_bytes((checksum + '  ' + archive.name + '\n').encode('utf-8'))
    print('[OK] ' + str(archive))


def activate(args):
    root = args.bundle.resolve()
    m = validate(root, pipeline=args.pipeline, revision=args.revision)
    if m['mode'] != 'singularity' or args.profile not in ('apptainer', 'singularity'):
        fail('Complete offline execution requires a container bundle and --profile apptainer or singularity; online Conda execution remains available.')
    if platform.system() == 'Linux' and platform.machine() != m['architecture']:
        fail('Bundle CPU architecture differs from this machine.')
    state = args.state.resolve()
    config = state / 'nfcore_offline.config'
    write_runtime_config(root, m, config, args.profile)
    env = bundle_env(root, True, m['nextflow_version'], state / 'nextflow-home')
    exports = {k: v for k, v in env.items() if k.startswith('NXF_') and k in (
        'NXF_HOME', 'NXF_DIST', 'NXF_PLUGINS_DIR', 'NXF_VER', 'NXF_OFFLINE',
        'NXF_DISABLE_CHECK_LATEST', 'NXF_PLUGINS_DEFAULT', 'NXF_SYNTAX_PARSER')}
    exports.update(NFCORE_RUN_TARGET=(root / m['workflow']).as_posix(),
                   NFCORE_NEXTFLOW=(root / 'bin/nextflow').as_posix(), NFCORE_OFFLINE_CONFIG=config.as_posix(),
                   NXF_SINGULARITY_CACHEDIR=(root / m['images']).as_posix(),
                   NXF_APPTAINER_CACHEDIR=(root / m['images']).as_posix())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', encoding='utf-8') as out:
        out.write('unset NXF_BIN NXF_JAR NXF_LAUNCHER NXF_PACK NXF_CLASSPATH NXF_ASSETS\n')
        out.write('unset NXF_PLUGINS_TEST_REPOSITORY NXF_PLUGINS_REGISTRY_URL\n')
        for key, value in exports.items():
            out.write('export ' + key + '=' + shlex.quote(value) + '\n')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    prepare_p = sub.add_parser('prepare', help='Prepare a pinned offline bundle on connected Linux')
    prepare_p.add_argument('--bundle', type=Path, required=True)
    prepare_p.add_argument('--pipeline', choices=['chipseq', 'atacseq'], default='chipseq')
    prepare_p.add_argument('--revision', required=True)
    prepare_p.add_argument('--nextflow-version', default='25.10.4')
    prepare_p.add_argument('--mode', choices=['singularity', 'code-only'], default='singularity')
    prepare_p.add_argument('--parallel-downloads', type=int, choices=range(1, 17), default=4)
    check_p = sub.add_parser('check', help='Check offline bundle integrity')
    check_p.add_argument('--bundle', type=Path, required=True)
    check_p.add_argument('--verify', action='store_true')
    check_p.add_argument('--runtime-check', action='store_true')
    check_p.add_argument('--state-dir', type=Path, help='Startup logs and runtime cache directory; defaults to a temporary directory')
    pack_p = sub.add_parser('pack', help='Create a portable tar archive and SHA256 checksum')
    pack_p.add_argument('--bundle', type=Path, required=True)
    pack_p.add_argument('--archive', type=Path, required=True)
    active_p = sub.add_parser('activate', help=argparse.SUPPRESS)
    active_p.add_argument('--bundle', type=Path, required=True)
    active_p.add_argument('--pipeline', required=True)
    active_p.add_argument('--revision', required=True)
    active_p.add_argument('--profile', required=True)
    active_p.add_argument('--state', type=Path, required=True)
    active_p.add_argument('--output', type=Path, required=True)
    return p


def main():
    args = parser().parse_args()
    try:
        {'prepare': prepare, 'check': check, 'pack': pack, 'activate': activate}[args.action](args)
    except (ValueError, OSError, subprocess.CalledProcessError, KeyError) as e:
        print('[ERROR] ' + str(e), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
