# Offline nf-core deployment with ChromaFlow

This guide applies to ChromaFlow 1.2.7, nf-core/chipseq 2.1.0, nf-core/atacseq 2.1.2 and Nextflow 25.10.4. DAP-seq uses the chipseq bundle. The preparation computer needs Linux/WSL, network access, Python, Java, nf-core/tools and Apptainer/Singularity. The server needs Linux, Python, Java and Apptainer/Singularity; it does not need the nf-core CLI or process-specific Conda environments.

## 1. Install ChromaFlow on the connected computer

Extract the ChromaFlow release and add its `bin` directory to PATH. Use a Linux path inside WSL:

```bash
cd /mnt/d/biosoft/ChromaFlow-1.2.7/ChromaFlow  # adjust to your extraction location
export CHROMAFLOW_HOME="$PWD"
export PATH="$CHROMAFLOW_HOME/bin:$PATH"
chromaflow version
chromaflow setup --mode prepare --components runtime --dry-run
chromaflow setup --mode prepare --components runtime
conda activate nfcore_env
java -version
python3 --version
nf-core --version
apptainer --version  # singularity --version is also acceptable
```

The package installer requires an existing Conda-compatible solver; it does not install Miniforge or the system container runtime. Existing environments are reused unchanged. If the required tools are already available, skip `setup` and use them. Windows `D:\...` paths cannot be passed directly to Bash `cd`; use `/mnt/d/...`.

## 2. Prepare complete bundles

```bash
mkdir -p "$HOME/nfcore-offline" "$HOME/nfcore-transfer"
chromaflow offline prepare --type chipseq \
  --revision 2.1.0 --bundle "$HOME/nfcore-offline/chipseq-2.1.0" --parallel-downloads 4
chromaflow offline prepare --type atacseq \
  --revision 2.1.2 --bundle "$HOME/nfcore-offline/atacseq-2.1.2" --parallel-downloads 4

for NAME in chipseq-2.1.0 atacseq-2.1.2; do
  chromaflow offline check --bundle "$HOME/nfcore-offline/$NAME" --verify --runtime-check
  chromaflow offline pack --bundle "$HOME/nfcore-offline/$NAME" \
    --archive "$HOME/nfcore-transfer/$NAME.tar"
done
```

ChromaFlow orchestrates workflow/configuration download, image retrieval, Nextflow initialization, pinned plugin installation, bundle inventory and hashing. Repeating `prepare` on the same incomplete bundle reuses completed downloads. Do not delete a partially completed bundle to troubleshoot a plugin error. `--mode code-only` excludes containers and cannot be used as a complete offline analysis bundle.

`offline check --verify` checks the file inventory and SHA256 hashes. `--runtime-check` starts the local workflow help with offline settings and checks plugins/runtime startup. It does not execute every container or analyze a dataset; run a small local dataset on the target server before general use.

If downloads encounter GitHub 403 responses, inspect authentication/rate-limit feedback. For a token you have already created, enter it privately in the preparation shell:

```bash
read -rsp 'GitHub token: ' GITHUB_TOKEN; printf '\n'
export GITHUB_TOKEN
# Repeat the same chromaflow offline prepare command.
unset GITHUB_TOKEN
```

A token cannot repair a TCP timeout to container registries. Fix connectivity/proxy routing on the preparation computer. The target server does not need this token, proxy or GitHub access. Lower `--parallel-downloads` to 1 or 2 if concurrent downloads are unstable.

## 3. Upload archives and checksums

Copy these files using your approved transfer method:

```text
ChromaFlow-1.2.7.tar.gz
ChromaFlow-1.2.7.tar.gz.sha256
chipseq-2.1.0.tar
chipseq-2.1.0.tar.sha256
atacseq-2.1.2.tar
atacseq-2.1.2.tar.sha256
```

For example, if SSH transfer is available, replace `USER` and `SERVER` with your own login and server address, and create `~/nfcore-upload` on the server first:

```bash
scp "$HOME/nfcore-transfer/"*.tar "$HOME/nfcore-transfer/"*.sha256 USER@SERVER:~/nfcore-upload/
```

The ChromaFlow release is copied separately from wherever you saved it. Complete uploads before checking/extracting. Checksums created by 1.2.7 use LF line endings on all platforms. If a legacy checksum file contains CRLF, remove only its CR characters before verifying; do not change the archive.

## 4. Shared bundle installation under /opt

Run on the server as an account with sudo. These commands assume the two bundle names above, which become the archive's top-level directories:

```bash
cd "$HOME/nfcore-upload"
sha256sum -c ChromaFlow-1.2.7.tar.gz.sha256
sha256sum -c chipseq-2.1.0.tar.sha256
sha256sum -c atacseq-2.1.2.tar.sha256
tar -tf chipseq-2.1.0.tar | head
tar -tf atacseq-2.1.2.tar | head

sudo install -d -m 0755 /opt/chromaflow/bundles

# Run once per NEW bundle; do not unpack over an installed version.
for NAME in chipseq-2.1.0 atacseq-2.1.2; do
  if sudo test -e "/opt/chromaflow/bundles/$NAME"; then
    printf '[SKIP] Already installed: %s\n' "$NAME"
  else
    sudo tar --no-same-owner -xf "$NAME.tar" -C /opt/chromaflow/bundles
    sudo chmod -R a+rX "/opt/chromaflow/bundles/$NAME"
  fi
done
```

Use a dedicated version directory for each bundle. Ordinary accounts need read permission for files and traversal permission on directories; they do not need write permission inside `/opt`. Files extracted by sudo above are administrator-owned. Do not run analyses with sudo. Per-user projects, Nextflow work, state and logs remain in the user's writable project, with startup checks using a temporary state directory by default.

Install ChromaFlow in each user's `~/biosoft/ChromaFlow-1.2.7` as shown in the README, or have the administrator install one readable shared code directory. A dedicated service user is unnecessary for interactive analyses. Each account can retain its own Conda environments.

## 5. Validate on the target server

Reuse existing Python and Java if available. If package sources are reachable during provisioning, a minimal environment can be created with:

```bash
chromaflow setup --mode offline --components runtime
conda activate ChromaFlow_runtime
```

On a fully disconnected server, that installer cannot fetch Conda packages. Use the existing working runtime, system packages provisioned by the administrator, a local package mirror, or environments prepared on a compatible Linux system. A WSL Ubuntu environment is not automatically certified to work on CentOS: CPU architecture and host-library compatibility need checking. Downstream environments have the same separate provisioning requirement; they are not in the nf-core tar files.

```bash
java -version
python3 --version
apptainer --version
chromaflow offline check --bundle /opt/chromaflow/bundles/chipseq-2.1.0 --verify --runtime-check
chromaflow offline check --bundle /opt/chromaflow/bundles/atacseq-2.1.2 --verify --runtime-check
```

Both `apptainer` and `singularity` may be entry points to the same Apptainer installation. Installing another runtime does not solve an inaccessible registry; complete local images avoid those downloads.

## 6. Configure and run a project

For an existing validated project, retain its reference paths and sample settings. Add the matching offline configuration:

```bash
NFCORE_PROFILE="singularity"
NFCORE_OFFLINE=true
NFCORE_OFFLINE_BUNDLE="/opt/chromaflow/bundles/chipseq-2.1.0"
NFCORE_CHIPSEQ_VERSION="2.1.0"
NFCORE_ATACSEQ_VERSION="2.1.2"
KEGG_ORGANISM=""  # no online KEGG request from downstream annotation
```

For a new project use `chromaflow init`, then edit its source FASTA/GFF3, sample records and resource limits. `OUTPUT_DIR="${PROJECT_DIR}"` gives the direct `results/`, `logs/`, `work/` layout; otherwise the default anchor is `outputs/`. Supply local organelle FASTA paths or leave `ORGANELLE_SOURCES=()`. Download reference annotations and motif databases before the disconnected run.

```bash
chromaflow run -c /path/to/chromaflow.config --type chipseq --dry-run
chromaflow run -c /path/to/chromaflow.config --type chipseq --yes

# The ATAC project needs its own samplesheet and suitable assay settings:
chromaflow run -c /path/to/atac.config --type atacseq --offline \
  --bundle /opt/chromaflow/bundles/atacseq-2.1.2 --dry-run
chromaflow run -c /path/to/atac.config --type atacseq --offline \
  --bundle /opt/chromaflow/bundles/atacseq-2.1.2 --yes
```

The selected bundle must match `--type` and the configured pipeline revision. `--bundle` overrides the config path for that invocation. Offline execution uses local workflow code, Nextflow and plugins, disables update/plugin downloads, and uses the transferred image cache. Extra user configuration or optional external resources can still introduce network dependencies; the actual sample run is the final deployment check.

## 7. Manage upgrades and diagnose failures

- **ChromaFlow-only update:** install a new software directory, check `version`, retain project configs and valid bundles. For 1.2.7, neither GFF reference regeneration nor another nf-core run is required to get the new logger/plotter.
- **Pipeline update:** create a new bundle directory, select explicit pipeline/Nextflow versions, validate and transfer it, then change the project bundle path and corresponding `NFCORE_*_VERSION`. Keep the previous bundle for reproducibility. Do not silently replace files inside an installed version.
- **Missing image or plugin:** repair/complete the bundle on the connected computer, verify, repack to a new archive filename and re-transfer. Manually copying a random Nextflow plugin cache does not establish version compatibility.
- **Checksum mismatch:** verify upload completion and compare the source archive; do not bypass verification.
- **Permission error:** confirm `/opt` is readable/traversable and the project is writable by the analysis account.
- **Downstream failure:** inspect `logs/downstream/<type>/*.log`, its final exit code and the active environment. A working nf-core bundle does not imply that MEME or ChIPseeker is installed.

Shared directories hold software only. Each project records its own config, derived values, normalized-annotation reports and logs. Preserve these alongside the results when archiving an analysis.
