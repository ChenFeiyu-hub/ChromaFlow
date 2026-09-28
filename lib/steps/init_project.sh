#!/usr/bin/env bash
# ============================================================================
# init_project.sh — 生成新项目骨架（由 `chromaflow init` 调用）
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME}/lib/common.sh"

PROJECT_DIR="${1:?Usage: chromaflow init <project_dir> [type]}"
TYPE="${2:-chipseq}"
TEMPLATES="${CHROMAFLOW_HOME}/templates"

[[ -e "${PROJECT_DIR}" ]] && die "Destination already exists: ${PROJECT_DIR}(use a new directory to avoid overwriting)"

mkdir -p "${PROJECT_DIR}"/{inputs/{raw,reference},outputs}
cp "${TEMPLATES}/chromaflow.config" "${PROJECT_DIR}/chromaflow.config"

# 样本表模板
if [[ -f "${TEMPLATES}/samples_${TYPE}.txt" ]]; then
    cp "${TEMPLATES}/samples_${TYPE}.txt" "${PROJECT_DIR}/samples_${TYPE}.txt"
fi

# 把绝对 PROJECT_DIR 写进新项目的 config
abs_dir="$(cd "${PROJECT_DIR}" && pwd)"
sed -i "s|^PROJECT_DIR=.*|PROJECT_DIR=\"${abs_dir}\"|" "${PROJECT_DIR}/chromaflow.config"

cat <<EOF

[DONE] Project created: ${abs_dir}

  ${abs_dir}/
  ├── chromaflow.config        ← Edit: species and reference source paths
  ├── samples_${TYPE}.txt       ← Edit: samples and FASTQ paths
  ├── inputs/
  │   ├── raw/                  ← Place FASTQ here, or configure another input path
  │   └── reference/            ← Place reference source files here
  └── outputs/                  ← Generated results, logs and work files

Next step:
  cd ${abs_dir}
  vim chromaflow.config          # Species and reference source paths
  vim samples_${TYPE}.txt         # Samples and FASTQ paths

  chromaflow doctor              # Check configuration and tools
  chromaflow reference           # Prepare reference files if needed
  chromaflow samplesheet --type ${TYPE}
  chromaflow run         --type ${TYPE}

Help:
  chromaflow help                # Command overview
  chromaflow help config         # Configuration reference
  chromaflow help samples        # Sample records and FASTQ pairing
EOF
