#!/usr/bin/env bash
# ============================================================================
# init_project.sh — 生成新项目骨架（由 `chromaflow init` 调用）
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME}/lib/common.sh"

PROJECT_DIR="${1:?用法: chromaflow init <project_dir> [type]}"
TYPE="${2:-chipseq}"
TEMPLATES="${CHROMAFLOW_HOME}/templates"

[[ -e "${PROJECT_DIR}" ]] && die "目标已存在: ${PROJECT_DIR}（避免覆盖，请换目录）"

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

[DONE] 已创建项目: ${abs_dir}

  ${abs_dir}/
  ├── chromaflow.config        ← 编辑：物种 + 源文件路径
  ├── samples_${TYPE}.txt       ← 编辑：样本与 fastq 路径
  ├── inputs/
  │   ├── raw/                  ← 放 fastq（或在 config 指向别处）
  │   └── reference/            ← 放参考基因组源文件
  └── outputs/                  ← 所有结果/日志/work（自动生成）

下一步:
  cd ${abs_dir}
  vim chromaflow.config          # 物种 + 源文件路径（派生文件不用填）
  vim samples_${TYPE}.txt         # 样本与 fastq 路径

  chromaflow doctor              # 自检
  chromaflow reference           # 可选：整理参考（已有基因组/GTF 可跳过）
  chromaflow samplesheet --type ${TYPE}
  chromaflow run         --type ${TYPE}

帮助:
  chromaflow help                # 命令总览
  chromaflow help config         # config 里该写什么
  chromaflow help samples        # 样本表格式与双端配对规则
EOF
