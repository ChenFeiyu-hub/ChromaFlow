#!/usr/bin/env bash
# ============================================================================
# 02_run_nfcore_conda.sh — 基于 conda profile 运行 nf-core pipeline
# ============================================================================
# ⚠ 备选方案：当 Apptainer 容器下载失败时使用
#   nf-core 官方推荐容器方式，conda 可能存在环境冲突风险，
#   但在网络受限环境下是可行的替代方案。
#
# 与 02_run_nfcore.sh 的区别：
#   - 使用 -profile conda 替代 -profile singularity
#   - 首次运行时 Nextflow 会自动为每个 process 创建独立 conda 环境
#   - 环境缓存在 NXF_CONDA_CACHEDIR 中，后续运行复用
#
# 用法: chromaflow run --profile conda --type [chipseq|atacseq|dapseq]
#   或在 chromaflow.config 里设 NFCORE_PROFILE=conda（`chromaflow run` 会自动
#   路由到本脚本，无需手动调用）。NFCORE_PROFILE=mamba 亦可。
# ============================================================================

set -euo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
PROFILE="$(resolve_nfcore_profile)"
[[ "$(nfcore_profile_kind "${PROFILE}")" == "conda" ]] || PROFILE="conda"

# ─────────────────────────────────────────────
# conda 环境缓存目录
# ─────────────────────────────────────────────
export NXF_CONDA_CACHEDIR="${OUTPUT_DIR}/.conda_envs"
mkdir -p "${NXF_CONDA_CACHEDIR}"

# ─────────────────────────────────────────────
# 确定 pipeline 和版本
# ─────────────────────────────────────────────
case "${PIPELINE}" in
    chipseq|dapseq)
        NFCORE_NAME="chipseq"
        NFCORE_VERSION="${NFCORE_CHIPSEQ_VERSION}"
        ;;
    atacseq)
        NFCORE_NAME="atacseq"
        NFCORE_VERSION="${NFCORE_ATACSEQ_VERSION}"
        ;;
    *)
        echo "[ERROR] 未知 pipeline: ${PIPELINE} (支持: chipseq, atacseq, dapseq)"
        exit 1
        ;;
esac

SAMPLESHEET="${SAMPLESHEET_DIR}/samplesheet_${PIPELINE}.csv"
OUTDIR="${RESULTS_DIR}/${PIPELINE}"

open_log "02_run_nfcore_conda" "${PIPELINE}"

log_section "Run nf-core/${NFCORE_NAME} v${NFCORE_VERSION} (${PIPELINE}) [profile=${PROFILE}]"

echo "  ⚠  使用 conda profile（Apptainer 备选方案）"
echo "  首次运行需为每个 process 创建 conda 环境，耗时较长（约 30-60 分钟）"
echo "  conda 环境缓存: ${NXF_CONDA_CACHEDIR}"
echo ""

# ─────────────────────────────────────────────
# 前置检查
# ─────────────────────────────────────────────
echo "[CHECK] 环境检查"

check_tool nextflow || exit 1
echo "[  OK ] nextflow $(nextflow -version 2>&1 | grep -oP 'version \K\S+')"

# 检查 conda/mamba
if command -v mamba &>/dev/null; then
    echo "[  OK ] mamba $(mamba --version 2>&1 | head -1)"
    # Nextflow 优先使用 mamba（更快）
    export NXF_CONDA_USE_MAMBA=true
elif command -v conda &>/dev/null; then
    echo "[  OK ] conda $(conda --version 2>&1)"
else
    echo "[ERROR] 未检测到 conda/mamba"; exit 1
fi

[[ ! -f "${SAMPLESHEET}" ]] && echo "[ERROR] samplesheet 不存在: ${SAMPLESHEET}" \
    && echo "  请先运行: chromaflow samplesheet --type ${PIPELINE}" && exit 1
[[ ! -f "${GENOME_COMBINED}" ]] && echo "[ERROR] 参考基因组不存在: ${GENOME_COMBINED}" \
    && echo "  请先运行: chromaflow reference" && exit 1
[[ ! -f "${GENOME_GTF}" ]] && echo "[ERROR] GTF 不存在: ${GENOME_GTF}" && exit 1

echo "[  OK ] samplesheet: ${SAMPLESHEET}"
echo "[  OK ] genome:      ${GENOME_COMBINED}"
echo "[  OK ] GTF:         ${GENOME_GTF}"
echo ""

# ─────────────────────────────────────────────
# 配置 conda channel 镜像（加速国内下载）
# ─────────────────────────────────────────────
echo "[CONDA] 检查 channel 配置"

# conda channel 镜像：NFCORE_CONDA_MIRROR=tuna(默认) | none
#   none = 不改动 channel，适合境外网络或已配好内网自建源的环境
if [[ "${NFCORE_CONDA_MIRROR:-tuna}" == "none" ]]; then
    echo "[SKIP] NFCORE_CONDA_MIRROR=none，使用系统默认 channel"
else
    CONDA_RC="${OUTPUT_DIR}/.condarc_nfcore"
    if [[ ! -f "${CONDA_RC}" ]]; then
        cat > "${CONDA_RC}" << 'EOF'
channels:
  - bioconda
  - conda-forge
  - defaults
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
show_channel_urls: true
EOF
        echo "[DONE] 已生成镜像配置: ${CONDA_RC}"
    else
        echo "[SKIP] 镜像配置已存在"
    fi
    export CONDARC="${CONDA_RC}"
    echo "[  OK ] 使用清华镜像加速 conda 下载（NFCORE_CONDA_MIRROR=none 可关闭）"
fi

echo ""

# ─────────────────────────────────────────────
# 构建运行命令
# ─────────────────────────────────────────────
# ─────────────────────────────────────────────
# 有效基因组大小 / 读长
# ─────────────────────────────────────────────
# 两条路都不需要用户手填：
#   A. 已有实测 MACS_GSIZE（00 用 khmer 测的，存在 derived.env）→ 直接传
#   B. 没有 → 只传 --read_length，由 nf-core 内部的 KHMER_UNIQUEKMERS 自己算
#      （nf-core/chipseq 2.x 起内置该逻辑，与我们 00 里用的是同一个工具）
GSIZE_ARGS=()
if [[ -n "${MACS_GSIZE}" ]]; then
    GSIZE_ARGS+=(--macs_gsize "${MACS_GSIZE}")
    echo "[  OK ] MACS_GSIZE: ${MACS_GSIZE} (来源: $(value_origin MACS_GSIZE)$(
        m="$(derived_get MACS_GSIZE_METHOD 2>/dev/null || true)"; [[ -n "${m}" ]] && echo ", ${m}")）"
elif [[ -n "${READ_LENGTH}" ]]; then
    echo "[INFO ] 无 MACS_GSIZE，交给 nf-core 用 --read_length=${READ_LENGTH} 自行计算"
else
    die "既无 MACS_GSIZE 也无 READ_LENGTH。先跑 chromaflow reference（自动实测），
     或 chromaflow samplesheet（实测读长），或在 config 里显式指定其一"
fi
if [[ -n "${READ_LENGTH}" ]]; then GSIZE_ARGS+=(--read_length "${READ_LENGTH}"); fi

COMMON_ARGS=(
    --input "${SAMPLESHEET}"
    --outdir "${OUTDIR}"
    --fasta "${GENOME_COMBINED}"
    --gtf "${GENOME_GTF}"
    --max_cpus "${MAX_CPUS}"
    --max_memory "${MAX_MEMORY}"
    --max_time "${MAX_TIME}"
    -profile "${PROFILE}"
    -r "${NFCORE_VERSION}"
    -w "${WORK_DIR}"
)
COMMON_ARGS+=(${GSIZE_ARGS[@]+"${GSIZE_ARGS[@]}"})

EXTRA_ARGS=()
case "${PIPELINE}" in
    chipseq|dapseq)
        EXTRA_ARGS+=(--narrow_peak)
        [[ -f "${BLACKLIST_BED}" ]] && EXTRA_ARGS+=(--blacklist "${BLACKLIST_BED}")
        ;;
    atacseq)
        [[ -f "${BLACKLIST_BED}" ]] && EXTRA_ARGS+=(--blacklist "${BLACKLIST_BED}")
        ;;
esac

echo "============================================"
echo " 运行命令"
echo "============================================"
echo ""
echo "nextflow run nf-core/${NFCORE_NAME} \\"
for arg in "${COMMON_ARGS[@]}"; do echo "  ${arg} \\"; done
for arg in "${EXTRA_ARGS[@]}"; do echo "  ${arg} \\"; done
echo ""
echo "环境变量:"
echo "  NXF_CONDA_CACHEDIR=${NXF_CONDA_CACHEDIR}"
echo "  NXF_CONDA_USE_MAMBA=${NXF_CONDA_USE_MAMBA:-false}"
echo "  CONDARC=${CONDARC}"
echo ""

# ─────────────────────────────────────────────
# 确认执行
# ─────────────────────────────────────────────
read -p "确认运行? [y/N] " confirm
if [[ "${confirm}" != "y" && "${confirm}" != "Y" ]]; then
    echo "[ABORT] 已取消"
    exit 0
fi

echo ""
echo "[START] $(date)"
echo ""

cd "${OUTPUT_DIR}"

nextflow run "nf-core/${NFCORE_NAME}" \
    "${COMMON_ARGS[@]}" \
    "${EXTRA_ARGS[@]}" \
    -resume

echo ""
echo "[DONE] nf-core/${NFCORE_NAME} (${PIPELINE}, conda mode) 完成 $(date)"
echo "  结果: ${OUTDIR}"
