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
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true

# 透传给 nextflow 的额外参数：--skip-spp 等快捷开关，或 -- 之后的任意原样参数
#   chromaflow run --skip-spp
#   chromaflow run -- --skip_preseq --skip_igv -resume
USER_ARGS=(); NO_CONDA_PATCH="${NO_CONDA_PATCH:-false}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-spp|--skip_spp) USER_ARGS+=(--skip_spp); shift ;;
        --skip-plot-profile)   USER_ARGS+=(--skip_plot_profile); shift ;;
        --skip-multiqc)        USER_ARGS+=(--skip_multiqc); shift ;;
        --no-conda-patch)      NO_CONDA_PATCH=true; shift ;;
        --) shift; while [[ $# -gt 0 ]]; do USER_ARGS+=("$1"); shift; done ;;
        *) shift ;;
    esac
done
# config 里也可以设 SKIP_SPP=true 长期关闭
if [[ "${SKIP_SPP:-false}" == "true" ]]; then
    case " ${USER_ARGS[*]} " in *" --skip_spp "*) ;; *) USER_ARGS+=(--skip_spp) ;; esac
fi
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
        echo "[ERROR] Unknown pipeline: ${PIPELINE} (Supported: chipseq, atacseq, dapseq)"
        exit 1
        ;;
esac

SAMPLESHEET="${SAMPLESHEET_DIR}/samplesheet_${PIPELINE}.csv"
OUTDIR="${CF_BATCH_NFCORE_OUT:-${RESULTS_DIR}/${PIPELINE}}"

open_log "02_run_nfcore_conda" "${PIPELINE}"

log_section "Run nf-core/${NFCORE_NAME} v${NFCORE_VERSION} (${PIPELINE}) [profile=${PROFILE}]"

source "${CHROMAFLOW_HOME}/lib/nfcore_runtime.sh"
configure_nfcore_runtime

echo "  ⚠  Using the online Conda profile"
echo "  The first run creates per-process Conda environments and may take additional time"
echo "  Conda environment cache: ${NXF_CONDA_CACHEDIR}"
echo ""

# ─────────────────────────────────────────────
# 前置检查
# ─────────────────────────────────────────────
echo "[CHECK] Environment check"

check_tool "${NFCORE_NEXTFLOW}" || exit 1
echo "[  OK ] nextflow $("${NFCORE_NEXTFLOW}" -version 2>&1 | grep -oP 'version \K\S+')"

# 检查 conda/mamba
if command -v mamba &>/dev/null; then
    echo "[  OK ] mamba $(mamba --version 2>&1 | head -1) —— enabled through conda.useMamba"
    # 注意: NXF_CONDA_USE_MAMBA 这个环境变量并不存在，nextflow 只认配置里的
    #       conda.useMamba；真正的开关写在 write_nfcore_conda_patch 生成的配置里。
elif command -v conda &>/dev/null; then
    echo "[  OK ] conda $(conda --version 2>&1)"
else
    echo "[ERROR] Conda/mamba was not detected"; exit 1
fi

[[ ! -f "${SAMPLESHEET}" ]] && echo "[ERROR] samplesheet does not exist: ${SAMPLESHEET}" \
    && echo "  Run first: chromaflow samplesheet --type ${PIPELINE}" && exit 1
[[ ! -f "${GENOME_COMBINED}" ]] && echo "[ERROR] Reference genome does not exist: ${GENOME_COMBINED}" \
    && echo "  Run first: chromaflow reference" && exit 1
[[ ! -f "${GENOME_GTF}" ]] && echo "[ERROR] GTF does not exist: ${GENOME_GTF}" && exit 1

echo "[  OK ] samplesheet: ${SAMPLESHEET}"
echo "[  OK ] genome:      ${GENOME_COMBINED}"
echo "[  OK ] GTF:         ${GENOME_GTF}"
echo ""

# ─────────────────────────────────────────────
# 配置 conda channel 镜像（加速国内下载）
# ─────────────────────────────────────────────
echo "[CONDA] Checking channel configuration"

# conda channel 镜像：NFCORE_CONDA_MIRROR=tuna(默认) | none
#   none = 不改动 channel，适合境外网络或已配好内网自建源的环境
#
# ⚠ channel 顺序应为 conda-forge > bioconda > defaults（bioconda 官方建议），
#   同时开 strict 优先级避免跨 channel 混装。
#
#   注意 nf-core 启动时那条 "There is a problem with your Conda configuration"
#   的 WARN：它执行 `conda config --show channels`，而 conda 会把系统级、
#   用户级 ~/.condarc 和 CONDARC 指向的文件【合并】起来看，所以即使这里写对了，
#   只要用户的 ~/.condarc 里 bioconda 在前，那条 WARN 依然会出现。
#   实际建环境用的是各 module 自带的 environment.yml（里面自己声明了 channels），
#   因此这条 WARN 多数情况下只是噪音，不是失败原因——真正的失败原因请看
#   "Failed to create Conda environment ... status : 120"（超时，见下方 createTimeout）。
CONDA_RC="${OUTPUT_DIR}/.condarc_nfcore"

# 镜像选择：tuna(默认) | bfsu | nju | sjtu | none | <自建源地址>
case "${NFCORE_CONDA_MIRROR:-tuna}" in
    tuna) MIRROR_HOST="https://mirrors.tuna.tsinghua.edu.cn/anaconda" ;;
    bfsu) MIRROR_HOST="https://mirrors.bfsu.edu.cn/anaconda" ;;
    nju)  MIRROR_HOST="https://mirror.nju.edu.cn/anaconda" ;;
    sjtu) MIRROR_HOST="https://mirror.sjtu.edu.cn/anaconda" ;;
    none) MIRROR_HOST="" ;;
    *)    MIRROR_HOST="${NFCORE_CONDA_MIRROR}" ;;
esac

# 基础部分：channel 顺序 + 下载韧性
# （慢链路上 conda 默认超时/重试太保守，断一次整个环境就失败：
#   curl: (18) end of response with N bytes missing）
cat > "${CONDA_RC}" << 'EOF'
channels:
  - conda-forge
  - bioconda
  - defaults
channel_priority: strict
show_channel_urls: true
remote_connect_timeout_secs: 30.0
remote_read_timeout_secs: 600.0
remote_max_retries: 5
remote_backoff_factor: 2
EOF

if [[ -z "${MIRROR_HOST}" ]]; then
    echo "[SKIP] NFCORE_CONDA_MIRROR=none; leaving channel URLs unchanged"
else
    cat >> "${CONDA_RC}" << EOF
default_channels:
  - ${MIRROR_HOST}/pkgs/main
  - ${MIRROR_HOST}/pkgs/r
custom_channels:
  conda-forge: ${MIRROR_HOST}/cloud
  bioconda: ${MIRROR_HOST}/cloud
EOF
    echo "[DONE] Mirror: ${NFCORE_CONDA_MIRROR:-tuna} (${MIRROR_HOST})"
    echo "       If downloads are slow or interrupted, try: NFCORE_CONDA_MIRROR=bfsu|nju|sjtu|none"
fi
export CONDARC="${CONDA_RC}"
echo "[  OK ] Channel priority conda-forge > bioconda > defaults (strict)"
echo "[  OK ] Download retries: up to 5; read timeout: 600s"

echo ""

# ─────────────────────────────────────────────
# 构建运行命令
# ─────────────────────────────────────────────
# ─────────────────────────────────────────────
# 已知依赖冲突的修补配置（conda 模式专属）
# ─────────────────────────────────────────────
# nf-core 的 conda spec 只锁工具、不锁依赖，解析时会拉到 numpy 2.x，
# 而 deeptools 3.5.x 仍用 numpy 2.0 已移除的 np.NAN → COMPUTEMATRIX 必崩。
# 这里生成一份 nextflow 配置为受影响 process 补上 numpy<2 约束。
# ── 1) 清理上次超时留下的半成品环境（否则会被当成「已建好」而复用）──
clean_partial_conda_envs "${NXF_CONDA_CACHEDIR}"

# ── 2) 生成 conda 模式补充配置 ──
#    · createTimeout：默认 20 分钟建不完 R/Bioconductor 环境 → status 120
#    · deeptools 依赖锁：规避 numpy 2.x 的 np.NAN 崩溃
PATCH_ARGS=()
PATCH_CFG="$(write_nfcore_conda_patch "${NFCORE_NAME}" \
             "$([[ "${NO_CONDA_PATCH}" == "true" ]] && echo false || echo true)" || true)"
if [[ -s "${PATCH_CFG:-}" ]]; then
    PATCH_ARGS+=(-c "${PATCH_CFG}")
    echo "[  OK ] Conda environment creation timeout: ${NFCORE_CONDA_TIMEOUT:-3 h}(increased from the default timeout)"
    if [[ "${NO_CONDA_PATCH}" == "true" ]]; then
        echo "[SKIP] --no-conda-patch: dependency overrides are disabled"
    else
        _dt="$(grep -oP 'bioconda::deeptools=?\K[0-9.]*' "${PATCH_CFG}" | head -1)"
        echo "[  OK ] Dependency overrides: deeptools${_dt:+=${_dt}} + numpy=1.26 (avoids the removed np.NAN alias)"
        [[ -z "${_dt}" ]] && log_warn "Could not read the pinned deepTools version from assets; leaving its version unconstrained"
    fi
    echo "         Configuration: ${PATCH_CFG}"
fi
echo ""

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
    echo "[  OK ] MACS_GSIZE: ${MACS_GSIZE} (Source: $(value_origin MACS_GSIZE)$(
        m="$(derived_get MACS_GSIZE_METHOD 2>/dev/null || true)"; [[ -n "${m}" ]] && echo ", ${m}"))"
    # read_length only drives nf-core's fallback size estimator. Its schema
    # excludes lengths such as 130 even when an explicit macs_gsize is supplied.
    [[ -z "${READ_LENGTH}" ]] || echo "[INFO ] READ_LENGTH=${READ_LENGTH} retained in ChromaFlow; nf-core uses the supplied MACS_GSIZE"
elif [[ -n "${READ_LENGTH}" ]]; then
    GSIZE_ARGS+=(--read_length "${READ_LENGTH}")
    echo "[INFO ] MACS_GSIZE is unset; nf-core will estimate it using --read_length=${READ_LENGTH}"
else
    die "Both MACS_GSIZE and READ_LENGTH are unset. Run chromaflow reference,
     or chromaflow samplesheet, or explicitly configure one of these values"
fi

COMMON_ARGS=(
    --input "${SAMPLESHEET}"
    --outdir "${OUTDIR}"
    --fasta "${GENOME_COMBINED}"
    --gtf "${GENOME_GTF}"
    --max_cpus "${MAX_CPUS}"
    --max_memory "${MAX_MEMORY}"
    --max_time "${MAX_TIME}"
    -profile "${PROFILE}"
    -w "${WORK_DIR}"
)
COMMON_ARGS+=(${GSIZE_ARGS[@]+"${GSIZE_ARGS[@]}"})
COMMON_ARGS+=(${PATCH_ARGS[@]+"${PATCH_ARGS[@]}"})
COMMON_ARGS+=(${NFCORE_VERSION_ARGS[@]+"${NFCORE_VERSION_ARGS[@]}"})

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
echo " Run command"
echo "============================================"
echo ""
print_nfcore_command "${COMMON_ARGS[@]}" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} ${USER_ARGS[@]+"${USER_ARGS[@]}"} -resume
echo ""
echo "Environment variables:"
echo "  NXF_CONDA_CACHEDIR=${NXF_CONDA_CACHEDIR}"
echo "  Environment solver: $(command -v mamba &>/dev/null && echo 'mamba (conda.useMamba)' || echo 'conda')"
echo "  CONDARC=${CONDARC}"
echo ""
if [[ "${CHROMAFLOW_DRY_RUN:-false}" == true ]]; then
    echo "[DRY-RUN] Checks completed; analysis was not started."
    exit 0
fi

# ─────────────────────────────────────────────
# 确认执行
# ─────────────────────────────────────────────
confirm=y
if [[ "${CHROMAFLOW_YES:-false}" != true ]]; then read -r -p "Start analysis? [y/N] " confirm; fi
if [[ "${confirm}" != "y" && "${confirm}" != "Y" ]]; then
    echo "[ABORT] Cancelled"
    exit 0
fi

echo ""
echo "[START] $(date)"
echo ""

cd "${OUTPUT_DIR}"

execute_nfcore \
    "${COMMON_ARGS[@]}" \
    ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} \
    ${USER_ARGS[@]+"${USER_ARGS[@]}"} \
    -resume

echo ""
echo "[DONE] nf-core/${NFCORE_NAME} (${PIPELINE}, conda mode) completed $(date)"
echo "  Results: ${OUTDIR}"
