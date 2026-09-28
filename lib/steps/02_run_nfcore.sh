#!/usr/bin/env bash
# ============================================================================
# 02_run_nfcore.sh — 运行 nf-core pipeline（容器模式）
# 用法: chromaflow run --type [chipseq|atacseq|dapseq] [--profile singularity]
# ============================================================================
# 【执行模式由 NFCORE_PROFILE 控制】
#   config 里的 NFCORE_PROFILE（或命令行 --profile）决定 nextflow 的 -profile，
#   同时决定 `chromaflow run` 走哪个脚本：
#     singularity (默认) / apptainer / docker / podman → 本脚本（容器模式）
#     conda / mamba                                    → 02_run_nfcore_conda.sh
#   优先级: --profile > config 的 NFCORE_PROFILE > singularity
#   也支持组合，如 NFCORE_PROFILE="singularity,test"。
#
# 容器模式下所有生信工具由 nextflow 自动拉取镜像，可重复性最好；
# 镜像缓存在 ${OUTPUT_DIR}/.singularity_cache（NXF_SINGULARITY_CACHEDIR）。
# ============================================================================
# 【重要提示：强制重跑与缓存清理】
# 本脚本默认携带 -resume 参数，以支持断点续传并节约计算资源。
# 若你更改了输入数据，或缓存发生损坏，需要彻底从零开始运行：
# 请在执行本脚本前，确保在当前工作目录下手动执行以下命令清空缓存：
#   rm -rf "${OUTPUT_DIR}"/work "${OUTPUT_DIR}"/.nextflow "${OUTPUT_DIR}"/.nextflow.log*
# （警告：这将清除所有已有的中间步骤运行记录，请谨慎操作）
# ============================================================================

set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true

# 透传给 nextflow 的额外参数：--skip-spp 等快捷开关，或 -- 之后的任意原样参数
#   chromaflow run --skip-spp
#   chromaflow run -- --skip_preseq --skip_igv -resume
USER_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-spp|--skip_spp) USER_ARGS+=(--skip_spp); shift ;;
        --) shift; while [[ $# -gt 0 ]]; do USER_ARGS+=("$1"); shift; done ;;
        *) shift ;;
    esac
done
# config 里也可以设 SKIP_SPP=true 长期关闭
if [[ "${SKIP_SPP:-false}" == "true" ]]; then
    case " ${USER_ARGS[*]} " in *" --skip_spp "*) ;; *) USER_ARGS+=(--skip_spp) ;; esac
fi
PROFILE="$(resolve_nfcore_profile)"
ENGINE="$(nfcore_container_engine "${PROFILE}")"

# ─────────────────────────────────────────────
# 确定 nf-core pipeline 和版本
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
        echo "[ERROR] Unknown pipeline: ${PIPELINE}"
        echo "  Supported: chipseq, atacseq, dapseq"
        exit 1
        ;;
esac

SAMPLESHEET="${SAMPLESHEET_DIR}/samplesheet_${PIPELINE}.csv"
OUTDIR="${RESULTS_DIR}/${PIPELINE}"

open_log "02_run_nfcore" "${PIPELINE}"

log_section "Run nf-core/${NFCORE_NAME} v${NFCORE_VERSION} (${PIPELINE}) [profile=${PROFILE}]"

source "${CHROMAFLOW_HOME}/lib/nfcore_runtime.sh"
configure_nfcore_runtime

# ─────────────────────────────────────────────
# 前置检查
# ─────────────────────────────────────────────
echo "[CHECK] Environment check"

check_tool "${NFCORE_NEXTFLOW}" || exit 1
echo "[  OK ] nextflow $("${NFCORE_NEXTFLOW}" -version 2>&1 | grep -oP 'version \K\S+')"

# 按 profile 检查对应的容器引擎（singularity 与 apptainer 可互相替代）
case "${ENGINE}" in
    singularity|apptainer)
        if command -v "${ENGINE}" &>/dev/null; then
            echo "[  OK ] ${ENGINE} $("${ENGINE}" --version)"
        else
            log_error "profile=${PROFILE} requires Apptainer/Singularity, but neither was detected"
            echo "  Switch to online Conda mode: chromaflow run --profile conda"
            echo "  Or set NFCORE_PROFILE=conda in chromaflow.config"
            exit 1
        fi
        echo "[  OK ] Image cache: ${NXF_SINGULARITY_CACHEDIR}"
        ;;
    docker|podman)
        check_tool "${ENGINE}" || { echo "  Use: chromaflow run --profile conda"; exit 1; }
        echo "[  OK ] ${ENGINE} $(${ENGINE} --version 2>&1 | head -1)"
        ;;
    *)
        log_warn "profile=${PROFILE} has no recognized container engine; passing it unchanged to Nextflow"
        ;;
esac


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
    echo "[  OK ] MACS_GSIZE: ${MACS_GSIZE} (Source: $(value_origin MACS_GSIZE)$(
        m="$(derived_get MACS_GSIZE_METHOD 2>/dev/null || true)"; [[ -n "${m}" ]] && echo ", ${m}")）"
elif [[ -n "${READ_LENGTH}" ]]; then
    echo "[INFO ] MACS_GSIZE is unset; nf-core will estimate it using --read_length=${READ_LENGTH}"
else
    die "Both MACS_GSIZE and READ_LENGTH are unset. Run chromaflow reference,
     or chromaflow samplesheet, or explicitly configure one of these values"
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
    -w "${WORK_DIR}"
)
COMMON_ARGS+=(${GSIZE_ARGS[@]+"${GSIZE_ARGS[@]}"})
COMMON_ARGS+=(${NFCORE_VERSION_ARGS[@]+"${NFCORE_VERSION_ARGS[@]}"})

# pipeline 特定参数
EXTRA_ARGS=()
case "${PIPELINE}" in
    chipseq|dapseq)
        EXTRA_ARGS+=(--narrow_peak)
        # 如有 blacklist
        [[ -f "${BLACKLIST_BED}" ]] && EXTRA_ARGS+=(--blacklist "${BLACKLIST_BED}")
        ;;
    atacseq)
        # atacseq 默认参数即可
        [[ -f "${BLACKLIST_BED}" ]] && EXTRA_ARGS+=(--blacklist "${BLACKLIST_BED}")
        ;;
esac

echo "============================================"
echo " Run command"
echo "============================================"
echo ""
print_nfcore_command "${COMMON_ARGS[@]}" ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} ${USER_ARGS[@]+"${USER_ARGS[@]}"} -resume
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

# 在 OUTPUT_DIR 下运行，使 work/.nextflow/日志 全部落在输出目录，代码目录零污染
cd "${OUTPUT_DIR}"

execute_nfcore \
    "${COMMON_ARGS[@]}" \
    ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"} \
    ${USER_ARGS[@]+"${USER_ARGS[@]}"} \
    -resume

echo ""
echo "[DONE] nf-core/${NFCORE_NAME} (${PIPELINE}) completed $(date)"
echo "  Results: ${OUTDIR}"
