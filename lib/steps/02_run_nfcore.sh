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
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"
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
        echo "[ERROR] 未知 pipeline: ${PIPELINE}"
        echo "  支持: chipseq, atacseq, dapseq"
        exit 1
        ;;
esac

SAMPLESHEET="${SAMPLESHEET_DIR}/samplesheet_${PIPELINE}.csv"
OUTDIR="${RESULTS_DIR}/${PIPELINE}"

open_log "02_run_nfcore" "${PIPELINE}"

log_section "Run nf-core/${NFCORE_NAME} v${NFCORE_VERSION} (${PIPELINE}) [profile=${PROFILE}]"

# ─────────────────────────────────────────────
# 前置检查
# ─────────────────────────────────────────────
echo "[CHECK] 环境检查"

check_tool nextflow || exit 1
echo "[  OK ] nextflow $(nextflow -version 2>&1 | grep -oP 'version \K\S+')"

# 按 profile 检查对应的容器引擎（singularity 与 apptainer 可互相替代）
case "${ENGINE}" in
    singularity|apptainer)
        if command -v apptainer &>/dev/null; then
            echo "[  OK ] apptainer $(apptainer --version)"
        elif command -v singularity &>/dev/null; then
            echo "[  OK ] singularity $(singularity --version)"
        else
            log_error "profile=${PROFILE} 需要 Apptainer/Singularity，但未检测到"
            echo "  改用 conda 模式: chromaflow run --profile conda"
            echo "  或在 chromaflow.config 里设 NFCORE_PROFILE=conda"
            exit 1
        fi
        echo "[  OK ] 镜像缓存: ${NXF_SINGULARITY_CACHEDIR}"
        ;;
    docker|podman)
        check_tool "${ENGINE}" || { echo "  改用: chromaflow run --profile conda"; exit 1; }
        echo "[  OK ] ${ENGINE} $(${ENGINE} --version 2>&1 | head -1)"
        ;;
    *)
        log_warn "profile=${PROFILE} 不含已知容器引擎，按原样交给 nextflow"
        ;;
esac


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
    -profile "${NFCORE_PROFILE}"
    -r "${NFCORE_VERSION}"
    -w "${WORK_DIR}"
)
COMMON_ARGS+=(${GSIZE_ARGS[@]+"${GSIZE_ARGS[@]}"})

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
echo " 运行命令"
echo "============================================"
echo ""
echo "nextflow run nf-core/${NFCORE_NAME} \\"
for arg in "${COMMON_ARGS[@]}"; do echo "  ${arg} \\"; done
for arg in "${EXTRA_ARGS[@]}"; do echo "  ${arg} \\"; done
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

# 在 OUTPUT_DIR 下运行，使 work/.nextflow/日志 全部落在输出目录，代码目录零污染
cd "${OUTPUT_DIR}"

nextflow run "nf-core/${NFCORE_NAME}" \
    "${COMMON_ARGS[@]}" \
    "${EXTRA_ARGS[@]}" \
    -resume

echo ""
echo "[DONE] nf-core/${NFCORE_NAME} (${PIPELINE}) 完成 $(date)"
echo "  结果: ${OUTDIR}"
