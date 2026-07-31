#!/usr/bin/env bash
# ============================================================================
# 01_make_samplesheet.sh — 生成 nf-core samplesheet
# ============================================================================
# 用法: chromaflow samplesheet [--type chipseq|atacseq|dapseq]
#                             [--samples FILE] [--pair-pattern '_good_{R}.fq.gz']
#                             [--single-end] [--list-fastq] [--strict]
#
# 样本表定位顺序: --samples > config 的 SAMPLES_FILE >
#                 <PROJECT_DIR|INPUT_DIR>/samples_<type>.{txt,tsv,csv} > samples.{txt,tsv,csv}
# 样本表格式（表头式 / 位置式 / 直接给 fastq_1,fastq_2）见 templates/samples_*.txt
# ============================================================================

set -euo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true

SAMPLES_ARG=""; PY_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples)       SAMPLES_ARG="$2"; shift 2 ;;
        --pair-pattern)  PY_ARGS+=(--pair-pattern "$2"); shift 2 ;;
        --single-end)    PY_ARGS+=(--single-end); shift ;;
        --list-fastq)    PY_ARGS+=(--list-fastq); shift ;;
        --strict)        PY_ARGS+=(--strict); shift ;;
        *) shift ;;
    esac
done

# DAP-seq 使用 chipseq pipeline
NFCORE_PIPELINE="${PIPELINE}"
[[ "${PIPELINE}" == "dapseq" ]] && NFCORE_PIPELINE="chipseq"

# config 里的全局配对设置（命令行优先，故放在 PY_ARGS 之前）
GLOBAL_ARGS=()
[[ -n "${FASTQ_PAIR_PATTERN:-}" ]] && GLOBAL_ARGS+=(--pair-pattern "${FASTQ_PAIR_PATTERN}")
[[ "${FASTQ_SINGLE_END:-false}" == "true" ]] && GLOBAL_ARGS+=(--single-end)

SAMPLES_TXT="$(resolve_samples_file "${PIPELINE}" "${SAMPLES_ARG}")" || true
OUTPUT_CSV="$(step_outdir "${SAMPLESHEET_DIR}")/samplesheet_${PIPELINE}.csv"

open_log "01_make_samplesheet" "${PIPELINE}"
log_section "Generate nf-core/${NFCORE_PIPELINE} samplesheet (${PIPELINE})"

echo "  样本表: ${SAMPLES_TXT}"
echo "  输出:   ${OUTPUT_CSV}"
echo ""

if [[ ! -f "${SAMPLES_TXT}" ]]; then
    log_error "找不到样本表: ${SAMPLES_TXT}"
    echo "  可用 --samples <file> 指定，或在 config 里设 SAMPLES_FILE，"
    echo "  或把样本表放在 ${PROJECT_DIR}/samples_${PIPELINE}.txt"
    echo "  模板: ${CHROMAFLOW_HOME}/templates/samples_${PIPELINE}.txt"
    exit 1
fi

mkdir -p "$(dirname "${OUTPUT_CSV}")"
python3 "${CHROMAFLOW_HOME}/lib/steps/make_samplesheet.py" \
    "${SAMPLES_TXT}" -o "${OUTPUT_CSV}" --pipeline "${NFCORE_PIPELINE}" \
    ${GLOBAL_ARGS[@]+"${GLOBAL_ARGS[@]}"} ${PY_ARGS[@]+"${PY_ARGS[@]}"}

# ── 从 fastq 实测读长 ──
# nf-core 只在缺 --macs_gsize 时用 read_length 推算基因组大小，但实测一次几乎
# 不花时间，且能让 00 的 k-mer 口径与真实数据一致。结果写进 derived.env，
# 不改 config。config 里显式写了 READ_LENGTH 就以你的为准。
if [[ "$(value_origin READ_LENGTH)" == "config" ]]; then
    echo ""
    echo "  读长: ${READ_LENGTH}（config 指定，不再实测）"
elif [[ -f "${OUTPUT_CSV}" ]]; then
    echo ""
    echo "[实测] 从 fastq 抽样测读长..."
    MAX_LEN=0
    while IFS=, read -r _s fq1 _rest; do
        [[ "${fq1}" == "fastq_1" || -z "${fq1}" ]] && continue
        L="$(detect_read_length "${fq1}" 2000 || true)"
        [[ -n "${L}" && "${L}" -gt "${MAX_LEN}" ]] && MAX_LEN="${L}"
    done < "${OUTPUT_CSV}"
    if [[ "${MAX_LEN}" -gt 0 ]]; then
        derived_set READ_LENGTH "${MAX_LEN}"
        echo "  [DONE] READ_LENGTH = ${MAX_LEN} bp（取抽样最大值；已记入 derived.env）"
        OLD_K="$(derived_get MACS_GSIZE_SIG 2>/dev/null || true)"
        if [[ -n "${OLD_K}" && "${OLD_K}" != *"k=${MAX_LEN}" ]]; then
            echo "  提示: 已有的 MACS_GSIZE 是按 ${OLD_K##*k=} 的 k 测的；"
            echo "        影响通常在 1% 内，若要严格对齐可跑 chromaflow reference --force"
        fi
    else
        log_warn "读长实测失败（fastq 不可读？）；02 会在需要时回退到 nf-core 自行计算"
    fi
fi

echo ""
echo "下一步: mamba activate nfcore_env && chromaflow run --type ${PIPELINE}"
