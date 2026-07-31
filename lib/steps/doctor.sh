#!/usr/bin/env bash
# ============================================================================
# doctor.sh — 运行前自检
# ============================================================================
# 用法: chromaflow doctor [-c <config>] [--type chipseq|atacseq|dapseq]
# 不修改任何文件；退出码 0=无致命问题, 1=有 FAIL。
# 逐项检查：软件完整性 / 配置 / 当前环境工具 / 参考 / 输入 / 各步产物就绪度。
# 对「用新版跑旧数据」特别有用：确认 OUTPUT_DIR 下能否找到 peaks/bigWig。
# ============================================================================
set -uo pipefail   # 故意不加 -e：即使某项检查失败也要继续跑完
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"

PIPELINE="${1:-chipseq}"

PASS=0; WARN=0; FAIL=0
ok()   { echo "  [ OK ] $*"; PASS=$((PASS+1)); }
info() { echo "  [INFO] $*"; }
skip() { echo "  [SKIP] $*"; }
warn() { echo "  [WARN] $*"; WARN=$((WARN+1)); }
bad()  { echo "  [FAIL] $*"; FAIL=$((FAIL+1)); }
have() { command -v "$1" &>/dev/null; }
# 检查一组工具，缺失则提示应激活的环境
need_tools() { # need_tools <env> <tool>...
    local env="$1"; shift
    local miss=()
    for t in "$@"; do have "$t" || miss+=("$t"); done
    if [[ ${#miss[@]} -eq 0 ]]; then ok "工具齐全 (${env}): $*"
    else warn "缺 ${miss[*]}（跑该步前 mamba activate ${env}）"; fi
}

echo ""
echo "============================================"
echo " chromaflow doctor  (pipeline: ${PIPELINE})"
echo " $(date)"
echo "============================================"

# ─────────────── 1. 软件完整性 ───────────────
echo ""; echo "[1] 软件本体"
ok "CHROMAFLOW_HOME = ${CHROMAFLOW_HOME}"
for f in lib/common.sh lib/help.sh bin/chromaflow \
         lib/steps/00_prepare_reference.sh lib/steps/normalize_gff.py \
         lib/steps/make_samplesheet.py lib/steps/05_peak_annotation.sh \
         lib/steps/06_signal_deeptools.sh; do
    [[ -f "${CHROMAFLOW_HOME}/${f}" ]] && ok "存在 ${f}" || bad "缺失 ${f}"
done

# ─────────────── 2. 配置 ───────────────
echo ""; echo "[2] 项目配置"
CONFIG_OK=0
if try_load_config "${CHROMAFLOW_CONFIG:-}"; then
    CONFIG_OK=1
    ok "已加载: ${CHROMAFLOW_CONFIG}"
    [[ -n "${SPECIES_SHORT:-}" ]] && ok "SPECIES_SHORT = ${SPECIES_SHORT}" || bad "SPECIES_SHORT 未设置"
    echo "         INPUT_DIR  = ${INPUT_DIR}"
    echo "         OUTPUT_DIR = ${OUTPUT_DIR}"
    echo "         REFERENCE  = ${REFERENCE_DIR}"
    # 可写性
    if mkdir -p "${OUTPUT_DIR}" 2>/dev/null && [[ -w "${OUTPUT_DIR}" ]]; then
        ok "OUTPUT_DIR 可写"
    else
        bad "OUTPUT_DIR 不可写: ${OUTPUT_DIR}"
    fi
else
    warn "未找到配置（-c <config> 或 ./chromaflow.config）；仅做软件/环境检查"
fi

# ─────────────── 3. 当前环境与工具 ───────────────
echo ""; echo "[3] 当前环境工具 (激活环境: ${CONDA_DEFAULT_ENV:-未知})"
echo "    （按步骤列出；只需在跑某步前，其对应环境的工具齐全即可）"
echo "  -- reference (00) --";   need_tools ChromaFlow_ref        samtools gffread bedtools
echo "  -- run (02) --"
if have nextflow; then ok "nextflow 就绪"; else warn "缺 nextflow（跑 02 前 mamba activate nfcore_env）"; fi

# 执行模式：NFCORE_PROFILE 同时决定 nextflow -profile 与 chromaflow run 的执行路径
_PROF="$(resolve_nfcore_profile)"; _KIND="$(nfcore_profile_kind "${_PROF}")"
info "执行模式 NFCORE_PROFILE=${_PROF} → $( [[ "${_KIND}" == conda ]] && echo 'conda 模式 (02_run_nfcore_conda.sh)' || echo '容器模式 (02_run_nfcore.sh)' )"
case "${_KIND}" in
    container)
        if have apptainer || have singularity || have docker || have podman; then
            ok "容器引擎就绪（镜像缓存: ${NXF_SINGULARITY_CACHEDIR}）"
        else
            warn "profile=${_PROF} 需容器引擎但未检测到；可改用 chromaflow run --profile conda"
        fi ;;
    conda)
        if have mamba || have conda; then ok "conda/mamba 就绪（nf-core 将自建各 process 环境）"
        else warn "profile=${_PROF} 需 conda/mamba 但未检测到"; fi ;;
    *)  info "profile=${_PROF} 非内置引擎名，将原样传给 nextflow" ;;
esac
echo "  -- idr (03) --";        need_tools ChromaFlow_idr        idr bedtools
echo "  -- motif (04) --";      need_tools ChromaFlow_motif      meme-chip bedtools
echo "  -- annotate (05, ChIPseeker/R) --"
if have Rscript; then
    if Rscript -e 'quit(status=!requireNamespace("ChIPseeker",quietly=TRUE))' 2>/dev/null; then
        ok "R + ChIPseeker 就绪 (05)"
    else warn "有 Rscript 但缺 ChIPseeker（mamba activate ChromaFlow_R）"; fi
else warn "缺 Rscript（05 需 ChromaFlow_R；或改用 annotate-py）"; fi
echo "  -- annotate-py (05 备选) --"; need_tools ChromaFlow_ref bedtools python3
echo "  -- candidates (06 整合) --"
need_tools ChromaFlow_downstream pyGenomeTracks
have python3 && ok "python3 可用 (读04 FIMO + 合并)" || warn "缺 python3"
echo "  -- signal (06) --";     need_tools ChromaFlow_downstream computeMatrix plotHeatmap plotProfile
echo "  -- tracks (07) --";     need_tools ChromaFlow_downstream pyGenomeTracks
echo "  -- go-plot (08) --"
have python3 && python3 -c "import matplotlib" 2>/dev/null \
    && ok "matplotlib 可用 (go-plot 出图)" || warn "无 matplotlib：go-plot 无法出图"

# 无配置则到此为止
if [[ ${CONFIG_OK} -eq 0 ]]; then
    echo ""; echo "[跳过 4-6：需要配置] 请提供 chromaflow.config 后重跑以检查参考/输入/产物"
    echo ""; echo "结果: ${PASS} OK / ${WARN} WARN / ${FAIL} FAIL"
    [[ ${FAIL} -eq 0 ]] && exit 0 || exit 1
fi

# ─────────────── 4. 参考基因组 ───────────────
echo ""; echo "[4] 参考基因组（路径由 config > manifest > 目录扫描 自动解析）"
chk_file() { [[ -f "$1" ]] && ok "$2: $1" || warn "$2 缺失: $1"; }
chk_file "${GENOME_COMBINED}" "基因组 FASTA"
chk_file "${GENOME_GTF}"      "GTF (02 nf-core 用)"
chk_file "${CHROM_SIZES}"     "chrom.sizes"
if [[ -f "${GENOME_GFF3}" ]]; then
    ok "GFF3 (05 注释用): ${GENOME_GFF3}"
    if grep -qP '\texon\t' "${GENOME_GFF3}" 2>/dev/null; then
        ok "GFF3 含 exon 行（UTR/exon/intron 判定可靠）"
    else
        warn "GFF3 无 exon 行 → 建议 chromaflow normalize-gff（可选步骤，无 AGAT 会用内置引擎）"
    fi
else
    warn "GFF3 缺失: ${GENOME_GFF3}（只做 peak calling 可忽略）"
fi

echo ""; echo "[4c] 实测量（自动测得，不写回 config）"
for _k in MACS_GSIZE READ_LENGTH GENOME_GC; do
    _v="$(eval "echo \"\${${_k}:-}\"")"
    [[ -z "${_v}" ]] && _v="$(derived_get "${_k}" 2>/dev/null || true)"
    if [[ -n "${_v}" ]]; then
        ok "${_k} = ${_v}  (来源: $(value_origin "${_k}"))"
    else
        case "${_k}" in
            MACS_GSIZE)  info "MACS_GSIZE 未测：跑 chromaflow reference 会用 khmer 实测；缺失时 02 交给 nf-core 计算" ;;
            READ_LENGTH) info "READ_LENGTH 未测：跑 chromaflow samplesheet 会从 fastq 抽样实测" ;;
            *)           info "${_k} 未测（仅供参考，不影响分析）" ;;
        esac
    fi
done

MF="$(reference_manifest_path)"
if [[ -f "${MF}" ]]; then
    info "产物登记表: ${MF}"
    sed 's/^/         /' "${MF}"
else
    info "尚无产物登记表（还没跑过 chromaflow reference / normalize-gff）"
fi

# ── 可选功能：缺了也不影响主流程，一律 INFO 不计入 WARN ──
echo ""; echo "[4b] 可选功能"
if [[ "${SKIP_ORGANELLE:-false}" == "true" || ${#ORGANELLE_SOURCES[@]} -eq 0 ]]; then
    skip "细胞器合并未启用（无叶绿体/线粒体基因组时属正常）"
else
    for _n in $(organelle_names); do
        _u=$(echo "${_n}" | tr '[:lower:]' '[:upper:]'); _f="$(eval echo "\${GENOME_${_u}:-}")"
        [[ -s "${_f}" ]] && ok "细胞器 ${_n}: ${_f}" \
            || info "细胞器 ${_n} 未就绪；chromaflow reference 会跳过它继续（--skip-organelle 可永久跳过）"
    done
fi
if have agat_convert_sp_gxf2gxf.pl; then
    ok "AGAT 可用 (normalize-gff --engine agat)"
else
    info "无 AGAT —— normalize-gff 将自动使用内置 Python 引擎（无需安装任何东西）"
    info "  需要 AGAT 时单独建环境: mamba create -n ChromaFlow_agat -c bioconda agat"
fi
if [[ -n "${ANNOTATION_INFO:-}" && -f "${ANNOTATION_INFO}" ]]; then ok "annotation_info (GO 用)"
else info "无 ANNOTATION_INFO：05 不产 gene2go（GO 富集受限，其余不受影响）"; fi
if [[ -n "${MEME_DB:-}" && -f "${MEME_DB}" ]]; then ok "MEME_DB (04 Tomtom 比对)"
else info "无 MEME_DB：04 跳过已知 motif 比对，de novo 发现照常"; fi

# ─────────────── 5. 输入 / samplesheet ───────────────
echo ""; echo "[5] 输入与 samplesheet"
SAMPLES_TXT="$(resolve_samples_file "${PIPELINE}" || true)"
[[ -f "${SAMPLES_TXT}" ]] && ok "样本表: ${SAMPLES_TXT}" \
    || warn "未找到样本表（找过 samples_${PIPELINE}.{txt,tsv,csv} / samples.*；仅跑下游可忽略）"
SS="${SAMPLESHEET_DIR}/samplesheet_${PIPELINE}.csv"
[[ -f "${SS}" ]] && ok "samplesheet: $(basename "${SS}")" || warn "无 samplesheet（跑 02 前需先 samplesheet）"

# ─────────────── 6. 各步产物就绪度 ───────────────
echo ""; echo "[6] 产物就绪度 (决定下游能否运行)"
PEAK_DIR=$(get_nfcore_peak_dir "${PIPELINE}")
n_np=$(find "${PEAK_DIR}" -maxdepth 1 -name "*_peaks.narrowPeak" 2>/dev/null | wc -l)
if [[ ${n_np} -gt 0 ]]; then ok "nf-core narrowPeak: ${n_np} 个 (${PEAK_DIR})"
else warn "未见 narrowPeak（02 未跑？或 OUTPUT_DIR 指错）: ${PEAK_DIR}"; fi

REC="${RESULTS_DIR}/downstream/${PIPELINE}/idr/recommended_peaks.bed"
CONS=$(get_consensus_peaks "${PIPELINE}")
if [[ -s "${REC}" ]]; then ok "recommended_peaks.bed 就绪 → 04/05/06 可跑 ($(wc -l < "${REC}") peaks)"
elif [[ -f "${CONS}" ]]; then warn "无 recommended（03 未跑）；将退回 consensus: $(basename "${CONS}")"
else warn "既无 recommended 也无 consensus → 04/05/06 暂无输入。先跑 idr，或检查 OUTPUT_DIR 是否指对了旧结果"; fi

BW_DIR="${RESULTS_DIR}/${PIPELINE}/bwa/merged_library/bigwig"
n_bw=$(find "${BW_DIR}" -maxdepth 1 \( -name "*.bigWig" -o -name "*.bw" \) 2>/dev/null | wc -l)
if [[ ${n_bw} -gt 0 ]]; then ok "bigWig: ${n_bw} 个 → 06 signal 可跑"
else warn "未见 bigWig（06 signal 需要）: ${BW_DIR}"; fi

# ─────────────── 汇总 ───────────────
echo ""
echo "============================================"
echo " 结果: ${PASS} OK / ${WARN} WARN / ${FAIL} FAIL"
if [[ ${FAIL} -eq 0 ]]; then
    echo " ✅ 无致命问题"
    [[ ${WARN} -gt 0 ]] && echo "    （WARN 多为「该步尚未运行」的正常提示，按流程推进即可）"
else
    echo " ❌ 有 ${FAIL} 项需处理，见上方 [FAIL]"
fi
echo "============================================"
[[ ${FAIL} -eq 0 ]] && exit 0 || exit 1
