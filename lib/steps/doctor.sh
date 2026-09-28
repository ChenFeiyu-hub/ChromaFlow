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
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"

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
    if [[ ${#miss[@]} -eq 0 ]]; then ok "Tools available (${env}): $*"
    else warn "Missing ${miss[*]} (activate ${env} before running that step)"; fi
}

echo ""
echo "============================================"
echo " chromaflow doctor  (pipeline: ${PIPELINE})"
echo " $(date)"
echo "============================================"

# ─────────────── 1. 软件完整性 ───────────────
echo ""; echo "[1] Software installation"
ok "CHROMAFLOW_HOME = ${CHROMAFLOW_HOME}"
for f in lib/common.sh lib/help.sh bin/chromaflow \
         lib/steps/00_prepare_reference.sh lib/steps/normalize_gff.py \
         lib/steps/make_samplesheet.py lib/steps/05_peak_annotation.sh \
         lib/steps/06_signal_deeptools.sh; do
    [[ -f "${CHROMAFLOW_HOME}/${f}" ]] && ok "Found ${f}" || bad "Missing ${f}"
done

# ─────────────── 2. 配置 ───────────────
echo ""; echo "[2] Project configuration"
CONFIG_OK=0
if try_load_config "${CHROMAFLOW_CONFIG:-}"; then
    CONFIG_OK=1
    ok "Loaded: ${CHROMAFLOW_CONFIG}"
    [[ -n "${SPECIES_SHORT:-}" ]] && ok "SPECIES_SHORT = ${SPECIES_SHORT}" || bad "SPECIES_SHORT unset"
    echo "         INPUT_DIR  = ${INPUT_DIR}"
    echo "         OUTPUT_DIR = ${OUTPUT_DIR}"
    echo "         REFERENCE  = ${REFERENCE_DIR}"
    # 可写性
    if mkdir -p "${OUTPUT_DIR}" 2>/dev/null && [[ -w "${OUTPUT_DIR}" ]]; then
        ok "OUTPUT_DIR is writable"
    else
        bad "OUTPUT_DIR is not writable: ${OUTPUT_DIR}"
    fi
else
    warn "Config not found (-c <config> or ./chromaflow.config); checking software/environment only"
fi

# ─────────────── 3. 当前环境与工具 ───────────────
echo ""; echo "[3] Tools on PATH (active environment: ${CONDA_DEFAULT_ENV:-unknown})"
echo "    (Each step only requires the tools in its own environment.)"
echo "  -- reference (00) --";   need_tools ChromaFlow_ref        samtools gffread bedtools
echo "  -- run (02) --"
if have nextflow; then ok "Nextflow is ready"; else warn "Nextflow not on PATH; offline bundles provide their own launcher"; fi

# 执行模式：NFCORE_PROFILE 同时决定 nextflow -profile 与 chromaflow run 的执行路径
_PROF="$(resolve_nfcore_profile)"; _KIND="$(nfcore_profile_kind "${_PROF}")"
info "Execution mode NFCORE_PROFILE=${_PROF} → $( [[ "${_KIND}" == conda ]] && echo 'Conda mode (02_run_nfcore_conda.sh)' || echo 'Container mode (02_run_nfcore.sh)' )"
case "${_KIND}" in
    container)
        if have apptainer || have singularity || have docker || have podman; then
            ok "Container engine is ready（Image cache: ${NXF_SINGULARITY_CACHEDIR}）"
        else
            warn "profile=${_PROF} requires a container engine; install/load it before running"
        fi ;;
    conda)
        if have mamba || have conda; then ok "Conda/mamba is ready (nf-core creates per-process environments)"
        else warn "profile=${_PROF} requires Conda/mamba, but neither was detected"; fi ;;
    *)  info "profile=${_PROF} is not a builtin engine profile; passing it unchanged to Nextflow" ;;
esac
echo "  -- idr (03) --";        need_tools ChromaFlow_idr        idr bedtools
echo "  -- motif (04) --";      need_tools ChromaFlow_motif      meme-chip bedtools
echo "  -- annotate (05, ChIPseeker/R) --"
if have Rscript; then
    if Rscript -e 'quit(status=!requireNamespace("ChIPseeker",quietly=TRUE))' 2>/dev/null; then
        ok "R + ChIPseeker are ready (05)"
    else warn "Rscript is available but ChIPseeker is missing (activate ChromaFlow_R)"; fi
else warn "Rscript is missing (activate ChromaFlow_R or use annotate-py)"; fi
# --peakset consensus 需要给 consensus 回填 signalValue（nf-core 的 consensus bed
# 第5列硬编码为 0）。bedtools 优先，缺了就用内置 python 回退。
if have bedtools; then ok "bedtools is available (consensus signalValue calculation)"
elif have python3;  then info "bedtools is unavailable; consensus signalValue will use the Python fallback"
else warn "Neither bedtools nor python3 is available for consensus signalValue calculation"; fi
echo "  -- annotate-py (step 05 alternative) --"; need_tools ChromaFlow_ref bedtools python3
echo "  -- candidates (step 06 integration) --"
need_tools ChromaFlow_downstream pyGenomeTracks
have python3 && ok "Python is available (step 04 FIMO integration)" || warn "Python is missing"
echo "  -- signal (06) --";     need_tools ChromaFlow_downstream computeMatrix plotHeatmap plotProfile
echo "  -- tracks (07) --";     need_tools ChromaFlow_downstream pyGenomeTracks
echo "  -- go-plot (08) --"
have python3 && python3 -c "import matplotlib" 2>/dev/null \
    && ok "matplotlib is available (go-plot)" || warn "matplotlib is missing; go-plot is unavailable"

# 无配置则到此为止
if [[ ${CONFIG_OK} -eq 0 ]]; then
    echo ""; echo "[SKIP 4-6] Supply chromaflow.config to check reference, inputs and outputs"
    echo ""; echo "Results: ${PASS} OK / ${WARN} WARN / ${FAIL} FAIL"
    [[ ${FAIL} -eq 0 ]] && exit 0 || exit 1
fi

# ─────────────── 4. 参考基因组 ───────────────
echo ""; echo "[4] Reference genome (path priority: config > manifest > directory discovery)"
chk_file() { [[ -f "$1" ]] && ok "$2: $1" || warn "$2 Missing: $1"; }
chk_file "${GENOME_COMBINED}" "Genome FASTA"
chk_file "${GENOME_GTF}"      "GTF (step 02 nf-core)"
chk_file "${CHROM_SIZES}"     "chrom.sizes"
if [[ -f "${GENOME_GFF3}" ]]; then
    ok "GFF3 (Step 05 annotation): ${GENOME_GFF3}"
    if grep -qP '\texon\t' "${GENOME_GFF3}" 2>/dev/null; then
        ok "GFF3 contains exon records for feature classification"
    else
        warn "GFF3 lacks exons; run chromaflow normalize-gff (builtin engine requires no AGAT)"
    fi
else
    warn "GFF3 Missing: ${GENOME_GFF3}(not needed for peak calling alone)"
fi

# ── 染色体命名一致性：FASTA 与 GFF3 必须用同一套名字 ──
# 这是换物种时最常见的失败（水稻 Chr1 / chr01 / 1 三种写法都存在）。
# 05 里虽然会拦，但那时 nf-core 已经跑了几小时，所以在这里提前查。
if [[ -s "${GENOME_COMBINED}" && -s "${GENOME_GFF3}" ]]; then
    _fa_names="$(grep '^>' "${GENOME_COMBINED}" 2>/dev/null | head -200 | sed 's/^>//; s/[ \t].*//' | sort -u)"
    _gff_names="$(awk -F'\t' '!/^#/ && NF>=8 {print $1}' "${GENOME_GFF3}" 2>/dev/null | sort -u | head -200)"
    if [[ -n "${_fa_names}" && -n "${_gff_names}" ]]; then
        _shared="$(comm -12 <(echo "${_fa_names}") <(echo "${_gff_names}") | wc -l)"
        _gff_total="$(echo "${_gff_names}" | wc -l)"
        if [[ "${_shared}" -eq 0 ]]; then
            bad "FASTA and GFF3 share no sequence names; annotation cannot match the genome"
            echo "         FASTA: $(echo "${_fa_names}" | head -3 | tr '\n' ' ')"
            echo "         GFF3 : $(echo "${_gff_names}" | head -3 | tr '\n' ' ')"
            echo "         Review reference versions and sequence identifiers before changing any names"
        elif [[ "${_shared}" -lt "${_gff_total}" ]]; then
            warn "GFF3 has $(( _gff_total - _shared ))/${_gff_total} sequence names absent from FASTA; review reference consistency"
        else
            ok "FASTA and GFF3 sequence names match (${_shared} sequences)"
        fi
    fi
fi

echo ""; echo "[4c] Measured values (stored separately from config)"
for _k in MACS_GSIZE READ_LENGTH GENOME_GC; do
    _v="$(eval "echo \"\${${_k}:-}\"")"
    [[ -z "${_v}" ]] && _v="$(derived_get "${_k}" 2>/dev/null || true)"
    if [[ -n "${_v}" ]]; then
        ok "${_k} = ${_v}  (Source: $(value_origin "${_k}"))"
    else
        case "${_k}" in
            MACS_GSIZE)  info "MACS_GSIZE is unset: run chromaflow reference, or let nf-core estimate it from read length" ;;
            READ_LENGTH) info "READ_LENGTH is unset: chromaflow samplesheet samples FASTQ read lengths" ;;
            *)           info "${_k} not measured (informational only)" ;;
        esac
    fi
done

MF="$(reference_manifest_path)"
if [[ -f "${MF}" ]]; then
    info "Reference manifest: ${MF}"
    sed 's/^/         /' "${MF}"
else
    info "No reference manifest yet (reference / normalize-gff has not run)"
fi

# ── 可选功能：缺了也不影响主流程，一律 INFO 不计入 WARN ──
echo ""; echo "[4b] Optional features"
if [[ "${SKIP_ORGANELLE:-false}" == "true" || ${#ORGANELLE_SOURCES[@]} -eq 0 ]]; then
    skip "Organelle merging is disabled (expected if no organelle references are supplied)"
else
    for _n in $(organelle_names); do
        _u=$(echo "${_n}" | tr '[:lower:]' '[:upper:]'); _f="$(eval echo "\${GENOME_${_u}:-}")"
        [[ -s "${_f}" ]] && ok "Organelle ${_n}: ${_f}" \
            || info "Organelle ${_n} is unavailable; reference skips it unless --strict is set"
    done
fi
if have agat_convert_sp_gxf2gxf.pl; then
    ok "AGAT is available (normalize-gff --engine agat)"
else
    info "AGAT is unavailable; normalization uses the builtin Python engine"
    info "  Install the optional AGAT environment if required: mamba create -n ChromaFlow_agat -c bioconda agat"
fi
if [[ -n "${ANNOTATION_INFO:-}" && -f "${ANNOTATION_INFO}" ]]; then ok "annotation_info (GO mappings)"
else info "ANNOTATION_INFO is unset; no local gene2go mappings will be generated"; fi
if [[ -n "${MEME_DB:-}" && -f "${MEME_DB}" ]]; then ok "MEME_DB (step 04 known-motif comparison)"
else info "MEME_DB is unset; step 04 performs de novo discovery without known-motif comparison"; fi

# ─────────────── 5. 输入 / samplesheet ───────────────
echo ""; echo "[5] Inputs and sample sheet"
SAMPLES_TXT="$(resolve_samples_file "${PIPELINE}" || true)"
[[ -f "${SAMPLES_TXT}" ]] && ok "Sample sheet: ${SAMPLES_TXT}" \
    || warn "Sample records not found (checked samples_${PIPELINE}.{txt,tsv,csv} / samples.*); not required for downstream-only runs"
SS="${SAMPLESHEET_DIR}/samplesheet_${PIPELINE}.csv"
[[ -f "${SS}" ]] && ok "samplesheet: $(basename "${SS}")" || warn "No sample sheet; run samplesheet before step 02"

# ─────────────── 6. 各步产物就绪度 ───────────────
echo ""; echo "[6] Output readiness for downstream analysis"
PEAK_DIR=$(get_nfcore_peak_dir "${PIPELINE}")
n_np=$(find "${PEAK_DIR}" -maxdepth 1 -name "*_peaks.narrowPeak" 2>/dev/null | wc -l)
if [[ ${n_np} -gt 0 ]]; then ok "nf-core narrowPeak: ${n_np} files (${PEAK_DIR})"
else warn "No narrowPeak files found; check step 02 and OUTPUT_DIR: ${PEAK_DIR}"; fi

REC="${RESULTS_DIR}/downstream/${PIPELINE}/idr/recommended_peaks.bed"
CONS=$(get_consensus_peaks "${PIPELINE}")
if [[ -s "${REC}" ]]; then ok "A recommended peak file is available; use an explicit --peakset for downstream analysis ($(wc -l < "${REC}") peaks)"
elif [[ -f "${CONS}" ]]; then warn "No recommended peaks; explicitly use --peakset consensus to select: $(basename "${CONS}")"
else warn "No recommended or consensus peaks; run IDR or check OUTPUT_DIR"; fi

BW_DIR="${RESULTS_DIR}/${PIPELINE}/bwa/merged_library/bigwig"
n_bw=$(find "${BW_DIR}" -maxdepth 1 \( -name "*.bigWig" -o -name "*.bw" \) 2>/dev/null | wc -l)
if [[ ${n_bw} -gt 0 ]]; then ok "bigWig: ${n_bw} files available for signal visualization"
else warn "No bigWig files found (required for signal visualization): ${BW_DIR}"; fi

# ─────────────── 汇总 ───────────────
echo ""
echo "============================================"
echo " Results: ${PASS} OK / ${WARN} WARN / ${FAIL} FAIL"
if [[ ${FAIL} -eq 0 ]]; then
    echo " ✅ No fatal issues found"
    [[ ${WARN} -gt 0 ]] && echo "    (Warnings may indicate a step or environment that is not available yet.)"
else
    echo " ❌ ${FAIL} checks need attention; see [FAIL] above"
fi
echo "============================================"
[[ ${FAIL} -eq 0 ]] && exit 0 || exit 1
