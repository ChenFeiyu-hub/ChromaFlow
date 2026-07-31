#!/usr/bin/env bash
# ============================================================================
# 03_idr_analysis.sh — IDR 重复一致性分析
# ============================================================================
# 支持: chipseq / atacseq / dapseq
# 环境: mamba activate ChromaFlow_idr
# 用法: chromaflow idr --type [chipseq|atacseq|dapseq]
# ============================================================================

set -euo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

PIPELINE="${1:-chipseq}"; shift || true
init_pipeline_dirs "${PIPELINE}"

# ---- 解析参数 ----
# --min-pairs N     : 要求 peak 在 ≥N 个两两配对中通过 IDR（默认 1 = 并集，旧行为）
# --exclude-rep S   : 排除某些重复（逗号分隔，匹配文件名子串，如 REP1 或 REP1,REP3）
MIN_PAIRS=1; EXCLUDE_REPS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --min-pairs)   MIN_PAIRS="$2"; shift 2 ;;
        --exclude-rep) EXCLUDE_REPS="$2"; shift 2 ;;
        *) shift ;;
    esac
done

open_log "03_idr" "${PIPELINE}"

log_section "Step 03: IDR Reproducibility Analysis (${PIPELINE})"

# ─────────────────────────────────────────────
# 路径配置（根据 pipeline 自动派生）
# ─────────────────────────────────────────────
NFCORE_PEAKS=$(get_nfcore_peak_dir "${PIPELINE}")
IDR_OUT="${RESULTS_DIR}/downstream/${PIPELINE}/idr"
mkdir -p "${IDR_OUT}"

echo "  Pipeline:    ${PIPELINE}"
echo "  Peak 目录:   ${NFCORE_PEAKS}"
echo "  输出目录:    ${IDR_OUT}"
echo ""

# ─────────────────────────────────────────────
# 检查工具和输入
# ─────────────────────────────────────────────
check_tool idr || { echo "请先: mamba activate ChromaFlow_idr"; exit 1; }

PEAK_FILES=($(find "${NFCORE_PEAKS}" -maxdepth 1 -name "*_peaks.narrowPeak" 2>/dev/null | sort))

# ATAC-seq 也可能产生 broadPeak
if [[ ${#PEAK_FILES[@]} -eq 0 ]] && [[ "${PIPELINE}" == "atacseq" ]]; then
    PEAK_FILES=($(find "$(get_nfcore_peak_dir "${PIPELINE}" broad_peak)" \
        -maxdepth 1 -name "*_peaks.broadPeak" 2>/dev/null | sort))
fi

N_REPS=${#PEAK_FILES[@]}

# ---- 排除指定重复（--exclude-rep）----
if [[ -n "${EXCLUDE_REPS}" ]]; then
    KEEP=(); DROPPED=()
    for f in "${PEAK_FILES[@]}"; do
        bn=$(basename "$f"); drop=false
        IFS=',' read -ra EXS <<< "${EXCLUDE_REPS}"
        for ex in "${EXS[@]}"; do
            ex="$(echo "$ex" | xargs)"
            [[ -n "$ex" && "$bn" == *"$ex"* ]] && drop=true
        done
        if $drop; then DROPPED+=("$bn"); else KEEP+=("$f"); fi
    done
    if [[ ${#DROPPED[@]} -gt 0 ]]; then
        echo "  [排除重复] ${DROPPED[*]}"
    else
        echo "  ⚠ --exclude-rep '${EXCLUDE_REPS}' 未匹配任何文件，保留全部"
    fi
    PEAK_FILES=("${KEEP[@]}")
    N_REPS=${#PEAK_FILES[@]}
fi

if [[ ${N_REPS} -lt 2 ]]; then
    echo "[ERROR] 在 ${NFCORE_PEAKS} 中未找到足够的 peak 文件 (找到 ${N_REPS}，至少需要 2)"
    echo "  请确认 nf-core pipeline 已成功运行"
    exit 1
fi

echo "  找到 ${N_REPS} 个重复的 peak 文件:"
for f in "${PEAK_FILES[@]}"; do
    n_peaks=$(wc -l < "$f")
    echo "    $(basename $f): ${n_peaks} peaks"
done
echo ""

# ─────────────────────────────────────────────
# STEP 0: QC 指标（如有）
# ─────────────────────────────────────────────
SPP_DIR="${RESULTS_DIR}/${PIPELINE}/bwa/merged_library/phantompeakqualtools"
if [[ -d "${SPP_DIR}" ]]; then
    echo "[INFO] phantompeakqualtools QC 指标:"
    echo "  (NSC > 1.05 理想, RSC > 1.0 理想)"
    echo ""
    for spp in "${SPP_DIR}"/*.spp.out; do
        [[ -f "$spp" ]] || continue
        sample=$(basename "$spp" .mLb.clN.sorted.bam.spp.out)
        nsc=$(awk '{print $9}' "$spp")
        rsc=$(awk '{print $10}' "$spp")
        quality=$(awk '{print $11}' "$spp")
        qual_label="PASS"
        [[ "$quality" == "0" ]] && qual_label="FAIL"
        [[ "$quality" == "-1" ]] && qual_label="LOW"
        printf "  %-25s NSC=%-10s RSC=%-10s %s\n" "$sample" "$nsc" "$rsc" "$qual_label"
    done
    echo ""
fi

# ─────────────────────────────────────────────
# STEP 1: 按 p-value 排序
# ─────────────────────────────────────────────
echo "[STEP 1] 按 p-value 排序 peak 文件"
echo ""

SORTED_PEAKS=()
for f in "${PEAK_FILES[@]}"; do
    ext="${f##*.}"
    base=$(basename "$f" ".${ext}")
    sorted="${IDR_OUT}/${base}.sorted.${ext}"
    sort -k8,8nr "$f" > "${sorted}"
    SORTED_PEAKS+=("${sorted}")
    echo "  排序完成: ${base}"
done
echo ""

# ─────────────────────────────────────────────
# STEP 2: 多阈值 IDR 分析
# ─────────────────────────────────────────────
echo "[STEP 2] 两两配对 IDR 分析"
echo ""

IDR_THRESHOLDS=(0.01 0.05 0.1)
echo "  评估阈值: ${IDR_THRESHOLDS[*]}"
echo ""

# 自动检测输入格式
INPUT_TYPE="narrowPeak"
if [[ "${PEAK_FILES[0]}" == *".broadPeak" ]]; then
    INPUT_TYPE="broadPeak"
fi

printf "  %-25s" "配对"
for t in "${IDR_THRESHOLDS[@]}"; do
    printf "IDR<%-8s" "$t"
done
printf "总peaks\n"
printf "  %s\n" "$(printf '─%.0s' {1..75})"

PAIR_COUNT=0
ALL_IDR_FILES=()

for (( i=0; i<N_REPS; i++ )); do
    for (( j=i+1; j<N_REPS; j++ )); do
        rep_i=$(basename "${SORTED_PEAKS[$i]}" | sed 's/.*_REP/REP/' | sed 's/_.*//')
        rep_j=$(basename "${SORTED_PEAKS[$j]}" | sed 's/.*_REP/REP/' | sed 's/_.*//')
        pair_name="${rep_i}_vs_${rep_j}"
        out_prefix="${IDR_OUT}/idr_${pair_name}"

        idr \
            --samples "${SORTED_PEAKS[$i]}" "${SORTED_PEAKS[$j]}" \
            --input-file-type "${INPUT_TYPE}" \
            --rank p.value \
            --output-file "${out_prefix}.txt" \
            --plot \
            --soft-idr-threshold 0.5 \
            --log-output-file "${out_prefix}.log" \
            2>/dev/null || true

        printf "  %-25s" "${pair_name}"
        n_total=0
        if [[ -f "${out_prefix}.txt" ]] && [[ -s "${out_prefix}.txt" ]]; then
            n_total=$(wc -l < "${out_prefix}.txt")
            for t in "${IDR_THRESHOLDS[@]}"; do
                log_t=$(python3 -c "import math; print(-math.log10(${t}))" 2>/dev/null || echo "2")
                # global IDR 在第 12 列(-log10 尺度)；$12 >= -log10(t) 即 global IDR <= t
                n_pass=$(awk -v lt="${log_t}" '$12 >= lt' "${out_prefix}.txt" 2>/dev/null | wc -l)
                printf "%-13s" "${n_pass}"
            done
            ALL_IDR_FILES+=("${out_prefix}.txt")
        else
            for t in "${IDR_THRESHOLDS[@]}"; do
                printf "%-13s" "N/A"
            done
        fi
        printf "%s\n" "${n_total}"

        PAIR_COUNT=$((PAIR_COUNT + 1))
    done
done

echo ""

# ─────────────────────────────────────────────
# STEP 3: 生成不同阈值的 peak 集合
# ─────────────────────────────────────────────
echo "[STEP 3] 生成 peak 集合"
echo ""

if [[ ${#ALL_IDR_FILES[@]} -eq 0 ]]; then
    echo "  [WARN] 无 IDR 输出文件，跳过"
else
    for t in "${IDR_THRESHOLDS[@]}"; do
        outbed="${IDR_OUT}/idr_peaks_threshold_${t}.bed"
        log_t=$(python3 -c "import math; print(-math.log10(${t}))" 2>/dev/null || echo "2")

        # IDR 输出列: 1chr 2start 3end 4name 5scaledIDR 6strand 7- 8signalValue
        #             9- 10summit 11localIDR 12globalIDR ...
        # 过滤: 第12列 globalIDR(-log10) >= -log10(t)  即 global IDR <= t
        # 保留: signalValue(第8列) 与 globalIDR(第12列)，供下游排序/筛选
        # 每个配对先各自合并，再统计每个区域被多少个配对支持（用于 --min-pairs）
        PAIRBEDS=(); n_raw=0
        for idrf in "${ALL_IDR_FILES[@]}"; do
            pb="${IDR_OUT}/.tmp_pair_$(basename "${idrf}" .txt)_${t}.bed"
            awk -v lt="${log_t}" 'BEGIN{OFS="\t"} $12 >= lt {print $1,$2,$3,$8,$12}' "${idrf}" \
                | sort -k1,1 -k2,2n \
                | bedtools merge -i - -c 4,5 -o max,max > "${pb}"
            n_raw=$(( n_raw + $(wc -l < "${pb}") ))
            [[ -s "${pb}" ]] && PAIRBEDS+=("${pb}")
        done

        if [[ ${#PAIRBEDS[@]} -gt 0 ]]; then
            # 合并所有配对，并数每个合并区域被几个配对覆盖
            cat "${PAIRBEDS[@]}" | sort -k1,1 -k2,2n \
                | bedtools merge -i - -c 4,5 -o max,max > "${IDR_OUT}/.tmp_union_${t}.bed"
            bedtools annotate -i "${IDR_OUT}/.tmp_union_${t}.bed" -files "${PAIRBEDS[@]}" -counts \
                | awk -v mp="${MIN_PAIRS}" 'BEGIN{OFS="\t"}{
                      np=0; for(i=6;i<=NF;i++) if($i+0>0) np++;
                      if(np>=mp) print $1,$2,$3,"IDR_peak_"(++n),$4,"+",$5,np }' \
                > "${outbed}"
            rm -f "${IDR_OUT}/.tmp_union_${t}.bed"
        else
            : > "${outbed}"
        fi
        rm -f "${PAIRBEDS[@]}" 2>/dev/null || true

        n_peaks=0
        [[ -s "${outbed}" ]] && n_peaks=$(wc -l < "${outbed}")
        echo "  IDR < ${t}: 过滤后 ${n_raw} 行 → 合并后满足 ≥${MIN_PAIRS} 配对: ${n_peaks} peaks → $(basename ${outbed})"
    done
    echo ""
    echo "  (输出列: chr start end name signalValue strand globalIDR n_pairs；n_pairs=支持该peak的配对数)"

    # 可重复性分层速览：各 peak 被多少个配对支持（帮助判断核心集大小）
    REC_T="${IDR_OUT}/idr_peaks_threshold_${IDR_THRESHOLD:-0.01}.bed"
    if [[ -s "${REC_T}" ]] && [[ ${#ALL_IDR_FILES[@]} -gt 1 ]]; then
        echo ""
        echo "  [可重复性分层] IDR<${IDR_THRESHOLD:-0.01} 的 peak 被几个配对支持:"
        awk '{c[$8]++} END{for(k in c) printf "    %s 个配对支持: %d peaks\n", k, c[k]}' "${REC_T}" | sort
        echo "    (用 --min-pairs 2 可只保留 ≥2 个配对支持的高置信核心集)"
    fi
fi

# ─────────────────────────────────────────────
# STEP 4: 选择推荐 peak 集合
# ─────────────────────────────────────────────
echo ""
echo "[STEP 4] 选择推荐 peak 集合"

RECOMMENDED=""

# 从严到松择优：0.01 → 0.05 → 0.1，取第一个 peak 数足够的
for t in 0.01 0.05 0.1; do
    candidate="${IDR_OUT}/idr_peaks_threshold_${t}.bed"
    if [[ -s "${candidate}" ]]; then
        n=$(wc -l < "${candidate}")
        if [[ ${n} -gt 10 ]]; then
            RECOMMENDED="${candidate}"
            echo "  使用 IDR < ${t} peak 集合 (${n} peaks，含 signalValue 富集值)"
            break
        fi
    fi
done

CONSENSUS=$(get_consensus_peaks "${PIPELINE}")
if [[ -z "${RECOMMENDED}" ]] || [[ ! -s "${RECOMMENDED}" ]]; then
    if [[ -f "${CONSENSUS}" ]]; then
        echo ""
        echo "  ⚠  IDR peak 数量不足 (<10)，退回 nf-core consensus peaks"
        # 回退时也把 MACS 富集挂回：与合并后的各重复 narrowPeak 求交，取最大 signalValue
        CONS_SIG="${IDR_OUT}/consensus_with_signal.bed"
        ALL_NP="${IDR_OUT}/.tmp_all_narrowpeak.bed"
        cat "${SORTED_PEAKS[@]}" > "${ALL_NP}"
        aw=$(awk 'NR==1{print NF; exit}' "${CONSENSUS}")   # consensus 列宽
        if [[ -n "${aw}" ]] && bedtools intersect -a "${CONSENSUS}" -b "${ALL_NP}" -wa -wb 2>/dev/null \
            | awk -v aw="${aw}" 'BEGIN{OFS="\t"}{
                  sig=$(aw+7);                              # narrowPeak signalValue = -b 第7列
                  key=$1"\t"$2"\t"$3"\t"$4;
                  if(sig>m[key]) m[key]=sig }
                  END{ for(k in m){ split(k,a,"\t"); print a[1],a[2],a[3],a[4],m[k],"+" } }' \
            | sort -k1,1 -k2,2n > "${CONS_SIG}" 2>/dev/null && [[ -s "${CONS_SIG}" ]]; then
            RECOMMENDED="${CONS_SIG}"
            echo "     consensus peaks (已回填 signalValue): $(wc -l < "${CONS_SIG}")"
        else
            RECOMMENDED="${CONSENSUS}"
            echo "     consensus peaks (无富集值): $(wc -l < "${CONSENSUS}")"
        fi
        rm -f "${ALL_NP}"
    fi
fi

if [[ -n "${RECOMMENDED}" ]] && [[ -f "${RECOMMENDED}" ]]; then
    ln -sf "$(readlink -f "${RECOMMENDED}")" "${IDR_OUT}/recommended_peaks.bed"
    echo ""
    echo "  ✅ 推荐 peak 集合: ${IDR_OUT}/recommended_peaks.bed"
    echo "     来源: $(basename ${RECOMMENDED})"
    echo "     数量: $(wc -l < "${RECOMMENDED}") peaks"
fi

# ─────────────────────────────────────────────
# STEP 5: 与 consensus 比较
# ─────────────────────────────────────────────
echo ""
echo "[STEP 5] 与 nf-core consensus peaks 比较"

if [[ -f "${CONSENSUS}" ]]; then
    n_consensus=$(wc -l < "${CONSENSUS}")
    echo "  nf-core consensus peaks: ${n_consensus}"

    for t in "${IDR_THRESHOLDS[@]}"; do
        idr_bed="${IDR_OUT}/idr_peaks_threshold_${t}.bed"
        if [[ -f "${idr_bed}" ]] && [[ -s "${idr_bed}" ]]; then
            n_idr=$(wc -l < "${idr_bed}")
            n_overlap=$(bedtools intersect -a "${idr_bed}" -b "${CONSENSUS}" -u 2>/dev/null | wc -l)
            echo "  IDR<${t}: ${n_overlap}/${n_idr} 与 consensus 重叠"
        fi
    done
else
    echo "  [SKIP] consensus 文件不存在: ${CONSENSUS}"
fi

# ─────────────────────────────────────────────
# 总结
# ─────────────────────────────────────────────
echo ""
echo "============================================"
echo " IDR 分析总结 (${PIPELINE})"
echo "============================================"
echo ""
echo "  重复数:    ${N_REPS}"
echo "  配对数:    ${PAIR_COUNT}"
echo "  输出目录:  ${IDR_OUT}"
echo ""
echo "[DONE] $(date)"
echo ""
echo "下一步: chromaflow motif --type ${PIPELINE}"