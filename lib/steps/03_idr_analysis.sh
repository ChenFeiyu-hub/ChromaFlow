#!/usr/bin/env bash
# ============================================================================
# 03_idr_analysis.sh — IDR 重复一致性分析
# ============================================================================
# 支持: chipseq / atacseq / dapseq
# 环境: conda activate ChromaFlow_idr
# 用法: chromaflow idr --type [chipseq|atacseq|dapseq]
# ============================================================================

set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
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
IDR_OUT="$(step_outdir "${RESULTS_DIR}/downstream/${PIPELINE}/idr")"
mkdir -p "${IDR_OUT}"

echo "  Pipeline:    ${PIPELINE}"
echo "  Peak directory:   ${NFCORE_PEAKS}"
echo "  Output directory:    ${IDR_OUT}"
echo ""

# ─────────────────────────────────────────────
# 检查工具和输入
# ─────────────────────────────────────────────
for tool in idr bedtools python3; do check_tool "$tool" || exit 1; done
printf 'running\n' > "${IDR_OUT}/idr.status"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "failed\n" > "${IDR_OUT}/idr.status"; fi; [[ -z "${SIGNAL_INPUT:-}" ]] || rm -f "${SIGNAL_INPUT}"' EXIT

mapfile -t PEAK_FILES < <(find "${NFCORE_PEAKS}" -maxdepth 1 -name "*_peaks.narrowPeak" | sort)

# ATAC-seq 也可能产生 broadPeak
if [[ ${#PEAK_FILES[@]} -eq 0 ]] && [[ "${PIPELINE}" == "atacseq" ]]; then
    mapfile -t PEAK_FILES < <(find "$(get_nfcore_peak_dir "${PIPELINE}" broad_peak)" \
        -maxdepth 1 -name "*_peaks.broadPeak" | sort)
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
        echo "  [EXCLUDED] ${DROPPED[*]}"
    else
        echo "  ⚠ --exclude-rep '${EXCLUDE_REPS}' did not match any files; keeping all inputs"
    fi
    PEAK_FILES=("${KEEP[@]}")
    N_REPS=${#PEAK_FILES[@]}
fi

if [[ ${N_REPS} -lt 2 ]]; then
    echo "[ERROR] Not enough peak files in ${NFCORE_PEAKS} (found ${N_REPS}; at least 2 required)"
    echo "  Confirm that the nf-core pipeline completed successfully"
    exit 1
fi

echo "  Found ${N_REPS} replicate peak files:"
[[ "${MIN_PAIRS}" =~ ^[1-9][0-9]*$ ]] || die '--min-pairs must be a positive integer'
[[ "${MIN_PAIRS}" -le $((N_REPS * (N_REPS - 1) / 2)) ]] || die '--min-pairs exceeds the number of replicate pairs'
python3 - "${PEAK_FILES[@]}" <<'PYEOF'
import pathlib, re, sys
groups, reps = set(), set()
for name in sys.argv[1:]:
    match = re.match(r'(.+)_REP([0-9]+)(?:[._]|$)', pathlib.Path(name).name)
    if not match:
        sys.exit('[ERROR] Cannot identify sample/replicate from peak filename: ' + name)
    groups.add(match[1])
    key = (match[1], int(match[2]))
    if key in reps:
        sys.exit('[ERROR] Multiple peak files for the same biological replicate: ' + str(key))
    reps.add(key)
if len(groups) != 1:
    sys.exit('[ERROR] IDR requires one biological condition per run; found: ' + ', '.join(sorted(groups)))
print('[GROUP] ' + next(iter(groups)))
PYEOF
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
    echo "[INFO] phantompeakqualtools QC metrics:"
    echo "  (NSC > 1.05 preferred, RSC > 1.0 preferred)"
    echo ""
    for spp in "${SPP_DIR}"/*.spp.out; do
        [[ -f "$spp" ]] || continue
        sample=$(basename "$spp" .mLb.clN.sorted.bam.spp.out)
        nsc=$(awk '{print $9}' "$spp")
        rsc=$(awk '{print $10}' "$spp")
        quality=$(awk '{print $11}' "$spp")
        case "$quality" in
            -2) qual_label=VERY_LOW ;; -1) qual_label=LOW ;; 0) qual_label=MEDIUM ;;
            1) qual_label=HIGH ;; 2) qual_label=VERY_HIGH ;; *) qual_label=UNKNOWN ;;
        esac
        printf "  %-25s NSC=%-10s RSC=%-10s %s\n" "$sample" "$nsc" "$rsc" "$qual_label"
    done
    echo ""
else
    # 跑了 --skip_spp（或该模块失败）时，用仍然存在的富集指标兜底。
    # FRiP = peak 内 reads 占比，是 ENCODE 的标准富集指标，与 NSC/RSC 回答
    # 的是同一个问题（这个样本到底有没有富集），只是算法不同。
    echo "[INFO] No phantompeakqualtools output (--skip_spp or step not run)"
    FRIP_FOUND=0
    while read -r f; do
        [[ -n "${f}" ]] || continue
        if [[ ${FRIP_FOUND} -eq 0 ]]; then
            echo "  FRiP (fraction of reads in peaks; interpret for your assay and reference):"
            FRIP_FOUND=1
        fi
        val=$(awk 'END{print $NF}' "${f}" 2>/dev/null)
        printf "    %-28s FRiP=%s\n" "$(basename "${f}" | sed 's/[._]FRiP.*//; s/\.mLb\.clN.*//')" "${val}"
    done < <(find "${RESULTS_DIR}/${PIPELINE}" -name '*FRiP*' -type f 2>/dev/null | sort)

    if [[ ${FRIP_FOUND} -eq 0 ]]; then
        echo "  No FRiP file found; inspect the following QC outputs:"
    else
        echo "  Other available QC outputs:"
    fi
    echo "    · deepTools plotFingerprint  —— enrichment profile"
    echo "    · deepTools plotProfile/Heatmap —— TSS enrichment profile"
    echo "    · MultiQC report: alignment, duplication and library complexity"
    echo "    · IDR pairwise table below: replicate reproducibility"
    echo ""
fi

# ─────────────────────────────────────────────
# STEP 1: 按 p-value 排序
# ─────────────────────────────────────────────
echo "[STEP 1] Sort peak files by p-value rank"
echo ""

SORTED_PEAKS=()
for f in "${PEAK_FILES[@]}"; do
    ext="${f##*.}"
    base=$(basename "$f" ".${ext}")
    sorted="${IDR_OUT}/${base}.sorted.${ext}"
    sort -k8,8nr "$f" > "${sorted}"
    SORTED_PEAKS+=("${sorted}")
    echo "  Sorted: ${base}"
done
echo ""

# ─────────────────────────────────────────────
# STEP 2: 多阈值 IDR 分析
# ─────────────────────────────────────────────
echo "[STEP 2] Pairwise IDR analysis"
echo ""

IDR_THRESHOLDS=(0.01 0.05 0.1)
echo "  Thresholds: ${IDR_THRESHOLDS[*]}"
echo ""

# 自动检测输入格式
INPUT_TYPE="narrowPeak"
GLOBAL_IDR_COL=12
if [[ "${PEAK_FILES[0]}" == *".broadPeak" ]]; then
    INPUT_TYPE="broadPeak"
    GLOBAL_IDR_COL=11
fi

printf "  %-25s" "Pair"
for t in "${IDR_THRESHOLDS[@]}"; do
    printf "IDR<%-8s" "$t"
done
printf "Total peaks\n"
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
            --log-output-file "${out_prefix}.log"

        printf "  %-25s" "${pair_name}"
        n_total=0
        if [[ -f "${out_prefix}.txt" ]] && [[ -s "${out_prefix}.txt" ]]; then
            n_total=$(wc -l < "${out_prefix}.txt")
            for t in "${IDR_THRESHOLDS[@]}"; do
                log_t=$(python3 -c "import math; print(-math.log10(${t}))" 2>/dev/null || echo "2")
                # global IDR 在第 12 列(-log10 尺度)；$12 >= -log10(t) 即 global IDR <= t
                n_pass=$(awk -v lt="${log_t}" -v gc="${GLOBAL_IDR_COL}" '$gc >= lt' "${out_prefix}.txt" | wc -l)
                printf "%-13s" "${n_pass}"
            done
            ALL_IDR_FILES+=("${out_prefix}.txt")
        else
            die "IDR returned no usable output for ${pair_name}; inspect ${out_prefix}.log"
        fi
        printf "%s\n" "${n_total}"

        PAIR_COUNT=$((PAIR_COUNT + 1))
    done
done

echo ""

# ─────────────────────────────────────────────
# STEP 3: 生成不同阈值的 peak 集合
# ─────────────────────────────────────────────
echo "[STEP 3] Generate peak sets"
echo ""

if [[ ${#ALL_IDR_FILES[@]} -eq 0 ]]; then
    die 'No nonempty IDR pair outputs; inspect the pair logs'
else
    SIGNAL_INPUT=$(mktemp "${IDR_OUT}/.replicate_signal_XXXXXX")
    cat "${PEAK_FILES[@]}" | LC_ALL=C sort -k1,1 -k2,2n > "${SIGNAL_INPUT}"
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
            awk -v lt="${log_t}" -v gc="${GLOBAL_IDR_COL}" 'BEGIN{OFS="\t"} $gc >= lt {print $1,$2,$3,0,$gc}' "${idrf}" \
                | LC_ALL=C sort -k1,1 -k2,2n \
                | bedtools merge -i - -c 4,5 -o max,max > "${pb}"
            n_raw=$(( n_raw + $(wc -l < "${pb}") ))
            [[ -s "${pb}" ]] && PAIRBEDS+=("${pb}")
        done

        if [[ ${#PAIRBEDS[@]} -gt 0 ]]; then
            # 合并所有配对，并数每个合并区域被几个配对覆盖
            cat "${PAIRBEDS[@]}" | LC_ALL=C sort -k1,1 -k2,2n \
                | bedtools merge -i - -c 4,5 -o max,max > "${IDR_OUT}/.tmp_union_${t}.bed"
            # annotate can emit intervals in internal-bin order; map requires coordinate order.
            bedtools annotate -i "${IDR_OUT}/.tmp_union_${t}.bed" -files "${PAIRBEDS[@]}" -counts \
                | LC_ALL=C sort -k1,1 -k2,2n \
                | awk -v mp="${MIN_PAIRS}" 'BEGIN{OFS="\t"}{
                      np=0; for(i=6;i<=NF;i++) if($i+0>0) np++;
                      if(np>=mp) print $1,$2,$3,"IDR_peak_"(++n),$4,"+",$5,np }' \
                > "${outbed}"
            rm -f "${IDR_OUT}/.tmp_union_${t}.bed"
            if [[ -s "${outbed}" ]]; then
                # Same definition as consensus: maximum overlapping MACS signalValue.
                # IDR's p-value rank column is not an enrichment score.
                bedtools map -a "${outbed}" -b "${SIGNAL_INPUT}" -c 7 -o max -null 0 \
                    | awk 'BEGIN{OFS="\t"} {$5=$9; NF=8; print}' > "${outbed}.signal.tmp"
                mv "${outbed}.signal.tmp" "${outbed}"
            fi
        else
            : > "${outbed}"
        fi
        rm -f "${PAIRBEDS[@]}" 2>/dev/null || true

        n_peaks=0
        [[ -s "${outbed}" ]] && n_peaks=$(wc -l < "${outbed}")
        echo "  IDR < ${t}: ${n_raw} filtered rows -> merged regions supported by >=${MIN_PAIRS} pairs: ${n_peaks} peaks → $(basename ${outbed})"
    done
    echo ""
    echo "  (Columns: chr start end name signalValue strand globalIDR n_pairs; n_pairs counts supporting pairs)"

    # 可重复性分层速览：各 peak 被多少个配对支持（帮助判断核心集大小）
    REC_T="${IDR_OUT}/idr_peaks_threshold_${IDR_THRESHOLD:-0.01}.bed"
    if [[ -s "${REC_T}" ]] && [[ ${#ALL_IDR_FILES[@]} -gt 1 ]]; then
        echo ""
        echo "  [SUPPORT] Number of pairs supporting IDR<=${IDR_THRESHOLD:-0.01} peaks:"
        awk '{c[$8]++} END{for(k in c) printf "    %s supporting pairs: %d peaks\n", k, c[k]}' "${REC_T}" | sort
        echo "    (Use --min-pairs 2 to require support from at least two pairs)"
    fi
fi

# ─────────────────────────────────────────────
# STEP 4: 选择推荐 peak 集合
# ─────────────────────────────────────────────
echo ""
echo "[STEP 4] Select a recommended peak set"

RECOMMENDED=""

# 从严到松择优：0.01 → 0.05 → 0.1，取第一个 peak 数足够的
for t in 0.01 0.05 0.1; do
    candidate="${IDR_OUT}/idr_peaks_threshold_${t}.bed"
    if [[ -s "${candidate}" ]]; then
        n=$(wc -l < "${candidate}")
        if [[ ${n} -gt 10 ]]; then
            RECOMMENDED="${candidate}"
            echo "  Using IDR <= ${t} (${n} peaks, with signalValue)"
            break
        fi
    fi
done

CONSENSUS=$(get_consensus_peaks "${PIPELINE}")
if [[ -z "${RECOMMENDED}" ]] || [[ ! -s "${RECOMMENDED}" ]]; then
    if [[ -f "${CONSENSUS}" ]]; then
        echo ""
        echo "  ⚠  Insufficient IDR peaks; recommendation falls back to nf-core consensus"
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
            echo "     consensus peaks (signalValue populated): $(wc -l < "${CONS_SIG}")"
        else
            RECOMMENDED="${CONSENSUS}"
            echo "     consensus peaks (no enrichment score): $(wc -l < "${CONSENSUS}")"
        fi
        rm -f "${ALL_NP}"
    fi
fi

if [[ -n "${RECOMMENDED}" ]] && [[ -f "${RECOMMENDED}" ]]; then
    ln -sf "$(readlink -f "${RECOMMENDED}")" "${IDR_OUT}/recommended_peaks.bed"
    echo ""
    echo "  ✅ Recommended peak set: ${IDR_OUT}/recommended_peaks.bed"
    echo "     Source: $(basename ${RECOMMENDED})"
    echo "     Count: $(wc -l < "${RECOMMENDED}") peaks"
fi

# ─────────────────────────────────────────────
# STEP 5: 与 consensus 比较
# ─────────────────────────────────────────────
echo ""
echo "[STEP 5] Compare with nf-core consensus peaks"

if [[ -f "${CONSENSUS}" ]]; then
    n_consensus=$(wc -l < "${CONSENSUS}")
    echo "  nf-core consensus peaks: ${n_consensus}"

    for t in "${IDR_THRESHOLDS[@]}"; do
        idr_bed="${IDR_OUT}/idr_peaks_threshold_${t}.bed"
        if [[ -f "${idr_bed}" ]] && [[ -s "${idr_bed}" ]]; then
            n_idr=$(wc -l < "${idr_bed}")
            n_overlap=$(bedtools intersect -a "${idr_bed}" -b "${CONSENSUS}" -u 2>/dev/null | wc -l)
            echo "  IDR<${t}: ${n_overlap}/${n_idr} overlap consensus"
        fi
    done
else
    echo "  [SKIP] Consensus file does not exist: ${CONSENSUS}"
fi

# ─────────────────────────────────────────────
# 总结
# ─────────────────────────────────────────────
echo ""
echo "============================================"
echo " IDR analysis summary (${PIPELINE})"
echo "============================================"
echo ""
echo "  Replicate count:    ${N_REPS}"
echo "  Pair count:    ${PAIR_COUNT}"
echo "  Output directory:  ${IDR_OUT}"
echo ""
echo "[DONE] $(date)"
printf 'complete\n' > "${IDR_OUT}/idr.status"
echo ""
echo "Next step: chromaflow motif --type ${PIPELINE}"
