#!/usr/bin/env bash
# ============================================================================
# 00_prepare_reference.sh — 准备参考基因组（可选步骤，做过一次就不必再跑）
# ============================================================================
# 数据流向：
#   GENOME_FASTA_SOURCE (.fa/.fa.gz) ─→ GENOME_NUCLEAR ─┐
#   细胞器 fasta ×N（可选，ORGANELLE_SOURCES）──────────→├─cat─→ combined.fa → faidx / chrom.sizes
#                                                        ┘
#   GENOME_GFF3_SOURCE (.gff3/.gz)   ─→ GFF3 ──gffread──→ GTF
#
# 环境: mamba activate ChromaFlow_ref     用法: chromaflow reference [选项]
#
# 选项:
#   --skip-organelle     跳过细胞器合并（无叶绿体/线粒体基因组的物种）
#   --skip-gtf           跳过 GFF3→GTF（只跑下游、不跑 nf-core 时）
#   --force              已存在的产物也重新生成（含重测基因组大小）
#   --strict             细胞器缺失时报错退出（默认只警告并跳过）
#   --no-measure         不实测有效基因组大小/GC（大基因组想省时间时）
#   --kmer N             实测 gsize 用的 k（默认取实测读长，未知时 150）
#   --write-config       把实测值以带标记的块写回 chromaflow.config（默认不写）
#
# 产物不写进 config：本步生成的文件会登记到
#   <OUTPUT_DIR>/.chromaflow/reference.manifest，后续步骤自动识别。
# ============================================================================

set -euo pipefail

source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

SKIP_ORG="${SKIP_ORGANELLE:-false}"; SKIP_GTF="false"; FORCE="false"; STRICT="false"
MEASURE="true"; WRITE_CONFIG="false"; KMER=""
# config 里显式写的 MACS_GSIZE 要尊重；从 derived.env 读回来的则可以重测
MACS_GSIZE_FROM_CONFIG=""
[[ "$(value_origin MACS_GSIZE)" == "config" ]] && MACS_GSIZE_FROM_CONFIG="${MACS_GSIZE}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-organelle|--no-organelle) SKIP_ORG="true"; shift ;;
        --skip-gtf)  SKIP_GTF="true"; shift ;;
        --force)     FORCE="true";   shift ;;
        --strict)    STRICT="true";  shift ;;
        --no-measure) MEASURE="false"; shift ;;
        --kmer)      KMER="$2"; shift 2 ;;
        --write-config) WRITE_CONFIG="true"; shift ;;
        *) shift ;;
    esac
done

[[ "${FORCE}" == "true" ]] && export FORCE_DECOMPRESS=true   # --force 时也重做解压
REF_OUT="$(derived_ref_dir)"
open_log "00_prepare_reference"
log_section "Step 00: Prepare Reference Genome (${SPECIES_NAME:-${SPECIES_SHORT}})"

echo "  物种:       ${SPECIES_NAME:-${SPECIES_SHORT}}"
echo "  基因组版本: ${GENOME_VERSION:-未标注}"
echo "  产物目录:   ${REF_OUT}"
echo ""

exists() { [[ -s "$1" && "${FORCE}" != "true" ]]; }

# =============================================
# STEP 1: 核基因组 FASTA
# =============================================
echo "[STEP 1] 核基因组 FASTA"
if exists "${GENOME_NUCLEAR}"; then
    echo "  [SKIP] 已存在: ${GENOME_NUCLEAR}"
elif [[ -s "${GENOME_FASTA_SOURCE:-}" ]]; then
    echo "  源文件: ${GENOME_FASTA_SOURCE} ($(detect_compression "${GENOME_FASTA_SOURCE}"))"
    # 压缩才解压（按魔数判断，不看扩展名）；已是纯文本则原地使用，不做多余拷贝
    PLAIN_FA="$(ensure_plain_file "${GENOME_FASTA_SOURCE}" "${REF_OUT}" \
                 "$(basename "${GENOME_NUCLEAR}")")" || die "核基因组解压失败"
    GENOME_NUCLEAR="${PLAIN_FA}"
    echo "  [DONE] ${GENOME_NUCLEAR} ($(du -h "${GENOME_NUCLEAR}" | cut -f1))"
elif [[ -s "${GENOME_COMBINED:-}" ]]; then
    echo "  [SKIP] 未给源 fasta，但已有可用基因组: ${GENOME_COMBINED}"
    GENOME_NUCLEAR="${GENOME_COMBINED}"
else
    die "找不到核基因组。请在 config 设 GENOME_FASTA_SOURCE=<你的 .fa / .fa.gz>"
fi
manifest_register GENOME_NUCLEAR "${GENOME_NUCLEAR}" >/dev/null

# =============================================
# STEP 2: GFF3 注释
# =============================================
echo ""
echo "[STEP 2] GFF3 注释"
GFF3_OUT="${GENOME_GFF3}"
if [[ -s "${GENOME_GFF3_SOURCE:-}" ]]; then
    COMP="$(detect_compression "${GENOME_GFF3_SOURCE}")"
    echo "  源文件: ${GENOME_GFF3_SOURCE} (${COMP})"
    if [[ "${COMP}" == "none" ]]; then
        # 已解压的 .gff3/.gff：直接用，不复制、不改名
        GFF3_OUT="${GENOME_GFF3_SOURCE}"
        echo "  [USE ] 已是纯文本，直接使用（无需解压）"
    else
        GFF3_OUT="$(ensure_plain_file "${GENOME_GFF3_SOURCE}" "${REF_OUT}")" || die "GFF3 解压失败"
        echo "  [DONE] ${GFF3_OUT}"
    fi
elif exists "${GFF3_OUT}"; then
    echo "  [SKIP] 未设 GENOME_GFF3_SOURCE，沿用已识别到的: ${GFF3_OUT}"
else
    log_warn "无 GFF3 源（config 未设 GENOME_GFF3_SOURCE）——注释类步骤(05/07)将不可用"
    echo "        提示: 只做 peak calling 可忽略；否则请补 .gff3（压缩与否都行）"
    GFF3_OUT=""
fi
if [[ -n "${GFF3_OUT}" && -s "${GFF3_OUT}" ]]; then
    GENOME_GFF3="${GFF3_OUT}"
    manifest_register GENOME_GFF3 "${GENOME_GFF3}" >/dev/null
fi

# =============================================
# STEP 3: 细胞器基因组（完全可选）
# =============================================
echo ""
echo "[STEP 3] 细胞器基因组（可选）"
ORGANELLE_FILES=()
if [[ "${SKIP_ORG}" == "true" ]]; then
    echo "  [SKIP] 已指定 --skip-organelle / SKIP_ORGANELLE=true —— 只用核基因组"
elif [[ ${#ORGANELLE_SOURCES[@]} -eq 0 ]]; then
    echo "  [SKIP] 未配置 ORGANELLE_SOURCES —— 只用核基因组"
    echo "         （该物种没有可用的叶绿体/线粒体基因组时，这是正常情况）"
else
    missing=0
    for entry in "${ORGANELLE_SOURCES[@]}"; do
        [[ -n "${entry}" ]] || continue
        org_name="${entry%%:*}"; org_source="${entry#*:}"
        org_upper=$(echo "${org_name}" | tr '[:lower:]' '[:upper:]')
        org_file="$(eval echo "\${GENOME_${org_upper}:-}")"
        if [[ -s "${org_file}" ]]; then
            # 细胞器 fasta 若是压缩的，先解压（合并时需要纯文本）
            org_plain="$(ensure_plain_file "${org_file}" "${REF_OUT}")" || org_plain="${org_file}"
            echo "  [  OK ] ${org_name}: ${org_plain}"
            ORGANELLE_FILES+=("${org_plain}")
            manifest_register "GENOME_${org_upper}" "${org_plain}" >/dev/null
        else
            missing=1
            log_warn "缺 ${org_name}: ${org_file}"
            if [[ "${org_source}" == NC_* || "${org_source}" == NW_* || "${org_source}" == N[A-Z]_* ]]; then
                echo "         下载: curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${org_source}&rettype=fasta&retmode=text' \\"
                echo "                 | sed '1s/>.*/>${org_name}/' > ${org_file}"
            else
                echo "         请把 ${org_source} 放到 ${org_file}"
            fi
        fi
    done
    if [[ ${missing} -eq 1 ]]; then
        if [[ "${STRICT}" == "true" ]]; then
            die "细胞器缺失且指定了 --strict"
        fi
        echo ""
        echo "  [继续] 缺失的细胞器已跳过，用现有序列继续（如需强制中断加 --strict）"
        echo "         完全不需要细胞器时可加 --skip-organelle，或把 ORGANELLE_SOURCES 置空"
    fi
fi

# =============================================
# STEP 4: 合并基因组
# =============================================
echo ""
echo "[STEP 4] 合并基因组 FASTA"
COMBINED="${GENOME_COMBINED}"
if [[ ${#ORGANELLE_FILES[@]} -eq 0 ]]; then
    # 无细胞器：直接用核基因组，不做多余拷贝
    COMBINED="${GENOME_NUCLEAR}"
    echo "  [DONE] 无细胞器，直接使用核基因组: ${COMBINED}"
elif exists "${COMBINED}"; then
    echo "  [SKIP] 已存在: ${COMBINED}（--force 可重建）"
else
    cat "${GENOME_NUCLEAR}" "${ORGANELLE_FILES[@]}" > "${COMBINED}"
    echo "  [DONE] 核基因组 + ${#ORGANELLE_FILES[@]} 个细胞器 → ${COMBINED} ($(du -h "${COMBINED}" | cut -f1))"
fi
GENOME_COMBINED="${COMBINED}"
manifest_register GENOME_COMBINED "${GENOME_COMBINED}" >/dev/null

# =============================================
# STEP 5: 索引 + chrom.sizes
# =============================================
echo ""
echo "[STEP 5] 索引与 chrom.sizes"
if [[ -f "${GENOME_COMBINED}.fai" && "${FORCE}" != "true" ]]; then
    echo "  [SKIP] faidx 已存在"
else
    check_tool samtools || die "需要 samtools: mamba activate ChromaFlow_ref"
    samtools faidx "${GENOME_COMBINED}"
    echo "  [DONE] samtools faidx"
fi
CHROM_SIZES=""            # 清空后由 ensure_chrom_sizes 重新生成/登记
CS="$(ensure_chrom_sizes || true)"
CHROM_SIZES="${CS}"       # ensure_chrom_sizes 在子 shell 中运行，需回写
[[ -s "${CS}" ]] && echo "  [DONE] chrom.sizes: ${CS}" || log_warn "chrom.sizes 生成失败"

# =============================================
# STEP 6: GFF3 → GTF（nf-core 需要）
# =============================================
echo ""
echo "[STEP 6] GFF3 → GTF"
if [[ "${SKIP_GTF}" == "true" ]]; then
    echo "  [SKIP] --skip-gtf"
elif [[ ! -s "${GENOME_GFF3}" ]]; then
    log_warn "无 GFF3，跳过 GTF 转换（02 nf-core 需要 GTF）"
else
    GTF_OUT="${GENOME_GTF}"
    [[ "${GTF_OUT}" == *.gtf ]] || GTF_OUT="${REF_OUT}/${SPECIES_SHORT}_genes.gtf"
    if exists "${GTF_OUT}"; then
        echo "  [SKIP] 已存在: ${GTF_OUT}"
        manifest_register GENOME_GTF "${GTF_OUT}" >/dev/null
    elif check_tool gffread; then
        gffread "${GENOME_GFF3}" -T -o "${GTF_OUT}"
        echo "  [DONE] ${GTF_OUT} ($(wc -l < "${GTF_OUT}") 行)"
        manifest_register GENOME_GTF "${GTF_OUT}" >/dev/null
    else
        log_warn "无 gffread，跳过 GTF（mamba activate ChromaFlow_ref）"
    fi
fi

# =============================================
# 统计
# =============================================
echo ""
echo "============================================"
echo " 参考基因组统计"
echo "============================================"
CS="${CHROM_SIZES:-}"
if [[ -s "${CS}" ]]; then
    echo ""
    echo "各序列长度："
    awk '{printf "  %-25s %12d bp\n", $1, $2}' "${CS}" | head -50
    total_bp=$(awk '{sum+=$2} END{print sum+0}' "${CS}")
    n_seqs=$(wc -l < "${CS}")
    org_pattern=""
    while read -r name; do
        [[ -n "${name}" ]] || continue
        org_pattern="${org_pattern:+${org_pattern}|}^${name}"
    done < <(organelle_names)
    if [[ -n "${org_pattern}" && ${#ORGANELLE_FILES[@]} -gt 0 ]]; then
        nuclear_bp=$(grep -v -E "${org_pattern}" "${CS}" | awk '{sum+=$2} END{print sum+0}')
    else
        nuclear_bp=${total_bp}
    fi
    organelle_bp=$((total_bp - nuclear_bp))
    echo ""
    echo "  总序列数:  ${n_seqs}"
    echo "  总长度:    ${total_bp} bp"
    echo "  核基因组:  ${nuclear_bp} bp"
    [[ ${organelle_bp} -gt 0 ]] && echo "  细胞器:    ${organelle_bp} bp"
fi

if [[ -s "${GENOME_GFF3}" ]]; then
    n_genes=$(grep -c $'\tgene\t' "${GENOME_GFF3}" 2>/dev/null || true)
    if [[ "${n_genes:-0}" -gt 0 ]]; then
        echo "  GFF3 基因数: ${n_genes}"
    else
        echo "  GFF3 无 gene 行 —— 建议跑 chromaflow normalize-gff 补全层级"
    fi
fi

echo ""
echo "============================================"
echo " STEP 7: 有效基因组大小 / GC（自动实测）"
echo "============================================"
# MACS 用 gsize 计算 Poisson 背景 λ = N·d/gsize。以前这里只打印一条让你手动
# 跑 unique-kmers.py 再回填 config 的提示；现在直接测，并写进 derived.env，
# 02 自动取用 —— 不改你的 chromaflow.config。
K_FOR_GSIZE="${KMER:-${READ_LENGTH:-150}}"
[[ -z "${KMER}" && -z "${READ_LENGTH:-}" ]] && \
    echo "  注: 读长尚未实测，k 暂用 150；01 生成样本表时会从 fastq 实测并记录"
GSIZE_SIG="$(file_signature "${GENOME_COMBINED}" "k=${K_FOR_GSIZE}")"
OLD_SIG="$(derived_get MACS_GSIZE_SIG 2>/dev/null || true)"

if [[ "${MEASURE}" == "false" ]]; then
    echo "  [SKIP] --no-measure"
elif [[ -n "${MACS_GSIZE_FROM_CONFIG}" ]]; then
    echo "  [SKIP] config 已显式指定 MACS_GSIZE=${MACS_GSIZE_FROM_CONFIG}，尊重你的设置"
    echo "         想让脚本实测：把 config 里那行注释掉即可"
elif [[ -n "$(derived_get MACS_GSIZE 2>/dev/null || true)" && "${OLD_SIG}" == "${GSIZE_SIG}" && "${FORCE}" != "true" ]]; then
    echo "  [SKIP] 已有实测值 $(derived_get MACS_GSIZE)（基因组与读长未变；--force 可重测）"
    MACS_GSIZE="$(derived_get MACS_GSIZE)"
else
    if command -v unique-kmers.py &>/dev/null; then
        echo "  方法: unique-kmers.py -k ${K_FOR_GSIZE}（唯一 k-mer 数，deepTools/nf-core 同款口径）"
        echo "  基因组 $(du -h "${GENOME_COMBINED}" | cut -f1)，大基因组可能需要几分钟..."
    else
        log_warn "未装 khmer，退到「非 N 碱基数」（偏大约 5-10%，作为上界可用）"
        echo "         精确做法: mamba install -c bioconda khmer，或让 02 交给 nf-core 内部的 KHMER_UNIQUEKMERS 计算"
    fi
    if GS_OUT="$(measure_effective_genome_size "${GENOME_COMBINED}" "${K_FOR_GSIZE}")"; then
        GS="${GS_OUT%% *}"; MEASURE_METHOD="${GS_OUT#* }"
        MACS_GSIZE="${GS}"
        derived_set MACS_GSIZE "${GS}"
        derived_set MACS_GSIZE_SIG "${GSIZE_SIG}"
        derived_set MACS_GSIZE_METHOD "${MEASURE_METHOD}"
        printf "  [DONE] MACS_GSIZE = %s  (%s)\n" "${GS}" "${MEASURE_METHOD}"
    else
        log_warn "有效基因组大小测量失败；02 会改为把 --read_length 交给 nf-core 自行计算"
    fi
fi

# GC%：便宜的一次扫描。注意它只作参考展示，04 的背景模型由
# fasta-get-markov 从真实序列直接估计，不依赖这个数。
if [[ "${MEASURE}" != "false" ]]; then
    if [[ -n "$(derived_get GENOME_GC 2>/dev/null || true)" && "${FORCE}" != "true" ]]; then
        echo "  [SKIP] GC = $(derived_get GENOME_GC)%（已实测）"
    elif GCV="$(measure_genome_gc "${GENOME_COMBINED}")" && [[ -n "${GCV}" ]]; then
        derived_set GENOME_GC "${GCV}"
        echo "  [DONE] 基因组 GC = ${GCV}%（仅供参考；MEME 背景模型另行由序列估计）"
    fi
fi

if [[ "${WRITE_CONFIG}" == "true" ]]; then
    write_measured_block_to_config
fi

echo ""
echo "============================================"
echo " 产物清单（已登记，无需写进 config）"
echo "============================================"
MF="$(reference_manifest_path)"
[[ -f "${MF}" ]] && sed 's/^/  /' "${MF}"
DF="$(derived_values_path)"
if [[ -f "${DF}" ]]; then
    echo ""
    echo "  实测值（derived.env）:"
    sed 's/^/    /' "${DF}"
fi
echo ""
echo "  manifest: ${MF}"
echo ""
echo "[DONE] 参考准备完成 $(date)"
echo ""
echo "下一步:"
[[ -s "${GENOME_GFF3}" ]] && echo "  chromaflow normalize-gff        # 可选：补全 GFF 层级，注释更准"
echo "  chromaflow samplesheet --type <chipseq|atacseq|dapseq>"
