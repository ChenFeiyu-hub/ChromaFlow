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
# 环境: conda activate ChromaFlow_ref     用法: chromaflow reference [选项]
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

source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI, or export CHROMAFLOW_HOME first}/lib/common.sh"
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

echo "  Species:       ${SPECIES_NAME:-${SPECIES_SHORT}}"
echo "  Genome version: ${GENOME_VERSION:-unspecified}"
echo "  Output directory:   ${REF_OUT}"
echo ""

exists() { [[ -s "$1" && "${FORCE}" != "true" ]]; }

# =============================================
# STEP 1: 核基因组 FASTA
# =============================================
echo "[STEP 1] Nuclear genome FASTA"
if exists "${GENOME_NUCLEAR}"; then
    echo "  [SKIP] Already exists: ${GENOME_NUCLEAR}"
elif [[ -s "${GENOME_FASTA_SOURCE:-}" ]]; then
    echo "  Source file: ${GENOME_FASTA_SOURCE} ($(detect_compression "${GENOME_FASTA_SOURCE}"))"
    # 压缩才解压（按魔数判断，不看扩展名）；已是纯文本则原地使用，不做多余拷贝
    PLAIN_FA="$(ensure_plain_file "${GENOME_FASTA_SOURCE}" "${REF_OUT}" \
                 "$(basename "${GENOME_NUCLEAR}")")" || die "Nuclear genome decompression failed"
    GENOME_NUCLEAR="${PLAIN_FA}"
    echo "  [DONE] ${GENOME_NUCLEAR} ($(du -h "${GENOME_NUCLEAR}" | cut -f1))"
elif [[ -s "${GENOME_COMBINED:-}" ]]; then
    echo "  [SKIP] No source FASTA specified; reusing existing genome: ${GENOME_COMBINED}"
    GENOME_NUCLEAR="${GENOME_COMBINED}"
else
    die "Nuclear genome not found. Set GENOME_FASTA_SOURCE=<your .fa / .fa.gz> in config"
fi
manifest_register GENOME_NUCLEAR "${GENOME_NUCLEAR}" >/dev/null

# =============================================
# STEP 2: GFF3 注释
# =============================================
echo ""
echo "[STEP 2] GFF3 annotation"
GFF3_OUT="${GENOME_GFF3}"
if [[ -s "${GFF3_OUT}" && -s "${GFF3_OUT}.report.json" ]]; then
    GFF_SOURCE_ARGS=()
    GFF_ORIGINAL="${GENOME_GFF3_RAW:-${GENOME_GFF3_SOURCE:-}}"
    [[ -n "${GFF_ORIGINAL}" ]] && GFF_SOURCE_ARGS+=(--source "${GFF_ORIGINAL}")
    [[ -n "${GFF_TYPE_MAP:-}" ]] && GFF_SOURCE_ARGS+=(--type-map "${GFF_TYPE_MAP}")
    python3 "${CHROMAFLOW_HOME}/lib/steps/prepare_gtf.py" "${GFF3_OUT}" unused \
        --check-source-only "${GFF_SOURCE_ARGS[@]}" || die "Normalized GFF does not match its source/report; rerun normalize-gff"
    echo "  [USE ] Validated normalized annotation: ${GFF3_OUT}"
elif [[ -s "${GENOME_GFF3_SOURCE:-}" ]]; then
    COMP="$(detect_compression "${GENOME_GFF3_SOURCE}")"
    echo "  Source file: ${GENOME_GFF3_SOURCE} (${COMP})"
    if [[ "${COMP}" == "none" ]]; then
        # 已解压的 .gff3/.gff：直接用，不复制、不改名
        GFF3_OUT="${GENOME_GFF3_SOURCE}"
        echo "  [USE ] Already plain text; using it directly without decompression"
    else
        GFF3_OUT="$(ensure_plain_file "${GENOME_GFF3_SOURCE}" "${REF_OUT}")" || die "GFF3 decompression failed"
        echo "  [DONE] ${GFF3_OUT}"
    fi
elif exists "${GFF3_OUT}"; then
    echo "  [SKIP] GENOME_GFF3_SOURCE is unset; reusing detected annotation: ${GFF3_OUT}"
else
    log_warn "No GFF3 source (GENOME_GFF3_SOURCE unset); annotation steps 05/07 will be unavailable"
    echo "        Peak calling alone does not need this; for annotation, provide a plain or compressed .gff3"
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
echo "[STEP 3] Organelle genomes (optional)"
ORGANELLE_FILES=()
if [[ "${SKIP_ORG}" == "true" ]]; then
    echo "  [SKIP] --skip-organelle / SKIP_ORGANELLE=true: using the nuclear genome only"
elif [[ ${#ORGANELLE_SOURCES[@]} -eq 0 ]]; then
    echo "  [SKIP] ORGANELLE_SOURCES is unset: using the nuclear genome only"
    echo "         (This is expected when no chloroplast/mitochondrial reference is available.)"
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
            log_warn "Missing ${org_name}: ${org_file}"
            if [[ "${org_source}" == NC_* || "${org_source}" == NW_* || "${org_source}" == N[A-Z]_* ]]; then
                echo "         Download: curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${org_source}&rettype=fasta&retmode=text' \\"
                echo "                 | sed '1s/>.*/>${org_name}/' > ${org_file}"
            else
                echo "         Place ${org_source} at ${org_file}"
            fi
        fi
    done
    if [[ ${missing} -eq 1 ]]; then
        if [[ "${STRICT}" == "true" ]]; then
            die "Organelle references are missing and --strict was specified"
        fi
        echo ""
        echo "  [CONTINUE] Missing organelles skipped; continuing with available sequences (use --strict to stop)"
        echo "         To omit organelles, use --skip-organelle or leave ORGANELLE_SOURCES empty"
    fi
fi

# =============================================
# STEP 4: 合并基因组
# =============================================
echo ""
echo "[STEP 4] Combine genome FASTA"
COMBINED="${GENOME_COMBINED}"
if [[ ${#ORGANELLE_FILES[@]} -eq 0 ]]; then
    # 无细胞器：直接用核基因组，不做多余拷贝
    COMBINED="${GENOME_NUCLEAR}"
    echo "  [DONE] No organelles; using the nuclear genome directly: ${COMBINED}"
elif exists "${COMBINED}"; then
    echo "  [SKIP] Already exists: ${COMBINED}(use --force to rebuild)"
else
    cat "${GENOME_NUCLEAR}" "${ORGANELLE_FILES[@]}" > "${COMBINED}"
    echo "  [DONE] Nuclear genome + ${#ORGANELLE_FILES[@]} organelles → ${COMBINED} ($(du -h "${COMBINED}" | cut -f1))"
fi
GENOME_COMBINED="${COMBINED}"
manifest_register GENOME_COMBINED "${GENOME_COMBINED}" >/dev/null

# =============================================
# STEP 5: 索引 + chrom.sizes
# =============================================
echo ""
echo "[STEP 5] Index and chrom.sizes"
if [[ -f "${GENOME_COMBINED}.fai" && "${FORCE}" != "true" ]]; then
    echo "  [SKIP] faidx Already exists"
else
    check_tool samtools || die "Requires samtools: conda activate ChromaFlow_ref"
    samtools faidx "${GENOME_COMBINED}"
    echo "  [DONE] samtools faidx"
fi
CHROM_SIZES=""            # 清空后由 ensure_chrom_sizes 重新生成/登记
CS="$(ensure_chrom_sizes || true)"
CHROM_SIZES="${CS}"       # ensure_chrom_sizes 在子 shell 中运行，需回写
[[ -s "${CS}" ]] && echo "  [DONE] chrom.sizes: ${CS}" || log_warn "Failed to generate chrom.sizes"

# =============================================
# STEP 6: GFF3 → GTF（nf-core 需要）
# =============================================
echo ""
echo "[STEP 6] GFF3 → GTF"
if [[ "${SKIP_GTF}" == "true" ]]; then
    echo "  [SKIP] --skip-gtf"
elif [[ ! -s "${GENOME_GFF3}" ]]; then
    log_warn "No GFF3; skipping GTF conversion (step 02 requires GTF)"
else
    GTF_OUT="${GENOME_GTF}"
    [[ "${GTF_OUT}" == *.gtf ]] || GTF_OUT="${REF_OUT}/${SPECIES_SHORT}_genes.gtf"
    GTF_OUT="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${GTF_OUT}")"
    if check_tool gffread; then
        GTF_ARGS=()
        [[ "${FORCE}" == "true" ]] && GTF_ARGS+=(--force)
        [[ -n "${GFF_TYPE_MAP:-}" ]] && GTF_ARGS+=(--type-map "${GFF_TYPE_MAP}")
        python3 "${CHROMAFLOW_HOME}/lib/steps/prepare_gtf.py" "${GENOME_GFF3}" "${GTF_OUT}" \
            "${GTF_ARGS[@]}" || die "GTF conversion/gene-assignment validation failed; previous GTF was preserved"
        echo "  [DONE] ${GTF_OUT} ($(wc -l < "${GTF_OUT}") lines)"
        manifest_register GENOME_GTF "${GTF_OUT}" >/dev/null
    else
        die "gffread is required to validate/generate GTF; activate reference tools or explicitly use --skip-gtf"
    fi
fi

# =============================================
# 统计
# =============================================
echo ""
echo "============================================"
echo " Reference genome statistics"
echo "============================================"
CS="${CHROM_SIZES:-}"
if [[ -s "${CS}" ]]; then
    echo ""
    echo "Sequence lengths:"
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
    echo "  Sequence count:  ${n_seqs}"
    echo "  Total length:    ${total_bp} bp"
    echo "  Nuclear genome:  ${nuclear_bp} bp"
    [[ ${organelle_bp} -gt 0 ]] && echo "  Organelles:    ${organelle_bp} bp"
fi

if [[ -s "${GENOME_GFF3}" ]]; then
    n_genes=$(grep -c $'\tgene\t' "${GENOME_GFF3}" 2>/dev/null || true)
    if [[ "${n_genes:-0}" -gt 0 ]]; then
        echo "  GFF3 gene count: ${n_genes}"
    else
        echo "  GFF3 has no ordinary gene records; inspect the normalization report for other valid gene types"
    fi
fi

echo ""
echo "============================================"
echo " STEP 7: Effective genome size / GC (automatically measured)"
echo "============================================"
# MACS 用 gsize 计算 Poisson 背景 λ = N·d/gsize。以前这里只打印一条让你手动
# 跑 unique-kmers.py 再回填 config 的提示；现在直接测，并写进 derived.env，
# 02 自动取用 —— 不改你的 chromaflow.config。
K_FOR_GSIZE="${KMER:-${READ_LENGTH:-150}}"
[[ -z "${KMER}" && -z "${READ_LENGTH:-}" ]] && \
    echo "  Read length is not measured yet; using k=150. Step 01 records read length from FASTQ."
GSIZE_SIG="$(file_signature "${GENOME_COMBINED}" "k=${K_FOR_GSIZE}")"
OLD_SIG="$(derived_get MACS_GSIZE_SIG 2>/dev/null || true)"

if [[ "${MEASURE}" == "false" ]]; then
    echo "  [SKIP] --no-measure"
elif [[ -n "${MACS_GSIZE_FROM_CONFIG}" ]]; then
    echo "  [SKIP] Using explicitly configured MACS_GSIZE=${MACS_GSIZE_FROM_CONFIG}"
    echo "         To enable automatic measurement, remove/comment that config assignment"
elif [[ -n "$(derived_get MACS_GSIZE 2>/dev/null || true)" && "${OLD_SIG}" == "${GSIZE_SIG}" && "${FORCE}" != "true" ]]; then
    echo "  [SKIP] Existing measured value $(derived_get MACS_GSIZE)(genome and read length unchanged; use --force to remeasure)"
    MACS_GSIZE="$(derived_get MACS_GSIZE)"
else
    if command -v unique-kmers.py &>/dev/null; then
        echo "  Method: unique-kmers.py -k ${K_FOR_GSIZE}(unique k-mers, matching the deepTools/nf-core method)"
        echo "  Genome size $(du -h "${GENOME_COMBINED}" | cut -f1); large genomes may take several minutes..."
    else
        log_warn "khmer is unavailable; using non-N base count as an upper-bound fallback"
        echo "         For unique-kmer estimation, install khmer or let nf-core run KHMER_UNIQUEKMERS in step 02"
    fi
    if GS_OUT="$(measure_effective_genome_size "${GENOME_COMBINED}" "${K_FOR_GSIZE}")"; then
        GS="${GS_OUT%% *}"; MEASURE_METHOD="${GS_OUT#* }"
        MACS_GSIZE="${GS}"
        derived_set MACS_GSIZE "${GS}"
        derived_set MACS_GSIZE_SIG "${GSIZE_SIG}"
        derived_set MACS_GSIZE_METHOD "${MEASURE_METHOD}"
        printf "  [DONE] MACS_GSIZE = %s  (%s)\n" "${GS}" "${MEASURE_METHOD}"
    else
        log_warn "Genome size measurement failed; step 02 will pass --read_length for nf-core estimation"
    fi
fi

# GC%：便宜的一次扫描。注意它只作参考展示，04 的背景模型由
# fasta-get-markov 从真实序列直接估计，不依赖这个数。
if [[ "${MEASURE}" != "false" ]]; then
    if [[ -n "$(derived_get GENOME_GC 2>/dev/null || true)" && "${FORCE}" != "true" ]]; then
        echo "  [SKIP] GC = $(derived_get GENOME_GC)%(measured)"
    elif GCV="$(measure_genome_gc "${GENOME_COMBINED}")" && [[ -n "${GCV}" ]]; then
        derived_set GENOME_GC "${GCV}"
        echo "  [DONE] Genome GC = ${GCV}%(for reference; MEME estimates its background model independently from sequence)"
    fi
fi

if [[ "${WRITE_CONFIG}" == "true" ]]; then
    write_measured_block_to_config
fi

echo ""
echo "============================================"
echo " Registered outputs (no manual config entries needed)"
echo "============================================"
MF="$(reference_manifest_path)"
[[ -f "${MF}" ]] && sed 's/^/  /' "${MF}"
DF="$(derived_values_path)"
if [[ -f "${DF}" ]]; then
    echo ""
    echo "  Measured values (derived.env):"
    sed 's/^/    /' "${DF}"
fi
echo ""
echo "  manifest: ${MF}"
echo ""
echo "[DONE] Reference preparation completed $(date)"
echo ""
echo "Next step:"
[[ -s "${GENOME_GFF3}" ]] && echo "  chromaflow normalize-gff        # Optional: conservative normalization and validation; then use reference to update GTF"
echo "  chromaflow samplesheet --type <chipseq|atacseq|dapseq>"
