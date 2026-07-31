#!/usr/bin/env bash
# ============================================================================
# 00b_normalize_gff.sh — 规范化 GFF3（可选步骤，注释更准；跳过不影响主流程）
# ============================================================================
# 为什么需要：Phytozome/部分 GFF 可能缺 exon 行、或 UTR 命名不统一，导致
# makeTxDbFromGFF 建的基因模型不完整 → ChIPseeker 的 UTR/exon/intron 判定错误
# （典型症状：注释饼图里 5'UTR/3'UTR 占比为 0）。
#
# 用法: chromaflow normalize-gff [--input raw.gff3] [--output clean.gff3]
#                               [--engine auto|agat|builtin|gffread] [--force]
#
# 引擎（AGAT 装不上时不必卡住）:
#   auto     (默认) 有 AGAT 就用 AGAT，没有就自动退到 builtin
#   agat     强制用 AGAT（agat_convert_sp_gxf2gxf.pl），最完整
#   builtin  纯 Python，无任何额外依赖：补 exon / 推 UTR / 合成 gene 层
#   gffread  用 gffread 清理层级（ChromaFlow_ref 里本来就有）
#
# 产物写入 manifest，后续步骤自动使用，不需要回头改 config 的 GENOME_GFF3。
# ============================================================================
set -euo pipefail
source "${CHROMAFLOW_HOME:?请通过 chromaflow CLI 运行，或先 export CHROMAFLOW_HOME}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

IN=""; OUT_GFF=""; ENGINE="${GFF_NORMALIZE_ENGINE:-auto}"; FORCE="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)  IN="$2"; shift 2 ;;
        --output) OUT_GFF="$2"; shift 2 ;;
        --engine) ENGINE="$2"; shift 2 ;;
        --force)  FORCE="true"; shift ;;
        *) shift ;;
    esac
done

# 默认输入：原始 GFF 优先（GENOME_GFF3_RAW / GENOME_GFF3_SOURCE），否则当前 GENOME_GFF3
: "${IN:=${GENOME_GFF3_RAW:-${GENOME_GFF3_SOURCE:-${GENOME_GFF3}}}}"

open_log "00b_normalize_gff"
log_section "Step 00b: GFF 规范化 (engine=${ENGINE})"

[[ -s "${IN}" ]] || die "找不到输入 GFF: ${IN}（用 --input 指定，或在 config 设 GENOME_GFF3_SOURCE）"

# 输入若是压缩的（按魔数判断，不看扩展名）先解压；已解压则原样使用
[[ "${FORCE}" == "true" ]] && export FORCE_DECOMPRESS=true
IN_COMP="$(detect_compression "${IN}")"
if [[ "${IN_COMP}" != "none" ]]; then
    IN="$(ensure_plain_file "${IN}" "$(derived_ref_dir)")" || die "输入 GFF 解压失败: ${IN}"
fi

# 默认输出：与输入同名 + .agat.gff3（保持历史命名，便于自动识别优先级）
if [[ -z "${OUT_GFF}" ]]; then
    base="$(basename "${IN}")"; base="${base%.gz}"; base="${base%.gff3}"; base="${base%.gff}"
    OUT_GFF="$(derived_ref_dir)/${base}.agat.gff3"
fi
[[ "${OUT_GFF}" == "${IN}" ]] && OUT_GFF="${OUT_GFF%.gff3}.normalized.gff3"

echo "  输入: ${IN}"
echo "  输出: ${OUT_GFF}"
echo ""

if [[ -s "${OUT_GFF}" && "${FORCE}" != "true" ]]; then
    echo "  [SKIP] 输出已存在（--force 可重做）"
else
    # ── 选择引擎 ──
    if [[ "${ENGINE}" == "auto" ]]; then
        if command -v agat_convert_sp_gxf2gxf.pl &>/dev/null; then ENGINE="agat"
        else
            ENGINE="builtin"
            log_warn "未检测到 AGAT，自动改用内置规范化器 (builtin)"
            echo "         如需 AGAT: mamba create -n ChromaFlow_agat -c bioconda agat"
            echo "         （单独建环境，避免与 ChromaFlow_ref 的依赖冲突）"
        fi
    fi

    case "${ENGINE}" in
        agat)
            check_tool agat_convert_sp_gxf2gxf.pl || die \
                "未找到 AGAT。可改用 --engine builtin（无需依赖），或单独建环境：mamba create -n ChromaFlow_agat -c bioconda agat"
            echo "[AGAT] agat_convert_sp_gxf2gxf.pl 规范化中..."
            agat_convert_sp_gxf2gxf.pl -g "${IN}" -o "${OUT_GFF}" 2>&1 | tail -15 || die "AGAT 规范化失败"
            ;;
        builtin)
            echo "[builtin] 纯 Python 规范化中（补 exon / 推 UTR / 合成 gene 层）..."
            python3 "${CHROMAFLOW_HOME}/lib/steps/normalize_gff.py" "${IN}" "${OUT_GFF}" \
                || die "内置规范化失败"
            ;;
        gffread)
            check_tool gffread || die "未找到 gffread（mamba activate ChromaFlow_ref）"
            echo "[gffread] 清理层级中..."
            gffread "${IN}" --keep-genes -F -O -o "${OUT_GFF}" || die "gffread 规范化失败"
            ;;
        *) die "未知 --engine: ${ENGINE}（可选 auto|agat|builtin|gffread）" ;;
    esac
fi

[[ -s "${OUT_GFF}" ]] || die "输出为空: ${OUT_GFF}"

# ── 规范化后自检 ──
echo ""
echo "[检查] feature 计数:"
awk -F'\t' '!/^#/ && NF>=3 {c[$3]++} END{for(f in c) printf "    %-20s %d\n", f, c[f]}' "${OUT_GFF}" | sort

count_feat() { awk -F'\t' -v t="$1" '$3==t{n++} END{print n+0}' "${OUT_GFF}"; }
has_exon=$(count_feat exon); has_5utr=$(count_feat five_prime_UTR); has_gene=$(count_feat gene)
echo ""
[[ "${has_gene}" -gt 0 ]] && echo "  ✓ 含 gene 行 (${has_gene})"        || echo "  ⚠ 无 gene 行"
[[ "${has_exon}" -gt 0 ]] && echo "  ✓ 含 exon 行 (${has_exon})"        || echo "  ⚠ 无 exon 行，请检查原始 GFF 结构"
[[ "${has_5utr}" -gt 0 ]] && echo "  ✓ 含 five_prime_UTR (${has_5utr})" || echo "  ⚠ 无 UTR 行（原始注释可能不含 UTR，或全部编码区）"

# ── 登记产物：后续步骤自动使用，无需修改 config ──
manifest_register GENOME_GFF3 "${OUT_GFF}"

echo ""
echo "============================================"
echo " GFF 规范化完成 (engine=${ENGINE})"
echo "============================================"
echo "  输出: ${OUT_GFF}"
echo "  已登记为当前 GENOME_GFF3，直接跑注释即可（不必改 config）："
echo "    chromaflow annotate --type <chipseq|dapseq>"
echo "[DONE] $(date)"
