#!/usr/bin/env bash
# Validate first, publish atomically, then register the normalized annotation.
set -euo pipefail
source "${CHROMAFLOW_HOME:?Run through the chromaflow CLI}/lib/common.sh"
load_config "${CHROMAFLOW_CONFIG:-}"
init_dirs

IN=""; OUT_GFF=""; ENGINE="${GFF_NORMALIZE_ENGINE:-auto}"
EXTRA=()
[[ -n "${GFF_TYPE_MAP:-}" ]] && EXTRA+=("--type-map" "${GFF_TYPE_MAP}")
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input|--output|--engine|--type-map)
            [[ $# -ge 2 && -n "$2" ]] || die "$1 is missing a value"
            case "$1" in
                --input) IN="$2" ;;
                --output) OUT_GFF="$2" ;;
                --engine) ENGINE="$2" ;;
                --type-map) EXTRA+=("--type-map" "$2") ;;
            esac
            shift 2 ;;
        --force|--repair-missing-parents) EXTRA+=("$1"); shift ;;
        *) die "Unknown normalize-gff option: $1" ;;
    esac
done
: "${IN:=${GENOME_GFF3_RAW:-${GENOME_GFF3_SOURCE:-${GENOME_GFF3}}}}"
open_log "00b_normalize_gff"
log_section "Step 00b: Conservative GFF3 normalization (engine=${ENGINE}; auto=builtin)"
[[ -s "${IN}" ]] || die "Input GFF not found: ${IN}"
if [[ -z "${OUT_GFF}" ]]; then
    base="$(basename "${IN}")"; base="${base%.gz}"; base="${base%.gff3}"; base="${base%.gff}"
    OUT_GFF="$(derived_ref_dir)/${base}.normalized.gff3"
fi
OUT_GFF="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${OUT_GFF}")"
# Read gzip by magic bytes, without a potentially stale decompression cache.
python3 "${CHROMAFLOW_HOME}/lib/steps/normalize_gff.py" "${IN}" "${OUT_GFF}" \
    --engine "${ENGINE}" "${EXTRA[@]}" || die "GFF3 normalization/validation failed; this result was not registered"
manifest_register GENOME_GFF3 "${OUT_GFF}"
echo "[DONE] ${OUT_GFF}"
echo "  Structure validation passed; GTF and TxDb compatibility still require separate checks."
