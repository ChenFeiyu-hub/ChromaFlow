#!/usr/bin/env bash
set -euo pipefail
source "${CHROMAFLOW_HOME:?}/lib/common.sh"

PIPELINE="${1:-chipseq}"; shift || true
ACTION="${1:-}"; shift || true
if [[ -n "${CHROMAFLOW_CONFIG:-}" ]]; then
    load_config "${CHROMAFLOW_CONFIG}"
fi
case "${PIPELINE}" in
    chipseq|dapseq) NAME=chipseq; REVISION="${NFCORE_CHIPSEQ_VERSION:-2.1.0}" ;;
    atacseq) NAME=atacseq; REVISION="${NFCORE_ATACSEQ_VERSION:-2.1.2}" ;;
    *) die "Unknown --type: ${PIPELINE}" ;;
esac
case "${ACTION}" in prepare|check|pack) ;; *) die "Usage: chromaflow offline prepare|check|pack --bundle <directory>; see chromaflow help offline" ;; esac
check_tool python3 || exit 1
ARGS=("${ACTION}")
if [[ "${ACTION}" == prepare ]]; then
    ARGS+=(--pipeline "${NAME}" --revision "${REVISION}")
fi
BUNDLE="${CHROMAFLOW_BUNDLE:-${NFCORE_OFFLINE_BUNDLE:-}}"
[[ -z "${BUNDLE}" ]] || ARGS+=(--bundle "${BUNDLE}")
exec python3 "${CHROMAFLOW_HOME}/lib/offline_bundle.py" "${ARGS[@]}" "$@"
