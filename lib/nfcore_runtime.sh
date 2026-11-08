#!/usr/bin/env bash
# Shared runtime selection and diagnostic logging for both nf-core runners.

configure_nfcore_runtime() {
    NFCORE_NEXTFLOW=nextflow
    NFCORE_RUN_TARGET="nf-core/${NFCORE_NAME}"
    NFCORE_VERSION_ARGS=(-r "${NFCORE_VERSION}")
    NFCORE_GLOBAL_ARGS=()
    local offline="${CHROMAFLOW_OFFLINE:-${NFCORE_OFFLINE:-${NXF_OFFLINE:-false}}}"
    case "${offline}" in true|false) ;; *) die "NFCORE_OFFLINE must be true or false" ;; esac
    if [[ "${offline}" == true ]]; then
        local bundle="${CHROMAFLOW_BUNDLE:-${NFCORE_OFFLINE_BUNDLE:-}}"
        [[ -n "${bundle}" ]] || die "Offline execution requires --bundle <directory> or NFCORE_OFFLINE_BUNDLE"
        check_tool python3 || return 1
        local envfile="${OUTPUT_DIR}/.chromaflow/offline-runtime.env"
        python3 "${CHROMAFLOW_HOME}/lib/offline_bundle.py" activate \
            --bundle "${bundle}" --pipeline "${NFCORE_NAME}" --revision "${NFCORE_VERSION}" \
            --profile "${PROFILE}" --state "${OUTPUT_DIR}/.chromaflow" --output "${envfile}" || return 1
        # Values are generated locally with shlex.quote, never copied from a bundle script.
        source "${envfile}"
        NFCORE_VERSION_ARGS=()
        # Exclude implicit user/project Nextflow configs; use only the fixed local configs.
        NFCORE_GLOBAL_ARGS=(-C "${NFCORE_RUN_TARGET}/nextflow.config,${NFCORE_OFFLINE_CONFIG}")
        echo "[OFFLINE] Local pipeline: ${NFCORE_RUN_TARGET}"
        echo "[OFFLINE] Plugins: ${NXF_PLUGINS_DIR}"
        echo "[OFFLINE] Nextflow ${NXF_VER}; Automatic version checks and plugin downloads are disabled"
    fi
    if [[ -n "${CF_BATCH_EXECUTOR_CONFIG:-}" ]]; then
        [[ -s "${CF_BATCH_EXECUTOR_CONFIG}" ]] || die 'Batch executor config is missing'
        # Nextflow ignores -c when -C is present. Keep the batch limits in the
        # explicit offline config list instead of silently dropping them.
        if [[ "${offline}" == true ]]; then
            NFCORE_GLOBAL_ARGS=(-C "${NFCORE_RUN_TARGET}/nextflow.config,${NFCORE_OFFLINE_CONFIG},${CF_BATCH_EXECUTOR_CONFIG}")
        else
            NFCORE_GLOBAL_ARGS+=(-c "${CF_BATCH_EXECUTOR_CONFIG}")
        fi
        export NXF_SYNTAX_PARSER=v1
    fi
    NFCORE_DEBUG_LOG="${LOG_DIR}/nextflow_${PIPELINE}_${TIMESTAMP}_$$.log"
    echo "[LOG] Nextflow debug log: ${NFCORE_DEBUG_LOG}"
}

print_nfcore_command() {
    printf '%q ' "${NFCORE_NEXTFLOW}" -log "${NFCORE_DEBUG_LOG}" \
        ${NFCORE_GLOBAL_ARGS[@]+"${NFCORE_GLOBAL_ARGS[@]}"} run "${NFCORE_RUN_TARGET}" "$@"
    printf '\n'
}

execute_nfcore() {
    local status=0
    "${NFCORE_NEXTFLOW}" -log "${NFCORE_DEBUG_LOG}" \
        ${NFCORE_GLOBAL_ARGS[@]+"${NFCORE_GLOBAL_ARGS[@]}"} run "${NFCORE_RUN_TARGET}" "$@" || status=$?
    if [[ ${status} -ne 0 ]]; then
        log_error "Nextflow exit code ${status}; provide this step log and ${NFCORE_DEBUG_LOG}"
    fi
    return "${status}"
}
