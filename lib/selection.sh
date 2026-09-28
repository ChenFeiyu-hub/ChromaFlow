# Shared candidates/lookup filter arguments. Caller loads common.sh first.
selection_defaults() {
    export CFSEL_MOTIFS="${MOTIF_IDS:-}" CFSEL_QVAL="${FIMO_QVAL:-}"
    export CFSEL_EVALUE="${MOTIF_EVALUE_MAX:-}" CFSEL_CENTRIMO="${CENTRIMO_EVALUE_MAX:-}"
}
selection_option() {
    [[ -n "${2:-}" ]] || die "Missing value for $1"
    case "$1" in
        --motif) export CFSEL_MOTIFS="$2" ;;
        --qval) export CFSEL_QVAL="$2" ;;
        --motif-evalue) export CFSEL_EVALUE="$2" ;;
        --centrimo-evalue) export CFSEL_CENTRIMO="$2" ;;
    esac
}
