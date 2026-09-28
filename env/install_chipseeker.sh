#!/usr/bin/env bash
# ============================================================================
# install_chipseeker.sh — 独立安装 ChromaFlow_R (ChIPseeker + 可视化)
# ============================================================================
# 用法:
#   bash install_chipseeker.sh                 # 方法A: 硬化 conda（默认，先试这个）
#   bash install_chipseeker.sh --method bioc   # 方法B: conda 装 R+基础设施, BiocManager 装 ChIPseeker
#   bash install_chipseeker.sh --cn            # 用清华镜像（国内网络）
#   bash install_chipseeker.sh --no-test       # 不跑安装自检
# ----------------------------------------------------------------------------
# 解决 genomeinfodbdata 失败的关键：
#   · conda-forge 频道排在 bioconda 前 + strict priority
#   · python / r-base 显式锁版本，避免解析漂移
#   · 数据包(genomeinfodbdata/go.db)若 conda 卡住，方法B 交给 BiocManager 从
#     Bioconductor 官方渠道装（自带 release 一致性与下载重试），绕开 bioconda post-link
# ============================================================================
set -euo pipefail

ENV_NAME="ChromaFlow_R"
METHOD="conda"          # conda | bioc
USE_CN=false
RUN_TEST=true
R_VER="4.4"
BIOC_VER="3.20"         # 与 r-base 4.4 匹配的 Bioconductor 发行
PY_VER="3.12"
CONDA=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --method) METHOD="$2"; shift 2 ;;
        --solver) CONDA="$2"; shift 2 ;;
        --cn)     USE_CN=true; shift ;;
        --no-test) RUN_TEST=false; shift ;;
        --name)   ENV_NAME="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# 选 conda 前端
if [[ -z "${CONDA}" ]]; then
    for c in mamba micromamba conda; do command -v "$c" &>/dev/null && { CONDA="$c"; break; }; done
fi
case "${CONDA}" in mamba|micromamba|conda) ;; *) echo '[ERROR] Invalid solver'; exit 2 ;; esac
command -v "${CONDA}" >/dev/null || { echo "[ERROR] Solver is unavailable: ${CONDA}"; exit 1; }
[[ -n "${CONDA}" ]] || { echo "[ERROR] mamba/conda was not found"; exit 1; }
echo "[INFO] Using ${CONDA}; method=${METHOD}; environment=${ENV_NAME}"

# 镜像
if [[ "${USE_CN}" == "true" ]]; then
    CF_CH="https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge"
    BC_CH="https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda"
    CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
    BIOC_MIRROR="https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
else
    CF_CH="conda-forge"; BC_CH="bioconda"
    CRAN="https://cloud.r-project.org"; BIOC_MIRROR="https://bioconductor.org"
fi

# ─────────────── 方法 A：硬化 conda（一把梭）───────────────
install_conda() {
    echo "[A] Conda installation (conda-forge first, strict channel priority, pinned versions)"
    # 关键：-c conda-forge 在 -c bioconda 之前；--strict-channel-priority
    ${CONDA} create -y -n "${ENV_NAME}" \
        --override-channels -c "${CF_CH}" -c "${BC_CH}" --strict-channel-priority \
        "python=${PY_VER}" "r-base=${R_VER}" \
        bioconductor-chipseeker \
        bioconductor-clusterprofiler \
        bioconductor-genomicfeatures \
        bioconductor-txdbmaker \
        bioconductor-rtracklayer \
        bioconductor-org.hs.eg.db \
        bioconductor-go.db \
        bioconductor-genomeinfodbdata \
        r-ggplot2 r-biocmanager \
        bedtools \
    || { echo "[A FAILED] Check dependency solving and genomeinfodbdata; alternative: bash $0 --method bioc"; return 1; }
}

# ─────────────── 方法 B：conda 装 R+基础设施，BiocManager 装 ChIPseeker ───────────────
install_bioc() {
    echo "[B] Install R/build tools with Conda; install ChIPseeker/data packages with BiocManager"
    # 从 bioconda 拿「代码类」重型包的二进制（这些一般不出问题），
    # 数据类包(genomeinfodbdata/go.db)与 ChIPseeker 由 BiocManager 处理。
    ${CONDA} create -y -n "${ENV_NAME}" \
        --override-channels -c "${CF_CH}" -c "${BC_CH}" --strict-channel-priority \
        "python=${PY_VER}" "r-base=${R_VER}" r-biocmanager \
        bioconductor-genomicranges bioconductor-rtracklayer \
        bioconductor-genomicfeatures bioconductor-txdbmaker \
        r-ggplot2 r-xml r-rcurl \
        c-compiler cxx-compiler make pkg-config \
        curl libxml2 openssl zlib bzip2 xz \
    || { echo "[B FAILED] Conda base environment installation failed; see errors above"; return 1; }

    echo "[B] Install ChIPseeker using BiocManager inside the environment"
    ${CONDA} run -n "${ENV_NAME}" Rscript -e "
      options(repos=c(CRAN='${CRAN}'), timeout=600)
      options(BioC_mirror='${BIOC_MIRROR}')
      if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')
      BiocManager::install(version='${BIOC_VER}', ask=FALSE, update=FALSE)
      BiocManager::install(c('ChIPseeker','clusterProfiler','GO.db','GenomeInfoDbData'),
                           ask=FALSE, update=FALSE)
    " || { echo "[B FAILED] BiocManager installation failed; see errors above"; return 1; }
}

# ─────────────── 执行 ───────────────
if "${CONDA}" env list 2>/dev/null | grep -qE "^\s*${ENV_NAME}\s"; then
    echo "[WARN] Environment ${ENV_NAME} already exists; kept unchanged. Use --name to select a new environment."
    exit 1
fi

case "${METHOD}" in
    conda) install_conda ;;
    bioc)  install_bioc ;;
    *) echo "Unknown --method: ${METHOD}"; exit 1 ;;
esac

# ─────────────── 自检 ───────────────
if [[ "${RUN_TEST}" == "true" ]]; then
    echo ""
    echo "[TEST] Run installation self-test..."
    HERE="$(cd "$(dirname "$0")" && pwd)"
    ${CONDA} run -n "${ENV_NAME}" Rscript "${HERE}/test_chipseeker.R" \
        && echo "[✅ PASS] ChIPseeker environment is ready" \
        || { echo "[❌ FAIL] Self-test failed; see stage-specific errors above"; exit 1; }
fi

echo ""
echo "Completed. Activate with: conda activate ${ENV_NAME}"
