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

while [[ $# -gt 0 ]]; do
    case "$1" in
        --method) METHOD="$2"; shift 2 ;;
        --cn)     USE_CN=true; shift ;;
        --no-test) RUN_TEST=false; shift ;;
        --name)   ENV_NAME="$2"; shift 2 ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

# 选 conda 前端
CONDA=""
for c in mamba micromamba conda; do command -v "$c" &>/dev/null && { CONDA="$c"; break; }; done
[[ -n "${CONDA}" ]] || { echo "[ERROR] 未找到 mamba/conda"; exit 1; }
echo "[INFO] 使用 ${CONDA}；方法=${METHOD}；环境=${ENV_NAME}"

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
    echo "[A] conda 一次性安装（conda-forge 优先 + strict priority + 锁版本）"
    # 关键：-c conda-forge 在 -c bioconda 之前；--strict-channel-priority
    ${CONDA} create -y -n "${ENV_NAME}" \
        -c "${CF_CH}" -c "${BC_CH}" --strict-channel-priority \
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
    || { echo "[A 失败] 多半卡在 genomeinfodbdata/solve。改用: bash $0 --method bioc"; return 1; }
}

# ─────────────── 方法 B：conda 装 R+基础设施，BiocManager 装 ChIPseeker ───────────────
install_bioc() {
    echo "[B] conda 只装 R + 编译基础设施，ChIPseeker 及数据包交给 BiocManager"
    # 从 bioconda 拿「代码类」重型包的二进制（这些一般不出问题），
    # 数据类包(genomeinfodbdata/go.db)与 ChIPseeker 由 BiocManager 处理。
    ${CONDA} create -y -n "${ENV_NAME}" \
        -c "${CF_CH}" -c "${BC_CH}" --strict-channel-priority \
        "python=${PY_VER}" "r-base=${R_VER}" r-biocmanager \
        bioconductor-genomicranges bioconductor-rtracklayer \
        bioconductor-genomicfeatures bioconductor-txdbmaker \
        r-ggplot2 r-xml r-rcurl \
        c-compiler cxx-compiler make pkg-config \
        curl libxml2 openssl zlib bzip2 xz \
    || { echo "[B 失败] conda 基础设施装不上，见上方报错"; return 1; }

    echo "[B] 在环境内用 BiocManager 装 ChIPseeker（数据包走 Bioconductor 官方渠道）"
    ${CONDA} run -n "${ENV_NAME}" Rscript -e "
      options(repos=c(CRAN='${CRAN}'), timeout=600)
      options(BioC_mirror='${BIOC_MIRROR}')
      if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')
      BiocManager::install(version='${BIOC_VER}', ask=FALSE, update=FALSE)
      BiocManager::install(c('ChIPseeker','clusterProfiler','GO.db','GenomeInfoDbData'),
                           ask=FALSE, update=FALSE)
    " || { echo "[B 失败] BiocManager 安装出错，见上方"; return 1; }
}

# ─────────────── 执行 ───────────────
if "${CONDA}" env list 2>/dev/null | grep -qE "^\s*${ENV_NAME}\s"; then
    echo "[WARN] 环境 ${ENV_NAME} 已存在。删除重装请先: ${CONDA} env remove -n ${ENV_NAME}"
    exit 1
fi

case "${METHOD}" in
    conda) install_conda ;;
    bioc)  install_bioc ;;
    *) echo "未知 --method: ${METHOD}"; exit 1 ;;
esac

# ─────────────── 自检 ───────────────
if [[ "${RUN_TEST}" == "true" ]]; then
    echo ""
    echo "[TEST] 运行安装自检..."
    HERE="$(cd "$(dirname "$0")" && pwd)"
    ${CONDA} run -n "${ENV_NAME}" Rscript "${HERE}/test_chipseeker.R" \
        && echo "[✅ PASS] ChIPseeker 环境可用" \
        || { echo "[❌ FAIL] 自检未过，见上方分阶段报错"; exit 1; }
fi

echo ""
echo "完成。使用: mamba activate ${ENV_NAME}"
