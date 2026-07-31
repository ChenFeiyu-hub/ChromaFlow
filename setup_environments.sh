#!/usr/bin/env bash
# ============================================================================
# setup_environments.sh
# ============================================================================
# 环境策略:
#   - nf-core pipelines: 通过 Apptainer 容器运行（不使用 conda profile）
#     备选: conda profile（网络受限时使用 02_run_nfcore_conda.sh）
#   - Nextflow + nf-core CLI: mamba 安装（仅运行引擎）
#   - 下游补充分析: mamba 安装（拆分为多个环境避免依赖冲突）
#
# 前置条件:
#   - mamba/miniforge 已安装
#   - Apptainer 已在系统级安装（HPC: module load apptainer）
#
# 经验总结（来自实际部署中的问题）:
#   - idr 最高支持 Python 3.10，且需要 numpy<1.24（np.int 废弃）
#   - MEME Suite 需独立环境 + 补装 bedtools/samtools
#   - R/Bioconductor 需显式安装 bioconductor-go.db、genomeinfodbdata
#     且必须锁定 python=3.12 防止 conda 解析出 3.14
#   - nf-core conda profile 下 deeptools 可能遇到 numpy 2.x 不兼容
#     (np.NAN 废弃)，需手动降级
#    pip install 'numpy<2.0'
#   - AGAT (Perl 生态) 与 ChromaFlow_ref 里的包经常解析冲突，故【单独成环境】
#     且默认不装：normalize-gff 没有 AGAT 时会自动改用内置 Python 引擎。
#
# 可选开关:
#   WITH_AGAT=true  bash setup_environments.sh   # 额外创建 ChromaFlow_agat
#   SKIP_R=true     bash setup_environments.sh   # 不装 R（改用 annotate-py）
# ============================================================================

set -euo pipefail

echo "============================================"
echo " ChromaFlow 环境安装"
echo " $(date)"
echo "============================================"

# ─────────────────────────────────────────────
# 检查前置条件
# ─────────────────────────────────────────────
if ! command -v mamba &>/dev/null; then
    echo "[ERROR] mamba 未安装"
    echo "  curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    echo "  bash Miniforge3-Linux-x86_64.sh"
    exit 1
fi

echo ""
echo "[CHECK] Apptainer / Singularity"
if command -v apptainer &>/dev/null; then
    echo "[  OK ] apptainer $(apptainer --version)"
elif command -v singularity &>/dev/null; then
    echo "[  OK ] singularity $(singularity --version)"
else
    echo "[WARN] 未检测到 Apptainer/Singularity"
    echo "  nf-core 将使用 conda profile 作为备选方案"
    echo "  安装 Apptainer: sudo apt install apptainer  或  module load apptainer"
fi

# ─────────────────────────────────────────────
# 辅助函数
# ─────────────────────────────────────────────
create_env() {
    local name="$1"; shift
    echo ""
    echo "[ENV] 创建 ${name}"
    if conda env list | grep -qw "^${name} "; then
        echo "[SKIP] ${name} 已存在"
        return 0
    fi
    mamba create -y -n "${name}" "$@"
    echo "[DONE] ${name}"
}

# =============================================================
# 环境 1: nfcore_env — Nextflow 运行引擎 + nf-core CLI
# =============================================================
create_env nfcore_env \
    -c bioconda -c conda-forge \
    nextflow \
    nf-core \
    openjdk=17

# =============================================================
# 环境 2: ChromaFlow_ref — 参考基因组处理 + 轻量下游工具
# =============================================================
# 用途: reference(00), samplesheet(01), annotate(05), spacing
# 包含 bedtools + samtools + python + matplotlib，是最常用的通用环境
# matplotlib-base 供 05 出图（基因组分布 / TSS 距离直方图）
# =============================================================
create_env ChromaFlow_ref \
    -c bioconda -c conda-forge \
    python=3.11 \
    samtools=1.21 \
    gffread \
    bedtools=2.31 \
    khmer \
    pigz \
    matplotlib-base

# =============================================================
# 环境 3: ChromaFlow_downstream — deepTools + HOMER 可视化
# =============================================================
create_env ChromaFlow_downstream \
    -c bioconda -c conda-forge \
    python=3.12 \
    deeptools=3.5.6 \
    homer \
    pygenometracks \
    csvtk \
    bedtools=2.31

# =============================================================
# 环境 3b: ChromaFlow_agat — AGAT（可选，默认不装）
# =============================================================
# 为什么单独成环境: AGAT 拉一整套 Perl/BioPerl 依赖，与 ChromaFlow_ref 中的
#   python/samtools/gffread 组合经常解析冲突（这是本项目实际踩过的坑）。
# 不装也完全没问题: chromaflow normalize-gff 会自动改用内置 Python 引擎
#   (lib/steps/normalize_gff.py)，补 exon / 推 UTR / 合成 gene 层。
# 需要 AGAT 的完整能力时:
#   WITH_AGAT=true bash setup_environments.sh
#   mamba activate ChromaFlow_agat && chromaflow normalize-gff --engine agat
# =============================================================
if [[ "${WITH_AGAT:-false}" == "true" ]]; then
    create_env ChromaFlow_agat -c bioconda -c conda-forge agat
else
    echo ""
    echo "[SKIP] 未创建 ChromaFlow_agat（AGAT 可选）"
    echo "       normalize-gff 将使用内置 Python 引擎；需要 AGAT 时:"
    echo "       WITH_AGAT=true bash setup_environments.sh"
fi

# =============================================================
# 环境 4: ChromaFlow_idr — IDR 分析
# =============================================================
# ⚠ 重要: idr 2.0.4.2 最高支持 Python 3.10
# ⚠ 安装后必须执行: pip install 'numpy<1.24'
#    (idr 使用了已废弃的 np.int，numpy ≥1.24 会报错)
# =============================================================
create_env ChromaFlow_idr \
    -c bioconda -c conda-forge \
    python=3.10 \
    idr \
    bedtools=2.31

echo ""
echo "[POST] ChromaFlow_idr: 降级 numpy 以兼容 idr"
mamba run -n ChromaFlow_idr pip install 'numpy<1.24' --quiet 2>/dev/null || true

# =============================================================
# 环境 5: ChromaFlow_motif — MEME Suite + bedtools
# =============================================================
# ⚠ 必须包含 bedtools 和 samtools（序列提取需要）
# =============================================================
create_env ChromaFlow_motif \
    -c bioconda -c conda-forge \
    python=3.12 \
    meme=5.5.7 \
    bedtools=2.31 \
    samtools=1.21

# =============================================================
# 环境 6: ChromaFlow_R — R/Bioconductor（可选，默认不安装）
# =============================================================
# 当前流程完全不依赖 R：注释(05)用 Python+bedtools，信号(06)用 deepTools。
# 此环境仅在需要 DiffBind 差异结合 / ChIPseeker 交叉验证时才需要。
# 历史上 genomeinfodbdata 常导致依赖解析失败，故默认跳过。
#
# 如需安装（已验证可行的配方）：
#   INSTALL_R=true bash setup_environments.sh
# 关键点：锁定 python=3.12（防解析出 3.14）；显式含 go.db + genomeinfodbdata；用镜像。
# =============================================================
# ChromaFlow_R 现在是 annotate(05, ChIPseeker) 的默认环境，默认安装。
# 无 R 需求（只用 annotate-py）时可 SKIP_R=true 跳过。
if [[ "${SKIP_R:-false}" == "true" ]]; then
    echo ""
    echo "[SKIP] 跳过 ChromaFlow_R。注意: chromaflow annotate 需要它；"
    echo "       无 R 时请改用 chromaflow annotate-py (ChromaFlow_ref)。"
    echo "       以后单独安装: bash env/install_chipseeker.sh [--method bioc] [--cn]"
else
    echo ""
    echo "[R] 安装 ChromaFlow_R (ChIPseeker + clusterProfiler)，调用专门安装器 + 自检"
    echo "    国内网络/BiocManager 混合法: bash env/install_chipseeker.sh --method bioc --cn"
    _HERE="$(cd "$(dirname "$0")" && pwd)"
    bash "${_HERE}/env/install_chipseeker.sh" \
        || echo "[R] 安装未成功，见上方；可改用 --method bioc 重试，或先用 annotate-py"
fi

# ─────────────────────────────────────────────
# Apptainer 缓存目录
# ─────────────────────────────────────────────
CACHE_DIR="${HOME}/.singularity_cache"
mkdir -p "${CACHE_DIR}"

# ─────────────────────────────────────────────
# nf-core 容器预下载指引
# ─────────────────────────────────────────────
echo ""
echo "============================================"
echo " nf-core 容器预下载指引 (Apptainer)"
echo "============================================"
echo ""
echo "  nf-core 通过 Apptainer 自动拉取容器镜像。"
echo "  首次运行会下载较多镜像 (~5-10 GB)，国内网络建议先预下载："
echo ""
echo "  mamba activate nfcore_env"
echo ""
echo "  # 方式 1: nf-core download（推荐，自动解析完整容器列表）"
echo "  nf-core download chipseq -r 2.1.0 \\"
echo "    --container-system singularity \\"
echo "    --container-cache-utilisation amend \\"
echo "    --compress none"
echo ""
echo "  nf-core download atacseq -r 2.1.2 \\"
echo "    --container-system singularity \\"
echo "    --container-cache-utilisation amend \\"
echo "    --compress none"
echo ""
echo "  # 方式 2: 网络受限时可使用 conda profile 作为备选"
echo "  在 chromaflow.config 设 NFCORE_PROFILE=conda，或直接运行 lib/steps/02_run_nfcore_conda.sh"
echo ""

# ─────────────────────────────────────────────
# ~/.bashrc 配置
# ─────────────────────────────────────────────
echo "============================================"
echo " 请将以下内容加入 ~/.bashrc"
echo "============================================"
echo ""
echo "  export APPTAINER_CACHEDIR=${CACHE_DIR}"
echo "  export NXF_SINGULARITY_CACHEDIR=${CACHE_DIR}"
echo ""

# ─────────────────────────────────────────────
# 验证
# ─────────────────────────────────────────────
echo "============================================"
echo " 验证步骤"
echo "============================================"
echo ""
echo "  mamba activate nfcore_env   && nextflow -version"
echo "  mamba activate ChromaFlow_ref && bedtools --version && python -c 'import matplotlib; print(matplotlib.__version__)'"
echo "  # 可选: mamba activate ChromaFlow_agat && agat_convert_sp_gxf2gxf.pl --help | head -3"
echo "  mamba activate ChromaFlow_downstream && deeptools --version"
echo "  mamba activate ChromaFlow_idr && idr --version"
echo "  mamba activate ChromaFlow_motif && meme --version && bedtools --version"
echo "  # 仅当 INSTALL_R=true 时: mamba activate ChromaFlow_R && R -e 'library(ChIPseeker)'"
echo ""

echo "============================================"
echo " 环境一览"
echo "============================================"
echo ""
echo "  nfcore_env              Nextflow + nf-core CLI"
echo "  ChromaFlow_ref          samtools, gffread, bedtools, python, matplotlib (通用)"
echo "  ChromaFlow_agat         AGAT (可选, WITH_AGAT=true 才创建; 不装则用内置引擎)"
echo "  ChromaFlow_downstream   deepTools, HOMER, pyGenomeTracks"
echo "  ChromaFlow_idr          IDR (Python 3.10, numpy<1.24)"
echo "  ChromaFlow_motif        MEME Suite + bedtools + samtools (motif, candidates)"
echo "  ChromaFlow_R            ChIPseeker + clusterProfiler + GO.db (annotate, 默认安装)"
echo ""

echo "============================================"
echo " 命令 → 环境 对照表"
echo "============================================"
echo ""
echo "  chromaflow reference    → ChromaFlow_ref"
echo "  chromaflow normalize-gff → ChromaFlow_ref (内置引擎) / ChromaFlow_agat (--engine agat)"
echo "  chromaflow samplesheet  → (任意，仅需 python3)"
echo "  chromaflow run          → nfcore_env (Apptainer / conda)"
echo "  chromaflow idr          → ChromaFlow_idr"
echo "  chromaflow motif        → ChromaFlow_motif"
echo "  chromaflow annotate     → ChromaFlow_R (ChIPseeker + clusterProfiler)"
echo "  chromaflow annotate-py  → ChromaFlow_ref (无 R 时的 Python 备选)"
echo "  chromaflow candidates   → ChromaFlow_downstream (读04 FIMO)"
echo "  chromaflow lookup       → ChromaFlow_downstream (先验基因反查)"
echo "  chromaflow signal       → ChromaFlow_downstream (可选 deepTools 聚合图)"
echo "  chromaflow tracks       → ChromaFlow_downstream (pyGenomeTracks)"
echo "  chromaflow go-prep      → ChromaFlow_ref"
echo "  chromaflow go-plot      → ChromaFlow_ref"
echo "  chromaflow spacing      → ChromaFlow_ref"
echo ""
echo "[DONE] $(date)"