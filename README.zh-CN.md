# ChromaFlow（中文文档）

[English README](README.md) · **中文**

ChIP-seq / DAP-seq / ATAC-seq 下游分析框架，基于 nf-core pipelines，补齐从 peak 到**候选靶基因**的完整证据链。

以**软件**方式组织：代码装一处，数据与结果由每个项目自己的 `chromaflow.config` 决定，代码目录零污染。

---

## 目录

- [设计理念](#设计理念软件--配置--运行-三层分离)
- [安装](#安装)
- [快速开始](#快速开始)
- [帮助系统](#帮助系统)
- [核心概念：peakset](#核心概念peakset决定下游用哪套-peak)
- [命令一览](#命令一览)
- [分析流程详解](#分析流程详解)
- [可选步骤与跳过方式](#可选步骤与跳过方式)
- [样本表与双端配对](#样本表与双端配对)
- [配置与文件自动识别](#配置与文件自动识别)
- [典型分析场景](#典型分析场景recipes)
- [结果文件导读](#结果文件导读)
- [环境一览](#环境一览)
- [质量评估与结果解读](#质量评估与结果解读)
- [适配其他物种](#适配其他物种)

---

## 设计理念：软件 / 配置 / 运行 三层分离

| 层 | 内容 | 位置 |
|---|---|---|
| ① 软件本体 | CLI、lib、steps、templates | **本仓库**（clone 一次，加进 PATH） |
| ② 项目配置 | 物种参数、输入/输出路径 | 每个数据集一份 `chromaflow.config` |
| ③ 运行 I/O | fastq、参考、结果、日志 | 由配置的 `OUTPUT_DIR` 决定 |

本仓库**只含代码**，不含任何 `results/`、`reference/`、`logs/`。换一个数据集 = 复制一份 config 改路径，代码目录完全不动。

```
ChromaFlow/
├── bin/chromaflow              # 统一命令行入口
├── lib/
│   ├── common.sh               # 配置解析 + 路径派生/自动识别 + peakset + 日志
│   ├── help.sh                 # chromaflow help [command]
│   └── steps/                  # 各分析步骤
│       ├── 00_prepare_reference.sh   00b_normalize_gff.sh
│       ├── normalize_gff.py          make_samplesheet.py
│       ├── 01_make_samplesheet.sh    02_run_nfcore.sh
│       ├── 03_idr_analysis.sh        04_meme_motif.sh
│       ├── 05_peak_annotation.R/.sh  06_integrate.sh
│       ├── 07_gene_tracks.sh         08_go_analysis.sh
│       ├── 09_lookup.sh              plot_tracks.sh  ← 纯绘图核心
│       └── support/
├── templates/                  # 项目模板（config + samples_*）
└── setup_environments.sh       # conda/mamba 环境安装
```

---

## 安装

```bash
git clone <repo> ~/software/ChromaFlow
export CHROMAFLOW_HOME=~/software/ChromaFlow
export PATH="$CHROMAFLOW_HOME/bin:$PATH"     # 建议写入 ~/.bashrc

bash $CHROMAFLOW_HOME/setup_environments.sh                 # 创建各 conda 环境
WITH_AGAT=true bash $CHROMAFLOW_HOME/setup_environments.sh # 可选：单独装 AGAT 环境
SKIP_R=true    bash $CHROMAFLOW_HOME/setup_environments.sh # 不装 R（改用 annotate-py）
bash $CHROMAFLOW_HOME/env/install_chipseeker.sh            # 单独装 R/ChIPseeker（较慢）
```

---

## 快速开始

```bash
# 1. 建项目（生成 chromaflow.config 模板）
chromaflow init /data/projects/MyChIP --type chipseq
cd /data/projects/MyChIP
vim chromaflow.config        # 只填物种 + 你手里已有的源文件路径
vim samples_chipseq.txt      # 样本与 fastq 路径

# 2. 参考准备（可选：已有基因组/GTF/chrom.sizes 可整步跳过）
mamba activate ChromaFlow_ref
chromaflow reference          # 00  基因组 / GTF / chrom.sizes
chromaflow normalize-gff      # 00b 规范化 GFF（可选，注释更准；无 AGAT 会用内置引擎）
chromaflow doctor             # 自检：文件、工具是否就位

# 3. 主流程
chromaflow samplesheet        # 01  生成 nf-core 样本表
mamba activate nfcore_env
chromaflow run                # 02  跑 nf-core/chipseq

# 4. 下游（核心）
mamba activate ChromaFlow_idr        && chromaflow idr        # 03 可重复 peak
mamba activate ChromaFlow_motif      && chromaflow motif      # 04 motif 发现
mamba activate ChromaFlow_R          && chromaflow annotate   # 05 peak 注释
mamba activate ChromaFlow_downstream && chromaflow candidates # 06 候选靶基因
```

全局参数（任何命令可用）：

| 参数 | 作用 |
|---|---|
| `-c, --config <file>` | 指定配置（默认 `$CHROMAFLOW_CONFIG` > `./chromaflow.config`） |
| `-t, --type <type>` | `chipseq`(默认) / `atacseq` / `dapseq` |
| `-o, --outdir <dir>` | 把**本步输出**重定向到指定目录 |
| `--peakset <spec>` | 下游用哪套 peak（见下节） |
| `--profile <name>` | 02 的执行模式：`singularity`(默认)/`apptainer`/`docker`/`podman`/`conda`/`mamba` |
| `-h, --help` | 当前命令的详细帮助 |

---

## 帮助系统

```bash
chromaflow help                 # 命令总览
chromaflow help <command>       # 单个命令的详细说明
chromaflow <command> --help     # 同上
chromaflow help config          # config 里到底该写什么
chromaflow help samples         # 样本表格式 + 双端配对规则
chromaflow help peakset         # peakset 语义
chromaflow help envs            # 每步需要哪个 conda 环境
chromaflow version
```

> 英文文档为主文档（[README.md](README.md)），本文件为中文版，内容同步。

---

## 核心概念：peakset（决定下游用哪套 peak）

下游所有分析（04/05/06/07/09）都要回答同一个问题：**用哪一批 peak？** 这由 `--peakset` 控制，**命令行指定、不写进 config**，方便随时对比。

| spec | 含义 | 输出目录后缀 |
|---|---|---|
| `idr`（默认） | IDR<0.01 的可重复 peak | 无 |
| `idr:0.05` / `idr:0.1` | 更宽松的 IDR 阈值 | `_idr0.05` |
| `consensus` | 跳过 IDR，用 nf-core consensus（自动回填 signalValue） | `_consensus` |
| `<文件路径>` | 任意自定义 peak 文件 | `_custom` |

不同 peakset 写入**带后缀的独立目录**（`annotation_consensus/`、`candidates_consensus/` …），可并存对比、互不覆盖。

> **重要规则**：要分析某套 peak，就用**同一个 `--peakset`** 依次跑完 `motif → annotate → candidates → lookup`。
> motif 是在特定 peak 序列上提取的；若 04 与 06 用了不同 peakset，FIMO 命中会对不上 peak，候选结果失真。

```bash
# 严格路线（可重复 peak）
chromaflow motif && chromaflow annotate && chromaflow candidates

# 宽松路线（重复性不足时妥协）
chromaflow motif      --peakset consensus
chromaflow annotate   --peakset consensus
chromaflow candidates --peakset consensus
```

---

## 命令一览

| 命令 | 步骤 | 环境 |
|---|---|---|
| `chromaflow init` | 建项目骨架 + config 模板 | — |
| `chromaflow reference` | 00 参考基因组准备（**可选**） | `ChromaFlow_ref` |
| `chromaflow normalize-gff` | 00b 规范化 GFF（**可选**） | `ChromaFlow_ref` / `ChromaFlow_agat` |
| `chromaflow samplesheet` | 01 生成 nf-core 样本表 | — |
| `chromaflow run` | 02 运行 nf-core 主流程 | `nfcore_env` |
| `chromaflow idr` | 03 IDR 重复一致性 | `ChromaFlow_idr` |
| `chromaflow motif` | 04 MEME-ChIP motif 发现 | `ChromaFlow_motif` |
| `chromaflow annotate` | 05 Peak 注释（ChIPseeker） | `ChromaFlow_R` |
| `chromaflow annotate-py` | 05 备选：纯 Python 注释 | `ChromaFlow_ref` |
| `chromaflow candidates` | 06 候选靶基因整合 | `ChromaFlow_downstream` |
| `chromaflow lookup` | 09 先验基因反查 | `ChromaFlow_downstream` |
| `chromaflow tracks` | 07 单基因浏览器波峰图 | `ChromaFlow_downstream` |
| `chromaflow signal` | (可选) deepTools 聚合信号图 | `ChromaFlow_downstream` |
| `chromaflow go-prep` / `go-plot` | 08 GO 富集导出 / 结果可视化 | `ChromaFlow_ref` |
| `chromaflow spacing` | (support) motif 间距分析 | `ChromaFlow_ref` |
| `chromaflow doctor` | 环境与文件自检 | 任意 |

---

## 分析流程详解

### 00 参考基因组准备（可选步骤）
整理基因组 FASTA、GTF、chrom.sizes；可选把细胞器基因组（叶绿体/线粒体）并入核基因组，便于识别细胞器来源的假 peak。

源文件**压缩与否都行，且按文件头魔数判断而不是看扩展名**：已解压的 `.fa`/`.gff3` 原地使用（不做多余拷贝），gzip/bgzip、bzip2、xz、zstd、zip 自动解压一次并在后续复用；「名字叫 `.gff3` 实际是 gzip」这种下载常见的情况也能正确识别，目录自动扫描也不会把压缩文件误选给需要纯文本的步骤。

**产物不写回 config**：本步生成的文件登记到 `<OUTPUT_DIR>/.chromaflow/reference.manifest`，后续步骤自动识别。已经有基因组/GTF/chrom.sizes 的话，直接放进 `REFERENCE_DIR` 并跳过本步即可。

细胞器缺失不再中断流程：缺哪个跳哪个（只 WARN），完全不需要时用 `--skip-organelle` 或 `ORGANELLE_SOURCES=()`。

```bash
chromaflow reference --skip-organelle   # 没有合适的叶绿体/线粒体基因组
chromaflow reference --skip-gtf         # 只跑下游，不需要 GTF
chromaflow reference --strict           # 反过来：细胞器缺失就报错
```

### 00b GFF 规范化（可选步骤）— 注释准确性的前提
Phytozome 等来源的 GFF 可能**缺 `exon` 行**或 UTR 命名不统一，导致 ChIPseeker 建的基因模型不完整、UTR/exon/intron 判定错误（典型症状：注释饼图里 5'UTR/3'UTR 占比为 0）。

```bash
chromaflow normalize-gff                    # 自动：有 AGAT 用 AGAT，没有就用内置引擎
chromaflow normalize-gff --engine builtin   # 强制内置（纯 Python，零依赖）
chromaflow normalize-gff --engine agat      # 强制 AGAT
```

AGAT 好用但 Perl 依赖链常与其它包冲突，因此**单独成环境、默认不装**。没有 AGAT 时，内置规范化器（`lib/steps/normalize_gff.py`）会：从 CDS 补 exon、由 exon−CDS 推 5'/3' UTR、合成缺失的 gene 层、排序输出。

产物会自动登记为当前 `GENOME_GFF3`，**不需要回头改 config**。整步跳过也没关系：05 会退回用原始 GFF3，只是提示可能缺 exon 层级。

### 实测量：基因组大小、读长、GC

这三个原本是 config 里手填的常数。改为自动实测，理由很简单：手填的数在基因组换版本、测序换批次之后没人会回头改。

| 参数 | 谁来测 | 怎么测 | 真正的用途 |
|---|---|---|---|
| `MACS_GSIZE` | `reference` (00) | khmer `unique-kmers.py -k <读长>`，就是以前日志里让你手动跑的那条命令 | MACS 的 Poisson 背景 λ = N·d/gsize |
| `READ_LENGTH` | `samplesheet` (01) | 每个 fastq 抽 2000 条 read 取最大值 | 只在缺 `--macs_gsize` 时被 nf-core 用来推算基因组大小 |
| `GENOME_GC` | `reference` (00) | 扫一遍 FASTA | 仅供参考，见下方说明 |

结果写进 `<OUTPUT_DIR>/.chromaflow/derived.env`，**不回写你的 config**：

```
config      = 你的声明   （永远优先；一旦写了某个值，该值就不再实测）
derived.env = 脚本的观测 （输入变了会重新测）
```

这与参考文件 manifest 是同一套分离逻辑。好处是 config 始终是你手写的、可以安心进 git 的声明，不会被脚本在背后改写；而实测值随 outputs 一起存在，删掉 `outputs/` 就自然失效重算。确实想把它们归档进 config 的话，`chromaflow reference --write-config` 会以带标记的块追加（幂等、自动备份、值带引号）。

实测结果按「基因组文件签名 + k」缓存，只测一次；`--force` 重测，`--no-measure` 跳过。

**没装 khmer 时**：退到非 N 碱基计数（上界，偏大约 5–10%）；此时 02 会改为不传 `--macs_gsize`、只传 `--read_length`，由 nf-core 自己算——nf-core/chipseq 2.x 内部的 `KHMER_UNIQUEKMERS` 用的正是同一个工具。两条路都不需要你手填。

**这几个数真的有用吗？**

- `MACS_GSIZE` 有用，但很宽容。它只是缩放 Poisson 背景，几个百分点的误差几乎不会改变 peak 列表——比如衣藻用占位的 1.0e8 对上实测的 103,854,640，差 3.7%，实际影响可以忽略；差一倍才要紧。
- `READ_LENGTH` 在已经给了 `MACS_GSIZE` 时是**完全不起作用**的：nf-core 只在缺 `--macs_gsize` 时用它推算基因组大小。仍然测它，是因为几乎不花时间，且能让 k-mer 口径与真实数据一致。
- `GC_CONTENT` 已**删除**。它从来没有被用到：04 的背景模型是 `fasta-get-markov` 从真实序列估计碱基与二核苷酸频率得到的，config 里那个值只出现在一行日志里——也就是说日志会一本正经地打印你手填的任何数字，哪怕跟基因组对不上。现在 04 改为从它实际用的背景模型里读回 GC 并显示。

### 01–02 样本表与 nf-core 主流程
`samplesheet` 依据样本记录表生成 nf-core 样本表（格式见[样本表与双端配对](#样本表与双端配对)）；`run` 调用 nf-core/chipseq（或 atacseq）。这一步产出 BAM、bigWig、MACS peak、consensus peak。

**执行模式由 `NFCORE_PROFILE` 控制**（命令行 `--profile` 优先级更高）。它既是 nextflow 的 `-profile`，也决定 `chromaflow run` 走哪条执行路径：

| 取值 | 模式 | 说明 |
|---|---|---|
| `singularity`（默认） | 容器 | 自动拉取镜像，可重复性最好 |
| `apptainer` | 容器 | 装的是 Apptainer 1.x 时用 |
| `docker` / `podman` | 容器 | 换成对应引擎 |
| `conda` / `mamba` | conda | 为每个 process 自建环境，首次约 30–60 分钟；适合拉不动镜像的受限网络 |

```bash
chromaflow run --profile conda      # 临时用一次
NFCORE_PROFILE=conda                # 或写进 config 作为该项目默认
```

`NFCORE_PROFILE="singularity,test"` 这类组合会原样传给 nextflow。镜像缓存在 `$OUTPUT_DIR/.singularity_cache`（同时导出为 `NXF_SINGULARITY_CACHEDIR` 和 `NXF_APPTAINER_CACHEDIR`）；conda 环境缓存在 `$OUTPUT_DIR/.conda_envs`。conda 模式默认把 channel 改成清华镜像，`NFCORE_CONDA_MIRROR=none` 可保留你自己的 channel 配置。

### 03 IDR 重复一致性
对每两个重复做 IDR，筛出**统计上可重复**的 peak。输出 8 列：
`chr start end name signalValue strand globalIDR n_pairs`

```bash
chromaflow idr                          # 默认：任一配对通过即保留（并集）
chromaflow idr --min-pairs 2            # 只保留 ≥2 个配对都支持的高置信核心集
chromaflow idr --exclude-rep REP1       # 排除某个重复后重跑
chromaflow idr --exclude-rep REP1,REP3  # 可排除多个
```

- **`--min-pairs N`**：`n_pairs` 记录每个 peak 被多少个两两配对支持。重复数 ≥3 时，`--min-pairs 2` 是比"丢掉某个重复"更稳妥的收紧方式——不扔数据，只提高证据要求。
- **`--exclude-rep`**：怀疑某个重复异常时用它实证检验影响（**先看日志里的 pairwise 表格再决定**，见[质量评估](#质量评估与结果解读)）。

运行结束会打印**可重复性分层**（各 peak 被几个配对支持），据此判断核心集大小。

### 04 MEME-ChIP motif 发现
在 peak summit±100bp 序列上做 de novo motif 发现 + Tomtom 已知库比对 + CentriMo 中心富集 + FIMO 全序列扫描。

- 使用**全部 peak 序列**（不截断 top-N）。`MEME_SEARCHSIZE=0` 表示不限制；序列很多时 MEME 较慢（运行时间约随 searchsize 平方增长），可设 `MEME_SEARCHSIZE=200000` 限速。
- 输出目录随 peakset 带后缀：`motif/`、`motif_consensus/` …
- MEME 编号出现**跳号**（如 MEME-1,3,4）是正常的：MEME-ChIP 会聚类并过滤冗余 motif，被过滤的不产生 FIMO 输出。

### 05 Peak 注释（ChIPseeker / R）
用 `annotatePeak` 把 peak 归属到基因并分类特征。

- `PROMOTER_UP` / `PROMOTER_DOWN` 就是 ChIPseeker 的 `tssRegion=c(-UP,+DOWN)`，**控制饼图里 "Promoter" 的范围**。
- ChIPseeker 按优先级 `Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic` **覆盖判定**。启动子窗口若含 TSS 下游，5'UTR 前端会被 Promoter 吸收——这是设计行为。
- 每个 peak 只归属**一个**基因（最近 TSS）；下游一律以注释表的 `geneId` 为准。

产出 `peak_annotation_full.tsv`、`genes.bed`、`background_universe.txt`、注释饼图、TSS 分布图等。

### 06 候选靶基因整合（04 motif × 05 注释）
把 FIMO 的 motif 命中（基因组坐标）映射回 peak，与注释的位置特征取交集：

> **候选 = 位置合适（`TARGET_FEATURES`） + peak 内含 FIMO motif**

```bash
chromaflow candidates
chromaflow candidates --features "Promoter,1st Intron,Distal Intergenic" --pval 1e-4
```

`TARGET_FEATURES` 可选：`Promoter` / `5' UTR` / `3' UTR` / `1st Exon` / `Other Exon` / `1st Intron` / `Other Intron` / `Downstream` / `Distal Intergenic`。想纳入首内含子或远端增强子时按需添加。

产出：`candidate_targets.tsv`（逐 peak：坐标 + 各 motif 计数 + 判定）、每个 motif 独立的 `motif_<名>_hits.tsv`（位置/链/p 值/距 TSS）与 `_genes.txt`、`motif_summary.tsv`（整合 CentriMo）、`high_confidence_genes.txt`、`go/`。

### 07 单基因浏览器波峰图
```bash
chromaflow tracks --genes GENE1,GENE2
chromaflow tracks --gene-list mylist.txt [--outdir DIR]
```

- **自动识别基因列**：带表头且含 `gene_id`/`Gene`/`Genes`/`gene_name`/`id` 时取该列并跳过表头，因此可直接喂 `query_binding.tsv`、`candidate_targets.tsv`、`candidate_genes.txt`；无表头则取第一列。
- **peak 标签带 signalValue**（如 `IDR_peak_1 (sig=6299)`），用定量富集判断强弱。
- **固定叠加 log2(ChIP/Input) 差值轨道**（deepTools：多重复先 `bigwigAverage` 合并，再 `bigwigCompare` 求 log2 比值）。独立 y 轴，正=富集(红)/负=低于背景(蓝)，**平坦即与背景无异**。

> 覆盖轨道的 y 值（如 0–6）是归一化测序深度，与 signalValue（MACS 富集倍数，可达数千）是**两套量纲**，不要用波峰高度判断富集强弱。

绘图核心是 `lib/steps/plot_tracks.sh`（顶部 STYLE 段集中控制配色/高度/字体/DPI）。`tracks` 与 `lookup --plot` 都调用它——**调美术只改这一处，主逻辑不动**。

### 08 GO 富集
`go-prep` 导出前景/背景基因列表（供 agriGO / g:Profiler 等外部工具），`go-plot` 把返回结果画成气泡图。

### 09 先验基因反查（lookup）— 独立证据分解
拿一份你关心的基因（如缺磷响应基因），**独立**从 05 注释 + 04 FIMO 判定（**不依赖 candidate**，因此 candidate 没筛出预期结果时仍可用）。

```bash
chromaflow lookup --gene-list phosphate_genes.txt
chromaflow lookup --gene-list phosphate_genes.txt --peakset consensus --plot
chromaflow lookup --gene-list phosphate_genes.txt --min-signal 50 --level candidate
```

产出 `query_binding.tsv` 是**证据分解表**，`tier` 列直接说明每个基因**为什么（不）是候选**：

| tier | 含义 |
|---|---|
| `candidate` | 位置合适 + 含 motif（最可信） |
| `in_target` | 有 peak 在合适位置，但**无 motif 支撑** |
| `motif_only` | 有 motif，但**位置不合适**（如远端） |
| `peak_only` | 有 peak，位置与 motif 都不满足 |
| `none` | 无 peak |

同时输出三档统计与 **ID 校验**——区分"基因存在但无 peak"（真实生物学结果）与"ID 在注释中不存在"（多半是 ID 格式/物种不匹配），后者写入 `id_not_found.txt`。

---

## 可选步骤与跳过方式

`reference` / `normalize-gff` / `signal` / `spacing` 都是可选的，跳过不影响主线分析。

| 情况 | 怎么办 | 后果 |
|---|---|---|
| 物种没有合适的线粒体/叶绿体基因组 | `chromaflow reference --skip-organelle`，或 `ORGANELLE_SOURCES=()` / `SKIP_ORGANELLE=true` | 只用核基因组，其余完全一致 |
| 部分细胞器文件缺失 | 什么都不用做，自动跳过并 WARN | 想严格中断就加 `--strict` |
| 基因组/GTF/chrom.sizes 已备好 | 整步跳过 `reference`，文件放进 `REFERENCE_DIR` | 自动识别 |
| AGAT 装不上 | 照常跑 `normalize-gff` | 自动改用内置 Python 引擎 |
| 不需要规范化 GFF | 跳过 `normalize-gff` | 05 用原始 GFF3，仅在缺 exon 行时提示 |
| 不需要 GTF（只跑下游） | `chromaflow reference --skip-gtf` | 只有 02 需要 GTF |
| 没有 R | 改用 `chromaflow annotate-py` | 输出契约一致，06–09 不受影响 |
| 没有 motif 库 | `MEME_DB` 留空 | 跳过 Tomtom，de novo 发现照常 |
| 没有 annotation_info | `ANNOTATION_INFO` 留空 | GO 导出受限，其余不受影响 |

---

## 样本表与双端配对

样本表自动定位顺序：`--samples <file>` > config 的 `SAMPLES_FILE` > `<PROJECT_DIR|INPUT_DIR>/samples_<type>.{txt,tsv,csv}` > `samples.{txt,tsv,csv}`。

三种写法（TAB / 逗号 / 空格分隔均可）：

```
# ① 表头式（推荐）：列顺序随意、可缺省、可加自定义列
sample_id   antibody   control_id   replicate   fastq_path
Dof_ChIP    crDOF      Input        1           /data/raw/Dof_rep1
Input       -          -            1           /data/raw/input

# ② 位置式（旧格式，仍兼容）
Dof_ChIP    crDOF      Input        1           /data/raw/Dof_rep1

# ③ 直接给文件（跳过自动配对）
sample_id   replicate  fastq_1                  fastq_2
Dof_ChIP    1          /data/raw/A_1.fq.gz      /data/raw/A_2.fq.gz
```

列名别名（不分大小写）：`sample|name`、`ab|target`、`control|input`、`rep`、`fastq|path|dir`、`r1|read1`、`r2|read2`。

**双端配对不再依赖固定后缀表**，而是识别文件名中**最右侧的双端记号**，因此天然支持各家命名：

```
_1/_2   _R1/_R2   .1/.2   -1/-2   _R1_001/_R2_001   _good_1/_good_2
_clean_1/_clean_2   _val_1/_val_2   _1P/_2P   sampleR1/sampleR2
```

多 lane 文件会为同一样本生成多行；找不到配对的文件按单端处理并 WARN；`sample_rep1` / `sample_rep2` 这类重复编号**不会**被误判为双端。

命名特别怪时可覆盖（`{R}` 处填 1/2）：

```bash
# 写在样本表顶部
#@ pair_pattern: _good_{R}.fq.gz
#@ single_end: true

# 或写在 chromaflow.config
FASTQ_PAIR_PATTERN="_R{R}_001.fastq.gz"
FASTQ_SINGLE_END=true

# 或命令行
chromaflow samplesheet --pair-pattern '_R{R}_001.fastq.gz'
chromaflow samplesheet --single-end
chromaflow samplesheet --list-fastq        # 只预览配对结果，不写文件
```

---

## 配置与文件自动识别

config 只声明**你手里本来就有的东西**，脚本产出的文件一律不写回 config。

```bash
SPECIES_SHORT="Cre"
PROJECT_DIR="/data/projects/CrDof_ChIP"
REFERENCE_DIR="${PROJECT_DIR}/inputs/reference"

GENOME_FASTA_SOURCE=".../genome.fa"         # 压缩与否都行（.gz/.bz2/.xz/.zst 自动识别）
GENOME_GFF3_SOURCE=".../genes.gff3"         # 同上
ANNOTATION_INFO=""                          # 可选
MEME_DB=""                                  # 可选
ORGANELLE_SOURCES=("chloroplast:NC_005353.1" "mitochondria:NC_001638.1")   # 可选
```

派生文件（`*_combined.fa`、`*.fai`、`chrom.sizes`、`*.gtf`、`*.agat.gff3`）解析优先级：

```
config 显式赋值  >  00/00b 登记的 manifest  >  REFERENCE_DIR 目录扫描  >  兜底命名
```

manifest 位于 `<OUTPUT_DIR>/.chromaflow/reference.manifest`。目录扫描会优先选规范化过的注释（`*.agat.gff3` 优先于普通 `*.gff3`）和已建索引的基因组，所以跑完 `normalize-gff` 后注释步骤自动用新文件，**不用改 config**。config 里显式写的路径永远优先，旧配置继续可用。

---

## 典型分析场景（recipes）

### 场景 A：标准分析，找直接靶基因
```bash
mamba activate ChromaFlow_idr        && chromaflow idr
mamba activate ChromaFlow_motif      && chromaflow motif
mamba activate ChromaFlow_R          && chromaflow annotate
mamba activate ChromaFlow_downstream && chromaflow candidates
chromaflow tracks --gene-list results/downstream/chipseq/candidates/high_confidence_genes.txt
```

### 场景 B：重复性不足，用 consensus 妥协
IDR 只剩几十个 peak 时，改用 nf-core consensus 全集（灵敏度↑、特异性↓）：
```bash
for step in motif annotate candidates; do chromaflow $step --peakset consensus; done
chromaflow lookup --peakset consensus --gene-list my_genes.txt
```

### 场景 C：验证先验基因是否被结合
```bash
# 两套 peakset 各跑一遍，对比"我的靶点在哪个可重复性层级出现"
chromaflow lookup --gene-list my_genes.txt                      # IDR
chromaflow lookup --gene-list my_genes.txt --peakset consensus  # consensus
# 对候选级基因出图（自动带差值轨道 + sig 标签）
chromaflow tracks --gene-list results/downstream/chipseq/lookup/candidate_genes.txt \
                  --outdir figures/my_targets
```

### 场景 D：收紧到高置信核心集
```bash
chromaflow idr --min-pairs 2       # 只留 ≥2 个配对支持的 peak
# 随后用同一 peakset 重跑 motif → annotate → candidates
```

### 场景 E：检验某个重复的影响
```bash
chromaflow idr                                        # 先看日志里的 pairwise 表格
chromaflow idr --exclude-rep REP1 -o results_noREP1   # 排除后对比 peak 数
```

### 场景 F：调整关注的基因区间
```bash
# 默认只认启动子；纳入首内含子与远端增强子：
chromaflow candidates --features "Promoter,5' UTR,1st Exon,1st Intron,Distal Intergenic"
```

---

## 结果文件导读

```
results/downstream/<type>/
├── idr/
│   ├── idr_peaks_threshold_0.01.bed   # 可重复 peak（signalValue / globalIDR / n_pairs）
│   └── recommended_peaks.bed          # 自动择优的推荐集
├── motif[_peakset]/
│   └── meme-chip_output/              # meme-chip.html、fimo_out_*/、centrimo_out/
├── annotation[_peakset]/
│   ├── peak_annotation_full.tsv       # 逐 peak：基因归属、特征、signalValue
│   ├── genes.bed                      # 基因坐标（tracks / lookup 依赖）
│   ├── background_universe.txt        # GO 背景
│   └── anno_pie.png / dist_to_tss.png / tss_avgprof.png
├── candidates[_peakset]/
│   ├── candidate_targets.tsv          # 逐 peak 候选判定（含各 motif 计数）
│   ├── motif_<名>_hits.tsv / _genes.txt
│   ├── motif_summary.tsv              # 各 motif 命中数 + CentriMo 中心富集
│   ├── high_confidence_genes.txt      # 候选基因（按强度）
│   └── go/                            # agriGO / g:Profiler 输入
├── lookup[_peakset]/
│   ├── query_binding.tsv              # 证据分解（tier 列）
│   └── candidate_genes.txt / bound_genes.txt / id_not_found.txt
└── tracks[_peakset]/                  # 单基因波峰图 + tracks.ini
```

---

## 环境一览

| 环境 | 用途 | 关键工具 |
|---|---|---|
| `nfcore_env` | 主流程 | nextflow、nf-core |
| `ChromaFlow_ref` | 参考准备、GFF 规范化、GO | samtools、gffread、bedtools |
| `ChromaFlow_agat` | *(可选)* AGAT 规范化 | AGAT（`WITH_AGAT=true` 才创建） |
| `ChromaFlow_idr` | IDR | idr、numpy<1.24 |
| `ChromaFlow_motif` | motif | MEME Suite、bedtools |
| `ChromaFlow_R` | ChIPseeker 注释 | R、ChIPseeker、clusterProfiler |
| `ChromaFlow_downstream` | 整合、反查、绘图 | **deepTools**、pyGenomeTracks、bedtools |

`ChromaFlow_agat` 刻意独立且默认不装：AGAT 的 Perl 依赖树与其它包经常解析冲突。只有确实要用 AGAT 时才 `WITH_AGAT=true bash setup_environments.sh`。

`chromaflow doctor` 可自检工具与参考文件是否就位；可选功能（细胞器、AGAT、GO 注释表）一律报 `[INFO]`/`[SKIP]`，不计入失败。

---

## 质量评估与结果解读

### 关键指标怎么读

| 指标 | 来源 | 含义 |
|---|---|---|
| NSC / RSC | phantompeakqualtools（03 打印） | 富集的链交叉相关特征；NSC>1.05、RSC>1.0 为理想 |
| pairwise 通过率 | 03 的配对表格 | 各重复之间的一致性 |
| `n_pairs` | 03 输出第 8 列 | 该 peak 被几个配对支持 |
| signalValue | MACS fold enrichment | **定量**富集倍数（peak 标签会显示） |
| globalIDR | IDR | 可重复性（-log10 尺度，越大越可重复） |
| CentriMo E 值 | 04 | motif 是否在 summit **中心富集**（直接结合的位置学证据） |

### 判断重复质量的正确顺序

1. **先看 pairwise 表格**，而不是单个重复的 peak 数或 RSC。peak 多的重复不一定一致性好——深度或背景更高时会 call 出大量无法复现的边缘 peak。
2. **把 NSC 与 Input 对比**。若各 ChIP 的 NSC 与 Input 相近，说明整体富集有限，此时 RSC 是否"PASS"参考价值很低。
3. **想收紧就用 `--min-pairs`，而不是删重复**。删重复会直接损失配对；`--min-pairs 2` 在保留全部数据的前提下提高证据要求。

### 结论强度的分层

证据链上每一层回答不同问题，**不要混为一谈**：

- **MACS / IDR** → 这里有没有可重复的富集（vs Input 的统计检验）
- **ChIPseeker** → peak 落在哪个基因的什么位置
- **FIMO** → peak 序列里有没有 motif（**序列匹配质量，不等于结合强度**）
- **CentriMo** → motif 是否中心富集（直接结合的佐证）
- **候选判定** → 位置证据 + 序列证据的交集

因此一个基因可以是"候选级"却在波峰图上信号很弱——候选看的是**位置与序列**，强度要另看 signalValue 与差值轨道。要从"候选"走到"直接功能靶点"，还需 RNA-seq 差异表达或 EMSA / ChIP-qPCR 等独立验证。

---

## 适配其他物种

1. 改 config 的 `SPECIES_NAME`/`SPECIES_SHORT`、源文件路径、`MACS_GSIZE`、`ORGANELLE_SOURCES`。没有合适细胞器基因组的物种直接 `ORGANELLE_SOURCES=()` 或加 `--skip-organelle`。
2. 跑一次 `chromaflow normalize-gff`，确认注释饼图里 UTR 比例正常。
3. `MEME_DB` 换成对应类群的 motif 库（如 JASPAR plants / vertebrates）；留空则跳过 Tomtom。
4. `KEGG_ORGANISM` 填 KEGG 物种码；若基因 ID 体系与 KEGG 不一致会自动跳过。
5. **基因 ID 格式**：lookup 只剥离尾部版本号（`.1`、`.v6.1`），兼容 `AT1G01010.1`、`ENSG00000141510.17`、`LOC_Os01g01010.1`、`Cre02.g075050_4532` 等；归一化后若仍对应多个基因会判为歧义而不猜测。跑 lookup 时留意 **ID 校验**那行，可及时发现列表格式不匹配。

---

## 依赖版本

nf-core/chipseq 2.1.0 · nf-core/atacseq 2.1.2 · MACS3 · IDR · MEME Suite 5.x · ChIPseeker 1.42 · deepTools 3.5 · pyGenomeTracks · AGAT · bedtools 2.31
