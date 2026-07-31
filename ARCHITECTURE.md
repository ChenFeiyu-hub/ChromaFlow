# ChromaFlow 架构

## 核心：解耦「软件 / 配置 / 运行」三层

| 层 | 内容 | 位置 | 谁修改 | 进 git 吗 |
|---|---|---|---|---|
| ① 软件本体 | CLI、lib、steps、templates | 这个仓库（clone 一次，加进 PATH） | 开发者 | **是**（就是本仓库） |
| ② 项目配置 | 物种参数、输入/输出路径 | 每个数据集一份 `chromaflow.config` | 用户 | 否（在你的数据目录里） |
| ③ 运行 I/O | fastq、参考、结果、日志、work | 由配置里的 `INPUT_DIR`/`OUTPUT_DIR` 决定 | 自动 | 否 |

旧版把三层塞进同一个 `PROJECT_DIR`，于是 `results/ logs/ work/` 长在代码仓库里，换数据集要改源码。
v2 让**代码目录永远只读、零生成物**：所有产物落在项目的 `OUTPUT_DIR` 下。

## GitHub main 目录（= 软件本体，只放代码）

```
ChromaFlow/                          # 仓库根 = 软件，不含任何具体数据/结果
├── README.md                        # 英文主文档
├── README.zh-CN.md                  # 中文文档
├── ARCHITECTURE.md
├── .gitignore                       # 兜底忽略误建的项目产物
├── setup_environments.sh            # conda/mamba 环境安装
├── bin/
│   └── chromaflow                   # 唯一命令行入口（自动探测 CHROMAFLOW_HOME）
├── lib/
│   ├── common.sh                    # 配置解析 + 路径派生/自动识别 + 日志（架构核心）
│   ├── help.sh                      # chromaflow help [command] 的帮助文本
│   └── steps/
│       ├── init_project.sh          # chromaflow init 的实现
│       ├── 00_prepare_reference.sh
│       ├── 00b_normalize_gff.sh
│       ├── normalize_gff.py         # 无 AGAT 时的内置 GFF 规范化器
│       ├── 01_make_samplesheet.sh
│       ├── make_samplesheet.py
│       ├── 02_run_nfcore.sh
│       ├── 02_run_nfcore_conda.sh
│       ├── 03_idr_analysis.sh
│       ├── 04_meme_motif.sh
│       ├── 05_peak_annotation.sh
│       ├── 06_signal_deeptools.sh
│       └── support/
│           └── motif_spacing_analysis.sh
├── templates/                       # 项目模板：init 时拷贝到用户数据目录
│   ├── chromaflow.config
│   ├── samples_chipseq.txt
│   ├── samples_atacseq.txt
│   └── samples_dapseq.txt
└── conf/
    └── nextflow.config              # Nextflow 通用配置（不硬编码 profile）
```

关键点：仓库里**没有** `reference/`、`results/`、`logs/`、`samplesheets/`、`00_data/`。
这些属于「某个数据集」，由 `chromaflow init` 在你选择的位置生成。

## 一个「项目」的样子（= ②③，放在你的数据盘，不进这个仓库）

```
/data/CrDof_DAPseq/                  # chromaflow init 生成
├── chromaflow.config                # 这个数据集的配置（物种 + 路径）
├── samples_dapseq.txt               # 样本与 fastq 路径
├── inputs/
│   ├── raw/                         # fastq（也可在 config 指向别处）
│   └── reference/                   # 参考基因组源文件
└── outputs/                         # 全部结果/日志/work/缓存（自动生成）
    ├── results/{dapseq,downstream/dapseq/{idr,motif,annotation,signal}}
    ├── logs/
    ├── samplesheets/
    ├── work/                        # nextflow work
    ├── .chromaflow/                 # reference.manifest（参考产物登记表）
    └── .singularity_cache/
```

## 使用

```bash
# 安装（一次）
git clone https://github.com/<you>/ChromaFlow.git ~/ChromaFlow
export PATH="$HOME/ChromaFlow/bin:$PATH"      # 建议写进 ~/.bashrc

# 新建项目（放哪都行）
chromaflow init /data/CrDof_DAPseq --type dapseq
cd /data/CrDof_DAPseq
# 编辑 chromaflow.config（物种 + 路径）与 samples_dapseq.txt

# 逐步运行（项目目录内自动发现 ./chromaflow.config）
chromaflow reference
chromaflow samplesheet --type dapseq
chromaflow run         --type dapseq
# 下游逐步（各步 conda 环境不同，不提供一键 all）：
mamba activate ChromaFlow_idr        && chromaflow idr      --type dapseq
mamba activate ChromaFlow_motif      && chromaflow motif    --type dapseq
mamba activate ChromaFlow_ref        && chromaflow annotate --type dapseq
mamba activate ChromaFlow_downstream && chromaflow signal   --type dapseq
```

配置解析优先级：`-c <config>` > `$CHROMAFLOW_CONFIG` > `./chromaflow.config`

## 关键设计点

- **两遍 source 配置**：让 `chromaflow.config` 里能写 `ANNOTATION_INFO="${REFERENCE_DIR}/..."` 这类引用派生路径的行（第一遍取锚点 → `derive_paths` → 第二遍正确展开）。
- **步骤自包含**：每个 step 自己 `source common.sh + load_config`，既能被 CLI 以子进程调用，也能 `CHROMAFLOW_HOME=... CHROMAFLOW_CONFIG=... bash lib/steps/06_...sh dapseq` 单独跑，便于调试。
- **nextflow 产物归位**：`02` 里 `cd "${OUTPUT_DIR}"` 且传 `-w "${WORK_DIR}"`，使 `work/`、`.nextflow/`、日志全部落在 `OUTPUT_DIR`，代码目录不被污染。
- **路径全部派生**：正文只用 `${RESULTS_DIR}` `${GENOME_GTF}` 等变量，语义与旧版一致，只是根从代码目录换成了 `OUTPUT_DIR`/`REFERENCE_DIR`。
- **产物登记而非回填 config**：`00`/`00b` 生成的派生文件（combined.fa、chrom.sizes、gtf、*.agat.gff3）写入 `<OUTPUT_DIR>/.chromaflow/reference.manifest`，由 `load_reference_manifest` 读回；再叠加 `autodetect_reference_files` 的目录扫描兜底。解析优先级：**config 显式 > manifest > 目录扫描 > 兜底命名**。用户因此不需要把「脚本的输出」抄进「脚本的输入配置」。
- **可选步骤可跳过**：`reference`（细胞器、GTF）、`normalize-gff`（AGAT）都是可选路径。缺件默认 WARN 并继续，`--strict` 才中断；AGAT 缺失时 `normalize-gff` 自动切换到内置 Python 引擎，因此 AGAT 的环境冲突不再阻塞主线。
- **样本表宽松解析**：`make_samplesheet.py` 支持表头式/位置式/显式 fastq 三种记录方式，双端配对靠「文件名最右侧的双端记号」推断而非固定后缀表，可用 `pair_pattern` 模板覆盖。
