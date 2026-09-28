# ChromaFlow 1.2.7 使用指南

[English README](README.md) · [离线部署完整命令](docs/OFFLINE.md) · [下游分析完整命令](docs/DOWNSTREAM.md) · [统一筛选与绘图](docs/SELECTION.md) · [版本说明](CHANGELOG.md)

ChromaFlow 支持 ChIP-seq、DAP-seq 与 ATAC-seq。原始测序流程分别使用 nf-core/chipseq 2.1.0 和 nf-core/atacseq 2.1.2；软件自有下游功能包括 IDR、motif、ChIPseeker 注释、候选基因整合、信号图、基因轨道与 GO 输入导出。

## 1.2.7 的主要变化

- candidates 与 lookup 共用 `--motif`、`--pval`、`--qval`、`--motif-evalue`、`--centrimo-evalue`，先筛选所有 peak 的证据，再选每个基因的最佳 peak。
- 自动生成 `selected_genes.txt`、`selected_evidence.tsv` 和 `selection.json`；`tracks --selection DIR` 直接使用同一筛选结果，校验输入版本和 peak 分支。
- `lookup --level candidate --plot` 仅绘制候选；默认 `--level all --plot` 绘制精确匹配的查询基因。未匹配 ID 仍保留在报告中。
- CentriMo 区分局部/中央模式与阈值状态；E 值有记录不再自动标成中央富集。新增过滤默认关闭。
- tracks 按转录本与 motif 数量自动调整轨道高度，保留邻近结构；peak 强度显示三位有效数字。
- 项目可设置 `DEFAULT_PEAKSET=consensus`，仍允许用 `--peakset` 显式覆盖。
- [统一筛选与绘图指南](docs/SELECTION.md)。已有参考、nf-core、注释和 motif 结果可直接复用；无需新建 conda 环境。

1.2.6 的精确 ID 匹配、y/N 确认、ChIPseeker 位置规则与完整启动子核验继续保留。

## 1.2.5 修复

- 修复 IDR 合并后 bedtools 输出顺序变化导致的 `map` 排序错误；保留 IDR 拟合和阈值。
- IDR 失败时禁止 GO 导出复用旧注释；修正把 R 缺失值 `NA` 计为基因归属冲突的提示。
- 新增 `tracks --motifs --flank 5000`，按输入列表逐基因展示信号、整套 peak、邻近基因结构和已有 FIMO 命中；导出坐标表及逐基因绘图状态。
- motif 轨道只展示已经扫描的输入窗口；不会自动扫描整个启动子。当前 `lookup --plot` 遵守 `--level`，默认 all 绘制精确匹配的全部查询基因。
- 复用已有参考文件、nf-core 结果及环境；无需重跑成功的 nf-core 流程。

## 1.2.4 引入的功能

- 各下游输出目录按需建立，不会仅因执行 annotation 就产生空的 `signal/`。
- TSS 图改成上下游距离区间条形图，标出百分比与峰数，并导出 `dist_to_tss.tsv`。横轴是峰的百分比，纵轴是距离区间。
- 下游终端输出自动写入 `<OUTPUT_DIR>/logs/downstream/<type>/`。每次执行一个独立日志，包含命令、版本、配置路径和退出码；失败信息也记录。
- 程序消息、帮助、图中文字改为英文；不修改用户路径、基因 ID 等原始信息。
- `signal`、`go-prep`、`annotate-py` 正确区分 IDR 与 consensus 输出。显式请求的 IDR 结果不存在时会报错。
- 安装脚本区分离线服务器运行、联网制包、在线 Conda 运行，并可只安装运行时或下游环境。
- IDR 错误不再被隐藏；补回真实 MACS `signalValue`，避免把 p-value 排名值当作峰强度。详情见版本说明。

## 软件升级

Linux/WSL 上可直接从本仓库安装：

```bash
mkdir -p "$HOME/biosoft"
git clone https://github.com/ChenFeiyu-hub/ChromaFlow.git "$HOME/biosoft/ChromaFlow"
export CHROMAFLOW_HOME="$HOME/biosoft/ChromaFlow"
export PATH="$CHROMAFLOW_HOME/bin:$PATH"
hash -r
chromaflow version
git -C "$CHROMAFLOW_HOME" rev-parse HEAD
```

记录最后输出的提交号便于复现。内网服务器可接收这份完整源码目录，或使用已有的软件发行包。
等待当前任务结束后再切换软件；升级到独立目录，避免覆盖正在执行的代码。
如果已拿到版本化的压缩包和校验文件，可用以下方式安装：

```bash
cd ~/nfcore-upload
sha256sum -c ChromaFlow-1.2.7.tar.gz.sha256
mkdir -p "$HOME/biosoft/ChromaFlow-1.2.7"
tar -xzf ChromaFlow-1.2.7.tar.gz \
  -C "$HOME/biosoft/ChromaFlow-1.2.7" --strip-components=1
export CHROMAFLOW_HOME="$HOME/biosoft/ChromaFlow-1.2.7"
export PATH="$CHROMAFLOW_HOME/bin:$PATH"
chromaflow version
```

确认后可把两行 `export` 写入 `~/.bashrc`。继续使用已经验证的项目 config，不要用新模板覆盖它。**本次升级不要求重新下载两个 nf-core bundle，也不要求重跑已成功的 nf-core、GFF3 规范化或参考构建。** 已有下游结果不因换版本自动更新；要应用评分/分析修正，应重跑对应下游步骤。仅更换 TSS 图可使用下面的独立重绘命令。

## 安装环境：选择适用模式

```bash
# 查看安装计划，不执行下载或修改
chromaflow setup --mode offline --components runtime --dry-run

# 离线容器运行端：只装 Python + Java；需要已有 Conda 包源可用
chromaflow setup --mode offline --components runtime

# 联网 Linux/WSL 制包端：装 Python + Java + nf-core/tools
chromaflow setup --mode prepare --components runtime

# 只有使用在线 Conda profile 时才选这一项
chromaflow setup --mode conda --components runtime

# 自有下游环境独立安装，可选 --skip-r、--with-agat
chromaflow setup --components downstream
```

三种 runtime 方案按机器用途选择，不需要依次全装。默认 `--mode offline --components all`；只运行 nf-core 容器时显式加 `--components runtime`。已有同名环境会保留，不会自动升级或删除。

**`offline` 指 nf-core 的执行方式，不表示 Conda 安装器能够凭空离线装包。** 内网服务器已具备 Java、Python、Apptainer 以及之前测试过的下游环境时，可以直接复用，无需运行安装器。容器包包含 nf-core 自身的分析工具；ChIPseeker、MEME、IDR、deepTools 等 ChromaFlow 下游工具仍使用独立环境。Apptainer/Singularity 需要系统管理员预先安装。

| 环境 | 用途 |
|---|---|
| `ChromaFlow_runtime` | 本地 bundle 运行时：Python、Java |
| `nfcore_env` | 联机制包或在线 Conda nf-core，依安装 mode 而定 |
| `ChromaFlow_ref` | 参考处理、样本表、Python 注释、GO 辅助、spacing |
| `ChromaFlow_idr` | IDR 与 bedtools |
| `ChromaFlow_motif` | MEME Suite、bedtools、samtools |
| `ChromaFlow_R` | ChIPseeker、TxDb、clusterProfiler |
| `ChromaFlow_downstream` | deepTools、pyGenomeTracks、候选整合与基因查询 |

## 两个离线包：由 ChromaFlow 统一管理

在联网的 WSL/Linux 上运行：

```bash
conda activate nfcore_env
chromaflow offline prepare --type chipseq --bundle "$HOME/nfcore-offline/chipseq-2.1.0"
chromaflow offline prepare --type atacseq --bundle "$HOME/nfcore-offline/atacseq-2.1.2"

for BUNDLE in chipseq-2.1.0 atacseq-2.1.2; do
  chromaflow offline check --bundle "$HOME/nfcore-offline/$BUNDLE" --verify --runtime-check
  chromaflow offline pack --bundle "$HOME/nfcore-offline/$BUNDLE" --archive "$HOME/$BUNDLE.tar"
done
```

不用另行手动执行 nf-core 下载、Nextflow 初始化、插件安装。将 tar 和对应 `.sha256` 上传到服务器，校验后分别解压到 `/opt/chromaflow/bundles/`。多账号共享只读软件包，每个账号的项目、work、缓存与日志各自存放。完整的上传后命令、权限、版本管理和故障处理见 [OFFLINE.md](docs/OFFLINE.md)。

WSL 中 Windows 的 `D:\biosoft` 对应 `/mnt/d/biosoft`；推荐在 Linux 的 `$HOME/nfcore-offline` 生成 bundle。制包端仍须能访问 GitHub 和容器仓库。服务器无网络不影响上传完成且经过验证的本地 bundle 运行。

## 新项目和前置流程

```bash
chromaflow init "$HOME/projects/MyChIP" --type chipseq
cd "$HOME/projects/MyChIP"
# 先编辑 chromaflow.config 与 samples_chipseq.txt
conda activate ChromaFlow_ref
chromaflow samplesheet --type chipseq --strict
chromaflow reference
chromaflow normalize-gff --engine builtin
chromaflow reference
chromaflow doctor --type chipseq
conda activate ChromaFlow_runtime
chromaflow run --type chipseq --offline \
  --bundle /opt/chromaflow/bundles/chipseq-2.1.0 --dry-run
chromaflow run --type chipseq --offline \
  --bundle /opt/chromaflow/bundles/chipseq-2.1.0 --yes
```

已有通过验证的参考文件时可以跳过参考构建。`normalize-gff` 之后第二次 `reference` 用于根据新注释更新 GTF。明确配置的 `GENOME_GFF3` 等路径优先级最高；旧项目若写死旧路径，需要有意识地修改，而不能仅依赖自动识别。结构校验与 GTF/TxDb 验收是不同检查，参见 [GFF_NORMALIZATION.md](GFF_NORMALIZATION.md)。

模板仅需填写源参考文件等已有输入；`READ_LENGTH` 和 `MACS_GSIZE` 默认不填，由前置步骤测量并记录在 `.chromaflow/derived.env`。先执行 samplesheet，使 k-mer 估计有实际读长可用。无细胞器时 `ORGANELLE_SOURCES=()`；完全离线时细胞器只能引用本地 FASTA。`KEGG_ORGANISM=""` 可关闭需要联网的 KEGG 分析。

## 分别运行 IDR 与 consensus

以下命令在项目目录执行；也可每条加 `-c /绝对路径/项目.config`。

```bash
conda activate ChromaFlow_idr
chromaflow idr --type chipseq --min-pairs 1

for PEAKSET in idr consensus; do
  conda activate ChromaFlow_R
  chromaflow annotate --type chipseq --peakset "$PEAKSET"
  conda activate ChromaFlow_motif
  chromaflow motif --type chipseq --peakset "$PEAKSET"
  conda activate ChromaFlow_downstream
  chromaflow candidates --type chipseq --peakset "$PEAKSET"
  chromaflow signal --type chipseq --peakset "$PEAKSET"
  conda activate ChromaFlow_ref
  chromaflow go-prep --type chipseq --peakset "$PEAKSET" --set promoter
done
```

IDR 默认集为 0.01，也可以明确指定 `idr:0.05`。IDR 需要同一生物条件下至少两个重复；多个条件不会混合配对。`--min-pairs 1` 表示至少一个配对支持，并不代表所有重复一致。三重复有三个配对，可视分析目的改成 2 或 3。IDR 脚本每次会重写自己的阈值文件，比较不同 min-pairs 方案时应使用独立输出根目录。

两个分支分别输出 `annotation/` 和 `annotation_consensus/`、`motif/` 和 `motif_consensus/`、`signal/` 和 `signal_consensus/` 等。IDR 是可选步骤：只分析 consensus 时可以跳过 `idr`，将各下游命令设为 `--peakset consensus`。也可在 config 中设置 `DEFAULT_PEAKSET=consensus`，省略反复指定；显式 `--peakset` 的优先级更高。模板默认仍为 idr，请求的 IDR 文件不存在或为空时会报错。更多命令见 [DOWNSTREAM.md](docs/DOWNSTREAM.md)。

## 从候选筛选到基因查询、统一绘图

项目配置分别控制注释窗口和候选类别，例如：

```bash
# 示例窗口；应按分析目的自行选择。
DEFAULT_PEAKSET=consensus
PROMOTER_UP=3000
PROMOTER_DOWN=0
TARGET_FEATURES="Promoter,5' UTR,1st Exon"
TARGET_MAX_TSS_DIST=0
FIMO_PVAL=1e-4
```

ChIPseeker 的类别优先级为 Promoter > 5′ UTR > 3′ UTR > Exon > Intron > Downstream > Intergenic。
启动子窗口是相对链方向的闭区间 `[-UP,+DOWN]`。`TARGET_FEATURES` 在注释后筛选类别，
不改变最近 TSS 基因的归属；距离补救规则也只作用于已分配的基因，0 表示关闭。
改变窗口需要重新 annotate。内置启动子核验另外报告整个窗口内的全部重叠 peak，
包括分配给其他基因的 peak，输出 `promoter_peak_overlaps.tsv`、
`promoter_peak_summary.tsv` 和参数记录；不会据此自动改写归属或新增候选。

同一分支的 annotate、motif 完成后：

```bash
conda activate ChromaFlow_downstream
chromaflow candidates --type chipseq --peakset consensus
chromaflow lookup --type chipseq --peakset consensus \
  --gene-list genes_of_interest.txt --level candidate \
  --outdir results/downstream/chipseq/lookup_selected
chromaflow tracks --type chipseq --peakset consensus \
  --selection results/downstream/chipseq/lookup_selected \
  --outdir results/downstream/chipseq/tracks_selected

# 可选：指定本次分析的某个 motif，一条命令完成查询与绘图。
chromaflow lookup --type chipseq --peakset consensus \
  --gene-list genes_of_interest.txt --level candidate --motif MEME-4 --plot \
  --outdir results/downstream/chipseq/lookup_MEME4
```

lookup 先摘录两边 ID 并请求 y/N 确认，严格按完整 ID 匹配，不删除版本或其他后缀。
非交互执行时 `--yes` 表示明确接受精确匹配。`--level all` 为默认值，保留全部查询行并绘制
精确匹配的基因；`in_target` 选择有合适位置证据的基因，`candidate` 还要求通过筛选的 motif。
`--plot` 遵守这个选择，未匹配 ID 始终留在 `id_check.tsv`。

| candidates / lookup 参数 | 对应 config | 含义 |
|---|---|---|
| `--motif MEME-4,MEME-7` | `MOTIF_IDS` | 本次运行的精确 motif 标签或原始 ID |
| `--pval 1e-5` | `FIMO_PVAL` | 位点 FIMO p 值上限；默认 1e-4 |
| `--qval 0.05` | `FIMO_QVAL` | 位点 FIMO q 值上限 |
| `--motif-evalue 0.05` | `MOTIF_EVALUE_MAX` | 模型发现 E 值上限 |
| `--centrimo-evalue 0.05` | `CENTRIMO_EVALUE_MAX` | 位置富集 E 值上限 |
| `--min-signal 2` | 命令默认 0 | peak 信号分数下限 |

新增可选 config 过滤默认留空，表示关闭。所有 peak 先经过同一筛选，lookup 再选择每个基因的
最佳证据；不会仅因另一个 peak 的信号更强而漏掉有指定 motif 的 peak。
MEME-4 不一定对应 `fimo_out_4`，不同分支的相同 MEME 编号也不保证是相同模型。
指定未扫描模型、请求缺失 q 值或 E 值时明确报错。放宽阈值不能找回原 FIMO 没导出的位点，
这些命令不会自动重扫描 DNA。

自动生成的 `selected_genes.txt`、`selected_evidence.tsv`、`motif_catalog.tsv` 和
`selection.json` 保存基因、完整证据、模型身份、参数和文件摘要。
`tracks --selection` 检查分支、启动子参数及输入摘要，再使用保存的选择绘图；
区域内全部 peak 和邻近基因结构仍保留。`candidates --plot` 也可直接画候选。
完整 `candidate_targets.tsv` 包含不通过的 peak，候选图应使用 selected 列表或 selection。
历史文件名 `high_confidence_genes.txt` 是计算候选列表的兼容别名，不代表已完成实验验证。

CentriMo 区分 local/central/unknown 和“通过阈值、超过阈值、未报告”，显示标签使用 E≤0.05；
只有显式设置 `--centrimo-evalue` 才影响候选筛选。局部富集不自动解释成中央富集。
密集基因/motif 轨道按内容自动增高，可用 `TRACK_GENE_HEIGHT`、`TRACK_MOTIF_HEIGHT`
设置为 `auto` 或正数厘米；peak 强度显示三位有效数字。
完整语义、输出文件和限制见 [统一筛选与绘图指南](docs/SELECTION.md)。

## 日志与 TSS 图

统一日志位置为：

```text
<OUTPUT_DIR>/logs/downstream/chipseq/<时间>_<命令>_<peakset>_<随机标识>.log
```

无需再手写 `tee`；stdout 和 stderr 都记录。配置错误发生在项目路径尚未确定之前时只输出到终端。中断运行可能只留下未结束日志；已有历史终端内容不能追溯补录。`-o` 只改变该步骤的结果目录，日志仍在项目统一目录。

新 TSS 图按注释链方向区分上下游，分母为有有效 TSS 距离的峰数。`0 bp` 单列，其他区间不再包括 0；TSV 写出精确的绝对距离上下界及缺失数。该图适合检查峰相对 TSS 的分布，不应直接解释成相对随机背景的富集显著性。

只重绘既有 ChIPseeker 结果、不重复注释：

```bash
conda activate ChromaFlow_R
Rscript "$CHROMAFLOW_HOME/lib/steps/plot_tss.R" \
  results/downstream/chipseq/annotation/peak_annotation_full.tsv \
  results/downstream/chipseq/annotation/tss_plot_v1.2.7
```

重绘命令写入新子目录，可与旧图比较。它使用 R 版表中的 `distanceToTSS` 列。当前保留经验证的 1.2.3 GFF 规范化实现和 ID/Name 映射验收逻辑，不增加衣藻专用 ID 截断规则。

## 验证与更多文档

在软件仓库根目录运行：

```bash
python3 -m unittest discover -s tests -v
```

1.2.7 本地验收通过 113 项测试，覆盖离线包、命令日志、分支隔离、GFF/TxDb、真实 bedtools
后处理、统一 motif 筛选、CentriMo 标签及 selection 来源检查；缺少 bedtools 或 Linux/Bash 时，
相关测试会标为跳过。另用真实 pyGenomeTracks 验证了正负链、邻近 peak 和密集转录本绘图。
这些检查不能代替用户自己的数据验收。程序终端输出、帮助及图中文字保持英文。

- [离线部署与共享 /opt](docs/OFFLINE.md)：下载、打包、上传、权限与版本管理。
- [下游流程说明](docs/DOWNSTREAM.md)：IDR、注释、启动子核验、signal、GO、spacing。
- [GFF 规范化与 TxDb 验收](GFF_NORMALIZATION.md)：通用记录保持与显式身份映射规则。
- [软件架构](ARCHITECTURE.md)、[版本记录](CHANGELOG.md)、`chromaflow help <command>`。
