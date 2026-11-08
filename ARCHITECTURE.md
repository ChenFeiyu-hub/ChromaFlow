# ChromaFlow 1.2.9 架构

## 软件、项目配置与分析数据分开

| 层 | 内容 | 位置 |
|---|---|---|
| 软件 | CLI、步骤代码、模板、文档、测试 | 本仓库，可安装为共享只读目录 |
| 项目配置 | 物种输入、资源上限、分析参数 | 每个项目自己的 `chromaflow.config` |
| 分析数据 | FASTQ、参考文件、结果、work、缓存、日志 | 项目 INPUT_DIR / OUTPUT_DIR / REFERENCE_DIR |
| 离线软件包 | 固定版本 nf-core、Nextflow、插件和容器 | 独立 bundle，可共享在 `/opt/chromaflow/bundles/` |

配置发现顺序：`-c FILE` > `$CHROMAFLOW_CONFIG` > 当前目录 `chromaflow.config`。
物种前缀用于文件命名；不会据此截断、转换或猜测基因 ID。
安装、离线部署和分析命令见 [README](README.zh-CN.md)。

## 代码入口

```text
bin/chromaflow                      CLI、命令路由、统一下游日志
lib/common.sh                      配置读取、路径派生、输入发现、peakset 解析
lib/help.sh                        英文命令帮助与软件版本
lib/nfcore_runtime.sh               容器/Conda 入口共用的 nf-core 运行配置
lib/batch.py                        批量计划、执行、环境选择、结果隔离和断点状态
lib/offline_bundle.py               下载调度、版本锁定、文件校验与打包
lib/selection.sh                    candidates / lookup 共享参数解析
lib/steps/                         reference、samplesheet、run、各下游步骤
lib/steps/support/
  annotation_policy.py             特征类别、位置规则、注释参数一致性检查
  promoter_audit.py                 全窗口 peak 重叠核验，不改变基因归属
  motif_selection.py               模型/位点筛选、CentriMo 解释、selection 摘要
  candidates.py / lookup.py         候选整合、精确 ID 查询和最佳证据选择
  fimo_tracks.py / track_layout.py  已有 FIMO 坐标转换、区域轨道高度
 templates/                        init 使用的项目配置及三种分析样本表
 env/ + setup_environments.sh       按用途构建 runtime 和独立下游环境
 release/                          运行代码兼容关系及发布验收记录
 tools/acceptance.py                统一服务器验收及反馈收集入口
 docs/                             当前离线部署、下游与筛选操作指南
 tests/                            回归及可用时运行的真实工具集成检查
```

内部步骤由 CLI 调用。直接调用内部脚本需要显式环境变量，且会绕过 CLI 的完整日志记录，
不作为普通用户的推荐操作方式。

## 配置与派生值

`load_config` 分两遍读取配置，中间派生项目路径，允许配置引用 `${REFERENCE_DIR}`。
参考路径的优先级是：config 显式值 > `.chromaflow/reference.manifest` > 目录识别 > 默认命名。
旧配置写死参考路径时，新产物不会自动覆盖这些显式值。

`samplesheet` 可测量读长，`reference` 可测量 MACS 有效基因组大小与 GC。
派生值登记在 `<OUTPUT_DIR>/.chromaflow/derived.env`，不反复写回用户 config。
基因组、GFF/GTF、读长及其他分析输入属于项目，不进入软件仓库。

`normalize-gff` 的 `auto` 固定使用内置保守引擎，并非按 AGAT 是否安装而切换。
原始记录身份与 Parent 图保持、GTF 转换验收、TxDb 显式 ID/Name 映射是分开的检查；
非法关系或无法安全导出的模型明确失败。规则见 [GFF_NORMALIZATION.md](GFF_NORMALIZATION.md)。

## nf-core 与环境

`run` 根据 profile 选择容器或 Conda 入口。完整离线运行需要本地容器 bundle；
`offline prepare/check/pack` 负责制包、校验及迁移，`code-only` 不等于完整可运行的离线包。
Java、Python、Apptainer 以及 ChromaFlow 自有下游环境需另行提供。
安装器的 offline / prepare / conda 模式见 [OFFLINE.md](docs/OFFLINE.md)。

## 下游分支与统一筛选

`DEFAULT_PEAKSET` 指定项目默认分支；显式 `--peakset` 优先，模板默认 idr。
consensus、不同 IDR 阈值使用各自的目录后缀。请求缺失或失败的 IDR 集不会自动退回 consensus。
`-o` 只改当前命令输出，不自动改变后续输入发现；`tracks --selection DIR` 则是显式衔接接口。

ChIPseeker 的分类与最近 TSS 归属保持原始语义。位置规则由 candidates 和 lookup 共用，
启动子核验独立列出整个配置窗口内的重叠 peak。lookup 仅匹配完整 ID，先展示示例并确认。

共享 motif 模块先过滤模型及所有位点，然后将位点映射到 peak。lookup 在过滤后按证据级别、
信号选择每个基因的最佳 peak；完整选择证据另行导出。CentriMo 局部/中央模式、模型发现 E 值、
FIMO 位点 p/q 值各自记录，不混作同一统计量。

selection 保存参数、输入/输出摘要和模型身份；tracks 校验来源后使用保存的基因和 motif，
仍将绘图区间内的全部分支 peak、邻近基因结构传给绘图工具。详见 [SELECTION.md](docs/SELECTION.md)。

## 输出与日志

```text
project/
  chromaflow.config
  samples_chipseq.txt
  inputs/
    raw/
    reference/
  outputs/
    results/chipseq/
    results/downstream/chipseq/{idr,annotation,motif,candidates,lookup,tracks,...}
    logs/downstream/chipseq/
    samplesheets/
    work/
    .chromaflow/
```

设置 `OUTPUT_DIR="${PROJECT_DIR}"` 可把 outputs 下这些目录直接放在项目根。
下游目录按需创建。CLI 每次运行记录英文终端 stdout/stderr、命令、版本、config 路径和退出码，
使用独立日志文件；尚未识别项目时的参数错误及进程中断存在文档所述的日志边界。

## 批量 DAP-seq（1.2.9）

`chromaflow batch` 有独立的子命令解析器。用户只维护共享 config 和 samples.tsv；
`lib/batch.py` 生成快照，明确调用环境与绝对 CLI 路径，再复用既有参考和下游命令。
按 layout 与读长分组进入 nf-core/chipseq 2.1.0；不同 analysis_id 使用不同 antibody
分组键，控制样本显式引用。下游从精确输出 manifest 取得本分析的 peaks。
`CF_BATCH_*` 为内部生成配置字段，普通单项目不设置时保留原有目录布局。

批量用 SHA256 追踪参考/输出、用输入和参数签名控制复用；断点状态与结果摘要分别记录。
失败或无法评估的分析不填零候选。详情与限制见 [BATCH.md](docs/BATCH.md)。
