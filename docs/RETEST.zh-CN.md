# 1.2.7：仅使用 chromaflow 命令的逐步验收

本指南使用通用示例路径；请替换为自己的项目。安装与离线运行步骤见
[README](../README.zh-CN.md) 和 [离线指南](OFFLINE.md)。
已有参考、成功的 nf-core 结果和环境可以复用，无需重做分析来更新软件文档。

## 核对入口与配置

```bash
chromaflow version
command -v chromaflow
cd "$HOME/projects/MyChIP"
chromaflow doctor --type chipseq
```

应确认当前版本为 1.2.7，且项目目录中只有一份实际使用的 `chromaflow.config`。
无需删除参考目录。显式配置路径优先于自动发现；保留 `.chromaflow` 的登记与实测结果。
启动子窗口由 `PROMOTER_UP/DOWN` 控制，候选类别由 `TARGET_FEATURES` 控制，
距离补救由 `TARGET_MAX_TSS_DIST` 控制。先确定这些参数，再做下游分析。

## 逐步执行 consensus 分支

以下以已有 nf-core 输出和有效参考为前提。无需先执行 IDR。

```bash
conda activate ChromaFlow_R
chromaflow annotate --type chipseq --peakset consensus
conda activate ChromaFlow_motif
chromaflow motif --type chipseq --peakset consensus
conda activate ChromaFlow_downstream
chromaflow candidates --type chipseq --peakset consensus
chromaflow lookup --type chipseq --peakset consensus \
  --gene-list genes_of_interest.txt --level candidate \
  --outdir results/downstream/chipseq/lookup_selected
chromaflow tracks --type chipseq --peakset consensus \
  --selection results/downstream/chipseq/lookup_selected \
  --outdir results/downstream/chipseq/tracks_selected
```

按步骤核对日志后再继续。若已有与当前参数相同且通过校验的 annotation/motif 结果，
可从 candidates 开始。`-o` 只改变当前输出；上游 annotate/motif 使用默认目录，
使后续命令能自动发现它们。lookup 使用自定义目录后由 `tracks --selection` 显式衔接。

## 每步核对哪些结果

| 命令 | 核对项 |
|---|---|
| annotate | `peak_annotation_full.tsv`、`annotation_policy.json`、全启动子重叠报告 |
| motif | MEME-ChIP 完成状态、实际扫描窗口、FIMO 原始结果 |
| candidates | 候选列表、motif 目录清单、CentriMo 模式/阈值、selection 参数 |
| lookup | 两边精确 ID 示例、y/N 确认、`id_check.tsv`、tier、筛选后的 selected 列表 |
| tracks | `plot_status.tsv`、完整区域 peak/邻近基因、保存的 motif 过滤结果 |

全部下游命令自动记录到 `<OUTPUT_DIR>/logs/downstream/chipseq/`。
检查退出码、可选步骤警告和缺失序列提示；退出码 0 不代表每一项可选图形都已生成。
lookup 的未匹配 ID 保留为未匹配，不自动删除后缀或猜测别名。

## 可选对照与完整列表图

需要比较 IDR 时，先运行 `chromaflow idr --type chipseq --min-pairs 1`，
然后使用 `--peakset idr` 重复相应下游步骤。不要交叉使用不同分支的注释和 motif。

完整列表查询/绘图使用 `lookup --level all --plot`；仅候选使用
`lookup --level candidate --plot`。直接 `tracks --gene-list FILE --motifs` 使用
FIMO_PVAL；需要指定模型或 q/E 值过滤时，先 candidates/lookup，再用 selection 绘图。
[统一筛选指南](SELECTION.md) 列出所有参数、文件摘要校验和未扫描区域的解释边界。

[启动子核验](DOWNSTREAM.md) 报告窗口内分配给其他基因的 peak，但不改写 ChIPseeker 的原始归属。
GFF/TxDb 验收另见 [GFF_NORMALIZATION.md](../GFF_NORMALIZATION.md)。
