# GFF3 规范化与验收说明（1.2.3 修复，适用于 1.2.7）

2026-09-23 TxDb 验收补丁：TxDb 的公开转录本/基因名称不一定等于 GFF 的 ID。
验收现通过注释中明确的 `ID`、`Name`、`transcript_id`、`gene_id` 建立映射，
并继续严格核对完整 gene Parent 集、染色体、链、TSS 和导出的 transcript 范围。
重复或歧义别名会停止，不按版本后缀、物种或坐标猜测。
`TxDb-identity-map.tsv`、`TxDb-gene-identity-map.tsv` 和 `TxDb-identity-summary.json`
记录对应关系与校验结果。正式 ChIPseeker 建库继续使用相同的默认命名方式，
上述映射可用于关联 GFF/GTF 原始 ID；不需要改写已通过校验的 GFF/GTF。

1.2.2 的内置规范化器把部分已有但未识别的父类型误当成缺失 gene，
可能改写 Parent、合成多余 gene，并缩短已有 gene/transcript 坐标。
1.2.3 修复的是通用关系处理逻辑，不根据物种、染色体名或 ID 前缀判断。
旧版本生成的 GFF/GTF 不会因安装新版自动修正，必须从原始 GFF 重新生成。

## 常规使用

```bash
chromaflow normalize-gff -c chromaflow.config \
  --input /data/reference/original.gff3.gz \
  --output /data/reference/genes.normalized.gff3 --engine builtin
chromaflow reference -c chromaflow.config --skip-organelle --no-measure
```

输出文件默认使用 `.normalized.gff3`；显式指定旧 `.agat.gff3` 名称仍可用。
`auto` 固定解析为 `builtin`，不再由系统是否安装 AGAT 决定算法。
`--force` 只强制重建，**不授权破坏原始关系或放宽结构校验**。
manifest 登记成功产物；配置中的显式 `GENOME_GFF3`、`GENOME_GTF` 仍优先，
所以旧项目应确认这些路径与本次输出一致。

## 保守规则与边界

- 原始 feature 的类型、ID、Parent（含多个父节点）、坐标、链、CDS phase 和属性全部保留。
  合法分段 feature 可以共享 ID；跨类型、序列、链或 Parent 集不一致的共享 ID 会报错，
  需要单独审查，不自动合并。`Derives_from` 等属性原样保留；本版图校验针对 Parent。
- 未知类型与缺失 ID 分开处理。已有未知父节点保留并报告 `unknown-parent-type`。
  角色表是随包提供的、有限的标准类型及别名映射，**不是完整 SO 本体推理器**。
- 真正缺失的 Parent、循环关系、非法字段/坐标、冲突 ID 导致失败；错误不会被静默丢行。
  无父节点的顶层 feature 合法；无父节点的 transcript 会警告，GTF 导出需要明确基因归属。
- 缺 exon 时，只用已知 CDS/UTR 区间并集补建；不把无子特征的 transcript 整段当 exon。
  有原生 exon、可信 CDS 和一致明确链方向时，才补终端未覆盖 UTR。
  CDS 内部间隙不归为 3′UTR；不改写已有泛化 `UTR` 或提供者类型名。
- 不缩短 gene/transcript，不自动插入 mRNA 层。合法的 gene→CDS、圆形跨边界、
  跨序列等模型以保留/报告为原则；不安全的推断跳过。
  **直接 gene→CDS/exon、无法解释的模型，以及部分特殊结构的自动 GTF 导出暂不支持**，
  会明确停止。这不是该原始 GFF 无效，也不是该物种不受支持；需单独兼容导出及映射。
- `agat`、`gffread` 规范化也走同一保持校验。外部工具若改写/删除原有记录，结果不会发布。
  本版没有“忽略身份保持检查”开关。优先使用 builtin。
- 注释行顺序保留，新增行附在 feature 区末尾；支持向前引用。保留注释、元信息和内嵌 FASTA。
  因新增行可能改变分组，移除位置性的 `###` 分隔符并计数报告。
  输入按魔数支持 plain/gzip/bgzip/bzip2/xz、单文件 zip；zstd 输入需要本机 `zstd`。

只有在审查后确定缺失父节点确为 gene 时，才显式启用：

```bash
chromaflow normalize-gff --input original.gff3 --output repaired.gff3 \
  --engine builtin --repair-missing-parents
```

该开关仅恢复被 transcript 的原 Parent 明确引用的缺失 ID；所有引用者必须具有一致序列、链，
且只有该单一父节点。恢复一个共享 gene，不为每个转录本另造 `.gene`；没有 Parent 证据时不猜。

提供者扩展类型可以用本地 JSON 配置，例如 `{"provider_gene_model":"gene"}`。
命令支持 `--type-map FILE`；若后续还需要 GTF，应在项目配置中统一设置 `GFF_TYPE_MAP="/absolute/types.json"`。
可用角色为 `gene, transcript, exon, cds, utr, utr5, utr3`，禁止覆盖标准映射。
角色扩展只改变 ChromaFlow 的解释，不保证 gffread/TxDb 也支持该类型；下游仍必须通过检查。

## 报告、缓存与失败保护

每个规范化 GFF 附带：

- `.report.json`：输入/输出 SHA256、实际引擎及版本、策略/规则哈希、各类型计数、警告及示例。
- `.changes.tsv`：每条新增 feature 的 ID、类型、Parent、坐标和依据。

只有源内容、实现、角色映射、选项、引擎及其版本、输出和变更表哈希都一致才复用缓存。
旧版产物缺少这些元数据会重建；不能以“文件存在”作为有效性证据。
临时结果验证通过后原子替换 GFF，再登记 manifest；失败保留已有正式 GFF。
多个同目标文件的写操作应串行运行；不提供并发写锁。

`reference` 按 GFF 内容及转换规则使 GTF 缓存失效，旧版无报告的 GTF 会重新生成。
发布前检查所有识别到的 transcript→gene 归属及 exon 区间，阻止丢失基因归属或凭空生成 exon。
GTF 有自己的 `.report.json`。无 gffread 时明确失败；不需要 GTF 可使用 `--skip-gtf`。
本版不改变 MACS_GSIZE 算法，也不要求重新下载 nf-core 两套镜像。

## 在目标服务器验收

```bash
# CHROMAFLOW_HOME 需指向当前安装的 ChromaFlow 根目录。
python3 "$CHROMAFLOW_HOME/lib/steps/validate_annotation.py" \
  /data/reference/original.gff3 /data/reference/genes.normalized.gff3 \
  --gtf /data/reference/genes.gtf --report-dir /data/project/gff-acceptance

# 激活提供 Rscript / GenomicFeatures / txdbmaker 的注释环境后：
python3 "$CHROMAFLOW_HOME/lib/steps/validate_annotation.py" \
  /data/reference/original.gff3 /data/reference/genes.normalized.gff3 \
  --gtf /data/reference/genes.gtf --txdb --report-dir /data/project/txdb-acceptance
```

第二条读取实际安装的 R 包，核对原 transcript 对应的基因及 TSS，保留建库警告、R 版本和差异表。
R 注释步骤也不再隐藏建 TxDb 警告。不要把 GFF 结构通过等同于 TxDb 已验收。
此验收工具不修改输入文件。标准源中若确有缺失 Parent，先单独审查修复；原始结构校验不绕过错误。

本地发布验证使用合成结构案例和真实离线 gffread 0.12.7，不冒充服务器完整数据验收。
目标服务器的 810 条真实 TE 关系、现装 R 包和完整分析结果仍应运行上述验收。

类型/关系依据：[GFF3 规范](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)、
[Sequence Ontology](https://github.com/The-Sequence-Ontology/SO-Ontologies)、
[txdbmaker 源码](https://github.com/Bioconductor/txdbmaker/blob/devel/R/makeTxDbFromGRanges.R)。
