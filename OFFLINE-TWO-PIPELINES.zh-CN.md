# ChIP/DAP 与 ATAC 两套 nf-core 离线包

本页适用于 ChromaFlow 1.2.7。两套 bundle 独立保存，不能混合内部文件。

| 用途 | Pipeline | 版本 | Nextflow | 通用 bundle 目录示例 |
|---|---|---|---|---|
| ChIP-seq / DAP-seq | nf-core/chipseq | 2.1.0 | 25.10.4 | `$HOME/nfcore-offline/chipseq-2.1.0` |
| ATAC-seq | nf-core/atacseq | 2.1.2 | 25.10.4 | `$HOME/nfcore-offline/atacseq-2.1.2` |

联网准备端：

```bash
conda activate nfcore_env
chromaflow offline prepare --type chipseq --bundle "$HOME/nfcore-offline/chipseq-2.1.0"
chromaflow offline prepare --type atacseq --bundle "$HOME/nfcore-offline/atacseq-2.1.2"
for NAME in chipseq-2.1.0 atacseq-2.1.2; do
  chromaflow offline check --bundle "$HOME/nfcore-offline/$NAME" --verify --runtime-check
  chromaflow offline pack --bundle "$HOME/nfcore-offline/$NAME" --archive "$HOME/$NAME.tar"
done
```

完整安装、上传、共享目录权限和服务器验证命令见 [docs/OFFLINE.md](docs/OFFLINE.md)。
该指南的 `USER`、`SERVER`、`$HOME` 和 `/opt/chromaflow/bundles` 是通用示例，需按自己的部署填写。

每个项目选择与 `--type` 对应的 bundle。DAP-seq 使用 chipseq 包。
软件会准备所需插件，不需要手工改动 bundle 内的文件；升级 pipeline 或 Nextflow 时创建新目录。
升级 ChromaFlow 本身可以继续复用经过验证的这两套现有包。

当前检查确认本地 workflow/插件能够启动，不代替代表性数据的实际运行。
历史个人电脑、服务器和项目路径不属于公共安装要求，已从此文移除。
