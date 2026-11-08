# ChromaFlow 离线部署指南入口

当前版本为 1.2.8。离线下载、插件准备、校验、打包、上传和运行的统一文档是：

- [README 中文指南](README.zh-CN.md)：安装、环境模式和项目命令。
- [完整离线部署](docs/OFFLINE.md)：联网 Linux/WSL 制包，离线服务器共享安装与验证。
- [两套 pipeline 的选择](OFFLINE-TWO-PIPELINES.zh-CN.md)：ChIP/DAP 与 ATAC 的版本对应。

使用 `chromaflow offline prepare/check/pack` 管理软件包，使用
`chromaflow run --offline --bundle DIR` 运行项目。
容器 bundle 包含 nf-core workflow、Nextflow、插件和容器，不包含参考基因组、
FASTQ、motif 数据库、Java/Apptainer 安装或 ChromaFlow 独立下游环境。

历史个人部署记录已从公开文档中移除；请把当前指南的通用示例替换为自己的路径。
不要将 `--mode code-only` 当作完整离线分析包；纯 Conda 环境迁移也不在这个打包功能范围内。
