# R 语言生信分析 — 快速指南

这个仓库保存用于生物信息学分析的 R 脚本、示例数据说明与运行说明，目标是保持分析可复现并便于复用。

主要内容（示例）：
- r-00_rawdata.R
- r-01_DEG_GSE51588.R
- r-02_intersect.R
- r-03_enrichment.R
- r-05_SMR_code.R
- r-05_SMR_eqtl.R
- r-05_plot.R
- r-05_single.R
- r-06_Locater_eqtl.R
- r-07_MR.R
- r-08_CTD.R
- r-09_scan.R
- r-11_TF.R
- r-12_DGIdb.R
- r-13_Single_cell.R
- scripts/ (辅助脚本)

快速开始

1) 克隆仓库

```bash
git clone https://github.com/15290104813/1.git
cd 1
```

2) 恢复 R 环境（推荐使用 renv）

在 R 中运行：

```r
install.packages('renv', repos = 'https://cloud.r-project.org')
# 若已有 renv.lock，请运行：renv::restore()
```

3) 运行示例脚本（示例）

```bash
# 最简单方式（可能需要根据脚本内部路径调整）：
Rscript r-00_rawdata.R
# 或运行其它脚本，例如：
Rscript r-01_DEG_GSE51588.R
```

使用 Docker（可选）：

```bash
docker build -t r-workflows:latest .
docker run --rm -v $(pwd):/work -w /work r-workflows:latest Rscript r-00_rawdata.R
```

建议和注意事项

- 为保证可复现性，请使用 renv 并提交 renv.lock（或把依赖写在 DESCRIPTION）。
- 不要把大原始数据放在仓库中，建议提供下载脚本或在 .gitignore 中排除大文件。
- 如果仓库使用了 Git LFS，请在目标环境启用 LFS。

许可证

MIT
