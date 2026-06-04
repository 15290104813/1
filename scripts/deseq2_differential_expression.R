# ============================================================================
# R语言生信分析：DESeq2差异表达分析
# ============================================================================
# 功能：基于计数矩阵进行RNA-seq差异表达分析
# 依赖包：DESeq2
# 输入：计数矩阵（csv） + 样本信息表（csv）
# 输出：差异表达结果 + 可视化图表
# ============================================================================

# 1. 加载必要的R包
# ============================================================================
library(DESeq2)

# 如未安装，使用以下命令安装
# if (!require("DESeq2", quietly = TRUE)) {
#   BiocManager::install("DESeq2")
# }

# 2. 设置工作目录和路径
# ============================================================================
# 修改为您的实际路径
setwd("./")
data_dir <- "data"
results_dir <- "results"

# 创建输出目录（如不存在）
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# 3. 读取数据
# ============================================================================
# 读取基因计数矩阵 (行为基因，列为样本)
countData <- as.matrix(read.csv(file.path(data_dir, "count_matrix.csv"), 
                                 row.names = 1))

# 读取样本信息 (必须包含condition列用于分组)
sampleInfo <- read.csv(file.path(data_dir, "sample_info.csv"), 
                        row.names = 1)

# 验证数据
cat("计数矩阵维度：", dim(countData), "\n")
cat("样本信息维度：", dim(sampleInfo), "\n")
cat("\n样本信息头部：\n")
print(head(sampleInfo))

# 4. 创建DESeq2分析对象
# ============================================================================
# 确保样本顺序一致
stopifnot(colnames(countData) == rownames(sampleInfo))

# 创建DESeqDataSet对象
# design参数指定分析模型（这里为简单的单因子设计）
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = sampleInfo,
  design = ~ condition
)

# 5. 数据质量控制（可选但推荐）
# ============================================================================
# 去除低表达基因（保留至少在10个样本中counts >= 10的基因）
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]
cat("过滤后保留基因数：", nrow(dds), "\n")

# 6. 运行DESeq2差异表达分析
# ============================================================================
dds <- DESeq(dds)

# 7. 提取分析结果
# ============================================================================
# 获取结果（默认对比：最后一个condition组 vs 第一个condition组）
res <- results(dds)

# 转换为数据框并添加基因名
res_df <- as.data.frame(res)
res_df$gene <- rownames(res)
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]

# 按p值排序
res_df <- res_df[order(res_df$pvalue, na.last = NA), ]

# 8. 筛选显著差异基因
# ============================================================================
# 标准筛选条件：|log2FoldChange| >= 1 且 padj < 0.05
threshold_lfc <- 1
threshold_padj <- 0.05

sig_genes <- subset(res_df, 
                    abs(log2FoldChange) >= threshold_lfc & 
                    padj < threshold_padj)

# 分别获取上调和下调基因
up_genes <- subset(sig_genes, log2FoldChange > 0)
down_genes <- subset(sig_genes, log2FoldChange < 0)

cat("\n========== 差异分析结果摘要 ==========\n")
cat("显著差异基因总数：", nrow(sig_genes), "\n")
cat("显著上调基因数：", nrow(up_genes), "\n")
cat("显著下调基因数：", nrow(down_genes), "\n")

# 9. 导出结果
# ============================================================================
# 导出全部结果
write.csv(res_df, 
          file = file.path(results_dir, "all_genes_deseq2_results.csv"),
          row.names = FALSE)

# 导出显著差异基因
write.csv(sig_genes,
          file = file.path(results_dir, "significant_genes.csv"),
          row.names = FALSE)

# 导出上调基因
write.csv(up_genes,
          file = file.path(results_dir, "upregulated_genes.csv"),
          row.names = FALSE)

# 导出下调基因
write.csv(down_genes,
          file = file.path(results_dir, "downregulated_genes.csv"),
          row.names = FALSE)

# 10. 基础可视化
# ============================================================================

# 10.1 绘制MA图 (M-A plot)
png(file = file.path(results_dir, "MA_plot.png"), 
    width = 800, height = 600)
plotMA(res, main = "MA Plot (DESeq2)", 
       ylim = c(-5, 5),
       colNonSig = "gray60",
       colSig = "red",
       colLine = "gray40")
dev.off()

# 10.2 绘制火山图 (Volcano plot)
png(file = file.path(results_dir, "volcano_plot.png"), 
    width = 800, height = 600)
with(res_df, plot(log2FoldChange, -log10(pvalue), 
                   main = "Volcano Plot",
                   xlab = "log2(FoldChange)",
                   ylab = "-log10(p-value)",
                   col = ifelse(abs(log2FoldChange) >= threshold_lfc & 
                               padj < threshold_padj, "red", "gray"),
                   pch = 20, cex = 0.8))
abline(v = c(-threshold_lfc, threshold_lfc), lty = 2, col = "blue")
abline(h = -log10(threshold_padj), lty = 2, col = "blue")
legend("topright", 
       legend = c("Significant", "Not significant"),
       col = c("red", "gray"),
       pch = 20)
dev.off()

# 10.3 绘制前10个显著差异基因的表达量
if (nrow(sig_genes) > 0) {
  png(file = file.path(results_dir, "top10_genes_counts.png"), 
      width = 800, height = 600)
  
  # 获取前10个差异基因
  top_genes <- head(sig_genes$gene, min(10, nrow(sig_genes)))
  
  # 绘制多图
  par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
  for (gene in top_genes) {
    plotCounts(dds, gene = gene, intgroup = "condition", main = gene)
  }
  dev.off()
}

# 10.4 绘制样本聚类热图
png(file = file.path(results_dir, "sample_clustering.png"), 
    width = 800, height = 600)

# 进行VST转换（用于可视化）
vsd <- vst(dds, blind = FALSE)

# 计算样本间距离并聚类
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# 简单的聚类图
par(mar = c(5, 5, 2, 2))
heatmap(sampleDistMatrix,
        main = "Sample-to-Sample Distance",
        labRow = colnames(assay(vsd)),
        labCol = colnames(assay(vsd)))
dev.off()

# 11. 输出分析总结
# ============================================================================
cat("\n========== 分析完成 ==========\n")
cat("结果已保存到：", results_dir, "\n")
cat("输出文件列表：\n")
cat("  1. all_genes_deseq2_results.csv - 全部基因结果\n")
cat("  2. significant_genes.csv - 显著差异基因\n")
cat("  3. upregulated_genes.csv - 上调基因\n")
cat("  4. downregulated_genes.csv - 下调基因\n")
cat("  5. MA_plot.png - MA图\n")
cat("  6. volcano_plot.png - 火山图\n")
cat("  7. top10_genes_counts.png - 前10差异基因表达量\n")
cat("  8. sample_clustering.png - 样本聚类\n")

cat("\n========== 结果解读说明 ==========\n")
cat("baseMean: 基因的平均表达量\n")
cat("log2FoldChange: 对数2倍数变化\n")
cat("lfcSE: log2FoldChange的标准误\n")
cat("stat: Wald检验统计量\n")
cat("pvalue: 原始p值\n")
cat("padj: 校正后p值(Benjamini-Hochberg方法)\n")

# ============================================================================
# 脚本结束
# ============================================================================
