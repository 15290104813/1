rm(list = ls())
setwd("//")
if (!dir.exists("05_single")) {dir.create("05_single")}
setwd("05_single")

#eqtl--------------------------------------------
if (!dir.exists("eqtl")) {dir.create("eqtl")}

#install.packages("/data2/marui/pipeline/package/TeachingDemos_2.13.tar.gz", repos = NULL, type = "source")
library(magick)
library(TeachingDemos)
source("/data2/marui/pipeline/SMR/plot_SMR.r") 

# ENSG00000124813——RUNX2
# ENSG00000204217——BMPR2

#ENSG00000124813——RUNX2------------------
SMRData = ReadSMRData('//04_SMR/plot/eqtl.ENSG00000124813.txt')

#绘图步骤可能会报错 报错为正常 可以绘制出内容即可
pdf("eqtl/01.RUNX2_locusPlot.pdf", width=13, height=9,family = "Times")
SMRLocusPlot(data=SMRData, smr_thresh=0.0001, heidi_thresh=0.05, plotWindow=200, max_anno_probe=16)
dev.off()

png("eqtl/01.RUNX2_locusPlot.png", width = 13, height = 9,res = 600,units = 'in',family = "Times")
SMRLocusPlot(data=SMRData, smr_thresh=0.0001, heidi_thresh=0.05, plotWindow=200, max_anno_probe=16)
dev.off()

# smr_thresh：SMR测试的全基因组显著性水平
# heidi_thresh：异质性测试的阈值，默认0.05
# plotWindow：围绕探针中心选择顺式eqtl进行绘图的窗口大小，默认2000个字节kb
# max_anno_probe：图形上现实的最大探针名称数量，默认16

pdf("eqtl/02.RUNX2_effectPlot.pdf", width=7, height=5.5,family = "Times")
SMREffectPlot(data=SMRData, trait_name="RUNX2") 
dev.off()

png("eqtl/02.RUNX2_effectPlot.png", width = 7, height = 5.5,res = 600,units = 'in',family = "Times")
SMREffectPlot(data=SMRData, trait_name="RUNX2") 
dev.off()




#ENSG00000204217——BMPR2------------------
SMRData = ReadSMRData('//04_SMR/plot/eqtl.ENSG00000204217.txt')

#绘图步骤可能会报错 报错为正常 可以绘制出内容即可
pdf("eqtl/01.BMPR2_locusPlot.pdf", width=13, height=9,family = "Times")
SMRLocusPlot(data=SMRData, smr_thresh=0.0001, heidi_thresh=0.05, plotWindow=200, max_anno_probe=16)
dev.off()

png("eqtl/01.BMPR2_locusPlot.png", width = 13, height = 9,res = 600,units = 'in',family = "Times")
SMRLocusPlot(data=SMRData, smr_thresh=0.0001, heidi_thresh=0.05, plotWindow=200, max_anno_probe=16)
dev.off()

# smr_thresh：SMR测试的全基因组显著性水平
# heidi_thresh：异质性测试的阈值，默认0.05
# plotWindow：围绕探针中心选择顺式eqtl进行绘图的窗口大小，默认2000个字节kb
# max_anno_probe：图形上现实的最大探针名称数量，默认16

pdf("eqtl/02.BMPR2_effectPlot.pdf", width=7, height=5.5,family = "Times")
SMREffectPlot(data=SMRData, trait_name="BMPR2") 
dev.off()

png("eqtl/02.BMPR2_effectPlot.png", width = 7, height = 5.5,res = 600,units = 'in',family = "Times")
SMREffectPlot(data=SMRData, trait_name="BMPR2") 
dev.off()
