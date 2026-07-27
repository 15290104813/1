rm(list = ls())
setwd("//")
if (!dir.exists("05_plot")) {dir.create("05_plot")}
setwd("05_plot")


#install.packages("/data2/marui/pipeline/package/CMplot_4.5.1.tar.gz", repos = NULL, type = "source")
library(dplyr)
library(tidyr)
library(CMplot)

#eqtl-----
if (!dir.exists("eqtl")) {dir.create("eqtl")}
setwd("eqtl")

rt=read.csv('../../04_SMR/results_eqtl/results_eqtl.csv')
rt=na.omit(rt)
rt=rt[rt$p_HEIDI>0.05,]


geneRT=read.csv('../../04_SMR/results_eqtl/SMR_eqtl.csv', header=T, sep=",", check.names=F)
geneName=as.vector(geneRT[,"SYMBOL"])


data=rt[,c("SYMBOL","ProbeChr","Probe_bp","fdr")]
colnames(data)=c("Gene","chr","bp","pvalue")


CMplot(data,  plot.type="m",
       LOG10=TRUE, threshold=0.05, chr.den.col=NULL,ylab="-log10(FDR)",
       highlight=geneName, highlight.text=geneName, highlight.cex=1.2, highlight.text.cex=1,
       file="pdf", file.output=TRUE, width=15, height=8.5, verbose=TRUE)

CMplot(data,  plot.type="m",
       LOG10=TRUE, threshold=0.05, chr.den.col=NULL,ylab="-log10(FDR)",
       highlight=geneName, highlight.text=geneName, highlight.cex=1.2, highlight.text.cex=1,
       file="png", file.output=TRUE, width=15, height=8.5, verbose=TRUE)


CMplot(data,  plot.type="c",
       LOG10=TRUE, threshold=0.05, chr.den.col=NULL,
       highlight=geneName, highlight.text=geneName, highlight.cex=1.2, highlight.text.cex=1,
       file="pdf", file.output=TRUE, width=7, height=7, verbose=TRUE)

CMplot(data,  plot.type="c",
       LOG10=TRUE, threshold=0.05, chr.den.col=NULL,
       highlight=geneName, highlight.text=geneName, highlight.cex=1.2, highlight.text.cex=1,
       file="png", file.output=TRUE, width=7, height=7, verbose=TRUE)
