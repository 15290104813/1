rm(list = ls())
setwd("/data2/marui/project/11.KXF2025011502/")
if (!dir.exists("04_SMR")) {dir.create("04_SMR")}
setwd("04_SMR")

# if (!dir.exists("besd_eqtl")) {dir.create("besd_eqtl")}

# if (!dir.exists("results_eqtl")) {dir.create("results_eqtl")}

# if (!dir.exists("plot")) {dir.create("plot")}

smr <- read.table("../04_SMR/results_eqtl/multiple.smr",header = T,sep = "\t")
gene_transform <- bitr(smr$probeID,
                       fromType = "ENSEMBL",
                       toType = c('SYMBOL'),
                       OrgDb = "org.Hs.eg.db")
colnames(gene_transform)
results <- merge(gene_transform,smr,by.x = 'ENSEMBL',by.y = 'probeID')
colnames(results)[1] <- 'probeID'

results$or_SMR <- exp(results$b_SMR)
results$lci_SMR <- exp(results$b_SMR - 1.96*results$se_SMR)
results$uci_SMR <- exp(results$b_SMR + 1.96*results$se_SMR)

pvalue <- results[,"p_SMR"]
fdr=p.adjust(pvalue, method="fdr")
results <- cbind(results, fdr)
write.csv(results,'../04_SMR/results_eqtl/results_eqtl.csv',row.names = F)

outTab <- results[((results$p_SMR<0.05)  & (results$p_HEIDI>0.05)) & (results$fdr < 0.05),]#
outTab <- outTab[order(outTab$p_SMR),]
outTab<- rstatix::filter(outTab,!duplicated(outTab$probeID))
write.csv(outTab, file="../04_SMR/results_eqtl/SMR_eqtl.csv", row.names=F)

outTab <- read.csv('../04_SMR/results_eqtl/SMR_eqtl.csv')
out_multi <- outTab
out_multi <- out_multi[order(out_multi$or_SMR,decreasing = F),]


hz <- paste(round(out_multi$or_SMR,5),
            "(",round(out_multi$lci_SMR,5),
            "-",round(out_multi$uci_SMR,5),")",sep = "")


tabletext <- cbind(c(NA,"eQTL",out_multi$SYMBOL),
                   #c(NA,"P value",round(out_multi$pval,5)),
                   c(NA,"P value",ifelse(out_multi$p_SMR<0.001,
                                         "< 0.001",
                                         round(out_multi$p_SMR,5))),
                   c(NA,"Odd Ratio(95% CI)",hz))
tabletext2 <- data.frame(tabletext)
write.csv(tabletext2,'../04_SMR/results_eqtl/tabletext.csv')

library(forestplot)

pdf(file = "../04_SMR/results_eqtl/01.forest.pdf", height = 10, width =10, onefile = F,family='Times')
forestplot(labeltext=tabletext,
           graph.pos=4,  
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="darkred", lines="darkblue", zero = "gray50"),
           mean=c(NA,NA,out_multi$or_SMR),
           lower=c(NA,NA,out_multi$lci_SMR), 
           upper=c(NA,NA,out_multi$uci_SMR), 
           boxsize=0.1,lwd.ci=3, 
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5, 
           colgap=unit(5,"mm"),  
           xticks = c(0.999,1,1.005,1.01,1.015),
           lwd.xaxis=2,           
           lineheight = unit(1.5,"cm"),
           graphwidth = unit(.5,"npc"),
           cex=0.9, fn.ci_norm = fpDrawCircleCI,
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Odd Ratio",
           grid = T,
           clip = c(0,2)) 
dev.off()


png(file = "../04_SMR/results_eqtl/01.forest.png", height = 10, width = 10, units='in',res=600,family='Times')
forestplot(labeltext=tabletext,
           graph.pos=4,  
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="darkred", lines="darkblue", zero = "gray50"),
           mean=c(NA,NA,out_multi$or_SMR),
           lower=c(NA,NA,out_multi$lci_SMR), 
           upper=c(NA,NA,out_multi$uci_SMR), 
           boxsize=0.1,lwd.ci=3,  
           ci.vertices.height = 0.08,ci.vertices=TRUE,
           zero=1,lwd.zero=0.5,    
           colgap=unit(5,"mm"),  
           xticks = c(0.999,1,1.005,1.01,1.015), 
           lwd.xaxis=2,         
           lineheight = unit(1.5,"cm"), 
           graphwidth = unit(.5,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Odd Ratio",
           grid = T,
           clip = c(0,2)) 
dev.off()
