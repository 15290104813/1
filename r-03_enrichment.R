rm(list = ls())
setwd("//")
if (!dir.exists("03_enrichment")) {dir.create("03_enrichment")}
setwd("03_enrichment")

library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(Rgraphviz)
library(ggnewscale)
library(DOSE)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggraph)
library(tidygraph)
gene <- read.csv('../02_intersect/02.Shared_genes.csv')
gene_transform <- bitr(gene$x,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")

ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
go_result <- data.frame(ego)
go_result2 <- go_result
go_result <- subset(go_result,p.adjust<0.05)  #2180
table(go_result$ONTOLOGY)
write.csv(go_result,file = "01.GO.csv",row.names = T)
# BP  CC  MF 
# 2363  113  144 

display_number = c(5, 5, 5)  
go_result <- go_result[order(go_result$p.adjust,decreasing = F),]
BP <- go_result[which(go_result$ONTOLOGY=='BP'),]
CC <- go_result[which(go_result$ONTOLOGY=='CC'),]
MF <- go_result[which(go_result$ONTOLOGY=='MF'),]
go_result_BP = as.data.frame(BP)[1:display_number[1], ]
go_result_CC = as.data.frame(CC)[1:display_number[2], ]
go_result_MF = as.data.frame(MF)[1:display_number[3], ]

go_enrich = data.frame(
  ID=c(go_result_BP$ID, go_result_CC$ID, go_result_MF$ID),                      
  Description=c(go_result_BP$Description,go_result_CC$Description,go_result_MF$Description),
  Count=c(go_result_BP$Count, go_result_CC$Count, go_result_MF$Count),
  type=factor(c(rep("Biological Process", display_number[1]), 
                rep("Cellular Component", display_number[2]),
                rep("Molecular Function", display_number[3])),
              levels=c("Biological Process", "Cellular Component","Molecular Function" )))

go_enrich <- na.omit(go_enrich)

go_enrich$type_order = factor(go_enrich$Description,levels=go_enrich$Description,ordered = T)


go_enrich$Description
# [1] "epithelial cell proliferation"               "regulation of epithelial cell proliferation"
# [3] "renal system development"                    "muscle system process"                      
# [5] "kidney development"                          "cytoplasmic vesicle lumen"                  
# [7] "vesicle lumen"                               "secretory granule lumen"                    
# [9] "collagen-containing extracellular matrix"    "apical part of cell"                        
# [11] "glycosaminoglycan binding"                   "heparin binding"                            
# [13] "sulfur compound binding"                     "growth factor activity"                     
# [15] "lipoprotein particle binding"     


library(ggpubr)
head <- rbind(go_result_BP,go_result_CC,go_result_MF)
head$GeneRatio <- head$Count / 861
colnames(head)
p0<-ggdotchart(head, x = "Description", y = "GeneRatio",
               dot.size = 'Count',
               color ='p.adjust',
               sorting = "descending",
               add = "segments",                          
               rotate = TRUE,
               ggtheme = theme_pubr(),                      
               ylab="GeneRatio",
               xlab='',
               title=""
)+scale_colour_gradient( high = "#F781BF",low = "#88c4e8")
p10<-p0+theme(legend.position = "right",
              panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
  theme(axis.title.x =element_text(size=14,family = "Times", face = "bold"),
        axis.text.x =element_text(size=12,family = "Times", face = "bold"),
        axis.title.y =element_text(size=14,family = "Times", face = "bold"),
        axis.text.y=element_text(size=12,family = "Times", face = "bold"),
        plot.title=element_text(size=14,family = "Times", face = "bold",hjust=0.5),
        legend.text = element_text(size = 14, family = "Times"),
        legend.title = element_text(size = 14, family = "Times",face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
pdf("01.GO.pdf", width = 9, height = 6,family = "Times")
p10
dev.off()
png("01.GO.png", width = 9, height = 6,res = 600,units = 'in',family = "Times")
p10
dev.off()



KEGG_database <- 'hsa'
kk <- enrichKEGG(gene = gene_transform$ENTREZID,organism =KEGG_database, 
                 keyType = "kegg",pAdjustMethod = "none",
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)
enrichKK<-DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
KEGG_results<-data.frame(enrichKK) 

KEGG_results2 <- KEGG_results
KEGG_results <- subset(KEGG_results,p.adjust<0.05)  #108

write.csv(KEGG_results,'02.KEGG.csv',row.names = T)


#棒棒糖图--------------
library(ggpubr)
KEGG_results <- KEGG_results[order(KEGG_results$p.adjust,decreasing = F),]
head <- KEGG_results[1:15,]
head$Description
# [1] "AGE-RAGE signaling pathway in diabetic complications"
# [2] "Lipid and atherosclerosis"                           
# [3] "PI3K-Akt signaling pathway"                          
# [4] "Fluid shear stress and atherosclerosis"              
# [5] "Cholesterol metabolism"                              
# [6] "Amoebiasis"                                          
# [7] "Malaria"                                             
# [8] "Cytoskeleton in muscle cells"                        
# [9] "Human papillomavirus infection"                      
# [10] "Type II diabetes mellitus"                           
# [11] "Cytokine-cytokine receptor interaction"              
# [12] "Proteoglycans in cancer"                             
# [13] "TGF-beta signaling pathway"                          
# [14] "Adipocytokine signaling pathway"                     
# [15] "Hedgehog signaling pathway" 

head$GeneRatio <- head$Count / 591
colnames(head)
p0<-ggdotchart(head, x = "Description", y = "GeneRatio",
               dot.size = 'Count',
               color ='p.adjust',
               sorting = "descending",
               add = "segments",                             
               rotate = TRUE,
               ggtheme = theme_pubr(),                       
               ylab="GeneRatio",
               xlab='',
               title=""
)+scale_colour_gradient( high = "#e8c559",low = "#f06152")
p10<-p0+theme(legend.position = "right",
              panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
  theme(axis.title.x =element_text(size=14,family = "Times", face = "bold"),
        axis.text.x =element_text(size=12,family = "Times", face = "bold"),
        axis.title.y =element_text(size=14,family = "Times", face = "bold"),
        axis.text.y=element_text(size=12,family = "Times", face = "bold"),
        plot.title=element_text(size=14,family = "Times", face = "bold",hjust=0.5),
        legend.text = element_text(size = 14, family = "Times"),
        legend.title = element_text(size = 14, family = "Times",face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p10
pdf("02.KEGG.pdf", width = 10, height = 6,family = "Times")
p10
dev.off()
png("02.KEGG.png", width = 10, height = 6,res = 600,units = 'in',family = "Times")
p10
dev.off()

