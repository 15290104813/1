rm(list = ls())
setwd("//")
if (! dir.exists("./02_intersect")){
  dir.create("./02_intersect")
}
setwd("./02_intersect")

gene1 <- read.csv('../00_rawdata/hypercholesterolemia.csv')    
gene2 <- read.csv('../01_DEG_GSE51588/02.DEG_sig.csv')  

library(VennDiagram)
venn_list<-list(gene1=gene1$Gene.Symbol,
                gene2=gene2$X
)
venn.plot <-venn.diagram(venn_list,  filename = NULL ,
                         force.unique = T,     
                         print.mode = c("percent", "raw"),
                         fill = c("#A6D854",'orange'), alpha = 0.6,
                         col = FALSE,
                         cex = 1.2,                                   
                         cat.fontface = "bold",                         
                         cat.col = c("darkgreen","#FF7F00"),  
                         category.names = c('HRGs','DEGs'),    
                         cat.pos = c(20,-20),
                         cat.cex = 1.3,                                    
                         cat.fontfamily = 'serif',
                         margin = 0.1,
                         scaled =FALSE
)

inter <- get.venn.partitions(venn_list)
pdf("01.Shared_genes.pdf")
grid.draw(venn.plot)
dev.off()

venn.diagram(venn_list, filename = '01.Shared_genes.png', imagetype = 'png',
             fill = c("#A6D854",'orange'), alpha = 0.6,
             col = FALSE,
             cat.fontface = "bold",
             cex = 1.2,
             force.unique = T,
             cat.col = c("darkgreen","#FF7F00"),
             category.names = c('HRGs','DEGs'),
             print.mode = c("percent", "raw"),
             fontfamily = 'serif',
             cat.pos = c(20,-20),
             cat.cex = 1.3,
             margin = 0.1,
             scaled =FALSE)
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
Candidate_genes<-inter[1,]$..values.. $`1`
write.csv(Candidate_genes,'02.Shared_genes.csv',quote = F,row.names = F)

library(clusterProfiler)
gene_transform <- bitr(Candidate_genes,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
write.csv(gene_transform,'03.gene_transform.csv',quote = F,row.names = F)
