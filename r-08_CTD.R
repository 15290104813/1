rm(list = ls())
setwd("//")
if (! dir.exists("./08_CTD")){
  dir.create("./08_CTD")
}
setwd("./08_CTD")

gene <- read.csv("RUNX2.csv")
colnames(gene)
gene1 <- subset(gene,Direct.Evidence != "" | Inference.Score>=100)
gene1$gene <- 'RUNX2'
write.csv(gene1,file = 'gene_result_RUNX2.csv',row.names = F)


gene <- read.csv("BMPR2.csv")
colnames(gene)
gene2 <- subset(gene,Direct.Evidence != "" | Inference.Score>=100)
gene2$gene <- 'BMPR2'
write.csv(gene2,file = 'gene_result_BMPR2.csv',row.names = F)


library(tidyverse)
library(magrittr)
library(tidygraph)
library(ggraph)
library(igraph)
library(MetBrewer)
library(ggforce)

df <- rbind(gene1,gene2)
df <- df[,c('gene','Disease.Name')]
df$Disease.Name <-gsub(' ','_',df$Disease.Name,fixed = T)

df1 <- df %>% 
  select(1) %>% dplyr::rename("Var"="gene") %>% 
  bind_rows(df  %>% 
              select(2) %>% dplyr::rename("Var"="Disease.Name")) %>% 
  separate(col=Var,into = c("gene","Type"), sep = " ") %>%
  mutate(Type=case_when(grepl(c("RUNX2"),gene) ~ "mRNA",
                        TRUE ~ 'Disease')) %>% 
  mutate(Type=case_when(grepl(c("BMPR2"),gene) ~ "mRNA",
                        TRUE ~ Type))

nodes <- df1 %>% distinct()

edges <- df %>% 
  separate(col=gene,into = c("gene"), sep = " ") %>% 
  separate(col=`Disease.Name`,into = c("Disease.Name"), sep = " ") %>% 
  set_colnames(c("from","to"))


graph <- graph_from_data_frame(edges,nodes,directed = FALSE)
V(graph)$degree <- igraph::degree(graph, mode = "all")
tidy_graph <- tidygraph::as_tbl_graph(graph)


p <- ggraph::ggraph(graph = tidy_graph, layout = "fr") + 
  geom_edge_link(aes(width = 0.1), color = "#bac4d0", alpha = 0.5) +
  geom_node_point(aes(size = sqrt(degree), color = Type,shape = Type),
                  alpha = 0.9,show.legend = FALSE) +
  geom_node_text(aes(label = name), #,filter = degree > median(degree)
                 repel = TRUE,
                 size = 3,
                 fontface = "bold",
                 color = "gray10",
                 box.padding = unit(0.35, "lines"),
                 point.padding = unit(0.3, "lines")) +
  scale_color_manual(values=c("#984EA3","#FF7F00"))+
  scale_edge_width(range = c(0.5, 1.5)) +
  scale_size(range = c(5, 10)) + 
  ggraph::theme_graph() +
  theme(legend.position = "none")
p
pdf("01.cluster.pdf", width = 11, height = 12,family = "Times")
print(p)
dev.off()

png("01.cluster.png", width = 11, height = 12,res = 600,units = 'in',family = "Times")
print(p)
dev.off()

