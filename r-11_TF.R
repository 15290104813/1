rm(list = ls())
setwd("//")
if (! dir.exists("./11_TF")){
  dir.create("./11_TF")
}
setwd("./11_TF")

library(tidyverse)
library(magrittr)
library(tidygraph)
library(ggraph)
library(igraph)
library(MetBrewer)
library(ggforce) 


df1 <- read.csv('TFs.csv') %>%
  select(1) %>% dplyr::rename("Var"="gene") %>% 
  bind_rows(read.csv('TFs.csv')  %>% 
              select(2) %>% dplyr::rename("Var"="TF")) %>% 
  separate(col=Var,into = c("gene","Type"), sep = " ") %>%
  mutate(Type=case_when(grepl("RUNX2",gene) ~ "mRNA",
                        TRUE ~ "TFs")) %>% 
  mutate(Type=case_when(grepl("BMPR2",gene) ~ "mRNA",
                        TRUE ~ Type))

nodes <- df1 %>% distinct()

edges <- read.csv('TFs.csv') %>% select(1,2) %>% 
  separate(col=gene,into = c("gene"), sep = " ") %>% 
  separate(col=TF,into = c("TF"), sep = " ") %>% 
  set_colnames(c("from","to"))


graph <- graph_from_data_frame(edges,nodes,directed = FALSE)
V(graph)$degree <- igraph::degree(graph, mode = "all")
tidy_graph <- tidygraph::as_tbl_graph(graph)


p <- ggraph::ggraph(graph = tidy_graph, layout = "fr") + 
  geom_edge_link(aes(width = 0.1), color = "#bac4d0", alpha = 0.5) +
  geom_node_point(aes(size = sqrt(degree), color = Type,shape = Type),
                  alpha = 0.7,show.legend = FALSE) +
  geom_node_text(aes(label = name), #, filter = degree > median(degree)
                 repel = TRUE,
                 size = 3,
                 fontface = "bold",
                 color = "gray10",
                 box.padding = unit(0.35, "lines"),
                 point.padding = unit(0.3, "lines")) +
  scale_color_manual(values=c("#F0027F","#F98400"))+
  scale_edge_width(range = c(0.5, 1.5)) +
  scale_size(range = c(5, 10)) + 
  ggraph::theme_graph() +
  theme(legend.position = "none")
p
pdf("01.cluster.pdf", width = 7, height = 6,family = "Times")
print(p)
dev.off()
png("01.cluster.png", width = 7, height = 6,res = 600,units = 'in',family = "Times")
print(p)
dev.off()
