rm(list = ls())
setwd('//')
if (! dir.exists("./01_DEG_GSE51588")){
  dir.create("./01_DEG_GSE51588")
}
setwd("./01_DEG_GSE51588")


library(limma)
df = read.csv('../00_rawdata/dat(GSE51588).csv',row.names = 1)
df.group = read.csv('../00_rawdata/group(GSE51588).csv')
df = df[df.group$sample]
table(df.group$group)
df.group$group = factor(df.group$group, levels = c("Control", "OA"))      
design.mat = cbind(Control = ifelse(df.group$group == "Control", 1, 0),   
                   OA = ifelse(df.group$group == "Control", 0, 1))
contrast.mat = makeContrasts(contrasts="OA-Control", levels=design.mat)  
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
DEG=na.omit(fit)
logFC_cutoff <- 0.5
DEG$change = as.factor(
  ifelse(DEG$P.Value <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')      
)
sig_diff <- subset(DEG,DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff)     
dim(DEG)     ## 
dim(sig_diff)  ##    
summary(sig_diff$change)
write.csv(DEG,file = "01.DEG_all.csv",quote = F,row.names = T)
write.csv(sig_diff,file = "02.DEG_sig.csv",quote = F,row.names = T)



library(dplyr)
library(ggfun)
library(grid)
library(ggplot2)
library(ggthemes)
library(Ipaper)
#library(scales)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(tidyverse)

DEG$symbol <- rownames(DEG)
volcano_plot <- ggplot(data = DEG) + 
  geom_point(aes(x = logFC, y = -log10(P.Value), 
                 color = logFC,
                 size = -log10(P.Value))) + 
  geom_text_repel(data =  DEG %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(change != "NOT") %>%
                    dplyr::arrange(desc(logFC)) %>%
                    dplyr::slice(1:5,) %>%
                    dplyr::filter(change == "UP"),
                  aes(x = logFC, y = -log10(P.Value), label = symbol),
                  box.padding = 1,
                  nudge_x = 3,
                  nudge_y = 0,
                  segment.curvature = 0,
                  segment.ncp = 3,
                  segment.angle = 0,
                  direction = "both",
                  hjust = "left",
                  max.overlaps = 200
  )+
  geom_text_repel(data =  DEG %>%
                    tidyr::drop_na() %>%
                    dplyr::filter(change != "NOT") %>%
                    dplyr::filter(change != "UP") %>%
                    dplyr::arrange(desc(-logFC)) %>%
                    dplyr::slice(1:5) %>%
                    dplyr::filter(change == "DOWN"),
                  aes(x = logFC, y = -log10(P.Value), label = symbol),
                  box.padding = 1,
                  nudge_x = -3,
                  nudge_y = 0,
                  segment.curvature = 0,
                  segment.ncp = 3,
                  segment.angle = 0,
                  direction = "both", 
                  hjust = "left",
                  max.overlaps = 200
  ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  scale_size(range = c(1,7)) + 
  ggtitle(label = "Volcano Plot") + 
  #xlim(c(-15, 15)) + 
  ylim(c(-0.5, 25)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5)
  ) + 
  #annotate(geom = "text", x = 4, y = 0, label = "adj.p = 0.05", size = 5) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")
    ), 
    xmin = (-logFC_cutoff)-5, 
    xmax = -logFC_cutoff,
    ymin = 24,
    ymax = 24
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "DOWN",
      gp = grid::gpar(col = "#74add1")
    ),
    xmin = (-logFC_cutoff)-5, 
    xmax = -logFC_cutoff,
    ymin = 24,
    ymax = 24
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#d73027")
    ), 
    xmin = logFC_cutoff+5, 
    xmax = logFC_cutoff,
    ymin = 24,
    ymax = 24
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "UP",
      gp = grid::gpar(col = "#d73027")
    ),
    xmin = logFC_cutoff+5, 
    xmax = logFC_cutoff,
    ymin = 24,
    ymax = 24
  ) 
volcano_plot
pdf("03.volcano_DEG.pdf", width = 8, height = 6,family = "Times")
volcano_plot
dev.off()
png("03.volcano_DEG.png", width = 8, height = 6,res = 600,units = 'in',family = "Times")
volcano_plot
dev.off()


library(ComplexHeatmap)
library(circlize)
sig_diff <- sig_diff[order(sig_diff$logFC),]
dat_rep<-sig_diff[rownames(sig_diff)%in%
                    rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                                   head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]
#dat_rep<-dat_rep[-c(21,22),]
rt <- df
rt <- rt[rownames(dat_rep),]
group <- df.group[order(df.group$group),]
rt <- rt[,group$sample]
#rt <- rt[order(rownames(rt)),]
x<-rt
mat <- t(scale(t(x)))#归一化
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2

pdf('04.heatmap.pdf',  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("OA" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(7, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("darkgreen", "white","orange"))(100))
dev.off()

png('04.heatmap.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("OA" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(7, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("darkgreen", "white","orange"))(100))
dev.off()



