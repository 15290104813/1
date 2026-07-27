#Single cell
rm(list = ls())
setwd('//')
if (! dir.exists("./12_Single_cell")){
  dir.create("./12_Single_cell")
}
setwd("./12_Single_cell")

#00_rawdata-----------------------------------
setwd("//12_Single_cell/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")

dir.create("GSE255460_RAW", showWarnings = FALSE)  
target_dir <- file.path(getwd(), "GSE255460_RAW")

# 解压.tar文件到指定的目标文件夹
untar("GSE255460_RAW.tar", exdir = target_dir)

# 获取文件夹中的文件列表
folders <- list.files(target_dir)

# 提取唯一的样本名称（即去除文件类型和后缀）
samples <- unique(sub("_.*", "", folders))

# 循环处理每个样本
for (sample in samples) {
  # 创建每个样本的新文件夹
  sample_dir <- file.path(target_dir, sample)
  dir.create(sample_dir, showWarnings = FALSE)
  
  # 获取属于当前样本的文件
  sample_files <- grep(paste0("^", sample, "_"), folders, value = TRUE)
  
  # 移动并重命名文件
  for (file in sample_files) {
    # 生成新文件名（去掉所有前缀）
    new_name <- sub("^.*_", "", file)
    
    # 如果文件名是 genes.tsv.gz，则重命名为 features.tsv.gz
    if (new_name == "genes.tsv.gz") {
      new_name <- "features.tsv.gz"
    }
    
    # 生成原始文件路径和新文件路径
    old_path <- file.path(target_dir, file)
    new_path <- file.path(sample_dir, new_name)
    
    # 移动并重命名文件
    file.rename(old_path, new_path)
  }
}

# 检查整理后的目录结构
list.files(target_dir, recursive = TRUE)


rm(list = ls())
library(GEOquery)
library(data.table)
library(tidyverse)
#library(lance)
library(magrittr)
library(Seurat)
folders=list.files('./GSE255460_RAW/')
folders
dir=paste0('./GSE255460_RAW/',folders)
dir

i=1
scList <- list()
for (i in c(1:19)) {
  count =  CreateSeuratObject(counts = Read10X(dir[i]),
                              min.cells = 3,    #最小基因数
                              min.features = 200)  #最小表达量
  scList[[i]] <- count
}
saveRDS(scList, file="scList.Rds")

list1 <- c()
for (i in c(2:19)) {
  sample = scList[[i]]
  list1 = append(list1,sample)
}

scRNA<- merge(scList[[1]],
              y = list1,
              add.cell.ids =folders)

saveRDS(scRNA, file="scRNA.Rds")

#group--
gset<-getGEO("GSE255460",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
a=gset[[1]]
pd<-pData(a)
table(pd$characteristics_ch1)
pd$title2 <- c(rep('Control', 3), rep('OA', 16))
group<-data.frame(sample=pd$geo_accession,group=pd$title2)
table(group$group)
group<-group[order(group$group),]       
write.csv(group,file = 'group(GSE255460).csv')


library(Seurat)
table(scRNA$orig.ident)
head(scRNA@meta.data)
length(colnames(scRNA)) #150706 
length(rownames(scRNA)) #34024  

library(tidyverse)
smps = Cells(scRNA) %>% str_sub(1,10)
table(smps)
grps = group
group.label = grps$group[match(smps,grps$sample)]
names(group.label) = Cells(scRNA)
scRNA = AddMetaData(scRNA,group.label,col.name = 'group')

Patients_ID <- tstrsplit(rownames(scRNA@meta.data), "_", fixed = TRUE)[[1]]
scRNA$Patients_ID<- Patients_ID
table(scRNA$Patients_ID)
table(scRNA$group)

###线粒体基因表达占比
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")  #正则表达式，表示以MT-开头；scRNA[["percent.mt"]]这种写法会在meta.data矩阵加上一列

###计算核糖体基因比例
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

#红细胞
scRNA[["percent.hb"]] <- PercentageFeatureSet(scRNA, pattern = "^HB-")

table(scRNA$group)
head(scRNA@meta.data)
system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))



#01_QC----------------------------------------
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./01_QC")){
  dir.create("./01_QC")
}
setwd("./01_QC")
load('../00_rawdata/01.scRNA_orig.Rdata')
### 质控前 ###
length(colnames(scRNA))  #150706 细胞数量
length(rownames(scRNA))  #34024 基因数量

# 设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")   
# 质控前小提琴图
library(patchwork)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       group.by = "group",
                       cols = c("#A6D854","#E78AC3"),
                       # type.by=group,
                       pt.size = 0,
                       features = plot.featrures[i]) + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3,cols = palettes(category = "random",22,show_col=F))
violin
pdf("01.vlnplot_before_qc.pdf", width = 8, height = 6,family = "Times")
violin
dev.off()
png("01.vlnplot_before_qc.png", width = 8, height = 6,res = 600,units = 'in',family = "Times")
violin
dev.off()


#设置质控标准
pctMT=5
minGene=200
maxGene=5000
maxUMI=20000

scRNA <- subset(scRNA, subset =
                  nCount_RNA <maxUMI &
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)
length(colnames(scRNA))   ##128298
length(rownames(scRNA))   ##34024

plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
# 质控后
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, 
                       group.by = "group",
                       cols = c("#A6D854","#E78AC3"),
                       pt.size = 0,
                       features = plot.featrures[i]) + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1,ncol = 3)  
violin
pdf("02.vlnplot_after_qc.pdf", width = 8, height = 6,family = "Times")
violin
dev.off()
png("02.vlnplot_after_qc.png", width = 8, height = 6,res = 600,units = 'in',family = "Times")
violin
dev.off()

system.time(save(scRNA, file = "02.scRNA_qc.Rdata"))  


#02_PCA-----------------------------------
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./02_PCA")){
  dir.create("./02_PCA")
}
setwd("./02_PCA")

load('../01_QC/02.scRNA_qc.Rdata')
scRNA <- JoinLayers(scRNA)  
library(Seurat)
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)

all.genes<-rownames(scRNA.norm)
#all.genes2 <- data.frame(all.genes)
scRNA.norm<-FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)
# Identify the 5 most highly variable genes
top10 <- head(VariableFeatures(scRNA.norm), 5)
#"MSMP"    "CCL20"   "SPP1"    "C2orf40" "CXCL8"
#trace(VariableFeaturePlot,edit=T)
VariableFeaturePlot2 <- function (object, cols = c("black","orange"), pt.size = 1, log = NULL, 
                                  selection.method = NULL, assay = NULL, raster = NULL, raster.dpi = c(512, 
                                                                                                       512)) 
{
  if (length(x = cols) != 2) {
    stop("'cols' must be of length 2")
  }
  hvf.info <- HVFInfo(object = object, assay = assay, method = selection.method, 
                      status = TRUE)
  status.col <- colnames(hvf.info)[grepl("variable", colnames(hvf.info))][[1]]
  var.status <- c("no", "yes")[unlist(hvf.info[[status.col]]) + 
                                 1]
  if (colnames(x = hvf.info)[3] == "dispersion.scaled") {
    hvf.info <- hvf.info[, c(1, 2)]
  }
  else if (colnames(x = hvf.info)[3] == "variance.expected") {
    hvf.info <- hvf.info[, c(1, 4)]
  }
  else {
    hvf.info <- hvf.info[, c(1, 3)]
  }
  axis.labels <- switch(EXPR = colnames(x = hvf.info)[2], variance.standardized = c("Average Expression", 
                                                                                    "Standardized Variance"), dispersion = c("Average Expression", 
                                                                                                                             "Dispersion"), residual_variance = c("Geometric Mean of Expression", 
                                                                                                                                                                  "Residual Variance"))
  log <- log %||% (any(c("variance.standardized", "residual_variance") %in% 
                         colnames(x = hvf.info)))
  plot <- SingleCorPlot(data = hvf.info, col.by = var.status, 
                        pt.size = pt.size, raster = raster, raster.dpi = raster.dpi)
  if (length(x = unique(x = var.status)) == 1) {
    switch(EXPR = var.status[1], yes = {
      cols <- cols[2]
      labels.legend <- "Variable"
    }, no = {
      cols <- cols[1]
      labels.legend <- "Non-variable"
    })
  }
  else {
    labels.legend <- c("Non-variable", "Variable")
  }
  plot <- plot + labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) + 
    scale_color_manual(labels = paste(labels.legend, "count:", 
                                      table(var.status)), values = cols)
  if (log) {
    plot <- plot + scale_x_log10()
  }
  return(plot)
}

pdf(file="01.feature_selection.pdf",width=7,height=4,family='Times',onefile=F)           
theme.set= theme(axis.title.x=element_blank(),
                 axis.title = element_text(size = 16, face = "bold", family = "Times"),
                 axis.text.x = element_text(size = 10,  family = "Times"),
                 axis.text.y = element_text(size = 14,  family = "Times"),
                 legend.position = 'top',legend.direction = 'vertical',
                 legend.text = element_text(size = 14, family = "Times"),
                 legend.title = element_text(size = 16,face='bold',family = "Times"),
                 text = element_text(family = "Times"))
plot1 <- VariableFeaturePlot2(object = scRNA.norm)+theme.set
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1+plot2)
dev.off()

png(file="01.feature_selection.png",width=7,height=4,family='Times',units='in',res=600)             
theme.set= theme(axis.title.x=element_blank(),
                 axis.title = element_text(size = 16, face = "bold", family = "Times"),
                 axis.text.x = element_text(size = 10,  family = "Times"),
                 axis.text.y = element_text(size = 14,  family = "Times"),
                 legend.position = 'top',legend.direction = 'vertical',
                 legend.text = element_text(size = 14, family = "Times"),
                 legend.title = element_text(size = 16,face='bold',family = "Times"),
                 text = element_text(family = "Times"))
plot1 <- VariableFeaturePlot2(object = scRNA.norm)+theme.set
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1+plot2)
dev.off()

#library(harmony)
scRNA.nor.sca<-ScaleData(scRNA.norm)
scRNA.norm.pca<-RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca),npcs = 50)

pdf("02.PCA.pdf",w=6,h=5,family = "Times")
DimPlot(scRNA.norm.pca, reduction = "pca",group.by="Patients_ID",cols = c("#E7298A","#66A61E","#E6AB02","#A6761D","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7",'#669900','#99CC66',"#BD6263","#8EA325","#A9D179","#84CAC0","#F5AE6B","#BCB8D3","#4387B5"))#,split.by = 'group'
dev.off()
png("02.PCA.png",w=6,h=5,family = "Times",units='in',res=600)
DimPlot(scRNA.norm.pca, reduction = "pca",group.by="Patients_ID",cols = c("#E7298A","#66A61E","#E6AB02","#A6761D","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7",'#669900','#99CC66',"#BD6263","#8EA325","#A9D179","#84CAC0","#F5AE6B","#BCB8D3","#4387B5"))#,split.by = 'group'
dev.off()

scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50)
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
plot_pca <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)
plot_pca
pdf("03.pca_cluster.pdf", width = 13, height = 8,family = "Times")
plot_pca
dev.off()
png("03.pca_cluster.png", width = 13, height = 8,res = 600,units = 'in',family = "Times")
plot_pca
dev.off()

plot_elbow <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot_elbow
pdf("04.pca_sd.pdf", width = 6, height = 5,family = "Times")
plot_elbow
dev.off()
png("04.pca_sd.png", width = 6, height = 5,res = 600,units = 'in',family = "Times")
plot_elbow
dev.off()

scRNA.norm.pca.c<-FindNeighbors(scRNA.norm.pca,dims = 1:30)
scRNA.norm.pca.c<-FindClusters(scRNA.norm.pca.c,resolution = 0.5)
UMAP<-RunUMAP(scRNA.norm.pca.c,dims = 1:30)
UMAP<-RunTSNE(UMAP,dims = 1:30)

col = c("#E41A1C","#377EB8","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7",
        "#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#1B9E77","#D95F02","#7570B3",
        "#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7"
)

pdf('05.UMAP.pdf',w=8,h=7,family = "Times")
DimPlot(UMAP,reduction = 'umap',label = T,cols = col)
dev.off()
png('05.UMAP.png',w=700,h=600,family = "Times")
DimPlot(UMAP,reduction = 'umap',label = T,cols = col)
dev.off()

pdf('06.UMAP.group.pdf',w=8,h=7,family = "Times")
DimPlot(UMAP,reduction = 'umap',group.by = 'group',cols = c("#FFFF33","#A65628"))
dev.off()
png('06.UMAP.group.png',w=700,h=600,family = "Times")
DimPlot(UMAP,reduction = 'umap',group.by = 'group',cols = c("#FFFF33","#A65628"))
dev.off()



system.time(save(UMAP, file = "UMAP.Rdata"))
save.image('all.Rdata')

# 03_singleR -----------------------------------------------------
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./03_singleR")){
  dir.create("./03_singleR")
}
setwd("./03_singleR")
load('../02_PCA/UMAP.Rdata')
UMAP2 <- UMAP
UMAP <- UMAP2
library(SingleR)
head(UMAP@meta.data)
table(UMAP$seurat_clusters)

#PMID: 38325908
gene <- c(
  'IFI16','IFI27',  	#preInfC：18
  'CXCL8','GPR183',   #InfC：15
  'COL27A1', 	        #preFC：14
  'SPP1','COL10A1',   #HTC：13
  'HSPA1A','HSPA6','DDIT3',   #HomC：10 11 
  'PRG4','ABI3BP',  	#preHTC：9
  'FRZB','CYTL1',     #EC：6 16
  'C11orf96','BMP2', 	#ProC：5 12
  'CHI3L1','CHI3L2',  #RegC：3 4 7 8 
  'MMP2','COL1A1',    #FC：1 17
  'CILP','OGN'       	#RepC：0 2 
)

theme.set = theme(
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

pdf(file="01.DotPlot.pdf",width=12,height=15,family='Times')
DotPlot(UMAP,features = gene)+theme.set+ coord_flip()+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()

png(file="01.DotPlot.png",width=12,height=15,family='Times',units='in',res=600)
DotPlot(UMAP,features = gene)+theme.set+ coord_flip()+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()

new.cluster.ids = c(
  "0" = "RepC",
  "1" = "FC",
  "2" = "RepC",
  "3" = "RegC",
  "4" = "RegC",
  "5" = "ProC",
  "6" = "EC",
  "7" = "RegC",
  "8" = "RegC",
  "9" = "preHTC",
  "10" = "HomC",
  "11" = "HomC",
  "12" = "ProC",
  "13" = "HTC",
  "14" = "preFC",
  "15" = "InfC",
  "16" = "EC",
  "17" = "FC",
  "18" = "preInfC"
)

#手动注释
UMAP<- RenameIdents(UMAP, new.cluster.ids)
UMAP$celltype = Idents(UMAP)
table(UMAP$celltype)
# RepC      FC    RegC    ProC      EC  preHTC    HomC     HTC   preFC    InfC preInfC 
# 31264   19336   35435   12209    9006    5533    8650    3019    1670    1567     609 
system.time(save(UMAP, file = "UMAP.Rdata")) 

pdf(file="01.DotPlot.pdf",width=12,height=6,family='Times')
DotPlot(UMAP, features = gene,group.by = 'celltype') + RotatedAxis()+
  labs(title="", y="Subtype", x = "",size=40)+
  scale_color_gradientn(colours = c('lightgrey','white',"#fddbc8","#FF7F00"))
dev.off()

png(file="01.DotPlot.png",width=12,height=6,family='Times',units='in',res=600)
DotPlot(UMAP, features = gene,group.by = 'celltype') + RotatedAxis()+
  labs(title="", y="Subtype", x = "",size=40)+
  scale_color_gradientn(colours = c('lightgrey','white',"#fddbc8","#FF7F00"))
dev.off()

library(colorspace)
library(IOBR)
col = c("#fbbab6","#e1c548","#f0e2a3","#5fa664","#4e79a6","#e98741","#8fc0dc","#967568","#f2d3ca","#eebd85","#82c785","#e3d1db","#74a893","#ac9141","#5ac6e9","#ebce8e","#e5c06e","#7587b1","#c7deef","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac",
        "#e0bc58","#64abc0","#fab37f","#e98741","#8fc0dc","#967568","#f2d3ca","#eebd85","#82c785","#edeaa4","#cdaa9f","#794976","#bcacd3","#889b5d","#4e9592","#dbad5f","#64ae79","#ac5092"
        
)

pdf(file="02.celltype.pdf",width=8,height=6,family='Times')
DimPlot(UMAP ,reduction = "umap",label.size = 4,label.color = 'black',label = T,cols = col,repel = T)
dev.off()

png(file="02.celltype.png",width=8,height=6,family='Times',units='in',res=600)
DimPlot(UMAP ,reduction = "umap",label.size = 4,label.color = 'black',label = T,cols = col,repel = T)
dev.off()

UMAP1 <- UMAP
UMAP1@meta.data$seurat_clusters
Idents(UMAP1) <- 'seurat_clusters'

p1 <- DimPlot(UMAP1,reduction = 'umap',label = T,cols = col,repel = T)
p2 <- DimPlot(UMAP ,reduction = "umap",label.size = 4,label.color = 'black',label = T,cols = col,repel = T)

pdf(file="03.celltype.pdf",width=14,height=6,family='Times')
print(p1+p2)
dev.off()

png(file="03.celltype.png",width=1400,height=600,family='Times')
print(p1+p2)
dev.off()

library(circlize)
library(Seurat)
library(tidyverse)
library(stringr)
library(gplots)
library(RColorBrewer)
source('../../../../pipeline/source/plot_circlize.R')

mycolors <-c("#f9766e","#4e79a6","#e1c548","#45337f","#e98741","#e1a4c6","#5fa664","#916ba6","#7587b1","#c7deef","#e97371")

cluster_colors<-mycolors
sce <- UMAP
Idents(sce)=sce$celltype
circ_data <- prepare_circlize_data(sce,scale=0.8)
group_colors<-rand_color(length(names(table(sce$group))))
rep_colors<-rand_color(length(names(table(sce$Patients_ID))))

###plot and save figures
pdf(file="06.celltype_circ.pdf",width=8,height=8,family='Times')
plot_circlize(circ_data,do.label = T, pt.size = 0.2,
              col.use = cluster_colors,bg.color = 'white',
              kde2d.n = 200,
              repel = T,label.cex = 1)+
  add_track(circ_data,group="group",colors=group_colors,track_num=2)+
  add_track(circ_data,group="Patients_ID",colors=rep_colors,track_num=3)
dev.off()

png(file="06.celltype_circ.png",width=8,height=8,family='Times',units='in',res=600)
plot_circlize(circ_data,do.label = T, pt.size = 0.2,
              col.use = cluster_colors,bg.color = 'white',
              kde2d.n = 200,
              repel = T,label.cex = 1)+
  add_track(circ_data,group="group",colors=group_colors,track_num=2)+
  add_track(circ_data,group="Patients_ID",colors=rep_colors,track_num=3)
dev.off()


##tumor
table(UMAP$group)
scRNA_CC<-subset(UMAP,group=='OA')
phe=scRNA_CC@meta.data

table(phe$celltype)
cell_type_CC <- as.data.frame(sort(table(phe$celltype)))
colnames(cell_type_CC)<-c('Celltype','OA')

##Healthy 
scRNA_Normal<-subset(UMAP,group=='Control')
phe=scRNA_Normal@meta.data
table(phe$celltype)
cell_type_Normal <- as.data.frame(sort(table(phe$celltype)))
colnames(cell_type_Normal)<-c('Celltype','Control')

cell_type<-merge(cell_type_Normal,cell_type_CC,by='Celltype')

i<-1
while(i<12){ 
  cell_type$All[i]<-sum(cell_type$CC[i],cell_type$control[i])
  i<-i+1
}

write.csv(cell_type,'cell_type.csv',row.names = F)

##换个形式
all.count <- data.frame(celltype=UMAP$celltype, sample=UMAP$group)
table(all.count$celltype)

library(reshape2)
plot.celltype = table(all.count$celltype, all.count$sample) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype = cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
colnames(plot.celltype)[2:5]
colnames(plot.celltype)[2:5] = c("Count_Control","Count_OA","Proportion_Control","Proportion_OA")
plot.celltype <- plot.celltype[order(plot.celltype$Count_Control,decreasing = T),]
write.csv(plot.celltype, "CellCount.csv", row.names = F,quote = F)

dat.celltype = plot.celltype[c(1,4:5)]

colnames(dat.celltype)[2:3]
colnames(dat.celltype)[2:3] = c("Control","OA")
dat.celltype[2:3] = 100 * dat.celltype[2:3]

dat.plot = melt(dat.celltype, id.vars = "Cell", variable.name = "Sample", value.name = "value")
dat.plot$Cell = factor(dat.plot$Cell, levels = plot.celltype$Cell)
library(IOBR)

color2 <- c('#ABD0F1','#E56F5E','#92B4C8','#F19685','#F6C957','#FFB77F','#FBE8D5','#EEA599',"#cab2d6")
dat.plot$label = round(dat.plot$value, 2) %>% paste0(.,"%")
a <- ggplot(data = dat.plot, mapping = aes(x = value, y = Cell)) + 
  geom_text(aes(x = value + 1, label = label), hjust = 0) +
  geom_bar(aes(fill = Cell), stat = "identity") + 
  facet_wrap(~Sample) +
  #scale_fill_manual(values = palettes(category = "random",22,show_col=F)) + 
  scale_fill_manual(values = col) + 
  xlim(c(0,65)) +
  xlab("Cell Fraction") + ylab("Cell Type") +
  theme_minimal(base_size = 16) + 
  theme(panel.border = element_rect(color = "grey60", fill = "transparent"), legend.position = "none",
        panel.grid = element_blank())+theme.set
a
pdf("07.CellProportion.pdf", width = 8, height = 4,family = "Times")
a
dev.off()
png("07.CellProportion.png", width = 8, height = 4,res = 600,units = 'in',family = "Times")
a
dev.off()

# PMID：38325908
# PMID：38325908
# PMID：38325908

# 04_hubgene------------------------------------
#Single cell
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./04_hubgene")){
  dir.create("./04_hubgene")
}
setwd("./04_hubgene")
#install.packages('tidydr')
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
#BiocManager::install('tidydr')
#library(tidydr)
library(stringr)
library(viridis)
library(ggsci)
library(paletteer)
library(scCustomize)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggpointdensity)
#install.packages("ggpointdensity")
# devtools::install_github(repo = "samuel-marsh/scCustomize")
# remotes::install_github(repo = "samuel-marsh/scCustomize")


load('../03_singleR/UMAP.Rdata')
hubgene <- read.csv('../../04_SMR/results_eqtl/SMR_eqtl.csv')
gene <- hubgene$SYMBOL


pal <- viridis(n = 10, option = "C")
# viridis提供了五种色带：
# viridis：option D，为默认色带，翠绿色；
# magma：option A，岩溶色；
# inferno：option B，火焰色；
# plasma：option C，血色；
# cividis：option E;
i=1
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = UMAP,
                                  colors_use = viridis_magma_dark_high,
                                  features = gene[i])+NoLegend()+NoAxes()+
    theme(panel.border = element_rect(fill = NA,color = "black",
                                      size=1.5,linetype = "solid"))
}
p<-wrap_plots(plots, ncol = 2);p
pdf("01.featureplot.pdf", width = 8, height = 4,family = "Times")
p
dev.off()
png("01.featureplot.png", width = 8, height = 4,res = 600,units = 'in',family = "Times")
p
dev.off()


for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = UMAP,
                                  colors_use = colorRampPalette(c("#3288BD", "white", "#D53E4F" ))(50),
                                  features = gene[i])+NoAxes()#+
  #theme(panel.border = element_rect(fill = NA,color = "black",
  #size=1.5,linetype = "solid"))
}
p<-wrap_plots(plots, ncol = 2);p
pdf("01.hubgene_UMAP.pdf", width = 8, height = 4,family = "Times")
p
dev.off()
png("01.hubgene_UMAP.png", width = 8, height = 4,res = 600,units = 'in',family = "Times")
p
dev.off()




sce.all.int <- UMAP
dat <- Embeddings(sce.all.int, reduction = "umap")
head(dat)

df <- FetchData(object =sce.all.int, vars = c("umap_1", "umap_2", gene), layer = "data")
head(df)

p1 <- ggplot(df, aes(x= umap_1, y= umap_2 )) +
  geom_point(data = df %>% filter(RUNX2 == 0), color = "#440154FF", size = 0.6) +
  ggpointdensity::geom_pointdensity(data = df %>% filter(RUNX2 > 0), size = 0.6) +
  viridis::scale_color_viridis() +
  theme_classic(base_size = 10) +
  labs(color= 'RUNX2')
p1

p2 <- ggplot(df, aes(x= umap_1, y= umap_2 )) +
  geom_point(data = df %>% filter(BMPR2 == 0), color = "#440154FF", size = 0.6) +
  ggpointdensity::geom_pointdensity(data = df %>% filter(BMPR2 > 0), size = 0.6) +
  viridis::scale_color_viridis() +
  theme_classic(base_size = 10) +
  labs(color= 'BMPR2')
p2
pdf("01.hubgene_FetchData.pdf", width = 10, height = 4,family = "Times")
p1+p2
dev.off()
png("01.hubgene_FetchData.png", width = 10, height = 4,res = 600,units = 'in',family = "Times")
p1+p2
dev.off()

p1 <- FeaturePlot(UMAP, features = gene,
                  cols = c("#FFFFB3","#E41A1C"),ncol = 2)
pdf("01.FeaturePlot_umap.pdf", width = 8, height = 4,family = "Times")
p1
dev.off()
png("01.FeaturePlot_umap.png", width = 8, height = 4,res = 600,units = 'in',family = "Times")
p1
dev.off()

p2 <- VlnPlot(UMAP, features = gene,group.by = "celltype",pt.size = 0,
              cols = c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7") ,
              ncol =2)
p2
pdf("02.hubgene_violin.pdf", width = 7, height = 4,family = "Times")
p2
dev.off()
png("02.hubgene_violin.png", width = 7, height = 4,res = 600,units = 'in',family = "Times")
p2
dev.off()

pdf(file = '03.Prognostic_gene.pdf',w=8,h=6,family = "Times")
DotPlot(UMAP, features = gene,group.by = 'celltype',
        cols = c("#377EB8","#E41A1C")) + RotatedAxis()
dev.off()
png(file = '03.Prognostic_gene.png',w=800,h=600,family = "Times")
DotPlot(UMAP, features = gene,group.by = 'celltype',
        cols = c("#377EB8","#E41A1C")) + RotatedAxis()
dev.off()

library(reshape2)
library(tidyverse)
# hubgene <- data.frame(symbol=gene$symbol)
# hubgene <- data.frame(symbol=gene$symbol)
vln.df=as.data.frame(UMAP[["RNA"]]$data[gene,])
vln.df$symbol=rownames(vln.df)

vln.df=melt(vln.df,id="symbol")
colnames(vln.df)[c(2,3)]=c("sample","exp")
head(vln.df)
anno=UMAP@meta.data[,c("group","celltype")]
anno$sample <- rownames(anno)
vln.df=inner_join(vln.df,anno,by="sample")

library(rstatix)
stat.test<-vln.df%>%
  group_by(symbol)%>%
  wilcox_test(exp ~ group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.001,"***",
                    ifelse(stat.test$p<0.05,"**",
                           ifelse(stat.test$p<0.05,"*",'ns')))
write.csv(stat.test,file = 'stat.res.group.csv',row.names = F,quote = F)
de.stat <- stat.test[which(!stat.test$p=='ns'),]

library(ggpubr)
library(ggplot2)
library(ggsci)
table(UMAP$group)
vln.df$group <- factor(vln.df$group,levels = c('OA','Control'))
# vln.df <- vln.df[vln.df$symbol%in%de.stat$symbol,]
exp_plot <- ggplot(vln.df,aes(x = symbol, y = exp, color = group)) +
  geom_violin(trim=F,color="black",aes(fill=group)) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4",'seagreen'), name = "Group")+
  scale_color_manual(values = c("#4682B4",'seagreen'))+
  labs(title="", x="", y = "Expression level",size=20) +
  
  stat_compare_means(data = vln.df,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.35) +
  # stat_pvalue_manual(stat.test,
  #                    y.position = c(4,5,4,5,4,5),
  #                    size = 3.2,
  #                    family = "Times",
  #                    label = "p",
  #                    #parse = T,
  #                    face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+facet_wrap(~symbol,ncol = 5,scales = 'free')
exp_plot

pdf("04.expression(group).pdf", width = 7, height = 5,family = "Times")
exp_plot
dev.off()
png("04.expression(group).png", width = 7, height = 5,res = 600,units = 'in',family = "Times")
exp_plot
dev.off()

my_comparisons <- list(c('OA','Control'))
table(vln.df$celltype)
vln.dat <- vln.df
celltype <- as.data.frame(table(vln.dat$celltype))[1]
#celltype <- celltype[-10,] %>% as.data.frame()
colnames(celltype) <- 'Var1'
plots <- list()
for (i in c(1:nrow(celltype))) {
  plot.dat <- vln.dat[which(vln.dat$celltype==celltype$Var1[i]),]
  position=c(max(plot.dat$exp)+0.5,max(plot.dat$exp)+1,max(plot.dat$exp+1.5))
  
  p <- ggplot(plot.dat,aes(x=group,y=exp,fill=group))+
    geom_violin(trim = F,color='black')+
    stat_boxplot(geom = 'errorbar',
                 width=0.1,
                 position = position_dodge(0.9))+
    geom_boxplot(width=0.4,
                 position = position_dodge(0.9),
                 outlier.shape = NA,fill='white')+
    scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
    labs(title=celltype$Var1[i], x="", y = "",size=20) +
    geom_signif(comparisons = my_comparisons,
                test = t.test,
                map_signif_level = T,
                y_position = position,textsize = 2.5)+
    # stat_compare_means(data = train_phenotype3,
    #                    mapping = aes(group = stage),
    #                    label ="p",
    #                    method = 'kruskal.test',
    #                    paired = F,label.y = 4.5,label.x = 1.3) +
    ggtitle(plot.dat$celltype[1])+
    ylim(x=c(0,max(position+0.3))) +
    theme_bw()+
    theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=10),
          axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10),
          axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
          axis.title.x=element_text(size=16,face="bold"),
          axis.title.y=element_text(size=16,face="bold"),
          legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
          legend.title = element_text(face = "bold", size = 12),
          legend.position = "top",
          panel.grid = element_line(color = 'gray',size = 0.2),
          panel.border = element_rect(colour = 'black',fill = NA,size = 1))+
    facet_wrap(~symbol,ncol = 11,scales = 'free')
  p
  plots[[i]] <- p
}
p
library(gridExtra)
pdf(file = '05.gene.solo.pdf',w=15,h=8,family = "Times")
grid.arrange(grobs = plots, ncol = 6)
dev.off()

png(file = '05.gene.solo.png',w=15,h=8,units = 'in',res = 300,family = "Times")
grid.arrange(grobs = plots, ncol = 6)
dev.off()

table(UMAP$celltype)

UMAP.deg <- subset(UMAP,celltype %in% c('RegC','EC','preFC'))
table(UMAP.deg$celltype)
system.time(save(UMAP.deg, file = "UMAP.deg.Rdata")) 


#05_cellchat------
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./05_cellchat")){
  dir.create("./05_cellchat")
}
setwd("./05_cellchat")

load('../03_singleR/UMAP.Rdata')
table(UMAP$group)

if (! dir.exists("./OA")){dir.create("./OA")}
if (! dir.exists("./Control")){dir.create("./Control")}

UMAP1 <- subset(UMAP, group == 'OA')
UMAP2 <- subset(UMAP,group == 'Control')

library(CellChat)
library(Seurat)
library(ggplot2)
#Case-
#读入seurat处理后的rds文件
scRNA<-UMAP1
# DefaultAssay(scRNA) = 'RNA'
# scRNA <- JoinLayers(scRNA)  
#data.input = as.data.frame(GetAssayData(subset(scRNA), slot='counts'))# raw count
data.input<-scRNA@assays$RNA$data # normalized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$group)
# res <- strsplit(rownames(meta), "_")
# row<- sapply(res,"[",3)
# row <- data.frame(row)
# rownames(meta) <- row$row
# meta <- meta[!duplicated(meta$row),]
# rownames(meta) <- meta$row

cell.use = rownames(meta)
data.input<-scRNA@assays$RNA$data
#data.input@Dimnames[[2]] <- meta
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)

##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息9
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents<-droplevels(cellchat@idents,exclude=setdiff(levels(cellchat@idents),unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

##基于配受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的
showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  ###看一下CellChatDB的基本结构
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use
##表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat)#
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)
#df.net<-na.omit(df.net)
write.csv(df.net,'./OA/01.df.net.csv',row.names = F,quote = F)
##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat,file = './OA/cellchat.rds')
#cellchat <- readRDS('cellchat.rds')
## 相互作用数目
pdf('./OA/03.number_of_interactions.pdf',w=6,h=5,family='Times')
netVisual_heatmap(cellchat)
dev.off()
png('./OA/03.number_of_interactions.png',w=500,h=400,family='Times')
netVisual_heatmap(cellchat)
dev.off()

#### circle
pdf('./OA/01.net_number.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
png('./OA/01.net_number.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('./OA/02.net_weight.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
png('./OA/02.net_weight.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

#需要指定受体细胞和配体细胞
cellchat@data.signaling
pdf('./OA/03.buble.pdf',w=15,h=12,family='Times')
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
png('./OA/03.buble.png',w=15,h=12,units = 'in',res = 300,family='Times')
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()

mat <- cellchat@net$count
pdf('./OA/04.single_circle.pdf',w=18,h=8,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

png('./OA/04.single_circle.png',w=18,h=8,units='in',res=600,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

pdf('./OA/04.single_circle1.pdf',w=9,h=5,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,5,9)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

png('./OA/04.single_circle1.png',w=9,h=5,units='in',res=600,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,5,9)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()


mat <- cellchat@net$weight
pdf('./OA/04.single_circle_weight.pdf',w=18,h=8,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

png('./OA/04.single_circle_weight.png',w=18,h=8,units='in',res=600,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

pdf('./OA/04.single_circle1_weight.pdf',w=9,h=5,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,5,9)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

png('./OA/04.single_circle1_weight.png',w=9,h=5,units='in',res=600,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,5,9)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

#Control--
#读入seurat处理后的rds文件
scRNA<-UMAP2
# DefaultAssay(scRNA) = 'RNA'
# scRNA <- JoinLayers(scRNA)  
#data.input = as.data.frame(GetAssayData(subset(scRNA), slot='counts'))# raw count
data.input<-scRNA@assays$RNA$data # normalized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
table(meta$group)
# res <- strsplit(rownames(meta), "_")
# row<- sapply(res,"[",3)
# row <- data.frame(row)
# rownames(meta) <- row$row
# meta <- meta[!duplicated(meta$row),]
# rownames(meta) <- meta$row

cell.use = rownames(meta)
data.input<-scRNA@assays$RNA$data
#data.input@Dimnames[[2]] <- meta
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
table(meta$celltype)
unique(meta$celltype)

##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype",)
cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息9
cellchat <- setIdent(cellchat, ident.use = "celltype") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
unique(cellchat@idents)
cellchat@idents<-droplevels(cellchat@idents,exclude=setdiff(levels(cellchat@idents),unique(cellchat@idents)))

groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

##基于配受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的
showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  ###看一下CellChatDB的基本结构
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use
##表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat)#
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)
#df.net<-na.omit(df.net)
write.csv(df.net,'./Control/01.df.net.csv',row.names = F,quote = F)
##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat,file = './Control/cellchat.rds')
#cellchat <- readRDS('./Control/cellchat.rds')
## 相互作用数目
pdf('./Control/03.number_of_interactions.pdf',w=6,h=5,family='Times')
netVisual_heatmap(cellchat)
dev.off()
png('./Control/03.number_of_interactions.png',w=500,h=400,family='Times')
netVisual_heatmap(cellchat)
dev.off()

#### circle
pdf('./Control/01.net_number.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
png('./Control/01.net_number.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('./Control/02.net_weight.pdf',w=7,h=9,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()
png('./Control/02.net_weight.png',w=7,h=9,units='in',res=600,family='Times')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")
dev.off()

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

#需要指定受体细胞和配体细胞
cellchat@data.signaling
pdf('./Control/03.buble.pdf',w=10,h=15,family='Times')
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()
png('./Control/03.buble.png',w=10,h=15,units = 'in',res = 300,family='Times')
netVisual_bubble(cellchat, remove.isolate = FALSE)
dev.off()

mat <- cellchat@net$count
pdf('./Control/04.single_circle.pdf',w=18,h=8,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()
png('./Control/04.single_circle.png',w=18,h=8,units='in',res=600,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

pdf('./Control/04.single_circle1.pdf',w=9,h=5,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,4,8)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

png('./Control/04.single_circle1.png',w=9,h=5,units='in',res=600,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,4,8)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

mat <- cellchat@net$weight
pdf('./Control/04.single_circle_weight.pdf',w=18,h=8,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

png('./Control/04.single_circle_weight.png',w=18,h=8,units='in',res=600,family='Times')
par(mfrow=c(2,6),xpd=T)
#par(mar=c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

pdf('./Control/04.single_circle1_weight.pdf',w=9,h=5,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,4,8)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()

png('./Control/04.single_circle1_weight.png',w=9,h=5,units='in',res=600,family='Times')
par(mfrow=c(1,3),xpd=T)
for (i in c(3,4,8)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()


#06_function--------------------------
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./06_function")){
  dir.create("./06_function")
}
setwd("./06_function")

load('../04_hubgene/UMAP.deg.Rdata')
UMAP <- UMAP.deg
table(UMAP$celltype)
library(ReactomeGSA)
#DefaultAssay(sce.car) = 'RNA'
#sce.car <- JoinLayers(sce.car)
gsva_result <- analyse_sc_clusters(UMAP, verbose = TRUE)
# colnames(gsva_result@results$Seurat$fold_changes) <- gsub("X", "Cluster ", colnames(gsva_result@results$Seurat$fold_changes))
head(gsva_result@results$Seurat$fold_changes)

pathway_expression <- pathways(gsva_result)
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
colnames(pathway_expression) <- gsub("X", "Cluster ", colnames(pathway_expression))
pathway_expression[1:3,]
# find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min
# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]
head(max_difference)
# library(tidyverse)
# p1 <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])
# p1$data %>% mutate(absmy = ifelse(expr>=0, "Z","Fy")) -> df
# df %>% ggplot(aes(cluster_id,   expr ,fill=absmy    ))+
#   geom_bar(stat='identity') + theme_bw()+
#   theme(axis.text.x = element_text(angle =45,hjust = .9,size = 10,vjust = 0.9))+
#   ggtitle("Alanine metabolism") + theme(legend.position="none")

#Additional parameters are directly passed to gplots heatmap.2 function
col2 <- colorRampPalette(c("#4682B4","white" ,"#CD3700"),alpha = TRUE)
pdf('01.deg_cells.heat.pdf',w=10,h=6,family = 'Times')
plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(7,18),truncate_names = T,col = col2)
dev.off()
png('01.deg_cells.heat.png',w=10,h=6,units = 'in',res = 600,family = 'Times')
plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(7,18),truncate_names = T,col = col2)
dev.off()

write.csv(max_difference,file = 'deg_cells.path.csv',row.names = F)







#07_subtype----------------
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./07_subtype")){dir.create("./07_subtype")}
setwd("./07_subtype")

library(Seurat)
library(dplyr)
library(tidyverse)
library(data.table)
library(readr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggsci)
library(scales)
library(SingleR)
#library(celldex)
library(grDevices)
library(IOBR)
library(RColorBrewer)
library(ComplexHeatmap)
library(rstatix)
library(ggpubr)
#install.packages('vctrs')

load('../04_hubgene/UMAP.deg.Rdata')
mycolor <-  c("#f9766e","#fbbab6","#e1c548","#f0e2a3","#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#4e79a6","#bac4d0","#45337f","#a199be","#aedd2f","#d7ee96",
              "#dc8e97","#e3d1db","#74a893","#ac9141","#5ac6e9","#ebce8e","#e5c06e","#7587b1","#c7deef","#e97371","#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac",
              "#e0bc58","#64abc0","#fab37f","#e98741","#8fc0dc","#967568","#f2d3ca","#eebd85","#82c785","#edeaa4","#cdaa9f","#794976","#bcacd3","#889b5d","#4e9592","#dbad5f","#64ae79","#ac5092"
              
)
scRNA <- UMAP.deg
table(scRNA$celltype)
celllist <- c('RegC','EC','preFC')
pcSelect=30

theme.set <- theme(
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16,  family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  title=element_text(size=24, color='black',hjust=0.5, face = "bold",family='Times'),
  text = element_text(family = "Times"))

j=1
for (j in 1:length(celllist)){
  cellname <- celllist[j]
  print(cellname)
  tryCatch({
    scRNA1 <- scRNA[ ,scRNA$celltype %in% c(cellname)]
    scRNA1 <- NormalizeData(scRNA1, normalization.method = "LogNormalize", scale.factor = 10000)
    scRNA1 <- FindVariableFeatures(scRNA1, selection.method = 'vst', nfeatures = 1000)
    scRNA1 <- ScaleData(scRNA1)
    scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(object = scRNA1))
    scRNA1 <- JackStraw(scRNA1, num.replicate = 100, dims = 50)
    scRNA1 <- ScoreJackStraw(scRNA1, dims = 1:50)
    
    p1 <- JackStrawPlot(scRNA1, dims = 1:50)
    pdf(paste0('01.subtype_',cellname,'_JackStrawPlot.pdf'),w=10,h=5,family='Times')
    print(p1)
    dev.off()
    png(paste0('01.subtype_',cellname,'_JackStrawPlot.png'),w=10,h=5,units = 'in',res = 300,family='Times')
    print(p1)
    dev.off()
    
    
    p2 <- ElbowPlot(scRNA1, ndims = 50)
    pdf(paste0('02.subtype_',cellname,'_ElbowPlot.pdf'),w=6,h=5,family='Times')
    print(p2)
    dev.off()
    png(paste0('02.subtype_',cellname,'_ElbowPlot.png'),w=6,h=5,units = 'in',res = 300,family='Times')
    print(p2)
    dev.off()
    
    system.time(save(scRNA1, file = paste0("scRNA.raw.",cellname,".Rdata")))
    #pcSelect = 19
    scRNA1 <- FindNeighbors(scRNA1, dims = 1:pcSelect)
    scRNA1 <- FindClusters(scRNA1, resolution = 0.4)
    scRNA1 <- RunUMAP(scRNA1,dims = 1:pcSelect)
    scRNA1 <- RunTSNE(scRNA1,dims = 1:pcSelect)
    
    p3 <- DimPlot(scRNA1, reduction = "tsne", group.by = "seurat_clusters", cols = mycolor,label = T)+
      ggtitle(paste0('Subtype of ',cellname))  + xlab('TSNE_1') + ylab('TSNE_2')
    pdf(paste0('03.subtype_',cellname,'_TSNE.pdf'),w=7,h=6,family='Times')
    print(p3)
    dev.off()
    png(paste0('03.subtype_',cellname,'_TSNE.png'),w=7,h=6,units = 'in',res = 300,family='Times')
    print(p3)
    dev.off()
    
    p4 <- DimPlot(scRNA1, reduction = "umap", group.by = "seurat_clusters", cols = mycolor,label = T)+
      ggtitle(paste0('Subtype of ',cellname))
    pdf(paste0('03.subtype_',cellname,'_UMAP.pdf'),w=7,h=6,family='Times')
    print(p4)
    dev.off()
    png(paste0('03.subtype_',cellname,'_UMAP.png'),w=7,h=6,units = 'in',res = 300,family='Times')
    print(p4)
    dev.off()
    
    system.time(save(scRNA1, file = paste0("scRNA.",cellname,".Rdata")))
    
  }, error = function(e){
    print(paste(j,'error'))
  },finally = {
    print("----------------------------------------------------------")
  })
}



#08_Trajectory---------------
rm(list = ls())
setwd("//12_Single_cell/")
if (! dir.exists("./08_Trajectory")){
  dir.create("./08_Trajectory")
}
setwd("./08_Trajectory")

#'RegC','EC','preFC'
# if (! dir.exists("./RegC")){dir.create("./RegC")}
# if (! dir.exists("./EC")){dir.create("./EC")}
# if (! dir.exists("./preFC")){dir.create("./preFC")}

#devtools::load_all('/data/nas1/luchunlin/project/YQQL0106-11/monocle')
library(monocle)
library(Seurat)
library(tidyverse)

##关键细胞1------------
load('../07_subtype/scRNA.RegC.Rdata')
Mono_tj = scRNA1
Mono_matrix = GetAssayData(Mono_tj, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann) = rownames(Mono_matrix)
Mono_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = Mono_tj@meta.data
cell.type = Idents(Mono_tj) %>% data.frame
sample_ann = cbind(sample_ann, cell.type)
colnames(sample_ann)[ncol(sample_ann)] = "Subtype"
Mono_pd = new("AnnotatedDataFrame", data =sample_ann)
Mono.cds = newCellDataSet(Mono_matrix, phenoData =Mono_pd, featureData=Mono_fd, expressionFamily=negbinomial.size())
rm(Mono_tj)
Mono.cds = estimateSizeFactors(Mono.cds)
Mono.cds = estimateDispersions(Mono.cds)
WGCNA::collectGarbage()
Mono.cds = reduceDimension(Mono.cds, max_components = 2, verbose = T, 
                           norm_method = "log"
)
save(Mono.cds, file = "./RegC/mycds_order_raw.RData")
Mono.cds = orderCells(Mono.cds)
save(Mono.cds, file = "./RegC/mycds_order.RData")
my_cds_subset=Mono.cds
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数
saveRDS(my_pseudotime_de,file = './RegC/my_pseudotime.de.rds')


load('./RegC/mycds_order.RData')
table(Mono.cds$State)
Mono.cds$seurat_clusters = factor(Mono.cds$seurat_clusters)
table(Mono.cds$seurat_clusters)
table(Mono.cds$celltype)
#Mono.cds = orderCells(Mono.cds,root_state = 4)
p = plot_cell_trajectory(Mono.cds,color_by="Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  #  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p
pdf('./RegC/01.PSD_SMC.pdf',w=6,h=6,family='Times')
print(p)
dev.off()
png('./RegC/01.PSD_SMC.png',w=6,h=6,units = 'in',res = 300,family='Times')
print(p)
dev.off()

p3 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))
p3
pdf('./RegC/03.State_Trajectory.pdf',w=6,h=6,family='Times')
print(p3)
dev.off()
png('./RegC/03.State_Trajectory.png',w=6,h=6,units = 'in',res = 300,family='Times')
print(p3)
dev.off()

p4 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))+
  facet_wrap(~group,nrow = 1)
p4
pdf('./RegC/04.State_Trajectory(group).pdf',w=12,h=6,family='Times')
print(p4)
dev.off()
png('./RegC/04.State_Trajectory(group).png',w=12,h=6,units = 'in',res = 300,family='Times')
print(p4)
dev.off()


p5 = plot_cell_trajectory(Mono.cds,color_by="seurat_clusters",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen',"#e5ce81","#f47720","#459943","#bdc3d2","#606f8a","#e8c559","#ea9c9d","#005496","#93cc82","#4d97cd","#ea9c9d")) +
  theme(text = element_text(size = 14))
#facet_wrap(~group,nrow = 1)
p5
pdf('./RegC/05.seurat_clusters_Trajectory.pdf',w=6,h=7,family='Times')
print(p5)
dev.off()
png('./RegC/05.seurat_clusters_Trajectory.png',w=6,h=7,units = 'in',res = 300,family='Times')
print(p5)
dev.off()


##hub基因表达-
hubgene <- read.csv('../../04_SMR/results_eqtl/SMR_eqtl.csv')
hubgene <- hubgene$SYMBOL

p4 <- plot_genes_in_pseudotime(Mono.cds[c(hubgene),],
                               color_by = "State",
                               ncol = 2)+
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
p4
pdf('./RegC/06.sdheat_hubgene.pdf',w=8,h=4,family='Times')
print(p4)
dev.off()
png('./RegC/06.sdheat_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p4)
dev.off()
# 横轴：拟时序（Pseudotime），这是一个根据细胞转录组数据的相似性计算出的、
# 代表细胞发育或分化进程的虚拟时间轴。
# 纵轴：基因表达量，通常以标准化后的表达值（如z-score或log转换值）表示。
# 图形解释：这个函数会绘制一个或多个基因的表达量随拟时序变化的曲线。
# 不同的颜色代表不同的细胞状态。通过观察这些曲线，你可以了解基因在不同发育阶段或细胞状态下的
# 表达模式。

p444 <- plot_genes_violin(Mono.cds[c(hubgene),],
                          color_by = "State",
                          ncol = 2)+
  scale_fill_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))+theme_bw()
p444
pdf('./RegC/07.violin_hubgene.pdf',w=8,h=4,family='Times')
print(p444)
dev.off()
png('./RegC/07.violin_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p444)
dev.off()
# 横轴：细胞状态（State），这是根据细胞的转录组数据聚类得到的分类。
# 纵轴：基因表达量，通常以实际或标准化的表达值表示，但在这个图中，
# 表达量的分布是以小提琴图的形式展示的。
# 图形解释：这个函数会为每个基因在每个细胞状态下绘制一个小提琴图，
# 展示该状态下所有细胞的基因表达分布。不同的颜色代表不同的状态。通过观察小提琴图的形状和位置，
# 你可以了解基因在不同状态下的表达分布和差异。

p555 <- plot_genes_jitter(Mono.cds[c(hubgene),],
                          color_by = "State",
                          ncol = 2)+
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
p555
pdf('./RegC/08.jitter_hubgene.pdf',w=8,h=4,family='Times')
print(p555)
dev.off()
png('./RegC/08.jitter_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p555)
dev.off()
# 横轴：细胞状态（State），与plot_genes_violin相同。
# 纵轴：基因表达量，通常以实际或标准化的表达值表示。
# 图形解释：这个函数会在每个状态下为每个基因绘制一个点图，
# 每个点代表一个细胞的基因表达量。为了使点图更清晰，
# 表达量值通常会添加一些随机扰动（jittering）。不同的颜色代表不同的状态。
# 通过观察这些点的分布和颜色，你可以了解基因在不同状态下的表达水平和变化。



##动力热图--
my_cds_subset=Mono.cds
# 拟时序数据和细胞位置在pData 中
# head(pData(my_cds_subset))
# # 这个differentialGeneTest会比较耗费时间，测试每个基因的拟时序表达
# my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数
# saveRDS(my_pseudotime_de,file = 'my_pseudotime.de(npc).rds')

my_pseudotime_de <- readRDS('RegC/my_pseudotime.de.rds')
head(my_pseudotime_de)
library(tidyverse)

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[hubgene,],# num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
pdf('./RegC/09.pseudotime_hubgene.pdf',w=4,h=3,family='Times')
print(my_pseudotime_cluster)
dev.off()
png('./RegC/09.pseudotime_hubgene.png',w=4,h=3,units = 'in',res = 300,family='Times')
print(my_pseudotime_cluster)
dev.off()



##关键细胞2--------------------
load('../07_subtype/scRNA.EC.Rdata')
Mono_tj = scRNA1
Mono_matrix = GetAssayData(Mono_tj, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann) = rownames(Mono_matrix)
Mono_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = Mono_tj@meta.data
cell.type = Idents(Mono_tj) %>% data.frame
sample_ann = cbind(sample_ann, cell.type)
colnames(sample_ann)[ncol(sample_ann)] = "Subtype"
Mono_pd = new("AnnotatedDataFrame", data =sample_ann)
Mono.cds = newCellDataSet(Mono_matrix, phenoData =Mono_pd, featureData=Mono_fd, expressionFamily=negbinomial.size())
rm(Mono_tj)
Mono.cds = estimateSizeFactors(Mono.cds)
Mono.cds = estimateDispersions(Mono.cds)
WGCNA::collectGarbage()
Mono.cds = reduceDimension(Mono.cds, max_components = 2, verbose = T, 
                           norm_method = "log"
)
save(Mono.cds, file = "./EC/mycds_order_raw.RData")
Mono.cds = orderCells(Mono.cds)
save(Mono.cds, file = "./EC/mycds_order.RData")
my_cds_subset=Mono.cds
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数
saveRDS(my_pseudotime_de,file = './EC/my_pseudotime.de.rds')

load('./EC/mycds_order.RData')
table(Mono.cds$State)
Mono.cds$seurat_clusters = factor(Mono.cds$seurat_clusters)
table(Mono.cds$seurat_clusters)
table(Mono.cds$celltype)
#Mono.cds = orderCells(Mono.cds,root_state = 4)
p = plot_cell_trajectory(Mono.cds,color_by="Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  #  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p
pdf('./EC/01.PSD_SMC.pdf',w=6,h=6,family='Times')
print(p)
dev.off()
png('./EC/01.PSD_SMC.png',w=6,h=6,units = 'in',res = 300,family='Times')
print(p)
dev.off()

p3 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))
p3
pdf('./EC/03.State_Trajectory.pdf',w=6,h=6,family='Times')
print(p3)
dev.off()
png('./EC/03.State_Trajectory.png',w=6,h=6,units = 'in',res = 300,family='Times')
print(p3)
dev.off()

p4 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))+
  facet_wrap(~group,nrow = 1)
p4
pdf('./EC/04.State_Trajectory(group).pdf',w=12,h=6,family='Times')
print(p4)
dev.off()
png('./EC/04.State_Trajectory(group).png',w=12,h=6,units = 'in',res = 300,family='Times')
print(p4)
dev.off()


p5 = plot_cell_trajectory(Mono.cds,color_by="seurat_clusters",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen',"#e5ce81","#f47720","#459943","#bdc3d2","#606f8a","#e8c559","#ea9c9d","#005496","#93cc82","#4d97cd","#ea9c9d")) +
  theme(text = element_text(size = 14))
#facet_wrap(~group,nrow = 1)
p5
pdf('./EC/05.seurat_clusters_Trajectory.pdf',w=6,h=7,family='Times')
print(p5)
dev.off()
png('./EC/05.seurat_clusters_Trajectory.png',w=6,h=7,units = 'in',res = 300,family='Times')
print(p5)
dev.off()



##hub基因表达-
hubgene <- read.csv('../../04_SMR/results_eqtl/SMR_eqtl.csv')
hubgene <- hubgene$SYMBOL
p4 <- plot_genes_in_pseudotime(Mono.cds[c(hubgene),],
                               color_by = "State",
                               ncol = 2)+
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
p4
pdf('./EC/06.sdheat_hubgene.pdf',w=8,h=4,family='Times')
print(p4)
dev.off()
png('./EC/06.sdheat_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p4)
dev.off()
# 横轴：拟时序（Pseudotime），这是一个根据细胞转录组数据的相似性计算出的、
# 代表细胞发育或分化进程的虚拟时间轴。
# 纵轴：基因表达量，通常以标准化后的表达值（如z-score或log转换值）表示。
# 图形解释：这个函数会绘制一个或多个基因的表达量随拟时序变化的曲线。
# 不同的颜色代表不同的细胞状态。通过观察这些曲线，你可以了解基因在不同发育阶段或细胞状态下的
# 表达模式。

p444 <- plot_genes_violin(Mono.cds[c(hubgene),],
                          color_by = "State",
                          ncol = 2)+
  scale_fill_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))+theme_bw()
p444
pdf('./EC/07.violin_hubgene.pdf',w=8,h=4,family='Times')
print(p444)
dev.off()
png('./EC/07.violin_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p444)
dev.off()
# 横轴：细胞状态（State），这是根据细胞的转录组数据聚类得到的分类。
# 纵轴：基因表达量，通常以实际或标准化的表达值表示，但在这个图中，
# 表达量的分布是以小提琴图的形式展示的。
# 图形解释：这个函数会为每个基因在每个细胞状态下绘制一个小提琴图，
# 展示该状态下所有细胞的基因表达分布。不同的颜色代表不同的状态。通过观察小提琴图的形状和位置，
# 你可以了解基因在不同状态下的表达分布和差异。

p555 <- plot_genes_jitter(Mono.cds[c(hubgene),],
                          color_by = "State",
                          ncol = 2)+
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
p555
pdf('./EC/08.jitter_hubgene.pdf',w=8,h=4,family='Times')
print(p555)
dev.off()
png('./EC/08.jitter_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p555)
dev.off()
# 横轴：细胞状态（State），与plot_genes_violin相同。
# 纵轴：基因表达量，通常以实际或标准化的表达值表示。
# 图形解释：这个函数会在每个状态下为每个基因绘制一个点图，
# 每个点代表一个细胞的基因表达量。为了使点图更清晰，
# 表达量值通常会添加一些随机扰动（jittering）。不同的颜色代表不同的状态。
# 通过观察这些点的分布和颜色，你可以了解基因在不同状态下的表达水平和变化。



##动力热图--
my_cds_subset=Mono.cds
# 拟时序数据和细胞位置在pData 中
# head(pData(my_cds_subset))
# # 这个differentialGeneTest会比较耗费时间，测试每个基因的拟时序表达
# my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数
# saveRDS(my_pseudotime_de,file = 'my_pseudotime.de(npc).rds')

my_pseudotime_de <- readRDS('EC/my_pseudotime.de.rds')
head(my_pseudotime_de)
library(tidyverse)

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[hubgene,],# num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
pdf('./EC/09.pseudotime_hubgene.pdf',w=4,h=3,family='Times')
print(my_pseudotime_cluster)
dev.off()
png('./EC/09.pseudotime_hubgene.png',w=4,h=3,units = 'in',res = 300,family='Times')
print(my_pseudotime_cluster)
dev.off()



##关键细胞3--------------------
load('../07_subtype/scRNA.preFC.Rdata')
Mono_tj = scRNA1
Mono_matrix = GetAssayData(Mono_tj, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann) = rownames(Mono_matrix)
Mono_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = Mono_tj@meta.data
cell.type = Idents(Mono_tj) %>% data.frame
sample_ann = cbind(sample_ann, cell.type)
colnames(sample_ann)[ncol(sample_ann)] = "Subtype"
Mono_pd = new("AnnotatedDataFrame", data =sample_ann)
Mono.cds = newCellDataSet(Mono_matrix, phenoData =Mono_pd, featureData=Mono_fd, expressionFamily=negbinomial.size())
rm(Mono_tj)
Mono.cds = estimateSizeFactors(Mono.cds)
Mono.cds = estimateDispersions(Mono.cds)
WGCNA::collectGarbage()
Mono.cds = reduceDimension(Mono.cds, max_components = 2, verbose = T, 
                           norm_method = "log"
)
save(Mono.cds, file = "./preFC/mycds_order_raw.RData")
Mono.cds = orderCells(Mono.cds)
save(Mono.cds, file = "./preFC/mycds_order.RData")
my_cds_subset=Mono.cds
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数
saveRDS(my_pseudotime_de,file = './preFC/my_pseudotime.de.rds')

load('./preFC/mycds_order.RData')
table(Mono.cds$State)
Mono.cds$seurat_clusters = factor(Mono.cds$seurat_clusters)
table(Mono.cds$seurat_clusters)
table(Mono.cds$celltype)
#Mono.cds = orderCells(Mono.cds,root_state = 4)
p = plot_cell_trajectory(Mono.cds,color_by="Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  #  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p
pdf('./preFC/01.PSD_SMC.pdf',w=6,h=6,family='Times')
print(p)
dev.off()
png('./preFC/01.PSD_SMC.png',w=6,h=6,units = 'in',res = 300,family='Times')
print(p)
dev.off()

p3 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))
p3
pdf('./preFC/03.State_Trajectory.pdf',w=6,h=6,family='Times')
print(p3)
dev.off()
png('./preFC/03.State_Trajectory.png',w=6,h=6,units = 'in',res = 300,family='Times')
print(p3)
dev.off()

p4 = plot_cell_trajectory(Mono.cds,color_by="State",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(text = element_text(size = 14))+
  facet_wrap(~group,nrow = 1)
p4
pdf('./preFC/04.State_Trajectory(group).pdf',w=12,h=6,family='Times')
print(p4)
dev.off()
png('./preFC/04.State_Trajectory(group).png',w=12,h=6,units = 'in',res = 300,family='Times')
print(p4)
dev.off()

p5 = plot_cell_trajectory(Mono.cds,color_by="seurat_clusters",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen',"#e5ce81","#f47720","#459943","#bdc3d2","#606f8a","#e8c559","#ea9c9d","#005496","#93cc82","#4d97cd","#ea9c9d")) +
  theme(text = element_text(size = 14))
#facet_wrap(~group,nrow = 1)
p5
pdf('./preFC/05.seurat_clusters_Trajectory.pdf',w=6,h=7,family='Times')
print(p5)
dev.off()
png('./preFC/05.seurat_clusters_Trajectory.png',w=6,h=7,units = 'in',res = 300,family='Times')
print(p5)
dev.off()


##hub基因表达-
hubgene <- read.csv('../../04_SMR/results_eqtl/SMR_eqtl.csv')
hubgene <- hubgene$SYMBOL

p4 <- plot_genes_in_pseudotime(Mono.cds[c(hubgene),],
                               color_by = "State",
                               ncol = 2)+
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
p4
pdf('./preFC/06.sdheat_hubgene.pdf',w=8,h=4,family='Times')
print(p4)
dev.off()
png('./preFC/06.sdheat_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p4)
dev.off()
# 横轴：拟时序（Pseudotime），这是一个根据细胞转录组数据的相似性计算出的、
# 代表细胞发育或分化进程的虚拟时间轴。
# 纵轴：基因表达量，通常以标准化后的表达值（如z-score或log转换值）表示。
# 图形解释：这个函数会绘制一个或多个基因的表达量随拟时序变化的曲线。
# 不同的颜色代表不同的细胞状态。通过观察这些曲线，你可以了解基因在不同发育阶段或细胞状态下的
# 表达模式。

p444 <- plot_genes_violin(Mono.cds[c(hubgene),],
                          color_by = "State",
                          ncol = 2)+
  scale_fill_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))+theme_bw()
p444
pdf('./preFC/07.violin_hubgene.pdf',w=8,h=4,family='Times')
print(p444)
dev.off()
png('./preFC/07.violin_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p444)
dev.off()
# 横轴：细胞状态（State），这是根据细胞的转录组数据聚类得到的分类。
# 纵轴：基因表达量，通常以实际或标准化的表达值表示，但在这个图中，
# 表达量的分布是以小提琴图的形式展示的。
# 图形解释：这个函数会为每个基因在每个细胞状态下绘制一个小提琴图，
# 展示该状态下所有细胞的基因表达分布。不同的颜色代表不同的状态。通过观察小提琴图的形状和位置，
# 你可以了解基因在不同状态下的表达分布和差异。

p555 <- plot_genes_jitter(Mono.cds[c(hubgene),],
                          color_by = "State",
                          ncol = 2)+
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple',"firebrick",'skyblue','pink','cyan','lightseagreen')) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
p555
pdf('./preFC/08.jitter_hubgene.pdf',w=8,h=4,family='Times')
print(p555)
dev.off()
png('./preFC/08.jitter_hubgene.png',w=8,h=4,units = 'in',res = 300,family='Times')
print(p555)
dev.off()
# 横轴：细胞状态（State），与plot_genes_violin相同。
# 纵轴：基因表达量，通常以实际或标准化的表达值表示。
# 图形解释：这个函数会在每个状态下为每个基因绘制一个点图，
# 每个点代表一个细胞的基因表达量。为了使点图更清晰，
# 表达量值通常会添加一些随机扰动（jittering）。不同的颜色代表不同的状态。
# 通过观察这些点的分布和颜色，你可以了解基因在不同状态下的表达水平和变化。




my_cds_subset=Mono.cds

my_pseudotime_de <- readRDS('./preFC/my_pseudotime.de.rds')
head(my_pseudotime_de)
library(tidyverse)

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[hubgene,],# num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)
pdf('./preFC/09.pseudotime_hubgene.pdf',w=4,h=3,family='Times')
print(my_pseudotime_cluster)
dev.off()
png('./preFC/09.pseudotime_hubgene.png',w=4,h=3,units = 'in',res = 300,family='Times')
print(my_pseudotime_cluster)
dev.off()

