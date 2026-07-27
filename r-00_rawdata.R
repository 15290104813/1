rm(list = ls())
setwd('//')
if (! dir.exists('./00_rawdata')) {
  dir.create('./00_rawdata')
}
setwd("./00_rawdata")

library(readr)
library(readxl)
library(tidyverse)
library(tibble)
library(dplyr) 
library(tidyr)
library(tinyarray)
library(Biobase)
library(AnnoProbe) #idmap
library(clusterProfiler)
library(GEOquery)

#GSE51588-----
gset<-getGEO("GSE51588",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
gset$GSE51588_series_matrix.txt.gz@annotation
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
gpl <- idmap(gpl = 'GPL13497')
probe2symobl<-gpl
colnames(probe2symobl)<-c('ID','symbol')
dat.va<-expr
dat.va$ID<-rownames(dat.va)
dat.va$ID<-as.character(dat.va$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat.va<-dat.va %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%   
  dplyr::select(symbol,everything())%>%   
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%  
  arrange(desc(rowMean))%>%     
  distinct(symbol,.keep_all = T)%>%    
  dplyr::select(-rowMean)%>%    
  tibble::column_to_rownames(colnames(.)[1])   
pd<-pData(a)
table(pd$characteristics_ch1)
group<-data.frame(sample=pd$geo_accession,group=pd$characteristics_ch1)
table(group$group)
group$group <- ifelse(group$group == "disease state: Normal",'Control','OA')
group<-group[order(group$group),]       
dat.va <- dat.va[,group$sample]

write.csv(dat.va,file = 'dat(GSE51588).csv')
write.csv(group,file = 'group(GSE51588).csv',row.names = F)




library(clusterProfiler)
library(plinkbinr)
library(VariantAnnotation)
library(gwasglue)
library(data.table)
library(tidyr)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
# ###tsv outcome
library(tidyverse)
library(tinyarray)

outcome<-read_tsv('//00_rawdata/GCST90038686_buildGRCh37.tsv')
outcome$n_cases <- 39515
outcome$n_controls <- 445083
outcome$variant_id[1:50]
# #如果SNP没有信息，可以转换。根据chr:pos查找RS号  https://blog.csdn.net/wendy_milk/article/details/121642843
# outcome$`chromosome:start` <- paste(outcome$chromosome,outcome$base_pair_location,sep=':')
# #match <- fread("/data2/marui/pipeline/MR/snp150_hg19.txt.gz",header=T,check.names=F,sep="\t")
# #system.time(save(match, file = "/data2/marui/pipeline/MR/match.Rdata"))
# load('/data2/marui/pipeline/MR/match.Rdata')
# need <- dplyr::left_join(outcome,match,by="chromosome:start")
# ##使用format_data()函数将该数据框转化成TwoSampleMR的格式
# save(need, file='need.RData')
# #load('need.RData')
# colnames(need)
# outcome_dat <- format_data(
#   dat=need,
#   type = "outcome",
#   snps = NULL,
#   snp_col = "name",
#   beta_col = "beta",
#   se_col = "standard_error",
#   eaf_col = "effect_allele_frequency",
#   effect_allele_col = "effect_allele",
#   other_allele_col = "other_allele",
#   pval_col = "p_value",
#   ncase_col = "n_cases",
#   ncontrol_col = "n_controls",
#   chr_col = "chromosome",
#   pos_col = "base_pair_location"
# )

colnames(outcome)
outcome_all <- format_data(
  dat=outcome,
  type = "outcome",
  snps = NULL,
  header = TRUE,
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  ncase_col = "n_cases",
  ncontrol_col = "n_controls",
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)
#outcome_all <- outcome_dat
outcome_all$id.outcome <- "Osteoarthritis"
outcome_all$outcome <- "Osteoarthritis || id:ebi-a-GCST90038686"
system.time(save(outcome_all, file = "outcome_all.Rdata"))

head(outcome_all)
library(data.table)
df=outcome_all
head(df)
df1=data.frame(SNP=df$SNP,  #SNP
               A1=df$effect_allele.outcome,     #效应等位基因
               A2=df$other_allele.outcome,     #参考等位基因
               freq=df$eaf.outcome,#效应等位基因发生频率
               b=df$beta.outcome,     #beta值
               se=df$se.outcome,  #标准误
               p=df$pval.outcome)     #p值
df1$n=39515+445083                  #样本量(疾病+对照)
head(df1)
df1 <- df1 %>% as_tibble() %>% separate_rows(SNP, sep = ",")
df1 <- df1[df1$SNP!="",]
df1 <- df1[!duplicated(df1$SNP),]
fwrite(df1,"00.outcome.ma",sep = "\t",quote = F,row.names = F)
