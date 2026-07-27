rm(list = ls())
setwd("//")
if(!dir.exists("06_Locater_eqtl")){dir.create("06_Locater_eqtl")}
setwd("06_Locater_eqtl")

library(plyr)
library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(VariantAnnotation)
library(gwasglue)
library(S4Vectors)
library(coloc)
# library(remotes)
# install_github("chr1swallace/coloc",build_vignettes=TRUE)
library(locuscomparer)
# devtools::install_github("boxiangliu/locuscomparer")

load('../00_rawdata/outcome_all.Rdata')
head(outcome_all)
outcome <- outcome_all
GWAS <- outcome


hub2 <- read.csv('../04_SMR/results_eqtl/SMR_eqtl.csv')


path = '//MR/eQTL_vcf/'

info <- hub2
eqtlID<-info$probeID
GWASID<-'Osteoarthritis'     
outcome_all=484598    #样本总数
outcome_case=39515    #疾病样本数
result<-data.frame()
for (i in c(1:length(eqtlID))) {
  vcf <- VariantAnnotation::readVcf(paste0(path,'eqtl-a-',eqtlID[i],'.vcf.gz'), "hg19")
  eqtl <- gwasvcf_to_TwoSampleMR(vcf)
  eqtl$maf<-ifelse(eqtl$eaf.exposure < 0.5, 
                   eqtl$eaf.exposure,
                   1 - eqtl$eaf.exposure)
  input <- merge(eqtl, GWAS, by.x="SNP", by.y="SNP", all=FALSE, suffixes=c("_eqtl","_GWAS"))
  input <- input %>% filter(SNP != ""& !duplicated(SNP) & SNP != ".")
  colnames(input)
  input <- subset(input, beta.outcome!= 'NA')
  input <- subset(input,effect_allele.outcome!= 'NA')
  
  input <- subset(input, maf!= 'NA')
  
  input <- subset(input,SNP == hub2$topSNP[i])
  
  input$VAR.exposure = input$se.exposure^2
  input$VAR.outcome = input$se.outcome^2
  
  TEST1 <- try(res <- coloc.abf(dataset1=list(snp=input$SNP, 
                                              pvalues=input$pval.exposure, 
                                              type="quant", 
                                              beta=input$beta.exposure, 
                                              varbeta=input$VAR.exposure,
                                              N=input$samplesize.exposure, 
                                              MAF=input$maf),
                                
                                dataset2=list(snp=input$SNP, 
                                              pvalues=input$pval.outcome, 
                                              type="cc", 
                                              beta=input$beta.outcome, 
                                              varbeta=input$VAR.outcome,
                                              s=outcome_case/outcome_all),
                                p1=1e-4,p2=1e-4,p12=1e-5))
  if ('try-error' %in% class(TEST1)){
    next
  }
  
  res_coloc<-data.frame(nsnp=res[["summary"]][["nsnps"]],
                        PPH0=res[["summary"]][["PP.H0.abf"]],
                        PPH1=res[["summary"]][["PP.H1.abf"]],
                        PPH2=res[["summary"]][["PP.H2.abf"]],
                        PPH3=res[["summary"]][["PP.H3.abf"]],
                        PPH4=res[["summary"]][["PP.H4.abf"]],
                        eqtl=eqtlID[i],
                        GWAS=GWASID)
  result<-rbind(res_coloc,result)
  
  gwas_fn <- input[,c('SNP','pval.outcome')] 
  colnames(gwas_fn)<-c('rsid','pval') 
  
  eqtl_fn <-input[c('SNP','pval.exposure')] 
  colnames(eqtl_fn)<-c('rsid','pval') 
  
  pdf(paste0(eqtlID[i],'-',GWASID,'.pdf'), width=10, height=6,family = "Times")
  TEST2<-try(print(locuscompare(in_fn1= gwas_fn, in_fn2= eqtl_fn, title1= GWASID, title2 =eqtlID[i])))
  if('try-error'%in% class(TEST2)) {
    dev.off()
    next}
  dev.off()

  png(paste0(eqtlID[i],'-',GWASID,'.png'), width=10, height=6,res = 600,units = 'in',family = "Times")
  TEST2<-try(print(locuscompare(in_fn1= gwas_fn, in_fn2= eqtl_fn, title1= GWASID, title2 =eqtlID[i])))
  if('try-error'%in% class(TEST2)) {
    dev.off()
    next}
  dev.off()

  print(paste0("------已经做完 ", i, " 个------"))
}

write.csv(result, file="coloc.csv",row.names=FALSE)

