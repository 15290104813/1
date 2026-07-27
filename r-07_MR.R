rm(list = ls())
setwd("//")
if (! dir.exists("./07_MR")){
  dir.create("./07_MR")
}
setwd("./07_MR")

library(clusterProfiler)
library(plinkbinr)
library(VariantAnnotation)
library(gwasglue)
library(data.table)
library(tidyr)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(org.Hs.eg.db)

# ENSG00000124813——RUNX2
# ENSG00000204217——BMPR2

gene <- c('RUNX2','BMPR2')
id<-bitr(gene,fromType = "SYMBOL",
         toType = c("ENSEMBL" ),
         OrgDb = 'org.Hs.eg.db')
id$ENSEMBL<-paste('eqtl','a',id$ENSEMBL,sep='-')
ids <- id$ENSEMBL
write.csv(id,'id.csv')

load('../00_rawdata/outcome_all.Rdata')

head(outcome_all)
dat.res <- list()
expo.res <- list()
results_all<-data.frame()
steiger_all<-data.frame()
clumping_all<-data.frame()
pelpvalue_all <- data.frame()
hetpvalue_all <- data.frame()
i=1
for (i in 1:length(ids))  {
  setwd("//07_MR/")
  
  expos_vcf = paste('//MR/eQTL_vcf/',ids[i],'.vcf.gz',sep='')
  vcf1 <- VariantAnnotation::readVcf(expos_vcf,"hg19")
  exposure_dat <- gwasvcf_to_TwoSampleMR(vcf1,type = "exposure")
  exposure_dat$id.exposure = ids[i]
  exposure_dat$exposure = paste(' || id:',ids[i],sep='')
  exp_data <- exposure_dat
  temp_dat <- exp_data[exp_data$pval.exposure<5e-08,]
  if (nrow(temp_dat) == 0) {
    next
  }

  temp_dat$id <- temp_dat$id.exposure
  temp_dat$rsid <- temp_dat$SNP
  temp_dat$pval <- temp_dat$pval.exposure

  exp_dat <- ld_clump(temp_dat,
                      plink_bin = get_plink_exe(),
                      bfile = '//SMR/00.g1000_eur/g1000_eur' ,
                      clump_kb =100, clump_r2 = 0.001)
  if (nrow(exp_dat) ==0 ) {
    next
  }
  ## 结局变量
  outcome_dat<- subset(outcome_all,outcome_all$SNP %in% exp_dat$SNP)
  if (nrow(outcome_dat) == 0) {
    next
  }
  dat <- harmonise_data(exposure_dat =exp_dat,outcome_dat = outcome_dat)
  clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
  clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
  clumping <- clumping[!is.na(clumping$F), ]

  dat = clumping[clumping$F >10,]

  if (nrow(dat) < 3 ) {
    next
  }
  expo.res[[i]] <- exp_dat
  dat.res[[i]] <- dat
  
  result <- generate_odds_ratios(mr_res = mr(dat))
  result2 <- subset(result,method == 'Inverse variance weighted')
  result2
  if (result$nsnp[1] < 3) {
    next
  }
  
  ple <- mr_pleiotropy_test(dat)
  het <- mr_heterogeneity(dat)
  #steiger_sl<-steiger_filtering(dat)
  steiger_sl2<-directionality_test(dat)
  
  results_all<-rbind(results_all,result)
  pelpvalue_all <- rbind(pelpvalue_all,ple)
  hetpvalue_all <- rbind(hetpvalue_all,het)
  clumping_all <- rbind(clumping_all,clumping)
  steiger_all <- rbind(steiger_all,steiger_sl2)
  
  
  dir.create(paste(id$SYMBOL[i],id$ENSEMBL[i],sep='_'))
  setwd(paste(id$SYMBOL[i],id$ENSEMBL[i],sep='_'))
  system.time(save(dat, file = "MR.Rdata"))
  
  expo.res[[i]] <- exp_dat
  dat.res[[i]] <- dat
  
  write.csv(dat,file = '01.dat.csv',row.names=F)
  write.csv(result,file = '02.result.csv')
  write.csv(ple,file = '03.pleiotropy.csv',row.names=F)
  write.csv(het,file = '04.heterogeneity.csv',row.names=F)
  write.csv(clumping,file = '05.clumping.csv',row.names=F)
  clumping = clumping[,c("SNP","id.exposure","F")]
  write.csv(steiger_sl2,file = '06.steiger.csv',row.names=F)
  
  #绘图
  library(patchwork)
  library(ggplot2)
  p1 <- mr_scatter_plot(dat,mr_results = mr(dat))[[1]] +
    scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 15)+
    theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),  
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size = 15))+
    theme(legend.background = element_blank(),
          legend.position = 'top',
          legend.direction = 'vertical',
          legend.text = element_text(size = 15)) +plot_layout(nrow = 1, byrow = FALSE)
  
  pdf('07.scatter.pdf',w=9,h=9,family='Times')
  print(p1)
  dev.off()
  
  png('07.scatter.png',w=9,h=9,units='in',res=600,family='Times')
  print(p1)
  dev.off()
  
  
  res_single <- mr_singlesnp(dat, all_method = c("mr_ivw"))
  p2<-mr_forest_plot(res_single)[[1]] +
    theme_classic(base_size = 15)+
    theme(panel.border = element_rect(size = 1.2,fill = 'transparent'),  
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'none',
          # legend.title = element_blank(),
          # legend.position = c(0.1,0.9),
          legend.text = element_text(size = 15))+
    theme(legend.background = element_blank(),
          legend.title = element_blank(),
          legend.position = 'none',
          legend.direction = 'vertical',
          plot.title=element_text(size=6),
          legend.text = element_text(size = 15)) + plot_layout(nrow = 1, byrow = FALSE)
  
  pdf('08.forest.pdf',w=12,h=7+0.3*unique(result2$nsnp),family='Times')
  print(p2)
  dev.off()
  
  png('08.forest.png',w=12,h=7+0.3*unique(result2$nsnp),units='in',res=600,family='Times')
  print(p2)
  dev.off()
  
  p3<-mr_funnel_plot(res_single)[[1]] +
    scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 20)+
    theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),  
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size = 18))+
    theme(legend.background = element_blank(),
          legend.position = 'top',
          legend.direction = 'vertical',
          legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
  
  pdf('09.funnel.pdf',w=8,h=9,family='Times')
  print(p3)
  dev.off()
  
  png('09.funnel.png',w=8,h=9,units='in',res=600,family='Times')
  print(p3)
  dev.off()
  
  
  
  single <- mr_leaveoneout(dat)
  p4<-mr_leaveoneout_plot(single)[[1]] +
    theme_classic(base_size = 15)+
    theme(panel.border = element_rect(size = 1.6,fill = 'transparent'),  
          axis.ticks = element_line(size = 1),
          legend.background = element_blank(),
          legend.position = 'none',
          legend.text = element_text(size = 18))+
    theme(legend.background = element_blank(),
          legend.title = element_blank(),
          legend.position = 'none',
          legend.direction = 'vertical',
          plot.title=element_text(size=8),
          legend.text = element_text(size = 18)) + plot_layout(nrow = 1, byrow = FALSE)
  
  pdf('10.leaveoneout.pdf',w=11,h=7+0.3*unique(result2$nsnp),family='Times')
  print(p4)
  dev.off()
  
  png('10.leaveoneout.png',w=11,h=7+0.3*unique(result2$nsnp),units='in',res=600,family='Times')
  print(p4)
  dev.off()
  
  library(MRPRESSO)
  mr_presso<-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
                       SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                       data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
  write.csv(mr_presso$`MR-PRESSO results`$`Global Test`,file = '13.mr_presso.csv',row.names=F)
  i = i+1
}

setwd("//07_MR/")
# save.image("MR.dat.RData")
expo.all <- do.call(rbind,expo.res)
dat.all <- do.call(rbind,dat.res)
dat <- split(dat,list(dat$id.exposure))
write.csv(expo.all,file="expo.all.csv")
write.csv(dat.all,file="dat.all.csv")
write.csv(steiger_all,'steiger_all.csv')
write.csv(clumping_all,'clumping_all.csv')
#write.csv(mr_presso_all,'mr_presso_all.csv')

results_all <- results_all[results_all$nsnp > 2,]
write.csv(results_all,'01.results_all_TrMGs.csv')
results<-results_all[results_all$method=='Inverse variance weighted',]
results<-results[results$pval<0.05,]
results$ENSEMBL<-results$id.exposure
results<-merge(results,id,by='ENSEMBL')
write.csv(results,'02.results.csv')

pelpvalue <- pelpvalue_all[pelpvalue_all$pval>0.05,]
pelpvalue$ENSEMBL <- pelpvalue$id.exposure
pelpvalue <- merge(pelpvalue,id,by='ENSEMBL')
write.csv(pelpvalue,'pleiotropy_all.csv')

gene <- intersect(results$SYMBOL,pelpvalue$SYMBOL)
write.csv(gene,'gene.csv')
