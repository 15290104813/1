rm(list = ls())
setwd("//")
if (! dir.exists("./09_scan")){
  dir.create("./09_scan")
}
setwd("./09_scan")

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

#此步骤需要填入令牌 孟德尔令牌需要到IEU OpenGWAS官网申请https://api.opengwas.io/?next=%2Fprofile%2F
#需要注册Github账号
Sys.setenv(OPENGWAS_JWT="令牌")
#options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
ieugwasr::get_opengwas_jwt()
ieugwasr::api_status()

#BMPR2-----------------
gwas <- read.csv('gwas_result_BMPR2.csv')
gwas <- subset(gwas,X!= '')
write.csv(gwas,'GWAS_BMPR2.csv')
colnames(gwas) <-c('disease','GWAS') 

gene <- 'BMPR2'
id<-bitr(gene,fromType = "SYMBOL",
         toType = c("ENSEMBL" ),
         OrgDb = 'org.Hs.eg.db')
id$ENSEMBL<-paste('eqtl','a',id$ENSEMBL,sep='-')
ids <- id$ENSEMBL

info <- data.frame(disease = '',expo_id = "",res_pval = '',or = '',or_lci95 = '',or_uci95 = '',
                   ple_pval = '')
info <- info[-1,]
i=1
for (i in 1:length(gwas$GWAS)) {
  tryCatch({
    setwd("//09_scan/")
    exposure <- extract_instruments(outcomes = ids,p1 = 5e-06,r2 = 0.001,clump = T,kb=10)
    clumped_exposure <- exposure
    outcome <-  extract_outcome_data(snps = clumped_exposure$SNP,outcomes = gwas$GWAS[i],proxies = TRUE)
    dat <- harmonise_data(exposure_dat = clumped_exposure,outcome_dat = outcome)
    result <- generate_odds_ratios(mr_res = mr(dat,method_list = c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_simple_mode","mr_weighted_mode")))
    result2 <- subset(result,method == 'Inverse variance weighted')
    ple <- mr_pleiotropy_test(dat)
    
    info <- rbind(info,data.frame(disease = gwas$disease[i],expo_id = gwas$GWAS[i],res_pval = result2$pval,
                                  or = result2$or,or_lci95 = result2$or_lci95,or_uci95 = result2$or_uci95,ple_pval = ple$pval))
    
    dir.create(paste(gwas$disease[i],gwas$GWAS[i],sep='_'))
    setwd(paste(gwas$disease[i],gwas$GWAS[i],sep='_'))
    system.time(save(dat, file = "MR.Rdata"))
    
    write.csv(dat,file = '01.dat.csv',row.names=F)
    write.csv(result,file = '02.result.csv')
    write.csv(ple,file = '03.pleiotropy.csv',row.names=F)
    
    het <- mr_heterogeneity(dat)
    het
    write.csv(het,file = '04.heterogeneity.csv',row.names=F)

    clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
    clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
    write.csv(clumping,file = '05.clumping.csv',row.names=F)

    steiger_sl2<-directionality_test(dat)
    write.csv(steiger_sl2,file = '06.steiger.csv',row.names=F)
  }, error = function(e){
    print('error')
  },finally = {
    print(gwas$GWAS[i])
  })
}

setwd("//09_scan/")
write.csv(info,'all_results_BMPR2.csv',row.names = F)

info2 <- subset(info,res_pval < 0.05 & ple_pval > 0.05)
info2$gene <- 'BMPR2' 
write.csv(info2,'hub_results_BMPR2.csv',row.names = F)



#RUNX2-----------------
gwas <- read.csv('gwas_result_RUNX2.csv')
gwas <- subset(gwas,X!= '')
write.csv(gwas,'GWAS_RUNX2.csv')
colnames(gwas) <-c('disease','GWAS') 

gene <- 'RUNX2'
id<-bitr(gene,fromType = "SYMBOL",
         toType = c("ENSEMBL" ),
         OrgDb = 'org.Hs.eg.db')
id$ENSEMBL<-paste('eqtl','a',id$ENSEMBL,sep='-')
ids <- id$ENSEMBL

info <- data.frame(disease = '',expo_id = "",res_pval = '',or = '',or_lci95 = '',or_uci95 = '',
                   ple_pval = '')
info <- info[-1,]
i=1
for (i in 1:length(gwas$GWAS)) {
  tryCatch({
    setwd("//09_scan/")
    exposure <- extract_instruments(outcomes = ids,p1 = 5e-06,r2 = 0.001,clump = T,kb=10)
    clumped_exposure <- exposure
    outcome <-  extract_outcome_data(snps = clumped_exposure$SNP,outcomes = gwas$GWAS[i],proxies = TRUE)
    dat <- harmonise_data(exposure_dat = clumped_exposure,outcome_dat = outcome)
    result <- generate_odds_ratios(mr_res = mr(dat,method_list = c("mr_egger_regression","mr_weighted_median","mr_ivw","mr_simple_mode","mr_weighted_mode")))
    result2 <- subset(result,method == 'Inverse variance weighted')
    ple <- mr_pleiotropy_test(dat)
    
    info <- rbind(info,data.frame(disease = gwas$disease[i],expo_id = gwas$GWAS[i],res_pval = result2$pval,
                                  or = result2$or,or_lci95 = result2$or_lci95,or_uci95 = result2$or_uci95,ple_pval = ple$pval))
    
    dir.create(paste(gwas$disease[i],gwas$GWAS[i],sep='_'))
    setwd(paste(gwas$disease[i],gwas$GWAS[i],sep='_'))
    system.time(save(dat, file = "MR.Rdata"))
    
    write.csv(dat,file = '01.dat.csv',row.names=F)
    write.csv(result,file = '02.result.csv')
    write.csv(ple,file = '03.pleiotropy.csv',row.names=F)
    
    het <- mr_heterogeneity(dat)
    het
    write.csv(het,file = '04.heterogeneity.csv',row.names=F)

    clumping = mutate(dat,R=get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure))
    clumping = mutate(clumping,F=(samplesize.exposure-2)*((R*R)/(1-R*R)))
    write.csv(clumping,file = '05.clumping.csv',row.names=F)

    steiger_sl2<-directionality_test(dat)
    write.csv(steiger_sl2,file = '06.steiger.csv',row.names=F)
  }, error = function(e){
    print('error')
  },finally = {
    print(gwas$GWAS[i])
  })
}

setwd("//09_scan/")
write.csv(info,'all_results_RUNX2.csv',row.names = F)

info2 <- subset(info,res_pval < 0.05 & ple_pval > 0.05)
info2$gene <- 'RUNX2' 
write.csv(info2,'hub_results_RUNX2.csv',row.names = F)


gene1 <- read.csv('hub_results_BMPR2.csv')
gene2 <- read.csv('hub_results_RUNX2.csv')
results <- rbind(gene1,gene2)

out_multi <- results
hz <- paste(round(out_multi$or,5),
            "(",round(out_multi$or_lci95,5),
            "-",round(out_multi$or_uci95,5),")",sep = "")


tabletext <- cbind(c(NA,"Gene",out_multi$gene),
                   c(NA,"Outcomes",out_multi$disease),
                   c(NA,"P value",ifelse(out_multi$res_pval<0.001,
                                         "< 0.001",
                                         round(out_multi$res_pval,5))),
                   c(NA,"Odd Ratio(95% CI)",hz))
tabletext2 <- data.frame(tabletext)
write.csv(tabletext2,'tabletext.csv')

library(forestplot)
pdf(file = "01.forest.pdf", height = 6, width =12, onefile = F,family='Times')
forestplot(labeltext=tabletext,
           graph.pos=5, 
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="darkred", lines="darkblue", zero = "gray50"),
           mean=c(NA,NA,out_multi$or),
           lower=c(NA,NA,out_multi$or_lci95), 
           upper=c(NA,NA,out_multi$or_uci95), 
           boxsize=0.15,lwd.ci=3, 
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,  
           colgap=unit(5,"mm"),   
           xticks = c(0.8,0.9,1,1.1,1.2), 
           lwd.xaxis=2,          
           lineheight = unit(1.2,"cm"),
           graphwidth = unit(.5,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI,
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Odd Ratio",
           grid = T,
           clip = c(0,3)) 
dev.off()


png(file = "01.forest.png", height = 6, width = 12, units='in',res=600,family='Times')
forestplot(labeltext=tabletext,
           graph.pos=5,  
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="darkred", lines="darkblue", zero = "gray50"),
           mean=c(NA,NA,out_multi$or),
           lower=c(NA,NA,out_multi$or_lci95), 
           upper=c(NA,NA,out_multi$or_uci95), 
           boxsize=0.15,lwd.ci=3,  
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=0.5,  
           colgap=unit(5,"mm"),   
           xticks = c(0.8,0.9,1,1.1,1.2), 
           lwd.xaxis=2,         
           lineheight = unit(1.2,"cm"), 
           graphwidth = unit(.5,"npc"),
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Odd Ratio",
           grid = T,
           clip = c(0,3)) 
dev.off()
