#APOE
#phenotype	num_cases	num_controls
#Vascular dementia	3624	475484
#Vascular dementia (multiple infarctions)	620	474702
#Vascular dementia (mixed)	531	474702
#Vascular dementia (other)	200	474702
#Vascular dementia (subcortical)	1031	474702
#Vascular dementia (sudden onset)	184	474702
#Vascular dementia (undefined)	1656	474702
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(data.table)
library(ggplot2)
library(MendelianRandomization)
library(ieugwasr)
library(plinkbinr)
library(metafor)
library(openxlsx)
info <- read.table("./血管性痴呆下载性状.txt", header = T, comment.char = "#",sep="\t")
#暴露工具变量
exposure <- read.xlsx("./brain-csf-plasma-pqtl.xlsx")
brain <- exposure[exposure$Tissue == "Brain",] #608个  616iv
csf <- exposure[exposure$Tissue == "CSF",] #214个 233iv
plasma <- exposure[exposure$Tissue == "Plasma",] #612个 616iv
# 暴露
exp_dat_clump <- csf[which(csf$Protein == "APOE"),]

i = 1

t1 <- Sys.time()
# 结局
#结局
outcome <- fread(paste0("./finngen_R12_", info$phenocode[i], ".gz"))
outcome$phenotype <- info$phenotype[i]
outcome$samplesize <- info$num_cases[i]+info$num_controls[i]
outcome <- as.data.frame(outcome)
out_dat <- format_data(dat = outcome,
                       type = "outcome",
                       snps = exp_dat_clump$SNP,
                       phenotype_col = "phenotype",
                       snp_col = "rsids",
                       beta_col = "beta",
                       se_col = "sebeta",
                       effect_allele_col = "alt",
                       other_allele_col = "ref",
                       pval_col = "pval",
                       chr_col = "#chrom",
                       pos_col = "pos",
                       eaf_col = "af_alt",
                       samplesize_col= "samplesize")
out_dat <- out_dat %>% subset(., !duplicated(SNP))

exp_data_clump <- subset(exp_dat_clump, SNP %in% exp_dat_clump$SNP)

exp_data_clump <- format_data(dat = exp_data_clump,
                              type = "exposure",
                              phenotype_col = "protein",
                              snp_col = "SNP",
                              beta_col = "beta",
                              se_col = "se",
                              effect_allele_col = "effect_allele",
                              other_allele_col = "other_allele",
                              pval_col = "p")

mydata <- harmonise_data(exposure_dat=exp_data_clump,outcome_dat=out_dat,action=2)
mydata <- mydata[which(mydata$mr_keep==TRUE),]

mydata_clump <- mydata
####################################
mydata_filteroutcome <- mydata_clump
mydata_outcome_n<-dim(mydata_filteroutcome)[1]
Nbd <- 10000 #nrow(mydata_filteroutcome)/0.05+1

if(mydata_outcome_n<=3){
  presso_pval<-"NA"
  mydata_presso<-mydata_filteroutcome
}else{
  presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                      SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                      data = mydata_filteroutcome, 
                      NbDistribution = Nbd,  SignifThreshold = 0.05)
  presso
  presso_pval<-presso$`MR-PRESSO results`$`Global Test`$Pvalue
  presso_snp<-presso$`MR-PRESSO results`$`Outlier Test`
  ###mydata_presso pvalue
  ###presso_pval>0.05，不需去除outliers
  if(presso_pval>=0.05){
    #mydata_presso
    mydata_presso<-mydata_filteroutcome
  }else{
    ###presso_pval<0.05，说明存在多效性，需要去除outliers
    #去除outliers
    ############################
    out_order<-order(presso_snp$Pvalue)
    out_order
    for(ii in 1:length(out_order)){
      snp_ii<-out_order[1:ii]
      mydata_presso<-mydata_filteroutcome[-snp_ii,]
      mydata_pres_n<-dim(mydata_presso)[1]
      if(mydata_pres_n<=3){
        presso_pval<-"NA"
      }else{
        presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                             SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                             data = mydata_presso, 
                             NbDistribution = Nbd,  SignifThreshold = 0.05)
        presso_pval<-presso2$`MR-PRESSO results`$`Global Test`$Pvalue
        #mydata_presso
        mydata_presso<-mydata_presso
      }
      print(ii);
      if(presso_pval>=0.05){
        break;##终止for循环
      }
    }
    ############################
  }
}


#filter
exp_n<-dim(exp_dat_clump)[1]
out_n<-dim(out_dat)[1]
mydata_har_n<-dim(mydata)[1]
mydata_clu_n<-dim(mydata_clump)[1]
mydata_outcome_n<-dim(mydata_filteroutcome)[1]
mydata_pres_n<-dim(mydata_presso)[1]
filter_data<-data.frame(exp_n,out_n,mydata_clu_n,mydata_har_n,mydata_outcome_n,mydata_pres_n,presso_pval)
filter_data

##########MR分析
set.seed(5201314)
res<-mr(mydata_presso, method_list=c("mr_egger_regression", "mr_ivw","mr_ivw_mre","mr_ivw_fe","mr_wald_ratio","mr_weighted_median","mr_two_sample_ml","mr_simple_mode","mr_weighted_mode"))#常见8种方法
res2<-mr(mydata_presso, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_weighted_mode", "mr_simple_mode"))####出散点图用
res_or<-generate_odds_ratios(res)
##########异质性、多效性、敏感性
het<-mr_heterogeneity(mydata_presso, method_list=c("mr_egger_regression", "mr_ivw"))#the mr_heterogeneity() function can take an argument to only perform heterogeneity tests using specified methods
plt<-mr_pleiotropy_test(mydata_presso)
#I2   第一种方法 MendelianRandomization包
#I2object <- MendelianRandomization::mr_input(bx = mydata_presso$beta.exposure, bxse = mydata_presso$se.exposure, by = mydata_presso$beta.outcome, byse = mydata_presso$se.outcome, snps = mydata_presso$SNP)
#I2 <- MendelianRandomization::mr_ivw(object = I2object, model = "fixed")
#I2@Fstat

#I2   第二种方法
res_single1 <- mr_singlesnp(mydata_presso, all_method = "mr_ivw")
res_single2 <- res_single1[grep("^rs", res_single1$SNP),]
res_meta <- metafor::rma(yi = res_single2$b, sei = res_single2$se, weights = 1/mydata_presso$se.outcome^2,data = res_single2, method = 'FE')
I2 <- as.data.frame(cbind(res_meta$I2, res_meta$H2))
colnames(I2) <- c("I2", "H2") #I2是百分比，后面那个是Heterogeneity test statistic (Cochran's Q)

res_loo<-mr_leaveoneout(mydata_presso)
##########可视化
####散点图
p1<-mr_scatter_plot(res,mydata_presso)
p1.2<-mr_scatter_plot(res2,mydata_presso)
####森林图
res_single<-mr_singlesnp(mydata_presso)
p2<-mr_forest_plot(res_single)
####漏斗图
p3<-mr_funnel_plot(res_single)
####敏感性图
p4<-mr_leaveoneout_plot(res_loo)


###############输出数据
out_name <- info$phenotype[i]
setwd("./APOE_result/")
dir.create(paste0("csf_APOE_", out_name))#创建文件夹
filepath <- paste0("./APOE_result/", "csf_APOE_", out_name,"/")

setwd(filepath)
ggsave(p1[[1]], file="mr_scatter_plot.pdf",height = 7,width = 7)
ggsave(p1.2[[1]], file="mr_scatter_plot2.pdf",height = 7,width = 7)
ggsave(p2[[1]], file="mr_forest_plot.pdf",height = 7,width = 7)
ggsave(p3[[1]], file="mr_funnel_plot.pdf",height = 7,width = 7)
ggsave(p4[[1]], file="mr_leaveoneout_plot.pdf",height = 7,width = 7)
write.table(mydata_presso,"mydata_presso.txt",sep="\t",quote = FALSE,row.names = FALSE)#IVs
write.table(res,"res.txt",sep="\t",quote = FALSE,row.names = FALSE)#res
write.table(res_or,"res_or.txt",sep="\t",quote = FALSE,row.names = FALSE)#res_or
write.table(het,"het.txt",sep="\t",quote = FALSE,row.names = FALSE)#heterogeneity
write.table(plt,"plt.txt",sep="\t",quote = FALSE,row.names = FALSE)#pleiotropy
write.table(filter_data,"filter_data.txt",sep="\t",quote = FALSE,row.names = FALSE)#data filter process
#write.table(steiger_snpi,"steiger_snpi.txt",sep="\t",quote = FALSE,row.names = FALSE)#steiger_filtering
#write.table(out,"out.txt",sep="\t",quote = FALSE,row.names = FALSE)#steiger_pval
#write.table(res_phenoscanner,"res_phenoscanner.txt",sep="\t",quote = FALSE,row.names = FALSE)#res_phenoscanner
write.table(I2,"I2.txt",sep="\t",quote = FALSE,row.names = FALSE)
t2 <- Sys.time()
message("The ", i, " Done in ", round(difftime(t2, t1, units = "mins"), 3), " minutes.")