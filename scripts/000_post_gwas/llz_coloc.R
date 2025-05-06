# install.packages('R.utils')

#1 下载、安装coloc
# if(!require("remotes"))
#   install.packages("remotes")
# install.packages("dplyr")
# install.packages("remote")

library(remotes)
library(tidyverse)
# install_github("chr1swallace/coloc",build_vignettes=TRUE)
library("coloc")
library(dplyr)
library(data.table)

disease_gwas_files = list.files("fuma/ukb/disease/", full.names = TRUE)
pheno_gwas_files = list.files("fuma/ukb/new_pheno/", full.names = TRUE)

for (i in seq_along(disease_gwas_files)) {
  for (j in seq_along(pheno_gwas_files)) {
    disease_name = str_split(basename(disease_gwas_files[i]), pattern = "_")[[1]][1]
    new_pheno_name = str_split(basename(pheno_gwas_files[j]), pattern = "_")[[1]][1]
    df1 = fread(disease_gwas_files[i], sep = "\t")
    df2 = fread(pheno_gwas_files[j], sep = "\t")
    df1 = na.omit(df1)
    df2 = na.omit(df2)
    
    df1$varbeta = df1$se * df1$se
    df2$varbeta = df2$se * df2$se
    
    df1 = df1[!(duplicated(df1$rsid)), ]
    df2 = df2[!(duplicated(df2$rsid)), ]
    
    
    results = coloc.abf(list(snp=df1$rsid, beta=df1$beta, varbeta=df1$varbeta, type="quant"), 
                        list(snp=df2$rsid, beta=df2$beta, varbeta=df2$varbeta, type="cc"))
    
    
    # 3.5 筛选共定位的位点
    #通常情况下，很多文献认为PPA > 0.95的位点是共定位位点，也有一些文献会放松要求到0.75。
    #这里假定后验概率大于0.95为共定位位点：
    filter_result=results$results %>% filter(SNP.PP.H4 > 0.25)
    write.csv(filter_result, paste0("fuma/ukb/coloc/", disease_name, "_", new_pheno_name, ".csv"))
  }
}
