library(tidyverse)
library(data.table)
library(tools)
library("PMA")
library(writexl)

source("scripts/utils/pfunc.R", encoding = "UTF-8")
source("scripts/utils/Func.R")


# 小分子表型
moleculars = ReadMicrophenomics(dataPath = "Sandbox/")
moleculars = reduce(moleculars, my_join)
moleculars = remove_col_by_na_ratio(moleculars, na_ratio_threshold = 0.5)


features = read.csv("data/sandbox_feature.csv", row.names = "buId", check.names = FALSE)
top_feature_name = read.csv("data/sandbox_feature_net_sort_by_wholenet.csv")
features = features[, colnames(features) %in% top_feature_name$feature[1:50]]

new_pheno_files = list.files("data/sandbox_new_pheno_without_regression", full.names = TRUE, pattern = ".csv")

library(parallel)
# 包装循环逻辑的函数
process_pheno_feature_molecular <- function(new_pheno_file, feature_idx, molecular_idx) {
  new_pheno = read.csv(new_pheno_file, row.names = "buId")
  feature = features[, feature_idx, drop=FALSE]
  molecular = moleculars[, molecular_idx, drop=FALSE]
  
  df = reduce(list(feature, new_pheno, molecular), .f = my_join)
  df = na.omit(df)
  df = data.frame(scale(df), check.names = FALSE)
  
  tryCatch({
    temp_results = my_mediation(df = df, sims=500)
  }, 
  error=function(e) {
    print(paste0("When get error, feature_idx=", feature_idx, 
                 "molecular_idx=", molecular_idx, "new_pheno_file=", 
                 new_pheno_file))
    return(NULL)
  })
}

# 生成任务列表
tasks <- expand.grid(new_pheno_files = new_pheno_files,
                     feature_idx = seq_along(features),
                     molecular_idx = seq_along(moleculars))
tasks$new_pheno_files = as.character(tasks$new_pheno_files)


# 使用 mclapply 并行处理任务
results_list <- mclapply(1:nrow(tasks), function(i) {
  print(i)
  task <- tasks[i, ]
  process_pheno_feature_molecular(task$new_pheno_files, task$feature_idx, task$molecular_idx)
}, mc.cores = 96)  # 使用除一个以外的所有核心

# 合并结果
results <- do.call(rbind, results_list)

check_path("data/mediation/sandbox_med_50_whole.csv")
write.csv(results, "data/mediation/sandbox_med_top50_whole.csv")