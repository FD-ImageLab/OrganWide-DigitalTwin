library(tidyverse)
library(data.table)
library(tools)
library("PMA")
library(writexl)

source("scripts/utils/pfunc.R", encoding = "UTF-8")
source("scripts/utils/Func.R")



feature_dict = read.csv("Sandbox/ukb_category/feature.csv")
features = read.csv("data/feature_top.csv", row.names = "eid", check.names = FALSE)
# features = features[c("20445-0.0", "102280-0.0", "34-0.0", "21003-0.0", "21022-0.0",
#                       "1588-0.0", "4429-0.0", "22602-0.0", "20109-0.0", "30510-0.0")]

# moleculars phenos
moleculars = read.csv("data/moleculars_60k.csv", row.names = "eid", check.names = FALSE)

new_pheno_files = list.files("data/ukb_new_pheno_without_regression", full.names = TRUE, pattern = ".csv")

new_pheno_id = rownames(na.omit(new_pheno))
# filter features，> 1000 intersection with organs
features = features[map_dbl(seq_along(features), ~sum(rownames(na.omit(features[, .x, drop=FALSE])) %in% new_pheno_id)) >= 1000]

library(parallel)

process_pheno_feature_molecular <- function(new_pheno_file, feature_idx, molecular_idx) {
  new_pheno = read.csv(new_pheno_file, row.names = "eid")
  feature = features[, feature_idx, drop=FALSE]
  molecular = moleculars[, molecular_idx, drop=FALSE]
  
  df = reduce(list(feature, new_pheno, molecular), .f = my_join)
  df = na.omit(df)
  df = df[1:min(nrow(df), 3000), ]
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

# generate tasks
tasks <- expand.grid(new_pheno_files = new_pheno_files,
                     feature_idx = seq_along(features)[1],
                     molecular_idx = seq_along(moleculars)[1])
tasks$new_pheno_files = as.character(tasks$new_pheno_files)


results_list <- mclapply(1:nrow(tasks), function(i) {
  task <- tasks[i, ]
  process_pheno_feature_molecular(task$new_pheno_files, task$feature_idx, task$molecular_idx)
}, mc.cores = 120) 

# merge results
results <- do.call(rbind, results_list)

write.csv(results, "data/mediation/top10_whole.csv")