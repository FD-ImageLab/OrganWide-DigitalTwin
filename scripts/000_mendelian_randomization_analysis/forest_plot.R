# install.packages("forestploter")
library(forestploter)
library(tidyverse)
# Source custom utility functions
source("scripts/utils/Func.R")
paths = list.files("data/MR/", pattern = "*_mr.csv", full.names = TRUE)

dt = map_dfr(paths, read.csv, row.names=1)
dt$exposure = NULL
dt = dplyr::rename(dt, exposure=id.exposure)
dt = dt[dt$method %in% c("MR Egger", "Weighted median", "Inverse variance weighted"), ]
order = unique(arrange(dt, pval)[, c("exposure", "outcome")])
dt$unique_order <- with(dt, match(paste(exposure, outcome), paste(order$exposure, order$outcome)))
dt <- dt[order(dt$unique_order),]
dt$unique_order <- NULL  # Optionally remove the helper column
dt = dt[1:36, ]
dt$pval = sprintf("%.2e", dt$pval)

# dt$id.exposure = ifelse(duplicated(dt$id.exposure), NA, dt$id.exposure)
# dt$id.outcome = ifelse(duplicated(dt$id.outcome), NA, dt$id.outcome)

dt$` ` = paste(rep(" ", 30), collapse = " ")
dt$`(95% CI)` <- sprintf("(%.2f to %.2f)", 
                         dt$b-dt$se, dt$b+dt$se)

pheno_dict = read.csv("scripts/pheno_dict.csv")
dt$exposure = rename_based_on_df(dt$exposure, nmapdf = pheno_dict, from = "fid", "description")
dt$outcome = paste0(dt$outcome, "-0.0")
dt$outcome = rename_based_on_df(dt$outcome, nmapdf = pheno_dict, from = "fid", "description")

# # set row 3n+1,3n+3""
# rows_to_na <- c(seq(1, nrow(dt), by = 3), seq(3, nrow(dt), by=3))
# dt$exposure[rows_to_na] = "" 
# dt$outcome[rows_to_na] = "" 


check_path("plot/MR/forest.png")
png("plot/MR/forest.png", width = 1200, height = 1200)

theme <- forest_theme(base_size = 10,
                      arrow_col = "red",

                      refline_col = "black",
                      footnote_col = "#636363",
                      footnote_fontface = "italic")
forest(dt[, c(1, 3, 9, 4, 10, 8)], theme = theme, 
       est = dt$b,
       lower = dt$b - dt$se,
       upper = dt$b + dt$se,
       ci_column = 3, ticks_digits = 2, ref_line = 0, )
dev.off()


