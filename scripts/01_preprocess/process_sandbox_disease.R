setwd("E:/Projects/indNet/")
library(tidyverse)
data = read.csv("正常队列扫描名单及异常病变记录-中山.csv", check.names = FALSE)
data = data %>%
  dplyr::select(`编号`, "AN-head":"AN-lung")
data = data[!duplicated(data$编号), ]
rownames(data) = NULL
data = data %>% column_to_rownames("编号")
data[data == "0"] = 0
data[data == "1"] = 0
data[data == ""] = 0
data[] <- lapply(data, function(x) ifelse(nchar(x) > 1, 1, x))
data = data %>% rownames_to_column("buId")
write.csv(data, "Sandbox/disease_sandbox.csv", row.names = FALSE)
