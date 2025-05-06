# Plot network figure
# 2022.05.02

# install.packages("ggraph")
# install.packages("tidyverse")
# install.packages("tidygraph")
# install.packages("RColorBrewer")
# install.packages("extrafont")

library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(scales)
# library(extrafont)
source("scripts/utils/pfunc.R", encoding = "UTF-8")
source("scripts/utils/Func.R")


############################## 函数定义 ##################################
# 得到所有人新表型的均值，将其转换为complexnet的数据, 数据必须为new_pheno/whole下的文件
get_complexnet_data_from_new_pheno = function(data=data) {
  data = data.frame(map(data, ~mean(abs(.x), na.rm=TRUE)))
  
  
  # 提取行名和列名
  row_names = color_map$organ[c(1:7)]
  col_names = color_map$organ[c(1:7)]
  
  # 创建新的数据框
  new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
  colnames(new_data) = unique(col_names)
  rownames(new_data) = unique(row_names)
  
  # 填充值
  for (i in 1:length(row_names)) {
    for (j in 1:length(col_names)) {
      if (i < j) {
        data_col = paste0("y_res_", row_names[i], "_", col_names[j])
        if (!(data_col %in% colnames(data))) {
          new_data[row_names[i], col_names[j]] = 0
        }
        else {
          new_data[row_names[i], col_names[j]] = data[, data_col]
          
        }
      }
      if (i > j) {
        data_col = paste0("y_res_", col_names[j], "_", row_names[i])
        if (!(data_col %in% colnames(data))) {
          new_data[row_names[i], col_names[j]] = 0
        }
        else {
          new_data[row_names[i], col_names[j]] = data[, data_col]
        }
      }
      if (i == j) {
        data_col = paste0("y_res_", col_names[j], "_", row_names[i])
        new_data[row_names[i], col_names[j]] = 0
        
      }
    }
  }
  output_cca_heatmap_matrix = new_data
  new_data = new_data * upper.tri(new_data)
  new_data[new_data == 0] = NA
  new_data = new_data %>% rownames_to_column("organ1")
  output_cca_edge = pivot_longer(new_data, cols = c(2:8))
  colnames(output_cca_edge) = c("organ1", "organ2", "average_y_res")
  output_cca_edge = na.omit(output_cca_edge)
  output_cca_edge$organ1 = rename_based_on_df(output_cca_edge$organ1, color_map, "organ", "short_organ")
  output_cca_edge$organ2 = rename_based_on_df(output_cca_edge$organ2, color_map, "organ", "short_organ")
  return(output_cca_edge)
}
##########################################################################


# data = read.csv("data/sandbox_new_pheno_without_regression/whole/sandbox_new_pheno.csv")
# 
# 
# # regress data
# regression_data = read_regression_data()
# 
# 
# #################################不同年龄的网络 ################################
# young_id = rownames(regression_data)[regression_data$年龄 >= 20 & regression_data$年龄 < 30]
# mid_id = rownames(regression_data)[regression_data$年龄 >= 30 & regression_data$年龄 < 45]
# old_id = rownames(regression_data)[regression_data$年龄 >= 45 & regression_data$年龄 <= 60]
# 
# 
# young_data = data[data$buId %in% young_id, ]
# rownames(young_data) = NULL
# young_data = young_data %>% column_to_rownames("buId")
# 
# 
# mid_data = data[data$buId %in% mid_id, ]
# rownames(mid_data) = NULL
# mid_data = mid_data %>% column_to_rownames("buId")
# 
# 
# old_data = data[data$buId %in% old_id, ]
# rownames(old_data) = NULL
# old_data = old_data %>% column_to_rownames("buId")
# 
# 
# 
# young_data = get_complexnet_data_from_new_pheno(young_data)
# mid_data = get_complexnet_data_from_new_pheno(mid_data)
# old_data = get_complexnet_data_from_new_pheno(old_data)
# 
# 
# check_path("data/report_on_aging_and_gender/young.csv")
# # write.csv(young_data, "data/report_on_aging_and_gender/young.csv", row.names = FALSE)
# # write.csv(mid_data, "data/report_on_aging_and_gender/mid.csv", row.names = FALSE)
# # write.csv(old_data, "data/report_on_aging_and_gender/old.csv", row.names = FALSE)
# check_path("plot/report_on_aging_and_gender/young.png")


## 使用沙箱导出的数据
young_data = read.csv("liumeng_240108/data/report_on_aging_and_gender/young.csv")
mid_data = read.csv("liumeng_240108/data/report_on_aging_and_gender/mid.csv")
old_data = read.csv("liumeng_240108/data/report_on_aging_and_gender/old.csv")

draw_organ_complexnet(young_data,
                      output_name = "plot/report_on_aging_and_gender/sandbox/young.png",
                      limits = c(0.67, 0.97), node_size = c(0.5, 1.8))
draw_organ_complexnet(mid_data,
                      output_name = "plot/report_on_aging_and_gender/sandbox/mid.png",
                      limits = c(0.67, 0.97), node_size = c(0.5, 1.8))
draw_organ_complexnet(old_data,
                      output_name = "plot/report_on_aging_and_gender/sandbox/old.png",
                      limits = c(0.67, 0.97), node_size = c(0.5, 1.8))


point_size=c()
o = c("B", "H", "Lu", "Li", "S", "Pa", "K")
for (i in seq_along(o)){
  point_size[i]=mean(young_data$average_y_res[which(young_data$organ1==o[i]|young_data$organ2==o[i])])
}
rader_data = data.frame(row.names = o, young = point_size)


point_size=c()
for (i in seq_along(o)){
  point_size[i]=mean(mid_data$average_y_res[which(mid_data$organ1==o[i]|mid_data$organ2==o[i])])
}
rader_data$mid = point_size


point_size=c()
for (i in seq_along(o)){
  point_size[i]=mean(old_data$average_y_res[which(old_data$organ1==o[i]|old_data$organ2==o[i])])
}
rader_data$old = point_size
rader_data = rader_data %>% rownames_to_column("organ")
rader_data$organ = rename_based_on_df(rader_data$organ, color_map, from = "short_organ", to = "organ")
check_path("data/RadarMap/sandbox/organ_change_by_age.csv")
write.csv(rader_data, "data/RadarMap/sandbox/organ_change_by_age.csv", row.names = FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #################################不同性别的网络 ################################
# male_id = rownames(regression_data)[regression_data$性别 == 1]
# female_id = rownames(regression_data)[regression_data$性别 == 0]
# 
# 
# male_data = data[data$buId %in% male_id, ]
# rownames(male_data) = NULL
# male_data = male_data %>% column_to_rownames("buId")
# 
# 
# female_data = data[data$buId %in% female_id, ]
# rownames(female_data) = NULL
# female_data = female_data %>% column_to_rownames("buId")
# 
# 
# 
# 
# 
# male_data = get_complexnet_data_from_new_pheno(male_data)
# female_data = get_complexnet_data_from_new_pheno(female_data)


# check_path("data/report_on_aging_and_gender/male.csv")
# write.csv(male_data, "data/report_on_aging_and_gender/sandbox/male.csv", row.names = FALSE)
# write.csv(female_data, "data/report_on_aging_and_gender/sandbox/female.csv", row.names = FALSE)



male_data = read.csv("liumeng_240114/male.csv")
female_data = read.csv("liumeng_240114/female.csv")

check_path("plot/report_on_aging_and_gender/male.png")
limits1 = min(c(male_data$average_y_res, female_data$average_y_res))
limits2 = max(c(male_data$average_y_res, female_data$average_y_res))
draw_organ_complexnet(male_data,
                      output_name = "plot/report_on_aging_and_gender/sandbox/male.png",
                      limits = c(0.67, 0.97), node_size = c(0.5, 1.8))
draw_organ_complexnet(female_data,
                      output_name = "plot/report_on_aging_and_gender/sandbox/female.png",
                      limits = c(0.67, 0.97), node_size = c(0.5, 1.8))


# 
# point_size=c()
# o = c("B", "H", "Lu", "Li", "S", "Pa", "K")
# for (i in seq_along(o)){
#   point_size[i]=mean(male_data$average_y_res[which(male_data$organ1==o[i]|male_data$organ2==o[i])])
# }
# rader_data = data.frame(row.names = o, male = point_size)
# 
# 
# point_size=c()
# for (i in seq_along(o)){
#   point_size[i]=mean(female_data$average_y_res[which(female_data$organ1==o[i]|female_data$organ2==o[i])])
# }
# rader_data$female = point_size
# 
# 
# write.csv("data/RaderMap/sandbox/organ_change_by_age.csv")






male_position = c("left_brain", "right_brain", "heart", "left_lung", "right_lung", 
                  "liver", "spleen", "pancreas", "left_kidney", "right_kidney", "prostate")
female_position = c("left_brain", "right_brain", "heart", "left_lung", "right_lung", 
                    "liver", "spleen", "pancreas", "left_kidney", "right_kidney", "uterus")

male_organs = unique(c(male_data$organ1, male_data$organ2))
df = data.frame()
for (i in male_organs){
  point_size=mean(male_data$average_y_res[which(male_data$organ1==i|male_data$organ2==i)], na.rm=TRUE)
  temp_df = data.frame(organ = i, mean = point_size)
  df = rbind(df, temp_df)
}
male_df = data.frame(organ = male_position, 
                     value = c(df$mean[df$organ == "B"], df$mean[df$organ == "B"], 
                               df$mean[df$organ == "H"], df$mean[df$organ == "Lu"], 
                               df$mean[df$organ == "Lu"], df$mean[df$organ == "Li"], 
                               df$mean[df$organ == "S"], df$mean[df$organ == "Pa"], 
                               df$mean[df$organ == "K"], df$mean[df$organ == "K"], 
                               0
                     ))
check_path("data/report_on_aging_and_gender/sandbox/male_humanbody.csv")
write.csv(male_df, "data/report_on_aging_and_gender/sandbox/male_humanbody.csv")



female_organs = unique(c(female_data$organ1, female_data$organ2))
df = data.frame()
for (i in female_organs){
  point_size=mean(female_data$average_y_res[which(female_data$organ1==i|female_data$organ2==i)], na.rm=TRUE)
  temp_df = data.frame(organ = i, mean = point_size)
  df = rbind(df, temp_df)
}
female_df = data.frame(organ = female_position, 
                       value = c(df$mean[df$organ == "B"], df$mean[df$organ == "B"], 
                                 df$mean[df$organ == "H"], df$mean[df$organ == "Lu"], 
                                 df$mean[df$organ == "Lu"], df$mean[df$organ == "Li"], 
                                 df$mean[df$organ == "S"], df$mean[df$organ == "Pa"], 
                                 df$mean[df$organ == "K"], df$mean[df$organ == "K"], 
                                 0
                       ))
write.csv(female_df, "data/report_on_aging_and_gender/sandbox/female_humanbody.csv")













##################################ukb########################################


data = read.csv("data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv")

basic_info = read.csv("data/basic_info.csv", check.names = FALSE, row.names = "eid")
covariates = read.csv("data/covariates.csv", check.names = FALSE, row.names = "eid")




#################################不同年龄的网络 ################################
# t1 = quantile(covariates$age, probs = 0.33)
# t2 = quantile(covariates$age, probs = 0.66)

              
young_id = rownames(covariates)[covariates$age >= 45 & covariates$age < 60]
mid_id = rownames(covariates)[covariates$age >= 60 & covariates$age < 70]
old_id = rownames(covariates)[covariates$age >= 70 & covariates$age <= 85]


young_data = data[data$eid %in% young_id, ]
rownames(young_data) = NULL
young_data = young_data %>% column_to_rownames("eid")


mid_data = data[data$eid %in% mid_id, ]
rownames(mid_data) = NULL
mid_data = mid_data %>% column_to_rownames("eid")


old_data = data[data$eid %in% old_id, ]
rownames(old_data) = NULL
old_data = old_data %>% column_to_rownames("eid")



young_data = get_complexnet_data_from_new_pheno(young_data)
mid_data = get_complexnet_data_from_new_pheno(mid_data)
old_data = get_complexnet_data_from_new_pheno(old_data)


check_path("data/report_on_aging_and_gender/ukb/young.csv")
write.csv(young_data, "data/report_on_aging_and_gender/ukb/young.csv", row.names = FALSE)
write.csv(mid_data, "data/report_on_aging_and_gender/ukb/mid.csv", row.names = FALSE)
write.csv(old_data, "data/report_on_aging_and_gender/ukb/old.csv", row.names = FALSE)
check_path("plot/report_on_aging_and_gender/ukb/young.png")
draw_organ_complexnet(young_data,
                      output_name = "plot/report_on_aging_and_gender/ukb/young.png",
                      limits = c(0.44, 0.81), node_size = c(0.45, 2))
draw_organ_complexnet(mid_data,
                      output_name = "plot/report_on_aging_and_gender/ukb/mid.png",
                      limits = c(0.44, 0.81), node_size = c(0.45, 2))
draw_organ_complexnet(old_data,
                      output_name = "plot/report_on_aging_and_gender/ukb/old.png",
                      limits = c(0.44, 0.81), node_size = c(0.45, 2))


point_size=c()
o = c("B", "H", "Lu", "Li", "S", "Pa", "K")
for (i in seq_along(o)){
  point_size[i]=mean(young_data$average_y_res[which(young_data$organ1==o[i]|young_data$organ2==o[i])])
}
rader_data = data.frame(row.names = o, young = point_size)


point_size=c()
for (i in seq_along(o)){
  point_size[i]=mean(mid_data$average_y_res[which(mid_data$organ1==o[i]|mid_data$organ2==o[i])])
}
rader_data$mid = point_size


point_size=c()
for (i in seq_along(o)){
  point_size[i]=mean(old_data$average_y_res[which(old_data$organ1==o[i]|old_data$organ2==o[i])])
}
rader_data$old = point_size
rader_data = rader_data %>% rownames_to_column("organ")
rader_data$organ = rename_based_on_df(rader_data$organ, color_map, from = "short_organ", to = "organ")
check_path("data/RadarMap/ukb/organ_change_by_age.csv")
write.csv(rader_data, "data/RadarMap/ukb/organ_change_by_age.csv", row.names = FALSE)









#################################不同性别的网络 ################################
male_id = rownames(covariates)[covariates$isMale == 1]
female_id = rownames(covariates)[covariates$isMale == 0]


male_data = data[data$eid %in% male_id, ]
rownames(male_data) = NULL
male_data = male_data %>% column_to_rownames("eid")


female_data = data[data$eid %in% female_id, ]
rownames(female_data) = NULL
female_data = female_data %>% column_to_rownames("eid")





male_data = get_complexnet_data_from_new_pheno(male_data)
female_data = get_complexnet_data_from_new_pheno(female_data)


check_path("data/report_on_aging_and_gender/ukb/male.csv")
write.csv(male_data, "data/report_on_aging_and_gender/ukb/male.csv", row.names = FALSE)
write.csv(female_data, "data/report_on_aging_and_gender/ukb/female.csv", row.names = FALSE)


check_path("plot/report_on_aging_and_gender/ukb/male.png")
draw_organ_complexnet(male_data,
                      output_name = "plot/report_on_aging_and_gender/ukb/male.png",
                      limits = c(0.44, 0.81), node_size = c(0.45, 2))
draw_organ_complexnet(female_data,
                      output_name = "plot/report_on_aging_and_gender/ukb/female.png",
                      limits = c(0.44, 0.81), node_size = c(0.45, 2))


# 
# point_size=c()
# o = c("B", "H", "Lu", "Li", "S", "Pa", "K")
# for (i in seq_along(o)){
#   point_size[i]=mean(male_data$average_y_res[which(male_data$organ1==o[i]|male_data$organ2==o[i])])
# }
# rader_data = data.frame(row.names = o, male = point_size)
# 
# 
# point_size=c()
# for (i in seq_along(o)){
#   point_size[i]=mean(female_data$average_y_res[which(female_data$organ1==o[i]|female_data$organ2==o[i])])
# }
# rader_data$female = point_size
# 
# 
# write.csv("data/RaderMap/organ_change_by_age.csv")


male_position = c("left_brain", "right_brain", "heart", "left_lung", "right_lung", 
                  "liver", "spleen", "pancreas", "left_kidney", "right_kidney", "prostate")
female_position = c("left_brain", "right_brain", "heart", "left_lung", "right_lung", 
                    "liver", "spleen", "pancreas", "left_kidney", "right_kidney", "uterus")

male_organs = unique(c(male_data$organ1, male_data$organ2))
df = data.frame()
for (i in male_organs){
  point_size=mean(male_data$average_y_res[which(male_data$organ1==i|male_data$organ2==i)], na.rm=TRUE)
  temp_df = data.frame(organ = i, mean = point_size)
  df = rbind(df, temp_df)
}
male_df = data.frame(organ = male_position, 
                     value = c(df$mean[df$organ == "B"], df$mean[df$organ == "B"], 
                               df$mean[df$organ == "H"], df$mean[df$organ == "Lu"], 
                               df$mean[df$organ == "Lu"], df$mean[df$organ == "Li"], 
                               df$mean[df$organ == "S"], df$mean[df$organ == "Pa"], 
                               df$mean[df$organ == "K"], df$mean[df$organ == "K"], 
                               NA
                               ))
write.csv(male_df, "data/report_on_aging_and_gender/ukb/male_humanbody.csv")



female_organs = unique(c(female_data$organ1, female_data$organ2))
df = data.frame()
for (i in female_organs){
  point_size=mean(female_data$average_y_res[which(female_data$organ1==i|female_data$organ2==i)], na.rm=TRUE)
  temp_df = data.frame(organ = i, mean = point_size)
  df = rbind(df, temp_df)
}
female_df = data.frame(organ = female_position, 
                     value = c(df$mean[df$organ == "B"], df$mean[df$organ == "B"], 
                               df$mean[df$organ == "H"], df$mean[df$organ == "Lu"], 
                               df$mean[df$organ == "Lu"], df$mean[df$organ == "Li"], 
                               df$mean[df$organ == "S"], df$mean[df$organ == "Pa"], 
                               df$mean[df$organ == "K"], df$mean[df$organ == "K"], 
                               NA
                     ))
write.csv(female_df, "data/report_on_aging_and_gender/ukb/female_humanbody.csv")

