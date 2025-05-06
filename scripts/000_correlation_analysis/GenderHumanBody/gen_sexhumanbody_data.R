# user-defined path
mainPath <- "E:/seafile/Seafile/SandBox/"  # local
# mainPath <- "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
codePath <- paste(mainPath,"/Code/R/HumanBody/",sep="")
dataPath <- paste(mainPath,"/Sandbox/",sep="")
resultPath = paste(mainPath, "/Data/HumanBody/", sep="")
# set path and load data
setwd(codePath)

source(paste(mainPath,"/Code/R/Toolbox/Func.R",sep=""), encoding = "UTF-8")
library(tidyverse)


read_data = function(path, part_name="SaveName") {
  ll = list.files(path, full.names = TRUE, pattern = "csv")
  data = lapply(ll, function(l) {
    temp_df = read.csv(l, fileEncoding="UTF-8-BOM", header=TRUE)
    organ = gen_part(df = temp_df, File_encoding = "UTF-8", File_fileEncoding = "UTF-8-BOM",
                     RefFile_fileEncoding = "UTF-8-BOM", part_name = part_name)
    return(organ)
  })
  data_names = lapply(ll, function(l) {
    temp_df = read.csv(l, fileEncoding="UTF-8-BOM", header=TRUE)
    organ_name = names(table(temp_df$Organ[temp_df$Organ > 0]))[which.max(table(temp_df$Organ[temp_df$Organ > 0]))]
    return(organ_name)
  })
  names(data) = data_names
  return(data)
}
organs = read_data(paste0(mainPath, "/Code/R/Manhattan/Category/Organ/"), part_name = "DispName")

temp = organs %>%
  map(~map_chr(.x, is.data.frame)) %>%
  map(as.logical)
organs = map2(organs, temp, ~.x[.y])



regression_data = read_regression_data(SEX_ = FALSE, BMI_ = TRUE, AGE_ = TRUE)
sex = read_regression_data(SEX_ = TRUE, BMI_ = FALSE, AGE_ = FALSE)
# 回归后的器官数据
organs = map(organs, ~map(.x, get_res_with_regression_data, regression_data = regression_data))


organs = organs %>%
  map(~map(.x, rownames_to_column, var="buId"))

join_by_buId = function(x, y) {inner_join(x, y, by="buId")}
reduce_join = function(l) {l %>% reduce(join_by_buId)}
organs = organs %>%
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")
organs = organs[order(match(names(organs), color_map$organ))]


read_left_right_data = function() {
  left_right_organs = list()
  Brain = organs$Brain
  left_brain = Brain[grep("_l", colnames(Brain))]
  right_brain = Brain[grep("_r", colnames(Brain))]
  left_right_organs$Brain = list(left_brain = left_brain, right_brain = right_brain)
  
  Heart = organs$Heart
  left_heart = Heart[c(grep("LA", colnames(Heart)), grep("LV", colnames(Heart)), 
                       grep("LPA", colnames(Heart)), grep("LPV", colnames(Heart)))]
  right_heart = Heart[c(grep("RA", colnames(Heart)), grep("RV", colnames(Heart)), 
                        grep("RPA", colnames(Heart)), grep("RPV", colnames(Heart)))]
  left_right_organs$Heart = list(left_heart = left_heart, right_heart = right_heart)
  
  Kidney = organs$Kidney
  left_kidney = Kidney[grep("L_", colnames(Kidney))]
  right_kidney = Kidney[grep("R_", colnames(Kidney))]
  left_right_organs$Kidney = list(left_kidney = left_kidney, right_kidney = right_kidney)
  
  Lung = organs$Lung
  left_lung = Lung[grep("L_", colnames(Lung))]
  right_lung = Lung[grep("R_", colnames(Lung))]
  left_right_organs$Lung = list(left_lung = left_lung, right_lung = right_lung)
  return(left_right_organs)
}
left_right_organs = read_left_right_data()


male_idx = rownames(sex)[sex$性别 == 1]
female_idx = rownames(sex)[sex$性别 == 0]

male_position = c("left_brain", "right_brain", "heart", "left_lung", "right_lung", 
             "liver", "spleen", "pancreas", "left_kidney", "right_kidney", "prostate")
female_position = c("left_brain", "right_brain", "heart", "left_lung", "right_lung", 
                  "liver", "spleen", "pancreas", "left_kidney", "right_kidney", "uterus")

WholeBody = list(left_brain = left_right_organs$Brain$left_brain, 
                 right_brain = left_right_organs$Brain$right_brain, 
                 heart = organs$Heart, 
                 left_lung = left_right_organs$Lung$left_lung, 
                 right_lung = left_right_organs$Lung$right_lung,
                 liver = organs$Liver, 
                 spleen = organs$Spleen, 
                 pancreas = organs$Pancreas, 
                 left_kidney = left_right_organs$Kidney$left_kidney, 
                 right_kidney = left_right_organs$Kidney$right_kidney, 
                 prostate = organs$Prostate, 
                 uterus = organs$Uterus)

f = function(df, idx) {return(df[rownames(df) %in% idx, ])}
MaleWholeBody = map(WholeBody, f, idx=male_idx)
FemaleWholeBody = map(WholeBody, f, idx=female_idx)


output_cancor_coef_df = data.frame()
for (i in c(1: length(MaleWholeBody))) {
  for (j in c(1: length(MaleWholeBody))) {
    i_data_name = names(MaleWholeBody)[i]
    j_data_name = names(MaleWholeBody)[j]
    i_data = MaleWholeBody[[i]]
    j_data = MaleWholeBody[[j]]
    if (i != j & nrow(na.omit(merge(i_data, j_data, by=0, check.names = FALSE))) != 0
        & (!(i==1 & j==2)) & (!(i==4 & j==5)) & (!(i==9 & j==10))
        & (!(i==2 & j==1)) & (!(i==5 & j==4)) & (!(i==10 & j==9))) {
      temp = get_cca_result(i_data, j_data)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=i_data_name, organ2=j_data_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}
male_heatmap_data = output_cancor_coef_df %>%
  group_by(organ1) %>%
  summarise(mean_r = mean(CC1, na.rm = TRUE))

male_heatmap_data = male_heatmap_data[order(match(male_heatmap_data$organ1, male_position)), ]

output_cancor_coef_df = data.frame()
for (i in c(1: length(FemaleWholeBody))) {
  for (j in c(1: length(FemaleWholeBody))) {
    i_data_name = names(FemaleWholeBody)[i]
    j_data_name = names(FemaleWholeBody)[j]
    i_data = FemaleWholeBody[[i]]
    j_data = FemaleWholeBody[[j]]
    if (i != j & nrow(na.omit(merge(i_data, j_data, by=0, check.names = FALSE))) != 0
        & (!(i==1 & j==2)) & (!(i==4 & j==5)) & (!(i==9 & j==10))
        & (!(i==2 & j==1)) & (!(i==5 & j==4)) & (!(i==10 & j==9))) {
      temp = get_cca_result(i_data, j_data)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=i_data_name, organ2=j_data_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}

female_heatmap_data = output_cancor_coef_df %>%
  group_by(organ1) %>%
  summarise(mean_r = mean(CC1, na.rm = TRUE))
female_heatmap_data = female_heatmap_data[order(match(female_heatmap_data$organ1, female_position)), ]

check_path(paste0(resultPath, "/Male.csv"))
check_path(paste0(resultPath, "/Female.csv"))
write.csv(male_heatmap_data, paste0(resultPath, "/Male.csv"), row.names = FALSE)
write.csv(female_heatmap_data, paste0(resultPath, "/Female.csv"), row.names = FALSE)