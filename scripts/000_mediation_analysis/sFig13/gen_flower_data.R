# user-defined path
mainPath <- "E:/Projects/indNet/"  # local

library(tidyverse)
library(VennDiagram)

source("scripts/utils/pfunc.R", encoding = "UTF-8")
source("scripts/utils/Func.R")

molecular_dict = read.csv("Sandbox/Category/molecular.csv")
feature_dict = read.csv("Sandbox/Category/feature_sly.csv")
# filter p <= p_threshold
p_threshold = ifelse(grepl("Seafile", mainPath), 0.05, 0.05)


#######################IFMF#########################################
organs_df = read.csv("data/mediation/sandbox_med_top50_sig.csv")

organs_df = organs_df[(organs_df$x %in% feature_dict$Short_name[feature_dict$Select == "TRUE"]), ]
organs_df$x = rename_based_on_df(organs_df$x, nmapdf = feature_dict, from = "Short_name", to = "DispName")

organs_df = organs_df[, c("x", "y", "mediator", "prop.mediated", "prop.mediated.p")]
organs_df$y = str_remove(organs_df$y, pattern = "y_res_")
organs_df = separate_rows(organs_df, y, sep = "_", convert = TRUE)
organs_df = organs_df[organs_df$mediator %in% molecular_dict$micro, ]


organs_df = organs_df %>%
  mutate(
    category_name = case_when(
      (grepl(pattern = "P1", mediator, ignore.case = FALSE)) & (grepl(pattern = "DCs", mediator, ignore.case = FALSE)) ~ "DCs", 
      (grepl(pattern = "P1", mediator, ignore.case = FALSE)) & (grepl(pattern = "mono", mediator, ignore.case = FALSE)) ~ "Mono", 
      (grepl(pattern = "P2", mediator, ignore.case = FALSE)) ~ "NK", 
      (grepl(pattern = "P3", mediator, ignore.case = FALSE)) ~ "Treg", 
      (grepl(pattern = "P4", mediator, ignore.case = FALSE)) ~ "T-cell", 
      (grepl(pattern = "P5", mediator, ignore.case = FALSE)) ~ "Th", 
      (grepl(pattern = "P6", mediator, ignore.case = FALSE)) ~ "B-cell", 
      TRUE ~ "IF&MF"
    )
  )

organs_df = organs_df %>%
  group_split(y)
names(organs_df) = c("Brain", "Heart", "Kidney", "Liver", "Lung", "Pancreas", "Spleen")
organ_names = c("Brain", "Heart", "Kidney", "Liver", "Lung", "Pancreas", "Spleen")


organs_df = organs_df %>%
  map(filter, category_name=="IF&MF") %>%
  map(filter, prop.mediated.p<=p_threshold)
organs = organs_df %>%
  map(~.x$mediator)
names(organs) = organ_names
organs = organs[order(match(names(organs), color_map$organ))]
organ_total_num = unlist(map(organs, ~length(unique(.x))))


inter = get.venn.partitions(organs)
inter = inter %>%
  dplyr::mutate(inter_organ_num = rowSums(.[1:7]))

# flower表格
inter_table = inter %>%
  filter(..count.. !=0) %>%
  select(-..set..)
inter_table = apply(inter_table, 2, as.character)
# check_path(paste0(resultPath, "/FACS_table.csv"))
# write.csv(inter_table, paste0(resultPath, "/FACS_table.csv"), row.names = FALSE)

# write.csv
flower_data = inter %>%
  filter(inter_organ_num %in% c(1))
center_data = inter %>%
  group_by(inter_organ_num) %>%
  summarise(total=sum(..count..))
data = data.frame(organ_names = c(colnames(flower_data)[7:1]), 
                  inter_count = as.numeric(flower_data$..count..))
data$inter_count = paste0(data$inter_count, "/", organ_total_num[7:1])
check_path(paste0("data/sfig13/", "IFMF", "_flower.csv"))
check_path(paste0("data/sfig13/", "IFMF", "_flower_center.csv"))
write.csv(data, file = paste0("data/sfig13/", "IFMF", "_flower.csv"), row.names = FALSE)
write.csv(center_data, file = paste0("data/sfig13/", "IFMF", "_flower_center.csv"), row.names = FALSE)






######################facs#################################
organs_df = read.csv("data/mediation/sandbox_med_top50_sig.csv")

organs_df = organs_df[(organs_df$x %in% feature_dict$Short_name[feature_dict$Select == "TRUE"]), ]
organs_df$x = rename_based_on_df(organs_df$x, nmapdf = feature_dict, from = "Short_name", to = "DispName")

organs_df = organs_df[, c("x", "y", "mediator", "prop.mediated", "prop.mediated.p")]
organs_df$y = str_remove(organs_df$y, pattern = "y_res_")
organs_df = separate_rows(organs_df, y, sep = "_", convert = TRUE)

organs_df = organs_df %>%
  mutate(
    category_name = case_when(
      (grepl(pattern = "P1", mediator, ignore.case = FALSE)) & (grepl(pattern = "DCs", mediator, ignore.case = FALSE)) ~ "DCs", 
      (grepl(pattern = "P1", mediator, ignore.case = FALSE)) & (grepl(pattern = "mono", mediator, ignore.case = FALSE)) ~ "Mono", 
      (grepl(pattern = "P2", mediator, ignore.case = FALSE)) ~ "NK", 
      (grepl(pattern = "P3", mediator, ignore.case = FALSE)) ~ "Treg", 
      (grepl(pattern = "P4", mediator, ignore.case = FALSE)) ~ "T-cell", 
      (grepl(pattern = "P5", mediator, ignore.case = FALSE)) ~ "Th", 
      (grepl(pattern = "P6", mediator, ignore.case = FALSE)) ~ "B-cell", 
      TRUE ~ "IF&MF"
    )
  )

organs_df = organs_df %>%
  group_split(y)
names(organs_df) = c("Brain", "Heart", "Kidney", "Liver", "Lung", "Pancreas", "Spleen")
organ_names = c("Brain", "Heart", "Kidney", "Liver", "Lung", "Pancreas", "Spleen")


organs_df = organs_df %>%
  map(filter, category_name=="DCs" | category_name=="Mono" | category_name=="NK" | 
        category_name=="Treg" | category_name=="T-cell" | category_name=="Th" |
        category_name=="B-cell") %>%
  map(filter, prop.mediated.p<=p_threshold)
organs = organs_df %>%
  map(~.x$mediator)
names(organs) = organ_names
organs = organs[order(match(names(organs), color_map$organ))]
organ_total_num = unlist(map(organs, ~length(unique(.x))))


inter = get.venn.partitions(organs)
inter = inter %>%
  dplyr::mutate(inter_organ_num = rowSums(.[1:7]))

# flower表格
inter_table = inter %>%
  filter(..count.. !=0) %>%
  select(-..set..)
inter_table = apply(inter_table, 2, as.character)
# check_path(paste0(resultPath, "/FACS_table.csv"))
# write.csv(inter_table, paste0(resultPath, "/FACS_table.csv"), row.names = FALSE)

# write.csv
flower_data = inter %>%
  filter(inter_organ_num %in% c(1))
center_data = inter %>%
  group_by(inter_organ_num) %>%
  summarise(total=sum(..count..))
data = data.frame(organ_names = c(colnames(flower_data)[7:1]), 
                  inter_count = as.numeric(flower_data$..count..))
data$inter_count = paste0(data$inter_count, "/", organ_total_num[7:1])
check_path(paste0("data/sfig13/", "FACS", "_flower.csv"))
check_path(paste0("data/sfig13/", "FACS", "_flower_center.csv"))
write.csv(data, file = paste0("data/sfig13/", "FACS", "_flower.csv"), row.names = FALSE)
write.csv(center_data, file = paste0("data/sfig13/", "FACS", "_flower_center.csv"), row.names = FALSE)


