library(tidyverse)

source("scripts/utils/pfunc.R", encoding = "UTF-8")
source("scripts/utils/Func.R")


#######################################organ####################################
data = read.csv("data/mediation/sandbox_med_top50_whole.csv")
data = data[data$x == "sex", ]

a = data

data = data %>%
  mutate(
    med_cate = case_when(
      (grepl(pattern = "P1", mediator, ignore.case = FALSE)) & (grepl(pattern = "DCs", mediator, ignore.case = TRUE)) ~ "DCs", 
      (grepl(pattern = "P1", mediator, ignore.case = FALSE)) & (grepl(pattern = "mono", mediator, ignore.case = TRUE)) ~ "Mono", 
      (grepl(pattern = "P2", mediator, ignore.case = FALSE)) ~ "NK", 
      (grepl(pattern = "P3", mediator, ignore.case = FALSE)) ~ "Treg", 
      (grepl(pattern = "P4", mediator, ignore.case = FALSE)) ~ "T-cell", 
      (grepl(pattern = "P5", mediator, ignore.case = FALSE)) ~ "Th", 
      (grepl(pattern = "P6", mediator, ignore.case = FALSE)) ~ "B-cell", 
      TRUE ~ "IF&MF"
    )
  )
output_df = expand.grid(organ_name=color_map$organ[1:7], category_name=unique(data$med_cate))
nmap(output_df, )




microphenotypes = ReadMicrophenomics()
microphenotypes$Serum = NULL
microphenotypes$IF = full_join(rownames_to_column(microphenotypes$IF, var = "ID"), 
                               rownames_to_column(microphenotypes$ELISA, var = "ID")) %>% column_to_rownames(var = "ID")
microphenotypes$ELISA = NULL
microphenotypes$DCs = microphenotypes$FACS[grepl("P1", colnames(microphenotypes$FACS), ignore.case = FALSE) & 
                                             grepl("DCs", colnames(microphenotypes$FACS), ignore.case = TRUE)]
microphenotypes$Mono = microphenotypes$FACS[grepl("P1", colnames(microphenotypes$FACS), ignore.case = FALSE) & 
                                              grepl("mono", colnames(microphenotypes$FACS), ignore.case = TRUE)]
microphenotypes$NK = microphenotypes$FACS[grepl("P2", colnames(microphenotypes$FACS), ignore.case = FALSE)]
microphenotypes$Treg = microphenotypes$FACS[grepl("P3", colnames(microphenotypes$FACS), ignore.case = FALSE)]
microphenotypes$`T-cell` = microphenotypes$FACS[grepl("P4", colnames(microphenotypes$FACS), ignore.case = FALSE)]
microphenotypes$Th = microphenotypes$FACS[grepl("P5", colnames(microphenotypes$FACS), ignore.case = FALSE)]
microphenotypes$`B-cell` = microphenotypes$FACS[grepl("P6", colnames(microphenotypes$FACS), ignore.case = FALSE)]
microphenotypes$FACS = NULL
microphenotypes = map(microphenotypes, get_res_with_regression_data, 
                      regression_data = regression_data)


for (i in c(1: length(organs))) {
  output_csv = data.frame()
  organ_name <- names(organs)[i]
  organ <- organs[[i]]
  for (j in c(1: length(microphenotypes))) {
    category <- microphenotypes[[j]]
    category_name <- names(microphenotypes)[j]
    for (k in c(1: length(category))) {
      i_data = organ
      j_data = category[, k, drop=FALSE]
      temp_i_j_data = merge_and_divide(i_data, j_data)
      i_data  = temp_i_j_data[[1]]
      j_data = temp_i_j_data[[2]]
      N = dim(i_data)[1]
      tryCatch({    # 删除方差为0的列
        i_data <- i_data[c(map_dfc(i_data, var, na.rm=TRUE) != 0)]
        j_data <- j_data[c(map_dfc(j_data, var, na.rm=TRUE) != 0)]}, 
        error = function(e) {})
      i_data_len = length(i_data)
      j_data_len = length(j_data)
      
      if (i_data_len == 1) {
        i_data["copy_col"] = i_data[[1]]
      } 
      if (j_data_len == 1) {
        j_data["copy_col"] = j_data[[1]]
      }
      tryCatch({cc_result <- PMA::CCA(i_data, j_data, typex = "standard", typez = "standard",
                                      penaltyx = 1, penaltyz = 1, trace = FALSE, 
                                      standardize = TRUE)
      u <- cc_result$u
      v <- cc_result$v
      cor_main_coef <- cc_result$cors
      cc1_organ <- as.matrix(i_data) %*% u
      cc2_organ <- as.matrix(j_data) %*% v
      cor_main_coef.p = cor.test(cc1_organ, cc2_organ)$p.value
      temp_csv <- data.frame(organ_name=organ_name, category_name=category_name,
                             feature = colnames(j_data)[1], 
                             r=cor_main_coef, 
                             p=cor_main_coef.p, 
                             N=N, row.names = NULL)
      output_csv <- rbind(output_csv, temp_csv)
      }, 
      error = function(e) {
      })
    }
  }
  output_csv$p.adjust = p.adjust(output_csv$p, method = "fdr")
  output_csv = output_csv %>%
    dplyr::select(-"N", 'N')
  check_path(paste0(resultPath, "/Results/Organ/", organ_name, "_microphenotypes_cca_nto1_mht_v1.csv"))
  check_path(paste0(mainPath, "/Data/VennDiagram/Results/Organ/", organ_name, "_microphenotypes_cca_nto1_mht_v1.csv"))
  file_name = paste0(resultPath, "/Results/Organ/", organ_name, "_microphenotypes_cca_nto1_mht_v1.csv")
  write.csv(output_csv, file = file_name, row.names = FALSE)
  write.csv(output_csv, file = paste0(mainPath, "/Data/VennDiagram/Results/Organ/", organ_name, "_microphenotypes_cca_nto1_mht_v1.csv"),
            row.names = FALSE)
}


positions = list.files(paste0(resultPath, "/Results/Organ/"), full.names = TRUE, pattern = ".csv")
positions = positions[unlist(map(color_map$organ, grep, positions))]

data = map_dfr(positions, read.csv) 
data = data %>% filter(category_name != "Serum")
data$category_name[data$category_name == "IF"] = "IF&MF"
data$category_name[data$category_name == "MF"] = "IF&MF"

# data = data %>% filter(category_name != "IF")
# data = data %>% filter(category_name != "MF")

temp = expand.grid(organ_name=unique(data$organ_name), category_name=unique(data$category_name))

data = data %>%
  filter(p <= p_threshold) %>%
  group_by(organ_name, category_name) %>%
  summarise(bar_start=min(r, na.rm = TRUE), bar_end=max(r, na.rm = TRUE), number=n())
df = temp %>% left_join(data, c("organ_name", "category_name"))
df[is.na(df)] = 0
df = rename(df, name=category_name)
df = df %>%
  arrange(organ_name, name)

cell_order = c("T-cell", "Th", "Treg", "B-cell", "IF&MF", "Mono", "Dcs", "NK")
df = df[order(match(df$name, cell_order)), ]
df = df[order(match(df$organ_name, color_map$organ)), ]




# # remove IF and MF
# df = df %>% 
#   filter(name != "IF") %>%
#   filter(name != "MF")
df1 = df[df$name == unique(df$name)[1:4], ]
df2 = df[df$name == unique(df$name)[5:8], ]


check_path(paste0(resultPath, "/OveralBar_data/organ1.csv"))
check_path(paste0(resultPath, "/OveralBar_data/organ2.csv"))
write.csv(df1, paste0(resultPath, "/OveralBar_data/organ1.csv"), row.names = FALSE)
write.csv(df2, paste0(resultPath, "/OveralBar_data/organ2.csv"), row.names = FALSE)
