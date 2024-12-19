# user-defined path
mainPath <- "E:/seafile/Seafile/SandBox/"  # local
# mainPath <- "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
codePath <- paste(mainPath,"/Code/R/Manhattan/",sep="")
dataPath <- paste(mainPath,"/Sandbox/",sep="")
resultPath = paste(mainPath, "/Data/Manhattan/", sep="")
# set path and load data
setwd(codePath)
library(tidyverse)
args <- commandArgs(TRUE)
eval(parse(text = args[1]))
source(paste(mainPath,"/Code/R/Toolbox/Func.R",sep=""), encoding = "UTF-8")
source(paste(mainPath,"/Code/R/Toolbox/read_sandbox_data.R",sep=""), encoding = "utf-8")


p_threshold = ifelse(grepl("Seafile", mainPath), 1, 0.05)


#######################################organ####################################
#  OrganHeatMap category
regression_data = read_regression_data()
organs = read_data(paste0(mainPath, "/Code/R/Category/Modality/"), 
                   part_name = "Prefix", class = "Modality")
temp = organs %>% 
  map(~map_chr(.x, is.data.frame)) %>%
  map(as.logical)
organs = map2(organs, temp, ~.x[.y])


# organ data after regression
organs = map(organs, ~map(.x, get_res_with_regression_data, regression_data = regression_data))

organs = organs %>% 
  map(~map(.x, rownames_to_column, var="buId"))

join_by_buId = function(x, y) {inner_join(x, y, by="buId")}
reduce_join = function(l) {l %>% reduce(join_by_buId)}
organs = organs %>% 
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")
organs = organs[order(match(names(organs), color_map$organ))]



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
      tryCatch({    # remove cols with square = 0
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
