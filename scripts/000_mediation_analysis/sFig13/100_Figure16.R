# user-defined path
# mainPath <- "E:/seafile/Seafile/SandBox/"  # local
#mainPath <- "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
# codePath <- paste(mainPath,"/Code/R/MediationBarPlot/",sep="")
# dataPath <- paste(mainPath,"/Sandbox/",sep="")
# resultPath = paste(mainPath, "/Data/MediationBarPlot/", sep="")
# set path and load data
# setwd(codePath)

getwd()
dataPath <- './Sandbox/'
codePath <- './scripts/000_mediation_analysis/sFig13/'

source("scripts/utils/Func.R", encoding = "UTF-8")
source("scripts/utils/read_sandbox_data.R", encoding = "UTF-8")
library(tidyverse)

p_threshold = 0.05

read_data = function(path) {
  ll = list.files(path, full.names = TRUE, pattern = "csv")
  data = lapply(ll, function(l) {
    temp_df = read.csv(l, fileEncoding="UTF-8-BOM", header=TRUE)
    organ = gen_part(df = temp_df, File_encoding = "UTF-8", File_fileEncoding = "UTF-8-BOM",
                     RefFile_fileEncoding = "UTF-8-BOM", part_name = "SaveName")
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
organs = read_data(paste0(dataPath, "Organ"))

VA = read.csv(paste0(dataPath, "Organ/heart_atrial_volume.csv"), header = TRUE, row.names = 1, check.names = FALSE)
VV = read.csv(paste0(dataPath, "Organ/heart_ventricular_volume.csv"), header = TRUE, row.names = 1, check.names = FALSE)
cardiac = merge(VA, VV, by=0)
rownames(cardiac) = cardiac$Row.names
cardiac$Row.names = NULL
cardiac_cate = read.csv(paste0(dataPath, "Category/Heart.csv"), encoding = "UTF-8")
colnames(cardiac_cate)[1] = "Name"
select_category = c("Function", "Structure")
cate_cardiac = divide_data_by_category_and_region_df(cardiac, cardiac_cate, select_category)
cate_cardiac = map(cate_cardiac, rename_based_on_df, cardiac_cate, "Name", "Short_name")
names(cate_cardiac)[2] = "Shape"
organs$Heart = append(organs$Heart, cate_cardiac)
organs = organs[order(match(names(organs), color_map$organ))]

read_regression_data = function() {
  BMI <- read.csv(file = paste0(dataPath, "Demography/DXA体成分-体成分高分辨率数字成像化系统.csv"),
                    sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
  BMI <- BMI[, 1, drop=FALSE]
  sex_age <- read.csv(file = paste0(dataPath, "Demography/基本信息.csv"),
                        sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
  sexual <- sex_age[, 1, drop=FALSE]
  sexual <- map_df(sexual, ~ ifelse(. == "男", 1, 0))
  rownames(sexual) <- rownames(sex_age)
  age <- sex_age[, 2, drop = FALSE]
  partial <- cbind(BMI, age, sexual)
  return(partial)
}
regression_data = read_regression_data()

temp = organs %>% 
  map(~map_chr(.x, is.data.frame)) %>%
  map(as.logical)
organs = map2(organs, temp, ~.x[.y])

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
microphenotypes = map(microphenotypes, get_res_with_regression_data,
                      regression_data = regression_data)


output_csv = data.frame(matrix(NA, nrow = 0, ncol = 5))
colnames(output_csv) = c("organ", "micro_type", "microphenomics", 
                         "group", "value")
for (i in seq_along(organs)) {
  organ_name = names(organs)[i]
  organ = organs[[i]]
  for (j in seq_along(microphenotypes)) {
    micro_type_name = names(microphenotypes)[j]
    micro_type = microphenotypes[[j]]
    
    for (k in seq_along(micro_type)) {
      micro_name = names(micro_type)[k]
      micro = micro_type[, k, drop=FALSE]
      
      idx = intersect(rownames(na.omit(organ)), rownames(na.omit(micro)))
      organ = organ[idx, , drop=FALSE]
      micro = micro[idx, , drop=FALSE]
      # 删除方差为0数据
      organ <- organ[c(map_dfc(organ, var, na.rm=TRUE) != 0)]
      micro <- micro[c(map_dfc(micro, var, na.rm=TRUE) != 0)]
      
      # get_cca_result = function(X, Y) {
      #   if (dim(X)[1] != 0 & dim(Y)[1] != 0) {
      #     
      #     if (ncol(X) == 1) {
      #       X["copy_col"] = X[[1]]
      #     } 
      #     if (ncol(Y) == 1) {
      #       Y["copy_col"] = Y[[1]]
      #     }
      #     
      #     cc_result <- PMA::CCA(X, Y, typex = "standard", typez = "standard", 
      #                           penaltyx = 1, penaltyz = 1, trace = FALSE, 
      #                           standardize = TRUE)
      #     cor_main_coef <- cc_result$cors
      #     return(cor_main_coef)
      #   }
      # }
      
      
      cca = get_cca_result(organ, micro)$cca
      cca.p = get_cca_result(organ, micro)$cca.p
      
      temp_output_csv = data.frame(organ = organ_name, 
                                   microphenomics = micro_name, 
                                   group = micro_type_name, 
                                   value = cca, p = cca.p)
      output_csv = rbind(output_csv, temp_output_csv)
    }
  }
}


f = function(df) {
  df %>% 
    filter(p <= p_threshold) %>%
    slice_max(value, n=10) %>%
    rename(individual=microphenomics) %>%
    rownames_to_column(var = "id") %>%
    right_join(data.frame(id=as.character(c(1:10))))
}

data = output_csv %>% 
  group_split(organ) %>%
  map(f)

names(data) = output_csv %>%
  group_split(organ) %>%
  map_chr(~unique(.x$organ))


func = function(df) {
  df$value[is.na(df$value)] = 0
  df$individual[is.na(df$individual)] = " "
  return(df)
}
data = map(data, func)
check_path(paste0(resultPath, "/Organ/"))
walk(data, ~write.csv(.x, paste0(resultPath, 
                                              "/Organ/", 
                                              .x[which(.x[["id"]] == 1), "organ"], 
                                              ".csv"), 
                                   row.names = FALSE, fileEncoding = "UTF-8"))

