library(dplyr)
library(purrr)
source("D:/Desktop/2/Code/Func.R")

# 读取FumaIDMapping.csv文件
mapping_data <- read.csv("D:/Desktop/2/FumaIDMapping.csv")
dict_data <- read.csv("D:/Desktop/2/pheno_dict.csv")

# 生成文件夹路径
setwd("D:/Desktop/2/fuma_results")
current_directory <- getwd()
folder_numbers <- list.files("./")
folder_paths <- paste0(current_directory, "/", folder_numbers)
print(folder_paths)

# 获取所有leadSNPs.txt文件的路径，只存储文件存在的路径
file_paths <- map(folder_paths, ~ {
  leadSNPs_path <- file.path(.x, "leadSNPs.txt")
  if (file.exists(leadSNPs_path)) leadSNPs_path else NULL
}) %>% compact()

# 更新folder_numbers为只包含实际存在的leadSNPs.txt文件的文件夹路径
folder_numbers <- map(file_paths, ~ basename(dirname(.)))
print(folder_numbers)

# 打印结果
print(file_paths)

# 通过purrr::map函数读取所有leadSNPs.txt文件并创建一个包含所有数据框的列表
data_frames <- map2(file_paths, folder_numbers, ~ read.table(file = .x, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
                      mutate(Type = as.character(.y)))

# 利用dplyr::bind_rows函数将所有数据框合并成一个大的数据框
merged_data <- bind_rows(data_frames)

# 将Type列移到第一列
merged_data <- merged_data %>%
  select(Type, everything())

# 映射Type列
merged_data$Type <- rename_based_on_df(merged_data$Type, mapping_data, from = "job_id", to = "eid")
merged_data$Type <- rename_based_on_df(merged_data$Type, dict_data, from = "fid", to = "description")

# 从 merged_data 中提取需要的列
output_data <- merged_data %>%
  select(Type, rsID, chr, pos)

# 指定输出文件路径
output_file_path <- file.path(current_directory, "output_data.csv")

# 将 output_data 保存为 CSV 文件
write.csv(output_data, file = output_file_path, row.names = FALSE)

# 打印输出文件路径
print(paste("Output data saved to:", output_file_path))



