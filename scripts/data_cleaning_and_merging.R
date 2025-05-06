folder_paths <- list.files("script/table", full.names = TRUE, pattern = "FUMA")
print(folder_paths)

# Get the paths of all leadSNPs.txt files, storing only the paths where the files exist.
file_paths <- map(folder_paths, ~ {
    leadSNPs_path <- file.path(.x, "leadSNPs.txt")
    if (file.exists(leadSNPs_path)) leadSNPs_path else NULL
}) %>% compact()

# Update folder_numbers to only include the folder paths where the leadSNPs.txt file actually exists.
folder_numbers <- map(file_paths, ~ basename(dirname(.)))
print(folder_numbers)

# print result
print(file_paths)

fuma_dict <- read.csv("fuma/fuma_results/fuma/ukb/FumaIDMapping.csv")

# Use the purrr::map function to read all leadSNPs.txt files and create a list containing all the data frames.
data_frames <- map2(file_paths, folder_numbers, ~ read.table(file = .x, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(Type = as.character(.y)))

# Use the dplyr::bind_rows function to combine all the data frames into one large data frame.
merged_data <- bind_rows(data_frames)

merged_data <- merged_data %>%
    dplyr::select(Type, everything())
merged_data <- merged_data %>%
    dplyr::select(-No)
merged_data$Type <- rename_based_on_df(merged_data$Type, nmapdf = fuma_dict, from = "job_id", to = "eid")
merged_data$Type <- rename_based_on_df(merged_data$Type, nmapdf = pheno_dict, from = "fid", to = "description")
merged_data <- dplyr::rename(merged_data, Phenotype = "Type")

merged_data$Phenotype <- rename_based_on_df(merged_data$Phenotype, nmapdf = pheno_dict, from = )

write.csv(merged_data, "C:/Users/liumeng/OneDrive/Desktop/leadsnp.csv", row.names = FALSE)




















folder_paths <- list.files("fuma/fuma_results/fuma/ukb", full.names = TRUE, pattern = "FUMA")
print(folder_paths)

# 获取所有leadSNPs.txt文件的路径，只存储文件存在的路径
file_paths <- map(folder_paths, ~ {
    leadSNPs_path <- file.path(.x, "magma.genes.out")
    if (file.exists(leadSNPs_path)) leadSNPs_path else NULL
}) %>% compact()

# 更新folder_numbers为只包含实际存在的leadSNPs.txt文件的文件夹路径
folder_numbers <- map(file_paths, ~ basename(dirname(.)))
print(folder_numbers)

# 打印结果
print(file_paths)

fuma_dict <- read.csv("fuma/fuma_results/fuma/ukb/FumaIDMapping.csv")

# 通过purrr::map函数读取所有leadSNPs.txt文件并创建一个包含所有数据框的列表
data_frames <- map2(file_paths, folder_numbers, ~ read.table(file = .x, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(Type = as.character(.y)))

# 利用dplyr::bind_rows函数将所有数据框合并成一个大的数据框
merged_data <- bind_rows(data_frames)

merged_data <- merged_data %>%
    dplyr::select(Type, everything())

merged_data$Type <- rename_based_on_df(merged_data$Type, nmapdf = fuma_dict, from = "job_id", to = "eid")
merged_data$Type <- rename_based_on_df(merged_data$Type, nmapdf = pheno_dict, from = "fid", to = "description")
merged_data <- dplyr::rename(merged_data, Phenotype = "Type")
write.csv(merged_data, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/magma.csv", row.names = FALSE)
merged_data <- merged_data[merged_data$P < 5 * 10e-8, ]
write.csv(merged_data, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/magma_sig.csv", row.names = FALSE)











files <- list.files("data/MR/", pattern = "mr.csv", full.names = TRUE)
data <- map_dfr(files, read.csv, row.names = 1)
data <- data[, 3:ncol(data)]
data <- data %>% select(-c("outcome", "exposure"))
data$id.outcome <- paste0(data$id.outcome, "-0.0")
data$id.exposure <- rename_based_on_df(data$id.exposure, nmapdf = pheno_dict, from = "fid", to = "description")
data$id.outcome <- rename_based_on_df(data$id.outcome, nmapdf = pheno_dict, from = "fid", to = "description")
data <- data %>%
    arrange(pval)
write.csv(data, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/mr.csv", row.names = FALSE)











folder_paths <- list.files("fuma/fuma_results/fuma/ukb", full.names = TRUE, pattern = "FUMA")
print(folder_paths)

# 获取所有leadSNPs.txt文件的路径，只存储文件存在的路径
file_paths <- map(folder_paths, ~ {
    leadSNPs_path <- file.path(.x, "gwascatalog.txt")
    if (file.exists(leadSNPs_path)) leadSNPs_path else NULL
}) %>% compact()

# 更新folder_numbers为只包含实际存在的leadSNPs.txt文件的文件夹路径
folder_numbers <- map(file_paths, ~ basename(dirname(.)))
print(folder_numbers)

# 打印结果
print(file_paths)

fuma_dict <- read.csv("fuma/fuma_results/fuma/ukb/FumaIDMapping.csv")

# 通过purrr::map函数读取所有leadSNPs.txt文件并创建一个包含所有数据框的列表
data_frames <- map2(file_paths, folder_numbers, ~ fread(file = .x, header = TRUE) %>%
    mutate(Type = as.character(.y)))

# 利用dplyr::bind_rows函数将所有数据框合并成一个大的数据框
merged_data <- bind_rows(data_frames)

merged_data <- merged_data %>%
    dplyr::select(Type, everything())

merged_data$Type <- rename_based_on_df(merged_data$Type, nmapdf = fuma_dict, from = "job_id", to = "eid")
merged_data$Type <- rename_based_on_df(merged_data$Type, nmapdf = pheno_dict, from = "fid", to = "description")
merged_data <- dplyr::rename(merged_data, Phenotype = "Type")
write.csv(merged_data, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/magma.csv", row.names = FALSE)
merged_data <- merged_data[merged_data$P < 5 * 10e-8, ]
write.csv(merged_data, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/magma_sig.csv", row.names = FALSE)




################## mediation ################################
pheno_dict <- read.csv("scripts/backup/pheno_dict.csv")
data <- read.csv("data/mediation/sandbox_med_top50_sig_with_serum.csv")
data$y <- rename_based_on_df(data$y, nmapdf = pheno_dict, from = "names", to = "description")
feature_dict <- read.csv("Sandbox/Category/feature_sly.csv")
data <- data[data$x %in% feature_dict$Short_name[feature_dict$Select == "TRUE"], ]
data$x <- rename_based_on_df(data$x, nmapdf = feature_dict, from = "Short_name", to = "DispName")
data <- data %>%
    arrange(prop.mediated.p)

molecular_dict <- read.csv("Sandbox/Category/molecular.csv")
data <- data[data$cell_category %in% unique(molecular_dict$category), ]

# 过滤每个y和m对应的最高的x
data <- data %>%
    group_by(y, mediator) %>%
    slice(which.max(abs(prop.mediated))) %>%
    ungroup() %>%
    arrange(prop.mediated.p)


data <- data %>%
    group_split(cell_category)
names(data) <- map_chr(data, ~ unique(.x$cell_category))

writexl::write_xlsx(data, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/mediation.xlsx")






################# r value ##########################
data <- read.csv("data/organ_pheno_organ_pheno/sandbox_r.csv")
data$P_bonf <- p.adjust(data$p, method = "bonferroni")
data <- data[data$P_bonf < 0.01, ]
data <- data[data$organ1_pheno != data$organ2_pheno, ]
data <- data %>% dplyr::arrange(desc(r))
data <- data %>% dplyr::arrange(P_bonf)
data <- data[data$organ1_pheno != data$organ2_pheno, ]
cross <- data[data$organ1 != data$organ2, ]
inter <- data[data$organ1 == data$organ2, ]
list <- list(`cross organ` = cross, `inter organ` = inter)
writexl::write_xlsx(list, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/sandbox_r.xlsx")


ukb_organ_dict <- read.csv("Sandbox/ukb_category/organ.csv")
data <- read.csv("data/organ_pheno_organ_pheno/ukb_r.csv")
data$P_bonf <- p.adjust(data$p, method = "bonferroni")
data <- data[data$P_bonf < 0.01, ]

data <- data %>% dplyr::arrange(desc(r))
data <- data %>% dplyr::arrange(P_bonf)
data <- data[data$organ1_pheno != data$organ2_pheno, ]
data$organ1_pheno <- rename_based_on_df(data$organ1_pheno, nmapdf = ukb_organ_dict, from = "UDI", to = "Label")
data$organ2_pheno <- rename_based_on_df(data$organ2_pheno, nmapdf = ukb_organ_dict, from = "UDI", to = "Label")



cross <- data[data$organ1 != data$organ2, ]
inter <- data[data$organ1 == data$organ2, ]
list <- list(`cross organ` = cross, `inter organ` = inter)
writexl::write_xlsx(list, "C:/Users/liumeng/OneDrive/Desktop/器官轴table/ukb_r.xlsx")
