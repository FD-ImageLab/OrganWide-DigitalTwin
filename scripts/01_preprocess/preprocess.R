# user-defined path
mainPath <- "E:/Projects/indNet/"  # local

# mainPath <- "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
# mainPath <- "D:/seafile/Seafile/SandBox/"  # xumei
# mainPath <- "/Users/apple/Seafile/Research/SandBox"  # wangcy
# mainPath <- "/public/sandbox/workdir/wangchy/SharedFolder"  # sandbox
codePath <- paste(mainPath,"/scripts/01_preprocess",sep="")
dataPath <- paste(mainPath,"/Sandbox/",sep="")

# set path and load data

library(readxl)
library(tidyverse)
library(mice)
library(Hmisc)
source(paste(mainPath,"/scripts/utils/Func.R",sep=""), encoding = "UTF-8")


# brain_vessel_position = list.files(paste0(dataPath, '/brain_vessel_ROI'), full.names = TRUE, recursive = TRUE, pattern = "xls")
# Names = basename(dirname(brain_vessel_position))
# 
# brain_vessel_position1 = brain_vessel_position[1:12]
# brain_vessel_position2 = brain_vessel_position[13:24]
# brain_vessel_position3 = brain_vessel_position[25:36]
# brain_vessel_position4 = brain_vessel_position[37:48]
# brain_vessel_position5 = brain_vessel_position[49:60]
# brain_vessel_position6 = brain_vessel_position[61:72]
# 
# data = map(brain_vessel_position, read_xls, skip = 1, col_names = TRUE, .name_repair = "minimal")
# f = function(df, Name) {
#   df$Name = Name
#   df = df %>%
#     dplyr::select(ID, Name, everything())
#   return(df)
# } 
# data = map2(data, Names, f)
# data = map(data, column_to_rownames, var="ID")
# prefixs = str_extract(basename(brain_vessel_position), ".*?(?=_params)")
# f = function(prefix_in_colnames, df) {
#   colnames(df)[2: ncol(df)] = paste0(prefix_in_colnames, "_", colnames(df)[2: ncol(df)])
#   return(df)
# }
# 
# data = map2(prefixs, data, f) %>%
#   map(rownames_to_column, var="ID")
# f = function(df) {
#   df$Name = NULL
#   return(df)
# }
# data = map(data, f)
# 
# 
# data1 = data[1:12]
# data2 = data[13:24]
# data3 = data[25:36]
# data4 = data[37:48]
# data5 = data[49:60]
# data6 = data[61:72]
# 
# data1 = data1 %>%
#   reduce(full_join)
# data2 = data2 %>%
#   reduce(full_join)
# data3 = data3 %>%
#   reduce(full_join)
# data4 = data4 %>%
#   reduce(full_join)
# data5 = data5 %>%
#   reduce(full_join)
# data6 = data6 %>%
#   reduce(full_join)
# 
# 
# data = list(data1, data2, data3, data4, data5, data6) %>%
#   # map(column_to_rownames, var="ID") %>%
#   reduce(full_join)
# 
# data = data[!(data$ID %in% data$ID[duplicated(data$ID)]), ]
# row.names(data) = NULL

# data = data %>%
#   column_to_rownames(var="ID")

#########################插补数据##############################
# write.csv(data, paste0(dataPath, "/Vessel/CerebralVessel.csv"), row.names = FALSE, fileEncoding = "UTF-8")










##############################NA process########################################



remove_col_by_na_ratio = function(df, na_ratio_threshold=0.5) {
  calc_na_ratio = function(vec) {
    na_ratio = sum(is.na(vec)) / length(vec)
    return(na_ratio)
  }
  na_ratio = map(df, calc_na_ratio)
  df = df[, na_ratio <= na_ratio_threshold]
  return(df)
}
remove_row_by_na_ratio = function(df, na_ratio_threshold=0.5) {
  calc_na_ratio = function(vec) {
    na_ratio = sum(is.na(vec)) / length(vec)
    return(na_ratio)
  }
  na_ratio = apply(X = df, MARGIN = 1, FUN = calc_na_ratio)
  df = df[na_ratio <= na_ratio_threshold, ]
  return(df)
}
# if all features remove 6% data, the number of observers will decrease sharply
# outliner_remove = function(df, ratio=0.03) {
#   rnames = row.names(df)
#   remove_two_side = function(vec, r=ratio) {
#     n = floor(length(vec) * r)
#     vec[order(vec)[c(1: n, (length(vec) - n + 1): length(vec))]] = NA
#     return(vec)
#   }
#   df = map_dfc(df, remove_two_side)
#   rownames(df) = rnames
#   return(df)
# }

outliner_remove = function(df, n_sd = 3) {
  for (col_name in names(df)) {
    # 计算均值和标准差
    mean_val <- mean(df[[col_name]], na.rm = TRUE)
    std_val <- sd(df[[col_name]], na.rm = TRUE)
    
    # 将超过3个标准差的值设为NA
    df[[col_name]][df[[col_name]] > mean_val + n_sd * std_val | df[[col_name]] < mean_val - n_sd * std_val] <- NA

  }
  return(df)
  
}


# pipeline
## this pipeline is used to clean data
# pipeline = function(read_position=read_position, write_position=write_position) {
#   library(tools)
#   if (file_ext(read_position) == "csv") {
#     data = read.csv(read_position, header = TRUE)
#     data = data[!duplicated(data[[1]]),]
#     rownames(data) = data[[1]]
#     data = data[-1]
#     data = scale(data, center = FALSE)
#     data = as.data.frame(data)
#   }
#   if (file_ext(read_position) == "txt") {
#     data = read.table(read_position, header = TRUE)
#     data = data[!duplicated(data[[1]]),]
#     rownames(data) = data[[1]]
#     data = data[-1]
#     data = scale(data, center = FALSE)
#     data = as.data.frame(data)
#   }
#   
#   mid_data = data %>% 
#     # remove_row_by_na_ratio() %>%
#     remove_col_by_na_ratio() %>%
#     outliner_remove()
#   midrname = row.names(mid_data)
#   f = function(df, rname) {row.names(df) = rname; return(df)}
#   final_data = mid_data %>%
#     mice(m = 1, method = "norm", maxit = 5, seed = 500) %>%
#     complete(action = 5) %>%
#     f(midrname) %>%
#     rownames_to_column(var="ID")
#   check_path(write_position)
#   # write.csv(final_data, file = write_position, row.names = FALSE, fileEncoding = "UTF-8")
#   return()
# }

# read_position = list.files(paste0(dataPath, "/Organ/"), pattern = "\\.csv", full.names = TRUE)
# write_postion = paste0("E:/Onedrive/桌面/test/", list.files(paste0(dataPath, "/test/"), pattern = "\\.csv"))
# walk2(read_position, write_postion, pipeline)




# 
# 
# ###先merge到一张大表后进行插补
# read_position = read_position[!(grepl("prostate", read_position) | grepl("uterus", read_position))]
# 
# 
# data = map(read_position, read.csv, header=TRUE)
# f = function(df) {
#   df = df[!duplicated(df[[1]]), ]
#   colnames(df)[1] = "ID"
#   return(df)
# }
# data = map(data, f)
# big_df = reduce(data, full_join, by="ID") %>%
#   column_to_rownames(var="ID") %>%
#   scale(center=FALSE) %>%
#   as.data.frame()
# 
# mid_data = big_df %>% 
#   remove_col_by_na_ratio() %>% 
#   outliner_remove() 
# midrname = row.names(mid_data) 
# 
# 
# 
# f = function(df, rname) {row.names(df) = rname; return(df)}
# final_data = mid_data %>%
#   mice(m = 1, method = "pmm", maxit = 1, seed = 500) %>%
#   complete(action = 1) %>%
#   f(midrname) %>%
#   rownames_to_column(var="ID")
# 
# as.numeric(map(big_df, calc_na_ratio))




############用99列来插补器官数据#############
tt1 = read.csv("scripts/01_preprocess/map-image_240101.csv", header = TRUE)
tt1 = tt1 %>% group_split(part)

for (st1 in tt1) {
  read_position = paste0(mainPath, "/", st1$path)
  write_postion = paste0("C:/Users/liumeng/OneDrive/Desktop/test1/", 
                         basename(dirname(read_position)), "/", 
                         basename(read_position))
  for (i in seq_along(read_position)) {
    # 待插补数据
    t1 = read.csv(read_position[i], header=TRUE, check.names = FALSE)
    t1 = t1[!duplicated(t1[[1]]),]
    colnames(t1)[1] = "ID"
    rownames(t1) = NULL
    t1 = t1 %>% column_to_rownames(var="ID") %>% remove_col_by_na_ratio()
    t1 = t1[names(unlist(map(t1, var, na.rm=TRUE)))[unlist(map(t1, var, na.rm=TRUE)) != 0]]
    
    # 去除cor=1的列
    remove_cor_equal_1_col = function(df) {
      rmat = rcorr(as.matrix(df))$r == 1
      cor_equal_1_indexA = which(rmat * upper.tri(rmat, diag = FALSE) == TRUE)
      if (is_empty(cor_equal_1_indexA)) {
        return(df)
      }
      cor_equal_1_indexB = lapply(cor_equal_1_indexA, FUN = function(x) {indexA2indexB(rmat, x)})
      cor_equal_1_index = unlist(map(cor_equal_1_indexB, ~.x[2]))
      df = df[-cor_equal_1_index]
      return(df)
    }
    t1 = remove_cor_equal_1_col(t1)
    
    t1cn = colnames(t1)
    colnames(t1) = paste0("need_impute_", colnames(t1))
    impute_t1cn = colnames(t1)
    # 用来插补的数据
    read_t2_position = read_position[-i]
    t2 = map(read_t2_position, read.csv, header=TRUE)
    f = function(df, position) {
      df = df[!duplicated(df[[1]]), ]
      colnames(df) = paste0(str_split(basename(position), "\\.")[[1]][1], "_", colnames(df))
      colnames(df)[1] = "ID"
      return(df)
    }
    t2 = map2(t2, read_t2_position, f)
    big_df = reduce(t2, full_join, by="ID") %>%
      column_to_rownames(var="ID") %>% 
      remove_col_by_na_ratio() 
    big_df = big_df[names(unlist(map(big_df, var, na.rm=TRUE)))[unlist(map(big_df, var, na.rm=TRUE)) != 0]]
    # big_df = remove_cor_equal_1_col(big_df)
    
    ## 选择插补的列
    a1 = full_join(t1 %>% rownames_to_column(var = "ID"), 
                   big_df %>% rownames_to_column(var="ID"), 
                   by="ID")
    a2 = a1[colnames(big_df)]
    f = function(vec) {return(sum(!is.na(vec)) / length(vec))}
    tcn = unlist(map(a2, f)) %>%
      sort(decreasing = TRUE)
    tcn = names(tcn[1: min(100, length(tcn))])
    ## 选择100列用于插补的数据
    selected_big_df = big_df[tcn]
    
    # 插补
    mid_df = full_join(t1 %>% rownames_to_column(var = "ID"), 
                       selected_big_df %>% rownames_to_column(var="ID"), 
                       by="ID") %>% 
      column_to_rownames(var="ID")
      # scale(center=FALSE) %>%
      # as.data.frame() 
    midrname = rownames(mid_df)
    temp = matrix(FALSE, nrow=nrow(mid_df), ncol=ncol(mid_df), 
                  dimnames = list(rownames(mid_df), colnames(mid_df)))
    temp[, 1:ncol(t1)] = TRUE
    f = function(df, rname) {row.names(df) = rname; return(df)}
    colnames(mid_df) = make.names(colnames(mid_df))
    temp_mid_df = mid_df %>% outliner_remove()
    
    final_df = temp_mid_df %>% 
      mice(m = 1, method = "cart", maxit = 1, seed = 500, remove_collinear=FALSE) %>% 
      mice::complete(action = 1) %>% 
      f(midrname)
    final_origin_df = final_df[make.names(impute_t1cn)]
    colnames(final_origin_df) = t1cn
    final_origin_df = final_origin_df %>% rownames_to_column(var="ID")
    
    m1 = apply(mid_df[make.names(impute_t1cn)], 1, function(x) {sum(is.na(x))})
    m2 = apply(final_origin_df, 1, function(x) {sum(is.na(x))})
    m1 = table(cut(m1, breaks = round(c(seq(0, 2, 1), seq(3, max(ncol(final_origin_df), 7), length=5))), 
                   right = FALSE, include.lowest = TRUE))
    m2 = table(cut(m2, breaks = round(c(seq(0, 2, 1), seq(3, max(ncol(final_origin_df), 7), length=5))), 
                   right = FALSE, include.lowest = TRUE))
    x = full_join(as.data.frame(m1), as.data.frame(m2), by="Var1") %>% 
      column_to_rownames(var="Var1") %>%
      t()
    rownames(x) = c("before_imputation", "after_imputation")
    
    check_path(paste0("C:/Users/liumeng/OneDrive/Desktop/test1/pics/", str_split(basename(read_position[i]), "\\.")[[1]][1], "_original_table.png"))
    png(paste0("C:/Users/liumeng/OneDrive/Desktop/test1/pics/", str_split(basename(read_position[i]), "\\.")[[1]][1], "_original_table.png"))
    barplot(x, beside = TRUE, col = c('orange', 'steelblue'), 
            legend.text = c("before_imputation", "after_imputation"), 
            ylab = "人数", xlab = "每行NA个数", 
            main = str_split(basename(read_position[i]), "\\.")[[1]][1])
    dev.off()
    
    m1 = apply(mid_df, 1, function(x) {sum(is.na(x))})
    m2 = apply(final_df, 1, function(x) {sum(is.na(x))})
    m1 = table(cut(m1, breaks = round(c(seq(0, 2, 1), seq(3, max(ncol(final_origin_df), 7), length=5))), 
                   right = FALSE, include.lowest = TRUE))
    m2 = table(cut(m2, breaks = round(c(seq(0, 2, 1), seq(3, max(ncol(final_origin_df), 7), length=5))), 
                   right = FALSE, include.lowest = TRUE))
    x = full_join(as.data.frame(m1), as.data.frame(m2), by="Var1") %>% 
      column_to_rownames(var="Var1") %>%
      t()
    rownames(x) = c("before_imputation", "after_imputation")
    check_path(paste0("C:/Users/liumeng/OneDrive/Desktop/test1/pics/", str_split(basename(read_position[i]), "\\.")[[1]][1], "_original_table_add100.png"))
    png(paste0("C:/Users/liumeng/OneDrive/Desktop/test1/pics/", str_split(basename(read_position[i]), "\\.")[[1]][1], "_original_table_add100.png"))
    barplot(x, beside = TRUE, col = c('orange', 'steelblue'), 
            legend.text = c("before_imputation", "after_imputation"), 
            ylab = "人数", xlab = "每行NA个数", 
            main = str_split(basename(read_position[i]), "\\.")[[1]][1])
    dev.off()

    # y = x %>% 
    #   rownames_to_column(var="ID") %>%
    #   pivot_longer(cols = c(before_impute, after_impute))
    # y = xtabs(value ~ NA_num + name, data=y)
    # ggplot(y, aes(x=NA_num, fill=name)) + geom_bar(position = "stack") 
    check_path(write_postion[i])
    write.csv(final_origin_df, write_postion[i], row.names = FALSE)
  }
}













# 
# # read_position = read_position[!(grepl("prostate", read_position) | grepl("uterus", read_position))]
# write_postion = paste0("E:/Onedrive/桌面/test/", basename(read_position))
# 
# for (i in seq_along(read_position)) {
#   t1 = read.csv(read_position[i], header=TRUE, check.names = FALSE)
#   t1 = t1[!duplicated(t1[[1]]),]
#   colnames(t1)[1] = "ID"
#   rownames(t1) = NULL
#   t1 = t1 %>% column_to_rownames(var="ID") %>% remove_col_by_na_ratio()
#   t1 = t1[names(unlist(map(t1, var, na.rm=TRUE)))[unlist(map(t1, var, na.rm=TRUE)) != 0]]
#   t1 = t1 %>% 
#     rownames_to_column(var="ID")
#   t1cn = colnames(t1)
#   
#   read_t2_position = read_position[-i]
#   t2 = map(read_t2_position, read.csv, header=TRUE)
#   f = function(df) {
#     df = df[!duplicated(df[[1]]), ]
#     colnames(df)[1] = "ID"
#     return(df)
#   }
#   t2 = map(t2, f)
#   big_df = reduce(t2, full_join, by="ID") %>%
#     column_to_rownames(var="ID") %>% 
#     remove_col_by_na_ratio() 
#   big_df = big_df[names(unlist(map(big_df, var, na.rm=TRUE)))[unlist(map(big_df, var, na.rm=TRUE)) != 0]]
#   
#   
#   rnt1 = rownames(na.omit(column_to_rownames(t1, var = "ID")))
#   df = big_df[rownames(big_df) %in% rnt1, ]
#   f = function(vec) {return(sum(is.na(vec)) / length(vec))}
#   tcn = unlist(map(df, f)) %>%
#     sort(decreasing = FALSE)
#   tcn = names(tcn[1: 100])
#   selected_big_df = big_df[tcn] %>% 
#     rownames_to_column(var="ID")
#   
#   
#   mid_df = full_join(t1, selected_big_df) %>% 
#     column_to_rownames(var="ID") %>%
#     scale(center=FALSE) %>%
#     as.data.frame() 
#     
#   midrname = rownames(mid_df)
#   temp = matrix(FALSE, nrow=nrow(mid_df), ncol=ncol(mid_df))
#   temp[, 1:ncol(t1)-1] = TRUE
#   f = function(df, rname) {row.names(df) = rname; return(df)}
#   colnames(mid_df) = make.names(colnames(mid_df))
#   
#   final_df = mid_df %>% 
#     mice(m = 1, method = "pmm", maxit = 1, seed = 500, where=is.na(.)&temp, ) %>% 
#     complete(action = 1) %>% 
#     f(midrname) %>% 
#     rownames_to_column(var="ID")
#   final_origin_df = final_df[make.names(t1cn)]
#   colnames(final_origin_df) = t1cn
#   write.csv(final_origin_df, write_postion[i], row.names = FALSE)
# }


# 
# ############用99列来插补血管数据#############
# read_position = list.files(paste0(dataPath, "/Organ/"), pattern = "\\.csv", full.names = TRUE)
# # read_position = read_position[!(grepl("prostate", read_position) | grepl("uterus", read_position))]
# read_position_vessel = list.files(paste0(dataPath, "/Vessel/"), pattern = "\\.csv", full.names = TRUE)
# read_position = unlist(list(read_position_vessel, read_position))
# 
# write_postion = paste0("E:/Onedrive/桌面/vessel/", basename(read_position))
# for (i in seq_along(read_position[1:32])) {
#   t1 = read.csv(read_position[i], header=TRUE, check.names = FALSE)
#   t1 = t1[!duplicated(t1[[1]]),]
#   colnames(t1)[1] = "ID"
#   rownames(t1) = NULL
#   t1 = t1 %>% column_to_rownames(var="ID") %>% remove_col_by_na_ratio()
#   t1 = t1[names(unlist(map(t1, var, na.rm=TRUE)))[unlist(map(t1, var, na.rm=TRUE)) != 0]]
#   t1 = t1 %>% 
#     rownames_to_column(var="ID")
#   t1cn = colnames(t1)
#   
#   read_t2_position = read_position[-i]
#   t2 = map(read_t2_position, read.csv, header=TRUE)
#   f = function(df) {
#     df = df[!duplicated(df[[1]]), ]
#     colnames(df)[1] = "ID"
#     return(df)
#   }
#   t2 = map(t2, f)
#   big_df = reduce(t2, full_join, by="ID") %>%
#     column_to_rownames(var="ID") %>% 
#     remove_col_by_na_ratio() 
#   big_df = big_df[names(unlist(map(big_df, var, na.rm=TRUE)))[unlist(map(big_df, var, na.rm=TRUE)) != 0]]
#   
#   
#   rnt1 = rownames(na.omit(column_to_rownames(t1, var = "ID")))
#   df = big_df[rownames(big_df) %in% rnt1, ]
#   f = function(vec) {return(sum(is.na(vec)) / length(vec))}
#   tcn = unlist(map(df, f)) %>%
#     sort(decreasing = FALSE)
#   tcn = names(tcn[1: 100])
#   selected_big_df = big_df[tcn] %>% 
#     rownames_to_column(var="ID")
#   
#   
#   mid_df = full_join(t1, selected_big_df) %>% 
#     column_to_rownames(var="ID") %>%
#     scale(center=FALSE) %>%
#     as.data.frame() 
#   
#   midrname = rownames(mid_df)
#   temp = matrix(FALSE, nrow=nrow(mid_df), ncol=ncol(mid_df))
#   temp[, 1:ncol(t1)-1] = TRUE
#   f = function(df, rname) {row.names(df) = rname; return(df)}
#   colnames(mid_df) = make.names(colnames(mid_df))
#   
#   final_df = mid_df %>% 
#     mice(m = 1, method = "pmm", maxit = 1, seed = 500, where=is.na(.)&temp, ) %>% 
#     complete(action = 1) %>% 
#     f(midrname) %>% 
#     rownames_to_column(var="ID")
#   final_origin_df = final_df[make.names(t1cn)]
#   colnames(final_origin_df) = t1cn
#   write.csv(final_origin_df, write_postion[i], row.names = FALSE)
# }

############用99列来插补微观数据#############


# read_position = list.files(paste0(dataPath, "/MicroPhenomics/"), pattern = "\\.csv", full.names = TRUE, recursive = TRUE)
# write_postion = paste0("E:/Onedrive/桌面/test1/", 
#                        basename(dirname(read_position[1])), "/", 
#                        basename(read_position))
# 
# library(tools)
# for (i in seq_along(read_position)) {
#   if (file_ext(read_position[i]) == "csv") {
#     t1 = read.csv(read_position[i], header=TRUE, check.names = FALSE)
#     
#   }
#   if (file_ext(read_position[i]) == "txt") {
#     t1 = read.table(read_position[i], header = TRUE, check.names = FALSE, sep = "\t")
#   }
#   t1 = t1[!duplicated(t1[[1]]),]
#   colnames(t1)[1] = "ID"
#   rownames(t1) = NULL
#   t1 = t1 %>% column_to_rownames(var="ID") %>% remove_col_by_na_ratio()
#   t1 = t1[names(unlist(map(t1, var, na.rm=TRUE)))[unlist(map(t1, var, na.rm=TRUE)) != 0]]
#   t1 = t1 %>% 
#     rownames_to_column(var="ID")
#   t1cn = colnames(t1)
#   
#   read_t2_position = read_position[-i]
#   temp_f = function(pos) {
#     if (file_ext(pos) == "csv") {
#       tt = read.csv(pos, header = TRUE, check.names = FALSE)
#     }
#     if (file_ext(pos) == "txt") {
#       tt = read.table(pos, header = TRUE, sep = "\t", check.names = FALSE)
#     }
#     return(tt)
#   }
#   t2 = map(read_t2_position, temp_f)
#   f = function(df) {
#     df = df[!duplicated(df[[1]]), ]
#     colnames(df)[1] = "ID"
#     return(df)
#   }
#   t2 = map(t2, f)
#   big_df = reduce(t2, full_join, by="ID") %>%
#     column_to_rownames(var="ID") %>% 
#     remove_col_by_na_ratio() 
#   big_df = big_df[names(unlist(map(big_df, var, na.rm=TRUE)))[unlist(map(big_df, var, na.rm=TRUE)) != 0]]
#   
#   
#   rnt1 = rownames(na.omit(column_to_rownames(t1, var = "ID")))
#   df = big_df[rownames(big_df) %in% rnt1, ]
#   f = function(vec) {return(sum(is.na(vec)) / length(vec))}
#   tcn = unlist(map(df, f)) %>%
#     sort(decreasing = FALSE)
#   tcn = names(tcn[1: 100])
#   selected_big_df = big_df[tcn] %>% 
#     rownames_to_column(var="ID")
#   
#   
#   mid_df = full_join(t1, selected_big_df) %>% 
#     column_to_rownames(var="ID") %>%
#     scale(center=FALSE) %>%
#     as.data.frame() 
#   
#   midrname = rownames(mid_df)
#   temp = matrix(FALSE, nrow=nrow(mid_df), ncol=ncol(mid_df))
#   temp[, 1:ncol(t1)-1] = TRUE
#   f = function(df, rname) {row.names(df) = rname; return(df)}
#   colnames(mid_df) = make.names(colnames(mid_df))
#   
#   final_df = mid_df %>% 
#     mice(m = 1, method = "pmm", maxit = 1, seed = 500, where=is.na(.)&temp, ) %>% 
#     complete(action = 1) %>% 
#     f(midrname) %>% 
#     rownames_to_column(var="ID")
#   final_origin_df = final_df[make.names(t1cn)]
#   colnames(final_origin_df) = t1cn
#   write.csv(final_origin_df, write_postion[i], row.names = FALSE)
# }




###############################合并心房和心室表格##############################
heart_atrial = read.csv(paste0(mainPath, "/Sandbox/Organ/heart_atrial_volume.csv"), 
                        header = TRUE, check.names = FALSE)

heart_ventricular = read.csv(paste0(mainPath, "Sandbox/Organ/heart_ventricular_volume.csv"), 
                             header = TRUE, check.names = FALSE)
heart = inner_join(heart_atrial, heart_ventricular, by="ID")
write.csv(heart, paste0(mainPath, "/Sandbox/Organ/heart_volume.csv"), row.names = FALSE)
